#include "AleOptimizer.hpp" 
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <IO/IO.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <util/Paths.hpp>
#include <maths/Random.hpp>
#include <search/SpeciesSPRSearch.hpp>
#include <search/SpeciesTransferSearch.hpp>
#include <search/DatedSpeciesTreeSearch.hpp>



static bool testAndSwap(size_t &hash1, size_t &hash2) {
  std::swap(hash1, hash2);
  return hash1 != hash2;
}

class HighwayFunction: public FunctionToOptimize {
public: 
  HighwayFunction(AleEvaluator &evaluator,
    const std::vector<Highway*> &highways): _highways(highways), _evaluator(evaluator) {}
  
  virtual double evaluate(Parameters &parameters) {
    return evaluatePrint(parameters, false); 
  }
  
  virtual double evaluatePrint(Parameters &parameters, bool print) {
    assert(parameters.dimensions() == _highways.size());
    parameters.ensurePositivity();
    for (unsigned int i = 0; i < _highways.size(); ++i)  {
      Highway highwayCopy = *_highways[i];
      highwayCopy.proba = parameters[i];
      _evaluator.addHighway(highwayCopy);
    }
    auto ll = _evaluator.computeLikelihood();
    if (print) {
      std::string out = FileSystem::joinPaths(_evaluator.getOutputDir(),
        std::string("transferll_") + std::string(_highways[0]->src->label) + 
        std::string("_") + std::string(_highways[0]->dest->label));
      _evaluator.savePerFamilyLikelihoodDiff(out);
      Logger::info << "Saving diff with p=" << parameters[0] 
        << " into " << out << std::endl;
    }
    for (auto highway: _highways) {
      (void)(highway);
      _evaluator.removeHighway();
    }
    parameters.setScore(ll);
    return ll;
  }
private:
  const std::vector<Highway *> &_highways;
  AleEvaluator &_evaluator;
};


static void getSpeciesToCatRec(corax_rnode_t *node,
    std::vector<unsigned int> &speciesToCat,
    unsigned int currentCat,
    const std::map<std::string, unsigned int> &labelToCat)
{
  std::string label = node->label;
  auto it = labelToCat.find(label);
  if (it != labelToCat.end()) {
    currentCat = it->second;
  }
  speciesToCat[node->node_index] = currentCat;
  if (node->left) {
    getSpeciesToCatRec(node->left, speciesToCat, currentCat, labelToCat);
    getSpeciesToCatRec(node->right, speciesToCat, currentCat, labelToCat);
  }
}

void  getSpeciesToCat(const PLLRootedTree &speciesTree,
    const std::string &speciesCategoryFile,
    std::vector<unsigned int> &speciesToCat,
    std::vector<std::string> &catToLabel)
{
  unsigned int N = speciesTree.getNodeNumber();
  speciesToCat = std::vector<unsigned int>(N, 0);
  catToLabel.push_back("all");
  std::ifstream is(speciesCategoryFile);
  if (!is) {
    return;
  }
  std::map<std::string, unsigned int> labelToCat;
  std::string line;
  auto allLabels = speciesTree.getLabels(false);
  while (std::getline(is, line)) {
    IO::removeSpaces(line);
    if (line.size() == 0 || labelToCat.find(line) != labelToCat.end()) {
      continue;
    }
    if (allLabels.find(line) == allLabels.end()) {
      Logger::error << "Warning, label " << line << " from " << 
        speciesCategoryFile << " is not in the species tree" << std::endl;
      continue;
    }
    unsigned int cat = labelToCat.size() + 1;
    labelToCat.insert({line, cat});
    catToLabel.push_back(line);
  }
  getSpeciesToCatRec(speciesTree.getRoot(),
      speciesToCat,
      0,
      labelToCat);
  
}

AleOptimizer::AleOptimizer(
    const std::string speciesTreeFile, 
    const Families &families, 
    const RecModelInfo &info,
    bool optimizeRates,
    bool optimizeVerbose,
    const std::string &speciesCategoryFile,
    const std::string &outputDir):
  _speciesTree(std::make_unique<SpeciesTree>(speciesTreeFile)),
  _geneTrees(families, false, true),
  _info(info),
  _outputDir(outputDir),
  _searchState(*_speciesTree,
      Paths::getSpeciesTreeFile(_outputDir, "inferred_species_tree.newick"),
      _geneTrees.getTrees().size()),
  _rootLikelihoods(_geneTrees.getTrees().size())
{
  Parameters startingRates;
  switch (info.model) {
  case RecModel::UndatedDL:
    startingRates = Parameters(0.2, 0.2);
    break;
  case RecModel::UndatedDTL:
    startingRates = Parameters(0.1, 0.2, 0.1);
    break;
  default:
    assert(false);
    break;
  }
  std::vector<unsigned int> speciesToCat;
  std::vector<std::string> catToLabel;
  getSpeciesToCat(_speciesTree->getTree(),
      speciesCategoryFile,
      speciesToCat,
      catToLabel);
  _modelRates = AleModelParameters(startingRates, 
      speciesToCat, 
      catToLabel,
      info);
  _speciesTree->addListener(this);
  ParallelContext::barrier();
  _evaluator = std::make_unique<AleEvaluator>(
      *_speciesTree, 
      _modelRates, 
      optimizeRates,
      optimizeVerbose,
      families, 
      _geneTrees,
      _outputDir);
  Logger::timed << "Initial ll=" << getEvaluator().computeLikelihood() 
    << std::endl;
  saveCurrentSpeciesTreeId("starting_species_tree.newick");
  saveCurrentSpeciesTreeId();
}
  
double AleOptimizer::optimizeModelRates(bool thorough)
{
  getEvaluator().optimizeModelRates(thorough);
  PerFamLL perFamLL;
  auto ll = getEvaluator().computeLikelihood(&perFamLL);
  _searchState.betterLikelihoodCallback(ll, perFamLL);
  return ll;
}

void AleOptimizer::optimize()
{
  size_t hash1 = 0;
  size_t hash2 = 0;
  unsigned int index = 0;
  PerFamLL initialPerFamLL;
  _searchState.bestLL = getEvaluator().computeLikelihood(&initialPerFamLL);
  _searchState.farFromPlausible = true;
  _searchState.betterTreeCallback(_searchState.bestLL, initialPerFamLL);
  /**
   *  Alternate transfer search and normal
   *  SPR search, until one does not find
   *  a better tree. Run each at least once.
   */
  rootSearch(3);
  do {
    if (index++ % 2 == 0) {
      transferSearch();
    } else {
      sprSearch(1);
    }
    if (!_searchState.farFromPlausible) {
      rootSearch(3);
    }
    hash1 = _speciesTree->getHash();
  }
  while(testAndSwap(hash1, hash2));
  rootSearch(-1);
}

double AleOptimizer::sprSearch(unsigned int radius)
{
  SpeciesSPRSearch::SPRSearch(*_speciesTree,
      getEvaluator(),
      _searchState,
      radius);
  Logger::timed << "After normal search: LL=" << _searchState.bestLL << std::endl;
  saveSupportTree();
  return _searchState.bestLL;
}

void AleOptimizer::saveSupportTree()
{
  auto outKH = Paths::getSpeciesTreeFile(_outputDir, 
        "species_tree_support_kh.newick");
  _searchState.saveSpeciesTreeKH(outKH);
  auto outBP = Paths::getSpeciesTreeFile(_outputDir, 
        "species_tree_support_bp.newick");
  _searchState.saveSpeciesTreeBP(outBP);
}

void AleOptimizer::onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *nodesToInvalidate)
{
  getEvaluator().onSpeciesTreeChange(nodesToInvalidate);
}


void AleOptimizer::saveSpeciesTree()
{
  saveCurrentSpeciesTreeId();
}

std::string AleOptimizer::saveCurrentSpeciesTreeId(std::string name, bool masterRankOnly)
{
  std::string res = Paths::getSpeciesTreeFile(_outputDir, name);
  if (_evaluator->isDated()) {
    _speciesTree->getDatedTree().rescaleBranchLengths();
  }
  saveCurrentSpeciesTreePath(res, masterRankOnly);
  if (!_rootLikelihoods.isEmpty()) {
    auto newick = _speciesTree->getTree().getNewickString();
    PLLRootedTree tree(newick, false); 
    _rootLikelihoods.fillTree(tree);
    auto out = Paths::getSpeciesTreeFile(_outputDir, 
        "species_tree_llr.newick");
    tree.save(out);
  }
  if (!_rootLikelihoods.isEmpty()) {
    auto newick = _speciesTree->getTree().getNewickString();
    PLLRootedTree tree(newick, false); 
    _rootLikelihoods.fillTreeBootstraps(tree);
    auto out = Paths::getSpeciesTreeFile(_outputDir, 
        "species_tree_root_rell.newick");
    tree.save(out);
  }
  return res;
}

void AleOptimizer::saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly)
{
  _speciesTree->saveToFile(str, masterRankOnly);
  if (masterRankOnly) {
    ParallelContext::barrier();
  }
}
  
double AleOptimizer::rootSearch(unsigned int maxDepth, bool thorough)
{
  _rootLikelihoods.reset();
  if (thorough) {
    _searchState.farFromPlausible = thorough;
  }
  SpeciesRootSearch::rootSearch(
      *_speciesTree,
      getEvaluator(),
      _searchState,
      maxDepth,
      &_rootLikelihoods);
  
  return _searchState.bestLL;
}
  
double AleOptimizer::transferSearch()
{
  SpeciesTransferSearch::transferSearch(
    *_speciesTree,
    getEvaluator(),
    _searchState);
  Logger::timed << "After normal search: LL=" 
    << _searchState.bestLL << std::endl;
  saveSupportTree();
  return _searchState.bestLL;
}
 

std::vector<std::string> getLines(const std::string path)
{
  std::ifstream is(path);
  std::string line;
  std::vector<std::string> res;
  while (std::getline(is, line)) {
    if (line.size() > 2) {
      res.push_back(line);
    }
  }
  return res;
}


static void saveStr(const std::string &str, const std::string &path) {
  std::ofstream os(path);
  os << str;
}

void AleOptimizer::reconcile(unsigned int samples)
{
  if (samples == 0) {
    return;
  }
  auto recDir = FileSystem::joinPaths(_outputDir, "reconciliations");
  FileSystem::mkdir(recDir, true);
  auto allRecDir = FileSystem::joinPaths(recDir, "all");
  FileSystem::mkdir(allRecDir, true);
  auto summariesDir = FileSystem::joinPaths(recDir, "summaries");
  FileSystem::mkdir(summariesDir, true);
  ParallelContext::barrier();
  auto &families = _geneTrees.getTrees();
  std::vector<std::string> summaryPerSpeciesEventCountsFiles;
  std::vector<std::string> summaryTransferFiles;
  for (unsigned int i = 0; i < families.size(); ++i) {
    std::vector<std::string> perSpeciesEventCountsFiles;
    std::vector<std::string> transferFiles;
    std::string geneTreesPath = FileSystem::joinPaths(allRecDir, families[i].name + std::string(".newick"));
    ParallelOfstream geneTreesOs(geneTreesPath, false);
    std::vector< std::shared_ptr<Scenario> > scenarios;
    _evaluator->sampleScenarios(i, samples, scenarios);
    assert(scenarios.size() == samples);
    for (unsigned int sample = 0; sample < samples; ++ sample) {
      auto out = FileSystem::joinPaths(allRecDir, 
          families[i].name + std::string("_") +  std::to_string(sample) + ".xml");
      auto eventCountsFile = FileSystem::joinPaths(allRecDir, 
          families[i].name + std::string("_eventcount_") +  std::to_string(sample) + ".txt");
      auto perSpeciesEventCountsFile = FileSystem::joinPaths(allRecDir, 
          families[i].name + std::string("_perspecies_eventcount_") +  std::to_string(sample) + ".txt");
      auto transferFile = FileSystem::joinPaths(allRecDir, 
          families[i].name + std::string("_transfers_") +  std::to_string(sample) + ".txt");
      perSpeciesEventCountsFiles.push_back(perSpeciesEventCountsFile);
      transferFiles.push_back(transferFile);
      auto &scenario = *scenarios[sample];
      scenario.saveReconciliation(out, ReconciliationFormat::RecPhyloXML, false);
      scenario.saveReconciliation(geneTreesOs, ReconciliationFormat::NewickEvents);
      scenario.saveEventsCounts(eventCountsFile, false);
      scenario.savePerSpeciesEventsCounts(perSpeciesEventCountsFile, false);
      scenario.saveTransfers(transferFile, false);
      geneTreesOs <<  "\n";
    }
    geneTreesOs.close();
    auto newicks = getLines(geneTreesPath);
    auto consensusPrefix = FileSystem::joinPaths(summariesDir, 
        families[i].name + "_consensus_");
    saveStr(PLLUnrootedTree::buildConsensusTree(newicks, 0.1), consensusPrefix + "mre.newick");
    //saveStr(PLLUnrootedTree::buildConsensusTree(newicks, 0.5), consensusPrefix + "50.newick");
    auto perSpeciesEventCountsFile = FileSystem::joinPaths(summariesDir, families[i].name + 
        std::string("_perspecies_eventcount.txt"));
    Scenario::mergePerSpeciesEventCounts(_speciesTree->getTree(),
        perSpeciesEventCountsFile,
        perSpeciesEventCountsFiles, false, true);
    summaryPerSpeciesEventCountsFiles.push_back(perSpeciesEventCountsFile);
    auto transferFile = FileSystem::joinPaths(summariesDir, families[i].name + 
        std::string("_transfers.txt"));
    Scenario::mergeTransfers(_speciesTree->getTree(),
        transferFile,
        transferFiles, false, true);
    summaryTransferFiles.push_back(transferFile);
  }
  auto totalPerSpeciesEventCountsFile = FileSystem::joinPaths(recDir, "perspecies_eventcount.txt");
  Scenario::mergePerSpeciesEventCounts(_speciesTree->getTree(),
      totalPerSpeciesEventCountsFile, 
      summaryPerSpeciesEventCountsFiles, 
      true, false);
  auto totalTransferFile = FileSystem::joinPaths(recDir, "transfers.txt");
  Scenario::mergeTransfers(_speciesTree->getTree(),
      totalTransferFile, 
      summaryTransferFiles, 
      true, false);
  ParallelContext::makeRandConsistent();
  Logger::timed << "Reconciliations output directory: " << recDir << std::endl;
}
 

void AleOptimizer::optimizeDates(bool thorough)
{
  if (!_info.isDated()) {
    return; 
  }
  auto scoredBackups = DatedSpeciesTreeSearch::optimizeDatesFromReconciliation(*_speciesTree, getEvaluator(), 200);
  Logger::timed << "Sorted dating likelihoods from fast datings:" << std::endl;
  for (auto &sb: scoredBackups) {
    Logger::info << sb.score << " ";
  }
  Logger::info << std::endl;
 

  auto bestLL = scoredBackups[0].score;
  _searchState.bestLL = bestLL;
  _speciesTree->getDatedTree().restore(scoredBackups[0].backup);
  DatedSpeciesTreeSearch::optimizeDates(*_speciesTree,
      getEvaluator(),
      _searchState,
      bestLL,
      true);

  saveCurrentSpeciesTreeId();
}



void AleOptimizer::randomizeRoot()
{
  auto &tree = _speciesTree->getDatedTree();
  unsigned int N = tree.getOrderedSpeciations().size();
  for (unsigned int i = 0; i < N; ++i) {
    auto direction = Random::getInt() % 4;
    if (SpeciesTreeOperator::canChangeRoot(*_speciesTree, direction)) {
      SpeciesTreeOperator::changeRoot(*_speciesTree, direction);
    }
  }
}

static Parameters testHighwayFast(AleEvaluator &evaluator,
    const Highway &highway,
    double startingProbability = 0.1)
{
  std::vector<Highway *> highways;
  auto copy = highway;
  highways.push_back(&copy);
  HighwayFunction f(evaluator, highways);
  Parameters parameters(1);
  parameters[0] = startingProbability;
  f.evaluate(parameters);
  return parameters;
}
static Parameters testHighway(AleEvaluator &evaluator,
    Highway &highway,
    double startingProbability = 0.01)
{
  std::vector<Highway *> highways;
  highways.push_back(&highway);
  HighwayFunction f(evaluator, highways);
  Parameters startingParameter(1);
  startingParameter[0] = startingProbability;
  OptimizationSettings settings;
  settings.minAlpha = 0.0001;
  settings.epsilon = 0.00001;
  Logger::info << "Unoptimized ll=" << f.evaluatePrint(startingParameter, true) << std::endl;
  auto parameters = DTLOptimizer::optimizeParameters(
      f, 
      startingParameter, 
      settings);
  return parameters;
}

static Parameters testHighways(AleEvaluator &evaluator,
    const std::vector<Highway *> &highways,
    const Parameters &startingProbabilities,
    bool optimize)
{
  assert(highways.size() == startingProbabilities.dimensions());
  HighwayFunction f(evaluator, highways);
  if (optimize) {
    OptimizationSettings settings;
    settings.lineSearchMinImprovement = 1.0;
    settings.minAlpha = 0.001;
    settings.epsilon = 0.000001;
    settings.verbose = true;
    return DTLOptimizer::optimizeParameters(
        f, 
        startingProbabilities, 
        settings);
  } else {
    auto parameters = startingProbabilities;
    f.evaluate(parameters);
    return parameters;
  }
}

void AleOptimizer::getCandidateHighways(std::vector<ScoredHighway> &scoredHighways, unsigned int toTest)
{
  unsigned int minTransfers = 1;
  MovesBlackList blacklist;
  std::vector<TransferMove> transferMoves;
  SpeciesTransferSearch::getSortedTransferList(*_speciesTree,
    getEvaluator(),
    minTransfers,
    blacklist, 
    transferMoves);
  
  for (const auto &transferMove: transferMoves) {
    auto prune = _speciesTree->getNode(transferMove.prune); 
    auto regraft = _speciesTree->getNode(transferMove.regraft);
    Highway highway(regraft, prune);
    scoredHighways.push_back(ScoredHighway(highway, 0.0));
    if (scoredHighways.size() >= toTest) {
      break;
    }
  }
}

bool isHighwayCompatible(Highway &highway,
    const RecModelInfo &info,
    const DatedTree &tree)
{
  auto from = highway.src;
  auto to = highway.dest;
  switch (info.transferConstraint) {
  case TransferConstaint::NONE:
    return true;
  case TransferConstaint::PARENTS:
    while (from) {
      if (to == from) {
        return false;
      }
      from = from->parent;
    }
    return true;
  case TransferConstaint::RELDATED:
    return tree.canTransferUnderRelDated(from->node_index,
        to->node_index);
  }
  assert(false);
  return false;    
}

void AleOptimizer::filterCandidateHighwaysFast(const std::vector<ScoredHighway> &highways, std::vector<ScoredHighway> &filteredHighways)
{
  double proba = 0.01;
  Logger::timed << "Filering " << highways.size() << " candidate highways using p=" << proba << std::endl;
  double initialLL = getEvaluator().computeLikelihood(); 
  Logger::timed << "initial ll=" << initialLL << std::endl;
  for (const auto &scoredHighway: highways) {
    auto highway = scoredHighway.highway;
    if (!isHighwayCompatible(highway, _info, _speciesTree->getDatedTree())) {
      Logger::info << "Incompatible highway " << highway.src->label << "->" << highway.dest->label << std::endl;
      continue;
    }
    auto parameters = testHighwayFast(*_evaluator, highway, proba);
    auto llDiff = parameters.getScore() - initialLL;
    if (llDiff > 0.01) {
      Logger::timed << "Accepting candidate: ";
      highway.proba = parameters[0];
      filteredHighways.push_back(ScoredHighway(highway, 
            -llDiff));
    } else {
      Logger::timed << "Rejecting candidate: ";
    }
    Logger::info << initialLL << " " << parameters.getScore() << std::endl;
    Logger::info << highway.src->label << "->" << highway.dest->label << std::endl;
    Logger::info << " ll diff = " << llDiff << std::endl; 
  }
  std::sort(filteredHighways.begin(), filteredHighways.end());
}


void AleOptimizer::selectBestHighways(const std::vector<ScoredHighway> &highways, 
    std::vector<ScoredHighway> &bestHighways)
{
  Logger::timed << "Looking for the best highways candidates among " << highways.size() << " candidates (slow round)" << std::endl;
  double initialLL = getEvaluator().computeLikelihood(); 
  Logger::info << "initial ll=" << initialLL << std::endl;
  _evaluator->saveSnapshotPerFamilyLL();
  for (const auto &scoredHighway: highways) {
    Highway highway(scoredHighway.highway);
    Logger::timed << "Optimizing candidate highway " << highway.src->label << "->" << highway.dest->label << std::endl;
    auto parameters = testHighway(*_evaluator, highway);
    highway.proba = parameters[0];
    bestHighways.push_back(ScoredHighway(highway, 
          parameters.getScore(),
          initialLL - parameters.getScore()));
    Logger::timed << "ph=" << parameters[0] 
      << " lldiff=" << initialLL - parameters.getScore() << std::endl;
  }
  std::sort(bestHighways.rbegin(), bestHighways.rend());
}

void AleOptimizer::addHighways(const std::vector<ScoredHighway> &candidateHighways,
    std::vector<ScoredHighway> &acceptedHighways)
{
  Logger::timed << "Trying to add all candidate highways simultaneously" << std::endl;
  std::vector<Highway> highways;
  std::vector<Highway *> highwaysPtr;
  Parameters startingProbabilities;
  for (const auto candidate: candidateHighways) {
    highways.push_back(candidate.highway);
    startingProbabilities.addValue(candidate.highway.proba);
  }
  for (auto &highway: highways) {
    highwaysPtr.push_back(&highway);
  }
  auto parameters = testHighways(*_evaluator, highwaysPtr, startingProbabilities, true);
  Logger::info << parameters << std::endl;
  for (unsigned int i = 0; i < candidateHighways.size(); ++i) {
    ScoredHighway sh(candidateHighways[i]);
    sh.highway.proba = parameters[i];
    acceptedHighways.push_back(sh);
    _evaluator->addHighway(sh.highway);
  }
  std::sort(acceptedHighways.rbegin(), acceptedHighways.rend());
}

void AleOptimizer::saveBestHighways(const std::vector<ScoredHighway> &scoredHighways,
      const std::string &output)
{
  ParallelOfstream os(output, true);
  Logger::info << "Outputing the " << scoredHighways.size() << 
    " highays into " << output << std::endl;
  for (const auto &scoredHighway: scoredHighways) {
    os << scoredHighway.highway.proba << ", ";
    os << scoredHighway.score << ", ";
    os << scoredHighway.scoreDiff << ", ";
    os << scoredHighway.highway.src->label << ",";
    os << scoredHighway.highway.dest->label << std::endl;
  }
}
  

