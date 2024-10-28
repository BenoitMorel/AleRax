#include "AleOptimizer.hpp"

#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <IO/ParallelOfstream.hpp>
#include <maths/Random.hpp>
#include <parallelization/ParallelContext.hpp>
#include <search/DatedSpeciesTreeSearch.hpp>
#include <search/SpeciesRootSearch.hpp>
#include <search/SpeciesSPRSearch.hpp>
#include <search/SpeciesTransferSearch.hpp>
#include <util/Paths.hpp>


static bool testAndSwap(size_t &hash1, size_t &hash2) {
  std::swap(hash1, hash2);
  return hash1 != hash2;
}

AleOptimizer::AleOptimizer(
    const std::string &speciesTreeFile,
    const Families &families,
    const RecModelInfo &info,
    const ModelParametrization &modelParametrization,
    const std::string &optimizationClassFile,
    const Parameters &startingRates,
    bool optimizeRates,
    bool optimizeVerbose,
    const std::string &outputDir):
  _state(speciesTreeFile),
  _families(families), // global families
  _geneTrees(families, false, true), // init local families
  _info(info),
  _rootLikelihoods(getLocalFamilyNumber()),
  _outputDir(outputDir),
  _checkpointDir(getCheckpointDir(outputDir)),
  _enableCheckpoints(true)
{
  // init the current state (from checkpoint or de novo)
  if (checkpointExists()) {
    Logger::info << "Loading checkpoint information..." << std::endl;
    loadCheckpoint();
  } else {
    for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
      _state.localFamilyNames.push_back(_geneTrees.getTrees()[i].name);
      _state.perLocalFamilyModelParams.push_back(AleModelParameters(startingRates,
          getSpeciesTree().getTree().getNodeNumber()));
    }
  }
  // init the species tree interface
  _speciesTreeSearchState = std::make_unique<SpeciesSearchState>(
      getSpeciesTree(),
      speciesTreeFile,
      getLocalFamilyNumber());
  _speciesTreeSearchState->addListener(this);
  getSpeciesTree().addListener(this);
  ParallelContext::barrier();
  // init the evaluator
  _evaluator = std::make_unique<AleEvaluator>(
      *this,
      getSpeciesTree(),
      getRecModelInfo(),
      modelParametrization,
      optimizationClassFile,
      getMixtureAlpha(),
      getModelParameters(),
      getTransferHighways(),
      optimizeRates,
      optimizeVerbose,
      families,
      _geneTrees,
      _outputDir);
  Logger::timed << "Initial ll=" << getEvaluator().computeLikelihood() << std::endl;
}

void AleOptimizer::randomizeRoot()
{
  assert(getCurrentStep() == AleStep::SpeciesTreeOpt);
  auto &tree = getSpeciesTree().getDatedTree();
  unsigned int N = tree.getOrderedSpeciations().size();
  for (unsigned int i = 0; i < N; ++i) {
    auto direction = Random::getInt() % 4;
    if (SpeciesTreeOperator::canChangeRoot(getSpeciesTree(), direction)) {
      SpeciesTreeOperator::changeRoot(getSpeciesTree(), direction);
    }
  }
}

void AleOptimizer::optimize()
{
  assert(getCurrentStep() == AleStep::SpeciesTreeOpt);
  size_t hash1 = 0;
  size_t hash2 = 0;
  unsigned int index = 0;
  PerFamLL initialPerFamLL;
  auto ll = getEvaluator().computeLikelihood(&initialPerFamLL);
  _speciesTreeSearchState->bestLL = ll;
  _speciesTreeSearchState->farFromPlausible = true;
  _speciesTreeSearchState->betterTreeCallback(ll, initialPerFamLL);
  // Alternate transfer search and normal SPR search,
  // until one cannot find a better tree.
  // Run each at least once.
  Logger::info << std::endl;
  rootSearch(3);
  do {
    if (index++ % 2 == 0) {
      Logger::info << std::endl;
      transferSearch();
    } else {
      Logger::info << std::endl;
      sprSearch(1);
    }
    if (!_speciesTreeSearchState->farFromPlausible) {
      Logger::info << std::endl;
      rootSearch(3);
    }
    hash1 = getSpeciesTree().getHash();
  }
  while(testAndSwap(hash1, hash2));
  // Final root search
  Logger::info << std::endl;
  rootSearch(-1);
}

void AleOptimizer::reroot()
{
  assert(getCurrentStep() == AleStep::SpeciesTreeOpt);
  PerFamLL initialPerFamLL;
  auto ll = getEvaluator().computeLikelihood(&initialPerFamLL);
  _speciesTreeSearchState->bestLL = ll;
  _speciesTreeSearchState->farFromPlausible = false;
  _speciesTreeSearchState->betterTreeCallback(ll, initialPerFamLL);
  // Optimize rates and run a local root search.
  // If the species tree is dated, run the second
  // local root search with thorough date optimization
  Logger::info << std::endl;
  optimizeModelRates(true);
  if (_info.isDated()) {
    Logger::info << std::endl;
    Logger::timed << "First root search, non-thorough dating" << std::endl;
    rootSearch(5, false);
    Logger::info << std::endl;
    Logger::timed << "Second root search, thorough dating" << std::endl;
    rootSearch(2, true);
  } else {
    Logger::info << std::endl;
    rootSearch(5);
  }
}

double AleOptimizer::rootSearch(unsigned int maxDepth, bool thorough)
{
  _rootLikelihoods.reset();
  if (thorough) {
    _speciesTreeSearchState->farFromPlausible = true;
  }
  SpeciesRootSearch::rootSearch(
      getSpeciesTree(),
      getEvaluator(),
      *_speciesTreeSearchState,
      maxDepth,
      &_rootLikelihoods);
  auto ll = _speciesTreeSearchState->bestLL;
  saveCheckpoint();
  saveSpeciesRootSupports();
  return ll;
}

double AleOptimizer::transferSearch()
{
  SpeciesTransferSearch::transferSearch(
      getSpeciesTree(),
      getEvaluator(),
      *_speciesTreeSearchState);
  auto ll = _speciesTreeSearchState->bestLL;
  Logger::timed << "After normal search: LL=" << ll << std::endl;
  saveCheckpoint();
  return ll;
}

double AleOptimizer::sprSearch(unsigned int radius)
{
  SpeciesSPRSearch::SPRSearch(
      getSpeciesTree(),
      getEvaluator(),
      *_speciesTreeSearchState,
      radius);
  auto ll = _speciesTreeSearchState->bestLL;
  Logger::timed << "After normal search: LL=" << ll << std::endl;
  saveCheckpoint();
  saveSpeciesBranchSupports();
  return ll;
}

double AleOptimizer::optimizeModelRates(bool thorough)
{
  getEvaluator().optimizeModelRates(thorough);
  PerFamLL perFamLL;
  auto ll = getEvaluator().computeLikelihood(&perFamLL);
  _speciesTreeSearchState->betterLikelihoodCallback(ll, perFamLL);
  saveCheckpoint();
  return ll;
}

void AleOptimizer::optimizeDates(bool thorough)
{
  assert(getCurrentStep() == AleStep::RelDating);
  if (!_info.isDated()) {
    return;
  }
  auto initialLL = _speciesTreeSearchState->bestLL;
  unsigned int datingsNumber = 500; // number of independent searches from random datings
  unsigned int backupsNumber = 20; // max number of best datings to continue with
  Logger::timed << "[Species search] Inferring speciation order for " << datingsNumber
                << " random datings" << std::endl;
  auto scoredBackups = DatedSpeciesTreeSearch::optimizeDatesFromReconciliation(
      getSpeciesTree(),
      getEvaluator(),
      datingsNumber,
      backupsNumber);
  Logger::timed << "[Species search] Best " << backupsNumber
                << " dating likelihoods from random datings:" << std::endl;
  for (const auto &sb: scoredBackups) {
    Logger::info << sb.score << " ";
  }
  Logger::info << std::endl;
  auto bestLL = scoredBackups[0].score;
  auto llDiff = bestLL - initialLL;
  if (llDiff > 0.0) {
    Logger::timed << "[Species search] Accept improvement from random datings: "
                  << "bestLL=" << bestLL << ", llDiff=" << llDiff << std::endl;
    _speciesTreeSearchState->bestLL = bestLL;
    getSpeciesTree().getDatedTree().restore(scoredBackups[0].backup);
  } else {
    Logger::timed << "[Species search] No improvement from random datings: "
                  << "bestLL=" << bestLL << ", llDiff=" << llDiff << std::endl;
  }
  Logger::timed << "[Species search] Optimizing speciation order with random dating perturbations" << std::endl;
  DatedSpeciesTreeSearch::optimizeDates(getSpeciesTree(),
      getEvaluator(),
      *_speciesTreeSearchState,
      thorough);
  saveCheckpoint();
  saveSpeciesTree();
  Logger::timed << "Species tree with ordered speciations: "
                << _speciesTreeSearchState->pathToBestSpeciesTree << std::endl;
}

void AleOptimizer::onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *nodesToInvalidate)
{
  getEvaluator().onSpeciesTreeChange(nodesToInvalidate);
}

void AleOptimizer::betterTreeCallback()
{
  saveCheckpoint();
}

void AleOptimizer::onBetterParametersFoundCallback()
{
  if (_enableCheckpoints) {
    saveCheckpoint();
  }
}

void AleOptimizer::saveSpeciesTree()
{
  if (_info.isDated()) {
    getSpeciesTree().getDatedTree().rescaleBranchLengths();
  }
  auto out = _speciesTreeSearchState->pathToBestSpeciesTree;
  getSpeciesTree().saveToFile(out, true);
  Logger::info << "save species tree to " << out << std::endl;
  ParallelContext::barrier();
}

void AleOptimizer::saveSpeciesBranchSupports()
{
  auto outKH = Paths::getSpeciesTreeFile(_outputDir,
      "species_tree_support_kh.newick");
  _speciesTreeSearchState->saveSpeciesTreeKH(outKH);
  //Logger::info << "save KH branch support tree to " << outKH << std::endl;
  auto outBP = Paths::getSpeciesTreeFile(_outputDir,
      "species_tree_support_bp.newick");
  _speciesTreeSearchState->saveSpeciesTreeBP(outBP);
  //Logger::info << "save BP branch support tree to " << outBP << std::endl;
}

void AleOptimizer::saveSpeciesRootSupports()
{
  assert(!_rootLikelihoods.isEmpty());
  auto newick = getSpeciesTree().getTree().getNewickString();
  PLLRootedTree treeLLR(newick, false);
  _rootLikelihoods.fillTree(treeLLR);
  auto outLLR = Paths::getSpeciesTreeFile(_outputDir,
      "species_tree_root_llr.newick");
  treeLLR.save(outLLR);
  //Logger::info << "save LLR root support tree to " << outLLR << std::endl;
  PLLRootedTree treeRELL(newick, false);
  _rootLikelihoods.fillTreeBootstraps(treeRELL);
  auto outRELL = Paths::getSpeciesTreeFile(_outputDir,
      "species_tree_root_rell.newick");
  treeRELL.save(outRELL);
  //Logger::info << "save RELL root support tree to " << outRELL << std::endl;
}

void AleOptimizer::saveRatesAndLL()
{
  // save per-family likelihoods
  auto perFamilyLikelihoodPath = FileSystem::joinPaths(_outputDir, "per_fam_likelihoods.txt");
  std::vector<unsigned int> localIndices;
  std::vector<double> localLikelihoods;
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    auto famIndex = _geneTrees.getTrees()[i].familyIndex;
    auto ll = _evaluator->computeFamilyLikelihood(i);
    localIndices.push_back(famIndex);
    localLikelihoods.push_back(ll);
  }
  ParallelContext::barrier();
  std::vector<unsigned int> indices;
  std::vector<double> likelihoods;
  ParallelContext::concatenateHetherogeneousUIntVectors(localIndices, indices);
  ParallelContext::concatenateHetherogeneousDoubleVectors(localLikelihoods, likelihoods);
  assert(indices.size() == _families.size());
  if (ParallelContext::getRank() == 0) {
    std::vector<ScoredFamily> familiesAndLLs;
    for (unsigned int i = 0; i < indices.size(); ++i) {
      familiesAndLLs.push_back(ScoredFamily(_families[indices[i]].name, likelihoods[i]));
    }
    std::sort(familiesAndLLs.begin(), familiesAndLLs.end());
    std::ofstream llOs(perFamilyLikelihoodPath);
    for (const auto &sf: familiesAndLLs) {
      llOs << sf.familyName << " " << sf.score << std::endl;
    }
    llOs.close();
  }
  Logger::info << "save per-family likelihoods to " << perFamilyLikelihoodPath << std::endl;
  // save the DTL rates
  auto parameterNames = Enums::parameterNames(_info.model);
  auto ratesDir = FileSystem::joinPaths(_outputDir, "model_parameters");
  FileSystem::mkdir(ratesDir, true);
  ParallelContext::barrier();
  if (_info.perFamilyRates) {
    for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
      auto familyRatesPath = FileSystem::joinPaths(ratesDir,
          _geneTrees.getTrees()[i].name + "_rates.txt");
      std::ofstream ratesOs(familyRatesPath);
      ratesOs << "# ";
      for (auto names: parameterNames) {
        ratesOs << names << " ";
      }
      ratesOs << std::endl;
      auto parameters = getModelParameters()[i].getParameters();
      assert(getModelParameters()[i].getParamTypeNumber() == parameterNames.size());
      for (unsigned int j = 0; j < getModelParameters()[i].getParamTypeNumber(); ++j) {
        ratesOs << parameters[j] << " ";
      }
      ratesOs.close();
    }
  } else {
    auto globalRatesPath = FileSystem::joinPaths(ratesDir, "model_parameters.txt");
    ParallelOfstream ratesOs(globalRatesPath, true);
    for (auto node: getSpeciesTree().getTree().getNodes()) {
      ratesOs << node->label;
      for (unsigned int rate = 0; rate < getModelParameters()[0].getParamTypeNumber(); ++rate) {
        ratesOs << " " << getModelParameters()[0].getParameter(node->node_index, rate);
      }
      ratesOs << "\n";
    }
    ratesOs.close();
  }
  Logger::info << "save model rates to " << ratesDir << std::endl;
  ParallelContext::barrier();
}

void AleOptimizer::saveGeneConsensusTree(const std::string &geneSampleFile,
    const std::string &outputFile)
{
  std::ifstream is(geneSampleFile);
  std::string line;
  std::vector<std::string> newicks;
  while (std::getline(is, line)) {
    if (line.size() > 2) {
      newicks.push_back(line);
    }
  }
  is.close();
  auto cons_newick = PLLRootedTree::buildConsensusTree(newicks, 0.50001);
  std::ofstream os(outputFile);
  os << cons_newick << std::endl;
  os.close();
}

void AleOptimizer::saveFamiliesTakingHighway(const Highway &highway,
    const VectorDouble &perFamilyTransfers,
    const std::string &directory)
{
  std::vector<unsigned int> localIndices;
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    auto famIndex = _geneTrees.getTrees()[i].familyIndex;
    localIndices.push_back(famIndex);
  }
  ParallelContext::barrier();
  std::vector<unsigned int> indices;
  std::vector<double> transfers;
  ParallelContext::concatenateHetherogeneousUIntVectors(localIndices, indices);
  ParallelContext::concatenateHetherogeneousDoubleVectors(perFamilyTransfers, transfers);
  assert(indices.size() == _families.size());
  if (ParallelContext::getRank() == 0) {
    std::vector<ScoredFamily> scoredFamilies;
    for (unsigned int i = 0; i < indices.size(); ++i) {
      scoredFamilies.push_back(ScoredFamily(_families[indices[i]].name, transfers[i]));
    }
    std::sort(scoredFamilies.rbegin(), scoredFamilies.rend());
    auto outputFile = FileSystem::joinPaths(directory,
        std::string("families_highway_") + highway.src->label + "_to_" + highway.dest->label + ".txt");
    std::ofstream os(outputFile);
    for (const auto &sf: scoredFamilies) {
      if (sf.score == 0.0) {
        break;
      }
      os << sf.score << " " << sf.familyName << std::endl;
    }
    os.close();
  }
  ParallelContext::barrier();
}

void AleOptimizer::reconcile(unsigned int samples)
{
  assert(getCurrentStep() == AleStep::Reconciliation);
  if (samples == 0) {
    return;
  }
  Logger::timed << "[Reconciliation] Sampling reconciled gene trees (" << samples
                << " samples per gene family)" << std::endl;
  // Initiating directories
  auto recDir = FileSystem::joinPaths(_outputDir, "reconciliations");
  FileSystem::mkdir(recDir, true);
  auto allRecDir = FileSystem::joinPaths(recDir, "all");
  FileSystem::mkdir(allRecDir, true);
  auto summariesDir = FileSystem::joinPaths(recDir, "summaries");
  FileSystem::mkdir(summariesDir, true);
  auto originsDir = FileSystem::joinPaths(recDir, "origins");
  FileSystem::mkdir(originsDir, true);
  auto highwaysOutputDir = FileSystem::joinPaths(_outputDir, "highways");
  auto highwayFamiliesOutputDir = FileSystem::joinPaths(highwaysOutputDir, "families");
  if (FileSystem::dirExists(highwaysOutputDir)) {
    FileSystem::mkdir(highwayFamiliesOutputDir, true);
  }
  ParallelContext::barrier();
  Logger::timed << "[Reconciliation] Inferring reconciliations" << std::endl;
  const auto &localFamilies = _geneTrees.getTrees();
  std::vector<std::string> summaryPerSpeciesEventCountsFiles;
  std::vector<std::string> summaryTransferFiles;
  std::vector< std::shared_ptr<Scenario> > allScenarios;
  MatrixDouble perHighwayPerFamTransfers;
  const auto &highways = _evaluator->getHighways();
  if (highways.size() && localFamilies.size()) {
    assert(FileSystem::dirExists(highwayFamiliesOutputDir));
    perHighwayPerFamTransfers = MatrixDouble(highways.size(), VectorDouble(localFamilies.size(), 0.0));
  }
  for (unsigned int i = 0; i < localFamilies.size(); ++i) {
    std::vector<std::string> perSpeciesEventCountsFiles;
    std::vector<std::string> transferFiles;
    std::vector< std::shared_ptr<Scenario> > scenarios;
    // Warning:
    // Using Random::getProba() in the sampling function makes
    // the random state inconsistent between the MPI ranks.
    // Call ParallelContext::makeRandConsistent() right after
    // all MPI ranks passed the loop
    _evaluator->sampleFamilyScenarios(i, samples, scenarios);
    allScenarios.insert(allScenarios.end(), scenarios.begin(), scenarios.end());
    assert(scenarios.size() == samples);
    // writing in the reconciliations/all/ dir
    auto geneTreesPath = FileSystem::joinPaths(allRecDir, localFamilies[i].name + std::string(".newick"));
    ParallelOfstream geneTreesOs(geneTreesPath, false);
    auto geneTreesAlePath = FileSystem::joinPaths(allRecDir, localFamilies[i].name + std::string(".rec_uml"));
    ParallelOfstream geneTreesAleOs(geneTreesAlePath, false);
    for (unsigned int sample = 0; sample < samples; ++sample) {
      auto out = FileSystem::joinPaths(allRecDir,
          localFamilies[i].name + std::string("_") + std::to_string(sample) + ".xml");
      auto eventCountsFile = FileSystem::joinPaths(allRecDir,
          localFamilies[i].name + std::string("_eventcount_") + std::to_string(sample) + ".txt");
      auto perSpeciesEventCountsFile = FileSystem::joinPaths(allRecDir,
          localFamilies[i].name + std::string("_perspecies_eventcount_") + std::to_string(sample) + ".txt");
      auto transferFile = FileSystem::joinPaths(allRecDir,
          localFamilies[i].name + std::string("_transfers_") + std::to_string(sample) + ".txt");
      perSpeciesEventCountsFiles.push_back(perSpeciesEventCountsFile);
      transferFiles.push_back(transferFile);
      auto &scenario = *scenarios[sample];
      scenario.saveReconciliation(out, ReconciliationFormat::RecPhyloXML, false);
      scenario.saveReconciliation(geneTreesOs, ReconciliationFormat::NewickEvents);
      geneTreesOs << "\n";
      scenario.saveReconciliation(geneTreesAleOs, ReconciliationFormat::ALE);
      scenario.saveEventsCounts(eventCountsFile, false);
      scenario.savePerSpeciesEventsCounts(perSpeciesEventCountsFile, false);
      scenario.saveTransfers(transferFile, false);
      for (unsigned int hi = 0; hi < highways.size(); ++hi) {
        perHighwayPerFamTransfers[hi][i] += scenario.countTransfer(highways[hi].src->label,
            highways[hi].dest->label);
      }
    }
    for (unsigned int hi = 0; hi < highways.size(); ++hi) {
      perHighwayPerFamTransfers[hi][i] /= static_cast<double>(samples);
    }
    geneTreesOs.close();
    geneTreesAleOs.close();
    // writing in the reconciliations/summaries/ dir
    auto consensusFile = FileSystem::joinPaths(summariesDir, localFamilies[i].name +
        std::string("_consensus_50.newick"));
    saveGeneConsensusTree(geneTreesPath,
        consensusFile);
    auto perSpeciesEventCountsFile = FileSystem::joinPaths(summariesDir, localFamilies[i].name +
        std::string("_perspecies_eventcount.txt"));
    Scenario::mergePerSpeciesEventCounts(getSpeciesTree().getTree(),
        perSpeciesEventCountsFile,
        perSpeciesEventCountsFiles, false, true);
    summaryPerSpeciesEventCountsFiles.push_back(perSpeciesEventCountsFile);
    auto transferFile = FileSystem::joinPaths(summariesDir, localFamilies[i].name +
        std::string("_transfers.txt"));
    Scenario::mergeTransfers(getSpeciesTree().getTree(),
        transferFile,
        transferFiles, false, true);
    summaryTransferFiles.push_back(transferFile);
  }
  ParallelContext::barrier();
  ParallelContext::makeRandConsistent();
  Logger::timed << "[Reconciliation] Exporting reconciliation summaries" << std::endl;
  // export total per-branch event counts
  auto totalPerSpeciesEventCountsFile = FileSystem::joinPaths(recDir, "perspecies_eventcount.txt");
  Scenario::mergePerSpeciesEventCounts(getSpeciesTree().getTree(),
      totalPerSpeciesEventCountsFile,
      summaryPerSpeciesEventCountsFiles,
      true, false);
  // export origins
  Scenario::saveOriginsGlobal(getSpeciesTree().getTree(), allScenarios, samples, originsDir);
  // export total pairwise transfer counts
  auto totalTransferFile = FileSystem::joinPaths(recDir, "transfers.txt");
  Scenario::mergeTransfers(getSpeciesTree().getTree(),
      totalTransferFile,
      summaryTransferFiles,
      true, false);
  // export highway information
  for (unsigned int hi = 0; hi < highways.size(); ++hi) {
    saveFamiliesTakingHighway(highways[hi], perHighwayPerFamTransfers[hi], highwayFamiliesOutputDir);
  }
  assert(ParallelContext::isRandConsistent());
  Logger::timed << "Reconciliation output directory: " << recDir << std::endl;
}

void AleOptimizer::saveBestHighways(const std::vector<ScoredHighway> &scoredHighways,
    const std::string &outputFile)
{
  Logger::info << "save highways to " << outputFile << std::endl;
  assert(scoredHighways.size());
  ParallelOfstream os(outputFile, true);
  for (const auto &scoredHighway: scoredHighways) {
    os << scoredHighway.highway.proba << ", ";
    os << scoredHighway.highway.src->label << ",";
    os << scoredHighway.highway.dest->label << ", ";
    os << scoredHighway.scoreDiff << std::endl;
  }
  os.close();
  ParallelContext::barrier();
}

void AleOptimizer::inferHighways(const std::string &highwayCandidateFile,
    unsigned int highwayCandidatesStep1,
    unsigned int highwayCandidatesStep2)
{
  // let's infer highways of transfers!
  assert(getCurrentStep() == AleStep::Highways);
  auto highwaysOutputDir = getHighwaysOutputDir();
  // Step 1: select initial candidates
  auto candidateHighwayOutput = FileSystem::joinPaths(highwaysOutputDir,
      "highway_best_candidates.txt");
  std::vector<ScoredHighway> candidateHighways;
  if (highwayCandidateFile.size()) {
    // the user sets the candidates
    auto highways = HighwayCandidateParser::parse(highwayCandidateFile,
        getSpeciesTree().getTree());
    Highways::getSortedCandidatesFromList(*this,
        highways,
        candidateHighways);
  } else {
    // automatically search for candidates
    Highways::getCandidateHighways(*this,
        candidateHighways,
        highwayCandidatesStep1);
  }
  if (!candidateHighways.size()) {
    Logger::timed << "No candidate highways found!" << std::endl;
    return;
  }
  FileSystem::mkdir(highwaysOutputDir, true);
  saveBestHighways(candidateHighways, candidateHighwayOutput);
  // Step 2: candidate filtering. We add each highway candidate individually, set a small highway
  // probability and keep the highway if the likelihood improves. We also sort the highways
  // by likelihood and keep the int(highwayCandidatesStep2) best of them
  std::vector<ScoredHighway> filteredHighways;
  Highways::filterCandidateHighways(*this,
      candidateHighways,
      filteredHighways,
      highwayCandidatesStep2);
  if (!filteredHighways.size()) {
    Logger::timed << "No candidate highways passed the filtering!" << std::endl;
    return;
  }
  // Step 3: optimize all the filtered highways together
  auto acceptedHighwayOutput = FileSystem::joinPaths(highwaysOutputDir,
      "highway_accepted_highways.txt");
  std::vector<ScoredHighway> acceptedHighways;
  Highways::optimizeAllHighways(*this,
      filteredHighways,
      acceptedHighways,
      true);
  assert(acceptedHighways.size()); // there must be some since we can get here
  saveBestHighways(acceptedHighways, acceptedHighwayOutput);
  Logger::timed << "Highway output directory: " << highwaysOutputDir << std::endl;
}

void AleOptimizer::saveCheckpoint()
{
  if (_info.isDated()) {
    // We rescale the branch lengths because relative dates will be
    // unserialized from the branch lengths
    getSpeciesTree().getDatedTree().rescaleBranchLengths();
  }
  _state.serialize(_checkpointDir);
}

void AleOptimizer::loadCheckpoint()
{
  assert(checkpointExists());
  std::vector<std::string> localFamilyNames;
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    localFamilyNames.push_back(_geneTrees.getTrees()[i].name);
  }
  _state.unserialize(_checkpointDir, localFamilyNames);
}

bool AleOptimizer::checkpointExists(const std::string &outputDir)
{
  return FileSystem::dirExists(getCheckpointDir(outputDir));
}


