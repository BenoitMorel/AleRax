#include "AleOptimizer.hpp"
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <iomanip>
#include <maths/Random.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <search/DatedSpeciesTreeSearch.hpp>
#include <search/SpeciesSPRSearch.hpp>
#include <search/SpeciesTransferSearch.hpp>
#include <sys/stat.h>
#include <sys/types.h>
#include <util/Paths.hpp>

static bool dirExists(const std::string &pathname) {
  struct stat info;
  if (stat(pathname.c_str(), &info) != 0) {
    return false;
  } else if (info.st_mode & S_IFDIR) {
    return true;
  } else {
    return false;
  }
}

static bool testAndSwap(size_t &hash1, size_t &hash2) {
  std::swap(hash1, hash2);
  return hash1 != hash2;
}

AleOptimizer::AleOptimizer(const std::string speciesTreeFile,
                           const Families &families, const RecModelInfo &info,
                           ModelParametrization modelParametrization,
                           const Parameters &startingRates, bool optimizeRates,
                           bool optimizeVerbose,
                           const std::string &optimizationClassFile,
                           const std::string &outputDir)
    : _state(speciesTreeFile), _families(families),
      _geneTrees(families, false, true), _info(info), _outputDir(outputDir),
      _checkpointDir(getCheckpointDir(outputDir)),
      _rootLikelihoods(_geneTrees.getTrees().size()), _enableCheckpoints(true) {

  if (checkpointExists()) {
    Logger::info << "Loading checkpoint information..." << std::endl;
    loadCheckpoint();
  } else {
    for (unsigned int i = 0; i < _geneTrees.getTrees().size(); ++i) {
      _state.perFamilyModelParameters.push_back(AleModelParameters(
          startingRates, getSpeciesTree().getTree().getNodeNumber()));
      _state.familyNames.push_back(_geneTrees.getTrees()[i].name);
    }
  }
  _speciesTreeSearchState = std::make_unique<SpeciesSearchState>(
      getSpeciesTree(),
      Paths::getSpeciesTreeFile(_outputDir, "inferred_species_tree.newick"),
      _geneTrees.getTrees().size());
  _speciesTreeSearchState->addListener(this);
  Logger::info << getSpeciesTree().getTree().getNewickString() << std::endl;
  getSpeciesTree().addListener(this);
  ParallelContext::barrier();
  _evaluator = std::make_unique<AleEvaluator>(
      *this, getSpeciesTree(), getRecModelInfo(), modelParametrization,
      getModelParameters(), optimizeRates, optimizeVerbose, families,
      _geneTrees, optimizationClassFile, _outputDir);
  Logger::timed << "Initial ll=" << getEvaluator().computeLikelihood()
                << std::endl;
  saveCurrentSpeciesTreeId("starting_species_tree.newick");
  saveCurrentSpeciesTreeId();
  FileSystem::mkdir(_checkpointDir, true);
  saveCheckpoint();
}

double AleOptimizer::optimizeModelRates(bool thorough) {
  getEvaluator().optimizeModelRates(thorough);
  PerFamLL perFamLL;
  auto ll = getEvaluator().computeLikelihood(&perFamLL);
  _speciesTreeSearchState->betterLikelihoodCallback(ll, perFamLL);
  saveRatesAndLL();
  saveCheckpoint();
  return ll;
}

void AleOptimizer::optimize() {
  size_t hash1 = 0;
  size_t hash2 = 0;
  unsigned int index = 0;
  PerFamLL initialPerFamLL;
  _speciesTreeSearchState->bestLL =
      getEvaluator().computeLikelihood(&initialPerFamLL);
  _speciesTreeSearchState->farFromPlausible = true;
  _speciesTreeSearchState->betterTreeCallback(_speciesTreeSearchState->bestLL,
                                              initialPerFamLL);
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
    if (!_speciesTreeSearchState->farFromPlausible) {
      rootSearch(3);
    }
    hash1 = getSpeciesTree().getHash();
  } while (testAndSwap(hash1, hash2));
  rootSearch(-1);
}

double AleOptimizer::sprSearch(unsigned int radius) {
  SpeciesSPRSearch::SPRSearch(getSpeciesTree(), getEvaluator(),
                              *_speciesTreeSearchState, radius);
  Logger::timed << "After normal search: LL=" << _speciesTreeSearchState->bestLL
                << std::endl;
  saveCheckpoint();
  saveSupportTree();
  return _speciesTreeSearchState->bestLL;
}

void AleOptimizer::saveSupportTree() {
  auto outKH =
      Paths::getSpeciesTreeFile(_outputDir, "species_tree_support_kh.newick");
  _speciesTreeSearchState->saveSpeciesTreeKH(outKH);
  Logger::info << "save support tree " << outKH << std::endl;
  auto outBP =
      Paths::getSpeciesTreeFile(_outputDir, "species_tree_support_bp.newick");
  _speciesTreeSearchState->saveSpeciesTreeBP(outBP);
}

void AleOptimizer::onSpeciesTreeChange(
    const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) {
  getEvaluator().onSpeciesTreeChange(nodesToInvalidate);
}

void AleOptimizer::betterTreeCallback() { saveCheckpoint(); }

void AleOptimizer::saveSpeciesTree() { saveCurrentSpeciesTreeId(); }

std::string AleOptimizer::saveCurrentSpeciesTreeId(std::string name,
                                                   bool masterRankOnly) {
  std::string res = Paths::getSpeciesTreeFile(_outputDir, name);
  if (_evaluator->isDated()) {
    getSpeciesTree().getDatedTree().rescaleBranchLengths();
  }
  saveCurrentSpeciesTreePath(res, masterRankOnly);
  if (!_rootLikelihoods.isEmpty()) {
    auto newick = getSpeciesTree().getTree().getNewickString();
    PLLRootedTree tree(newick, false);
    _rootLikelihoods.fillTree(tree);
    auto out = Paths::getSpeciesTreeFile(_outputDir, "species_tree_llr.newick");
    tree.save(out);
  }
  if (!_rootLikelihoods.isEmpty()) {
    auto newick = getSpeciesTree().getTree().getNewickString();
    PLLRootedTree tree(newick, false);
    _rootLikelihoods.fillTreeBootstraps(tree);
    auto out =
        Paths::getSpeciesTreeFile(_outputDir, "species_tree_root_rell.newick");
    tree.save(out);
  }
  return res;
}

void AleOptimizer::saveCurrentSpeciesTreePath(const std::string &str,
                                              bool masterRankOnly) {
  getSpeciesTree().saveToFile(str, masterRankOnly);
  if (masterRankOnly) {
    ParallelContext::barrier();
  }
}

double AleOptimizer::rootSearch(unsigned int maxDepth, bool thorough) {
  _rootLikelihoods.reset();
  if (thorough) {
    _speciesTreeSearchState->farFromPlausible = thorough;
  }
  SpeciesRootSearch::rootSearch(getSpeciesTree(), getEvaluator(),
                                *_speciesTreeSearchState, maxDepth,
                                &_rootLikelihoods);
  saveCheckpoint();
  return _speciesTreeSearchState->bestLL;
}

double AleOptimizer::transferSearch() {
  SpeciesTransferSearch::transferSearch(getSpeciesTree(), getEvaluator(),
                                        *_speciesTreeSearchState);
  Logger::timed << "After normal search: LL=" << _speciesTreeSearchState->bestLL
                << std::endl;
  saveSupportTree();
  return _speciesTreeSearchState->bestLL;
}

std::vector<std::string> getLines(const std::string path) {
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

void AleOptimizer::saveRatesAndLL() {

  // save per-family likelihood
  auto perFamilyLikelihoodPath =
      FileSystem::joinPaths(_outputDir, "per_fam_likelihoods.txt");
  std::vector<unsigned int> indices;
  std::vector<double> likelihoods;
  for (unsigned int i = 0; i < _geneTrees.getTrees().size(); ++i) {
    auto famIndex = _geneTrees.getTrees()[i].familyIndex;
    auto ll = _evaluator->computeFamilyLikelihood(i);
    indices.push_back(famIndex);
    likelihoods.push_back(ll);
  }
  std::vector<unsigned int> allIndices;
  std::vector<double> allLikelihoods;
  ParallelContext::concatenateHetherogeneousDoubleVectors(likelihoods,
                                                          allLikelihoods);
  ParallelContext::concatenateHetherogeneousUIntVectors(indices, allIndices);
  assert(allLikelihoods.size() == allIndices.size());
  std::vector<std::pair<double, std::string>> likelihoodAndFamilies;
  for (unsigned int i = 0; i < allLikelihoods.size(); ++i) {
    likelihoodAndFamilies.push_back(
        {allLikelihoods[i], _families[allIndices[i]].name});
  }
  std::sort(likelihoodAndFamilies.begin(), likelihoodAndFamilies.end());
  ParallelOfstream llOs(perFamilyLikelihoodPath);
  for (auto p : likelihoodAndFamilies) {
    llOs << p.second << " " << p.first << std::endl;
  }
  // save the DTL rates
  auto parameterNames = Enums::parameterNames(_info.model);
  auto ratesDir = FileSystem::joinPaths(_outputDir, "model_parameters");
  FileSystem::mkdir(ratesDir, true);
  ParallelContext::barrier();
  if (_info.perFamilyRates) {
    for (unsigned int i = 0; i < _geneTrees.getTrees().size(); ++i) {
      auto &geneTree = _geneTrees.getTrees()[i];
      auto family = geneTree.name;
      auto ratesPath = FileSystem::joinPaths(ratesDir, family + "_rates.txt");
      std::ofstream ratesOs(ratesPath);
      ratesOs << "# ";
      for (auto names : parameterNames) {
        ratesOs << names << " ";
      }
      ratesOs << std::endl;
      auto parameters = getModelParameters()[i].getParameters();
      assert(getModelParameters()[i].getParamTypeNumber() ==
             parameterNames.size());
      for (unsigned int j = 0; j < getModelParameters()[i].getParamTypeNumber();
           ++j) {
        ratesOs << parameters[j] << " ";
      }
    }
  } else {
    auto globalRatesPath =
        FileSystem::joinPaths(ratesDir, "model_parameters.txt");
    ParallelOfstream ratesOs(globalRatesPath, true);
    if (getModelParameters().size() >
        0) { // avoid segfault if #cores > #families
      for (auto node : getSpeciesTree().getTree().getNodes()) {
        ratesOs << node->label;
        for (unsigned int rate = 0;
             rate < getModelParameters()[0].getParamTypeNumber(); ++rate) {
          ratesOs << " "
                  << getModelParameters()[0].getParameter(node->node_index,
                                                          rate);
        }
        ratesOs << "\n";
      }
    }
  }
}

static void saveFamiliesTakingHighway(const Highway &highway,
                                      const Families &families,
                                      const PerCoreGeneTrees &geneTrees,
                                      const VectorDouble &perFamilyTransfers,
                                      const std::string &directory) {
  struct ScoredFamily {
    ScoredFamily(double score, const std::string &familyName)
        : score(score), familyName(familyName) {}
    double score;
    std::string familyName;
    bool operator<(const ScoredFamily &other) const {
      return score < other.score;
    }
  };

  std::vector<unsigned int> localIndices(geneTrees.getTrees().size());
  for (unsigned int i = 0; i < geneTrees.getTrees().size(); ++i) {
    localIndices[i] = geneTrees.getTrees()[i].familyIndex;
  }
  std::vector<unsigned int> globalIndices;
  std::vector<double> globalTransfers;
  ParallelContext::concatenateHetherogeneousUIntVectors(localIndices,
                                                        globalIndices);
  ParallelContext::concatenateHetherogeneousDoubleVectors(perFamilyTransfers,
                                                          globalTransfers);
  assert(globalIndices.size() == globalTransfers.size());

  // only master rank will save the file
  if (ParallelContext::getRank() == 0) {
    std::vector<ScoredFamily> scoredFamilies;
    for (unsigned int i = 0; i < globalIndices.size(); ++i) {
      scoredFamilies.push_back(
          ScoredFamily(globalTransfers[i], families[globalIndices[i]].name));
    }
    std::sort(scoredFamilies.rbegin(), scoredFamilies.rend());
    std::string outputFileName = FileSystem::joinPaths(
        directory, std::string("families_highway_") + highway.src->label +
                       "_to_" + highway.dest->label + ".txt");
    std::ofstream os(outputFileName);
    for (auto sf : scoredFamilies) {
      if (sf.score == 0.0) {
        break;
      }
      os << sf.score << " " << sf.familyName << std::endl;
    }
  }
  ParallelContext::barrier();
}

void AleOptimizer::reconcile(unsigned int samples) {
  if (samples == 0) {
    return;
  }
  Logger::timed << "Sampling reconciliations..." << std::endl;
  auto recDir = FileSystem::joinPaths(_outputDir, "reconciliations");
  FileSystem::mkdir(recDir, true);
  auto allRecDir = FileSystem::joinPaths(recDir, "all");
  FileSystem::mkdir(allRecDir, true);
  auto summariesDir = FileSystem::joinPaths(recDir, "summaries");
  FileSystem::mkdir(summariesDir, true);

  auto highwaysOutputDir = FileSystem::joinPaths(_outputDir, "highways");
  auto highwayFamiliesOutputDir =
      FileSystem::joinPaths(highwaysOutputDir, "families");
  FileSystem::mkdir(highwayFamiliesOutputDir, true);
  ParallelContext::barrier();
  auto &localFamilies = _geneTrees.getTrees();
  std::vector<std::string> summaryPerSpeciesEventCountsFiles;
  std::vector<std::string> summaryTransferFiles;
  std::vector<std::shared_ptr<Scenario>> allScenarios;
  Logger::timed << "Exporting reconciliations..." << std::endl;
  const auto &highways = _evaluator->getHighways();
  MatrixDouble perHighwayPerFamTransfers;
  if (highways.size() && localFamilies.size()) {
    perHighwayPerFamTransfers =
        MatrixDouble(highways.size(), VectorDouble(localFamilies.size(), 0.0));
  }
  for (unsigned int i = 0; i < localFamilies.size(); ++i) {
    std::vector<std::string> perSpeciesEventCountsFiles;
    std::vector<std::string> transferFiles;
    std::string geneTreesPath = FileSystem::joinPaths(
        allRecDir, localFamilies[i].name + std::string(".newick"));
    std::string geneTreesAlePath = FileSystem::joinPaths(
        allRecDir, localFamilies[i].name + std::string(".rec_uml"));
    ParallelOfstream geneTreesOs(geneTreesPath, false);
    std::vector<std::shared_ptr<Scenario>> scenarios;
    _evaluator->sampleScenarios(i, samples, scenarios);
    allScenarios.insert(allScenarios.end(), scenarios.begin(), scenarios.end());
    assert(scenarios.size() == samples);

    ParallelOfstream geneTreesAleOs(geneTreesAlePath, false);
    for (unsigned int sample = 0; sample < samples; ++sample) {
      auto out = FileSystem::joinPaths(
          allRecDir, localFamilies[i].name + std::string("_") +
                         std::to_string(sample) + ".xml");
      auto eventCountsFile = FileSystem::joinPaths(
          allRecDir, localFamilies[i].name + std::string("_eventcount_") +
                         std::to_string(sample) + ".txt");
      auto perSpeciesEventCountsFile = FileSystem::joinPaths(
          allRecDir, localFamilies[i].name +
                         std::string("_perspecies_eventcount_") +
                         std::to_string(sample) + ".txt");
      auto transferFile = FileSystem::joinPaths(
          allRecDir, localFamilies[i].name + std::string("_transfers_") +
                         std::to_string(sample) + ".txt");
      perSpeciesEventCountsFiles.push_back(perSpeciesEventCountsFile);
      transferFiles.push_back(transferFile);
      auto &scenario = *scenarios[sample];
      scenario.saveReconciliation(out, ReconciliationFormat::RecPhyloXML,
                                  false);
      scenario.saveReconciliation(geneTreesOs,
                                  ReconciliationFormat::NewickEvents);
      scenario.saveReconciliation(geneTreesAleOs, ReconciliationFormat::ALE);
      scenario.saveEventsCounts(eventCountsFile, false);
      scenario.savePerSpeciesEventsCounts(perSpeciesEventCountsFile, false);
      scenario.saveTransfers(transferFile, false);
      geneTreesOs << "\n";
      for (unsigned int hi = 0; hi < highways.size(); ++hi) {
        perHighwayPerFamTransfers[hi][i] += scenario.countTransfer(
            highways[hi].src->label, highways[hi].dest->label);
      }
    }
    for (unsigned int hi = 0; hi < highways.size(); ++hi) {
      perHighwayPerFamTransfers[hi][i] /= static_cast<double>(samples);
    }

    geneTreesOs.close();
    auto newicks = getLines(geneTreesPath);
    auto consensusPrefix = FileSystem::joinPaths(
        summariesDir, localFamilies[i].name + "_consensus_");

    saveStr(PLLRootedTree::buildConsensusTree(newicks, 0.50001),
            consensusPrefix + "50.newick");
    auto perSpeciesEventCountsFile = FileSystem::joinPaths(
        summariesDir,
        localFamilies[i].name + std::string("_perspecies_eventcount.txt"));
    Scenario::mergePerSpeciesEventCounts(
        getSpeciesTree().getTree(), perSpeciesEventCountsFile,
        perSpeciesEventCountsFiles, false, true);
    summaryPerSpeciesEventCountsFiles.push_back(perSpeciesEventCountsFile);
    auto transferFile = FileSystem::joinPaths(
        summariesDir, localFamilies[i].name + std::string("_transfers.txt"));
    Scenario::mergeTransfers(getSpeciesTree().getTree(), transferFile,
                             transferFiles, false, true);
    summaryTransferFiles.push_back(transferFile);
  }
  ParallelContext::barrier();
  Logger::timed << "Exporting reconciliation summaries..." << std::endl;
  auto totalPerSpeciesEventCountsFile =
      FileSystem::joinPaths(recDir, "perspecies_eventcount.txt");
  Scenario::mergePerSpeciesEventCounts(
      getSpeciesTree().getTree(), totalPerSpeciesEventCountsFile,
      summaryPerSpeciesEventCountsFiles, true, false);
  auto originsDir = FileSystem::joinPaths(recDir, "origins");
  FileSystem::mkdir(originsDir, true);
  Scenario::saveOriginsGlobal(getSpeciesTree().getTree(), allScenarios, samples,
                              originsDir);
  auto totalTransferFile = FileSystem::joinPaths(recDir, "transfers.txt");
  Scenario::mergeTransfers(getSpeciesTree().getTree(), totalTransferFile,
                           summaryTransferFiles, true, false);
  for (unsigned int hi = 0; hi < highways.size(); ++hi) {
    saveFamiliesTakingHighway(highways[hi], _families, _geneTrees,
                              perHighwayPerFamTransfers[hi],
                              highwayFamiliesOutputDir);
  }
  ParallelContext::makeRandConsistent();
  Logger::timed << "Reconciliations output directory: " << recDir << std::endl;
}

void AleOptimizer::optimizeDates(bool) {
  if (!_info.isDated()) {
    return;
  }
  auto scoredBackups = DatedSpeciesTreeSearch::optimizeDatesFromReconciliation(
      getSpeciesTree(), getEvaluator(), 500, 20);
  Logger::timed << "Sorted dating likelihoods from fast datings:" << std::endl;
  for (auto &sb : scoredBackups) {
    Logger::info << sb.score << " ";
  }
  Logger::info << std::endl;

  auto bestLL = scoredBackups[0].score;
  _speciesTreeSearchState->bestLL = bestLL;
  getSpeciesTree().getDatedTree().restore(scoredBackups[0].backup);
  DatedSpeciesTreeSearch::optimizeDates(getSpeciesTree(), getEvaluator(),
                                        *_speciesTreeSearchState, bestLL, true);
  saveCheckpoint();
  saveCurrentSpeciesTreeId();
}

void AleOptimizer::randomizeRoot() {
  auto &tree = getSpeciesTree().getDatedTree();
  unsigned int N = tree.getOrderedSpeciations().size();
  for (unsigned int i = 0; i < N; ++i) {
    auto direction = Random::getInt() % 4;
    if (SpeciesTreeOperator::canChangeRoot(getSpeciesTree(), direction)) {
      SpeciesTreeOperator::changeRoot(getSpeciesTree(), direction);
    }
  }
}

void AleOptimizer::saveBestHighways(
    const std::vector<ScoredHighway> &scoredHighways,
    const std::string &output) {
  ParallelOfstream os(output, true);
  Logger::info << "Outputing the highays into " << output << std::endl;
  for (const auto &scoredHighway : scoredHighways) {
    if (scoredHighway.highway.proba >= 0.000001) {
      os << scoredHighway.highway.proba << ", ";
      os << scoredHighway.highway.src->label << ",";
      os << scoredHighway.highway.dest->label << ", ";
      os << scoredHighway.scoreDiff << std::endl;
    }
  }
}

bool cmpHighwayByProbability(const ScoredHighway &a, const ScoredHighway &b) {
  return a.highway.proba < b.highway.proba;
}

std::string AleOptimizer::getHighwaysOutputDir() const {
  return FileSystem::joinPaths(_outputDir, "highways");
}

void AleOptimizer::saveCheckpoint() const { _state.serialize(_checkpointDir); }

void AleOptimizer::loadCheckpoint() {
  assert(checkpointExists());
  std::vector<std::string> perFamilyNames;
  for (auto family : _geneTrees.getTrees()) {
    perFamilyNames.push_back(family.name);
  }
  _state.unserialize(_checkpointDir, perFamilyNames);
}

bool AleOptimizer::checkpointExists(const std::string &outputDir) {
  return dirExists(getCheckpointDir(outputDir));
}

void AleOptimizer::onBetterParametersFoundCallback() {
  if (_enableCheckpoints) {
    saveCheckpoint();
  }
}
