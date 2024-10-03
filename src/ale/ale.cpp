#include "AleArguments.hpp"
#include "AleOptimizer.hpp"
#include "Highways.hpp"
#include "TrimFamilies.hpp"
#include <IO/Families.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/FileSystem.hpp>
#include <IO/HighwayCandidateParser.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <ccp/ConditionalClades.hpp>
#include <cstdio>
#include <numeric>
#include <parallelization/ParallelContext.hpp>
#include <routines/Routines.hpp>
#include <routines/SlavesMain.hpp>
#include <util/Paths.hpp>
#include <util/RecModelInfo.hpp>
#include <util/enums.hpp>

const char *version = "AleRax v1.2.0";

/**
 *  Check the validity of each gene family:
 *  - check that the gene tree files exist
 *  - if set, check that the mapping files exist
 */
void filterInvalidFamilies(Families &families) {
  Logger::timed << "Filtering families" << std::endl;
  Families validFamilies;
  for (const auto &family : families) {
    // check that the gene tree distribution file exists
    std::ifstream is(family.startingGeneTree);
    if (!is || is.peek() == std::ifstream::traits_type::eof()) {
      Logger::error << "Can't open input gene trees for family " << family.name
                    << std::endl;
      continue;
    }
    if (family.mappingFile.size()) {
      // check that the mapping file exists if it is set
      std::ifstream isMap(family.mappingFile);
      if (!isMap || isMap.peek() == std::ifstream::traits_type::eof()) {
        Logger::error << "Can't open the mapping file for family "
                      << family.name << std::endl;
        continue;
      }
    }
    validFamilies.push_back(family);
  }
  families = validFamilies;
}

/**
 * Callback for getSortedIndices
 */
bool compare(unsigned int a, unsigned int b,
             const std::vector<unsigned int> &data) {
  return data[a] < data[b];
}

/**
 *  Return the indices of the input vector, sorted in descending order
 */
std::vector<unsigned int>
getSortedIndices(const std::vector<unsigned int> &values) {
  std::vector<unsigned int> indices(values.size());
  std::iota(std::begin(indices), std::end(indices), 0);
  std::sort(
      std::rbegin(indices), std::rend(indices),
      std::bind(compare, std::placeholders::_1, std::placeholders::_2, values));
  return indices;
}

/**
 *  Read the gene tree distribution files (or the .ale files), compute
 *  the corresponding CCP objects, and serialize them to a binary file
 *  Each MPI rank will run this function on a subset of the families
 *
 *  @param families The gene families to process
 *  @param ccpRooting Mode to decide how to deal with rooted gene trees
 *  @param sampleFrequency Number used to subsample the gene trees (if set to 2,
 * we only use half)
 *  @param ccpDimensionFile Output file with the per-family ccp sizes
 *
 */
unsigned int generateCCPs(const Families &families, CCPRooting ccpRooting,
                          unsigned int sampleFrequency,
                          const std::string ccpDimensionFile) {
  ParallelContext::barrier();
  Logger::timed << "Generating ccp files..." << std::endl;
  auto N = families.size();

  std::vector<unsigned int> familyIndices;
  std::vector<unsigned int> treeSizes;
  std::vector<unsigned int> ccpSizes;
  unsigned int numTrees = 0;
  for (auto i = ParallelContext::getBegin(N); i < ParallelContext::getEnd(N);
       i++) {
    // Read and seralize the ccps
    ConditionalClades ccp(families[i].startingGeneTree,
                          families[i].likelihoodFile, ccpRooting,
                          sampleFrequency);
    if (!ccp.isValid()) {
      Logger::error << "Invalid CCP for family " << families[i].name
                    << std::endl;
      Logger::error << "Aborting" << std::endl;
      ParallelContext::abort(1);
    }
    ccp.serialize(families[i].ccp);
    // record the dimensions for subsequent sorting
    familyIndices.push_back(i);
    treeSizes.push_back(ccp.getLeafNumber());
    ccpSizes.push_back(ccp.getCladesNumber());
    numTrees += ccp.getInputTreesNumber();
  }
  ParallelContext::barrier();
  // Gather the family ccp dimensions from all MPI ranks
  ParallelContext::concatenateHetherogeneousUIntVectors(familyIndices,
                                                        familyIndices);
  ParallelContext::concatenateHetherogeneousUIntVectors(treeSizes, treeSizes);
  ParallelContext::concatenateHetherogeneousUIntVectors(ccpSizes, ccpSizes);
  ParallelContext::sumUInt(numTrees);
  // Output the families and there cpp sizes, from the largest to the smallest
  ParallelOfstream os(ccpDimensionFile);
  assert(familyIndices.size() == families.size());
  auto sortedIndices = getSortedIndices(ccpSizes);
  for (unsigned int i = 0; i < familyIndices.size(); ++i) {
    auto j = sortedIndices[i];
    os << families[familyIndices[j]].name << ",";
    os << treeSizes[j] << "," << ccpSizes[j] << std::endl;
  }

  ParallelContext::barrier();
  return numTrees;
}

/**
 *  Remove the serialized CCP files to free some disk space
 */
void cleanupCCPs(Families &families) {
  ParallelContext::barrier();
  Logger::timed << "Cleaning up ccp files..." << std::endl;
  for (const auto &family : families) {
    std::remove(family.ccp.c_str());
  }
  ParallelContext::barrier();
}

/**
 *  Check that a family conditional clade probability object
 *  is valid: check that the mapping gene <-> species is bijective
 */
bool checkCCP(FamilyInfo &family, ConditionalClades &ccp,
              std::unordered_set<std::string> &allSpecies) {
  GeneSpeciesMapping mapping;
  mapping.fill(family.mappingFile, family.startingGeneTree);
  auto m = mapping.getMap();
  for (auto pair : ccp.getCidToLeaves()) {
    auto gene = pair.second;
    auto it = m.find(gene);
    if (it == m.end()) {
      Logger::error << "Gene " << gene << " not mapped to any species "
                    << "in family " << family.name << std::endl;
      return false;
    }
    auto species = it->second;
    if (allSpecies.find(species) == allSpecies.end()) {
      Logger::error << "Gene " << gene << " from family " << family.name
                    << " is mapped to species " << species
                    << " but this species is not in the species tree"
                    << std::endl;
      return false;
    }
  }
  return true;
}

/**
 *  Run checkCCP on all gene families
 */
void checkCCPAndSpeciesTree(Families &families,
                            const std::string &speciesTreePath) {
  Logger::timed << "Checking that ccp and mappings are valid..." << std::endl;
  PLLRootedTree speciesTree(speciesTreePath);
  auto labels = speciesTree.getLabels(true);
  auto N = families.size();
  bool ok = true;
  for (auto i = ParallelContext::getBegin(N); i < ParallelContext::getEnd(N);
       i++) {
    ConditionalClades ccp;
    ccp.unserialize(families[i].ccp);
    ok &= checkCCP(families[i], ccp, labels);
  }
  ParallelContext::barrier();
  if (!ok) {
    ParallelContext::abort(4);
  }
}

/**
 *  Exclude families that are too large or with too much uncertainty
 *
 *  @param families The set of families to trim
 *  @param minSpecies Species that cover less than minSpecies species will be
 * trimmed
 *  @param trimRatio Proportion of families to trim, starting from the largest
 * ones
 *  @param maxCladeSplitRatio Trim families with (ccp size) / (gene tree size) <
 * maxCladeSplitRatio
 */
void trimFamilies(Families &families, int minSpecies, double trimRatio,
                  double maxCladeSplitRatio) {
  Logger::timed << "Families: " << families.size() << std::endl;
  if (minSpecies != -1) {
    Logger::timed << "Triming families covering less than " << minSpecies
                  << " species " << std::endl;
    TrimFamilies::trimMinSpeciesCoverage(families, minSpecies);
    Logger::timed << "Remaining families: " << families.size() << std::endl;
  }
  if (trimRatio > 0.0) {
    Logger::timed << "Trimming families with too many clades (keeping "
                  << (1.0 - trimRatio) * 100.0 << "\% of the families) "
                  << std::endl;
    TrimFamilies::trimHighCladesNumber(families, (1.0 - trimRatio));
    Logger::timed << "Remaining families: " << families.size() << std::endl;
  }
  if (maxCladeSplitRatio > 0) {
    Logger::timed << "Triming families with a ratio clades/nodes > "
                  << maxCladeSplitRatio << std::endl;
    TrimFamilies::trimCladeSplitRatio(families, maxCladeSplitRatio);
    Logger::timed << "Remaining families: " << families.size() << std::endl;
  }
}

/**
 *  @brief Generate or read the initial species tree
 *
 *  Can be either:
 *  - random
 *  - use-defined
 *  - computed with MiniNJ
 *
 *  @brief args The program arguments
 *  @brief families The gene family descriptions
 */
void initStartingSpeciesTree(AleArguments &args, Families &families) {
  Logger::timed << "Initializing starting species tree..." << std::endl;
  auto startingSpeciesTree =
      Paths::getSpeciesTreeFile(args.output, "starting_species_tree.newick");
  std::unique_ptr<PLLRootedTree> speciesTree(nullptr);
  if (args.speciesTreeAlgorithm == SpeciesTreeAlgorithm::User) {
    unsigned int canRead = 1;
    if (ParallelContext::getRank() == 0) {
      try {
        SpeciesTree reader(args.speciesTree);
      } catch (const std::exception &e) {
        Logger::info << "Error while trying to parse the species tree:"
                     << std::endl;
        Logger::info << e.what() << std::endl;
        canRead = 0;
      }
    }
    ParallelContext::broadcastUInt(0, canRead);
    if (!canRead) {
      ParallelContext::abort(153);
    }
    // add labels to internal nodes
    PLLRootedTree::labelRootedTree(args.speciesTree, startingSpeciesTree);
  } else {
    Routines::computeInitialSpeciesTree(families, args.output,
                                        args.speciesTreeAlgorithm)
        ->save(startingSpeciesTree);
  }
  ParallelContext::barrier();
  args.speciesTree =
      Paths::getSpeciesTreeFile(args.output, "inferred_species_tree.newick");
  if (ParallelContext::getRank() == 0) {
    SpeciesTree copy(startingSpeciesTree);
    copy.getTree().save(args.speciesTree);
  }
  ParallelContext::barrier();
  Logger::timed << "Finished starting species tree initialization" << std::endl;
}

/**
 *  If --per-species-rate is set, we generate the file that lists
 *  the species categories (here, each species has its own category)
 *  that share the same DTL parameters
 */
void generatePerSpeciesRateFile(const std::string &perSpeciesRatesFile,
                                const std::string &speciesTreePath,
                                const RecModelInfo &info) {
  PLLRootedTree speciesTree(speciesTreePath);
  ParallelOfstream os(perSpeciesRatesFile);
  for (auto label : speciesTree.getLabels(false)) {
    os << label << " ";
    for (auto p : info.getParamTypes()) {
      os << p;
    }
    os << "\n";
  }
  os.close();
  ParallelContext::barrier();
}

/**
 *  Create AleRax top directories and set ccpDir
 */
void initAleRaxDirectories(const AleArguments &args, const std::string ccpDir) {
  FileSystem::mkdir(args.output, true);
  FileSystem::mkdir(args.output + "/species_trees", true);
  FileSystem::mkdir(ccpDir, true);
  Logger::initFileOutput(FileSystem::joinPaths(args.output, "alerax"));
}

/**
 *  Print AleRax version, the command line, and a summary of
 *  the parameters
 */
void printInitialMessage(const AleArguments &args) {
  Logger::timed << version << std::endl;
  args.printCommand();
  args.printSummary();
}

/**
 *  Compute all information about the different gene families
 *  Filter the families according to the different filtering
 *  and trimming options
 */
void filterFamilies(const AleArguments &args, Families &families) {
  if (!args.skipFamilyFiltering) {
    filterInvalidFamilies(families);
  }
  trimFamilies(families, args.minCoveredSpecies, args.trimFamilyRatio,
               args.maxCladeSplitRatio);
  if (families.size() == 0) {
    Logger::info << "No valid family, aborting" << std::endl;
    ParallelContext::abort(0);
  }
}

RecModelInfo buildRecModelInfo(const AleArguments &args) {
  return RecModelInfo(
      ArgumentsHelper::strToRecModel(args.reconciliationModelStr), args.recOpt,
      (args.modelParametrization ==
       ModelParametrization::PER_FAMILY), // per family rates
      args.gammaCategories, args.originationStrategy, args.pruneSpeciesTree,
      false, // rooted gene tree
      false, // force gene tree root
      false, // mad rooting
      -1.0,  // branch length threshold
      args.transferConstraint,
      false, // no dup
      args.noTL, args.fractionMissingFile, args.memorySavings);
}

Parameters buildStartingRates(const AleArguments &args,
                              const RecModelInfo &info) {
  Parameters res(info.modelFreeParameters());
  double average = 0.0;
  switch (info.model) {
  case RecModel::UndatedDL:
    res[0] = args.d;
    res[1] = args.l;
    average = (args.d + args.l) / 2.0;
    break;
  case RecModel::UndatedDTL:
    res[0] = args.d;
    res[1] = args.l;
    res[2] = args.t;
    average = (args.d + args.l + args.t) / 3.0;
    break;
  default:
    assert(false);
  }
  if (args.originationStrategy == OriginationStrategy::OPTIMIZE) {
    // The initial value of the origination probabilities doesn't matter
    // because they will be normalized in the likelihood computatation.
    // We set it to the average of the other parameters such that all parameters
    // are in the same range (to facilitate gradient optimization)
    res[res.dimensions() - 1] = average;
  }
  return res;
}

void runSpeciesTreeSearch(const AleArguments &args,
                          AleOptimizer &speciesTreeOptimizer) {
  if (args.randomSpeciesRoot) {
    Logger::timed << "Random root position!" << std::endl;
    speciesTreeOptimizer.randomizeRoot();
  }
  // Species tree search or root search
  switch (args.speciesSearchStrategy) {
  case SpeciesSearchStrategy::HYBRID:
    speciesTreeOptimizer.optimize();
    break;
  case SpeciesSearchStrategy::REROOT:
    speciesTreeOptimizer.optimizeModelRates(true);
    Logger::timed << "First root search, non thorough" << std::endl;
    speciesTreeOptimizer.rootSearch(5, false);
    Logger::timed << "Second root search, thorough" << std::endl;
    speciesTreeOptimizer.rootSearch(2, true);
    break;
  case SpeciesSearchStrategy::SKIP:
    break;
  default:
    assert(false); // not implemented yet
    break;
  }
}

void runDateOptimization(const AleArguments &args,
                         AleOptimizer &speciesTreeOptimizer) {
  if (!args.inferSpeciationOrders) {
    return;
  }
  speciesTreeOptimizer.optimizeDates(true);
  speciesTreeOptimizer.getEvaluator().computeLikelihood();
}

void runTransferHighwayInference(const AleArguments &args,
                                 AleOptimizer &speciesTreeOptimizer,
                                 size_t sample_size) {
  if (!args.highways) {
    return;
  }
  Logger::info << "Highway BIC sample size = " << sample_size << std::endl;
  auto highwaysOutputDir = FileSystem::joinPaths(args.output, "highways");
  FileSystem::mkdir(highwaysOutputDir, true);
  // let's infer highways of transfers!
  auto highwayOutput =
      FileSystem::joinPaths(highwaysOutputDir, "highway_best_candidates.txt");
  std::vector<ScoredHighway> candidateHighways;
  // initial candidates
  if (args.highwayCandidateFile.size()) {
    // the user sets the candidates
    auto highways = HighwayCandidateParser::parse(
        args.highwayCandidateFile,
        speciesTreeOptimizer.getSpeciesTree().getTree());

    candidateHighways =
        Highways::getSortedCandidatesFromList(speciesTreeOptimizer, highways);
  } else {
    // automatically search for candidates
    Highways::getCandidateHighways(speciesTreeOptimizer, candidateHighways,
                                   args.highwayCandidatesStep1);
  }
  // first filtering step: we add each highway candidate individually, set a
  // small highway probability, and keep the highway if the likelihood improves.
  // We also sort the highways per likelihood
  std::vector<ScoredHighway> filteredHighways;
  Highways::filterCandidateHighwaysFast(speciesTreeOptimizer, candidateHighways,
                                        filteredHighways, sample_size);
  filteredHighways.resize(
      std::min(filteredHighways.size(), size_t(args.highwayCandidatesStep2)));
  // now optimize all highways together
  auto acceptedHighwayOutput =
      FileSystem::joinPaths(highwaysOutputDir, "highway_accepted_highways.txt");
  std::vector<ScoredHighway> acceptedHighways;
  Highways::optimizeAllHighways(speciesTreeOptimizer, filteredHighways,
                                acceptedHighways, true);
  speciesTreeOptimizer.saveBestHighways(acceptedHighways,
                                        acceptedHighwayOutput);
}

void initStartingSpeciesTreeFromCheckpoint(AleArguments &args) {
  auto checkpointPath = AleOptimizer::getCheckpointDir(args.output);
  auto newick = AleState::readSpeciesTreeNewick(checkpointPath);
  args.speciesTree =
      Paths::getSpeciesTreeFile(args.output, "inferred_species_tree.newick");
  ParallelOfstream os(args.speciesTree);
  os << newick;
  os.close();
  ParallelContext::barrier();
}

/**
 * Main function of AleRax once the arguments have been parsed
 *
 * @param args The program arguments
 */
void run(AleArguments &args) {
  Random::setSeed(static_cast<unsigned int>(args.seed));
  bool checkpointDetected = AleOptimizer::checkpointExists(args.output);
  if (checkpointDetected) {
    Logger::info << "Checkpoint detected" << std::endl;
  }

  auto ccpDir = FileSystem::joinPaths(args.output, "ccps");
  if (!checkpointDetected) {
    initAleRaxDirectories(args, ccpDir);
  }
  printInitialMessage(args);
  auto families = FamiliesFileParser::parseFamiliesFile(args.families);
  auto ccpDimensionFile = FileSystem::joinPaths(args.output, "ccpdim.txt");
  for (auto &family : families) {
    family.ccp = FileSystem::joinPaths(ccpDir, family.name + ".ccp");
  }
  unsigned int sample_size = 1000 * families.size();
  sample_size = generateCCPs(families, args.ccpRooting, args.sampleFrequency,
                             ccpDimensionFile);
  filterFamilies(args, families);
  if (!checkpointDetected) {
    initStartingSpeciesTree(args, families);
    checkCCPAndSpeciesTree(families, args.speciesTree);
  } else {
    initStartingSpeciesTreeFromCheckpoint(args);
  }
  auto info = buildRecModelInfo(args);
  auto startingRates = buildStartingRates(args, info);
  /*
  if (args.perSpeciesRates) {
    generatePerSpeciesRateFile(args.optimizationClassFile, args.speciesTree,
  info);
  }
  */
  std::string coverageFile(
      FileSystem::joinPaths(args.output, "fractionMissing.txt"));
  std::string fractionMissingFile(
      FileSystem::joinPaths(args.output, "perSpeciesCoverage.txt"));
  Family::printStats(families, args.speciesTree, coverageFile,
                     fractionMissingFile);
  AleOptimizer speciesTreeOptimizer(args.speciesTree, families, info,
                                    args.modelParametrization, startingRates,
                                    !args.fixRates, args.verboseOptRates,
                                    args.optimizationClassFile, args.output);

  if (!checkpointDetected) {
    speciesTreeOptimizer.setCurrentStep(AleStep::SpeciesTreeOpt);
    speciesTreeOptimizer.saveCheckpoint();
  }

  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::SpeciesTreeOpt) {
    runSpeciesTreeSearch(args, speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::ModelRateOpt);
    speciesTreeOptimizer.saveCheckpoint();
  }

  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::ModelRateOpt) {
    speciesTreeOptimizer.optimizeModelRates(false);
    speciesTreeOptimizer.setCurrentStep(AleStep::RelDating);
    speciesTreeOptimizer.saveCheckpoint();
  }
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::RelDating) {
    runDateOptimization(args, speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::Highways);
    speciesTreeOptimizer.saveCheckpoint();
  }
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::Highways) {
    runTransferHighwayInference(args, speciesTreeOptimizer, sample_size);
    speciesTreeOptimizer.setCurrentStep(AleStep::ModelRateOpt2);
    speciesTreeOptimizer.saveCheckpoint();
  }

  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::ModelRateOpt2) {
    if (!args.skipThoroughRates) {
      speciesTreeOptimizer.optimizeModelRates(true);
    }
    speciesTreeOptimizer.setCurrentStep(AleStep::Reconciliation);
    speciesTreeOptimizer.saveCheckpoint();
  }

  Logger::timed << "Sampling reconciled gene trees... (" << args.geneTreeSamples
                << " samples)" << std::endl;
  speciesTreeOptimizer.reconcile(args.geneTreeSamples);
  speciesTreeOptimizer.saveSpeciesTree();
  speciesTreeOptimizer.saveRatesAndLL();
  Logger::timed << "Final log likelihood: ll="
                << speciesTreeOptimizer.getEvaluator().computeLikelihood()
                << std::endl;
  if (args.cleanupCCP) {
    cleanupCCPs(families);
  }
  Logger::timed << "End of the execution" << std::endl;
}

/**
 *  This main will be called by the master process and is the entry
 *  point to AleRax. It inits the parallel context, the logger,
 *  reads the arguments, runs AleRax and closes everything.
 *  If comm is set to nullptr, MPI_Init will be called to initiate the
 *  MPI context
 */
int alerax_main(int argc, char **argv, void *comm) {
  ParallelContext::init(comm);
  Logger::init();
  Logger::timed << version << std::endl;
  AleArguments args(argc, argv);
  args.checkValid();
  run(args);
  Logger::close();
  ParallelContext::finalize();
  return 0;
}

/**
 *  AleRax "internal" main (can be called  by a scheduler with
 *  an already existing MPI communicator, for instance to spawn
 *  a slave process
 */
int internal_main(int argc, char **argv, void *comm) {
  if (SlavesMain::isSlave(argc, argv)) {
    int slaveComm = -1;
    return static_scheduled_main(argc, argv, &slaveComm);
  } else {
    return alerax_main(argc, argv, comm);
  }
}

int main(int argc, char **argv) {
#ifdef WITH_MPI
  // the null communicator signals that we need to call MPI_Init
  return internal_main(argc, argv, nullptr);
#else
  // this fake communicator signals that we're not using MPI
  int noMPIComm = -1;
  return internal_main(argc, argv, &noMPIComm);
#endif
}
