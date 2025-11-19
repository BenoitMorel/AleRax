#include <DistanceMethods/MiniNJ.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <IO/ParallelOfstream.hpp>
#include <ccp/ConditionalClades.hpp>
#include <parallelization/ParallelContext.hpp>
#include <util/Paths.hpp>
#include <util/RecModelInfo.hpp>
#include <util/enums.hpp>

#include "AleArguments.hpp"
#include "AleOptimizer.hpp"
#include "TrimFamilies.hpp"

const char *version = "AleRax v1.4.0";

/**
 *  Create AleRax top directories
 */
void initAleRaxDirectories(const AleArguments &args,
                           const std::string &ccpDir) {
  FileSystem::mkdir(args.output, true);
  FileSystem::mkdir(args.output + "/species_trees", true);
  FileSystem::mkdir(ccpDir, true);
  ParallelContext::barrier();
}

/**
 *  Print AleRax version, the command line and a summary of
 *  the parameters
 */
void printInitialMessage(const AleArguments &args) {
  Logger::initFileOutput(FileSystem::joinPaths(args.output, "alerax"));
  Logger::timed << version << std::endl;
  args.printCommand();
  args.printSummary();
}

/**
 *  Check that there are enough gene families
 */
void checkEnoughFamilies(const Families &families) {
  if (families.size() == 0) {
    Logger::info << "\nError: No valid families, aborting\n" << std::endl;
    ParallelContext::abort(0);
  }
  if (families.size() < ParallelContext::getSize()) {
    Logger::info << "\nError: More MPI ranks set than valid families, "
                 << "aborting\n"
                 << std::endl;
    ParallelContext::abort(0);
  }
}

/**
 *  Check the validity of each gene family:
 *  - check that the gene tree files exist
 *  - if set, check that the mapping files exist
 */
void filterInvalidFamilies(const AleArguments &args, Families &families) {
  if (args.skipFamilyFiltering) {
    return;
  }
  Logger::timed << "Checking families..." << std::endl;
  Families validFamilies;
  for (const auto &family : families) {
    // check that the gene tree distribution file exists
    std::ifstream isTrees(family.startingGeneTree);
    if (!isTrees || isTrees.peek() == std::ifstream::traits_type::eof()) {
      Logger::info << "Excluding family " << family.name
                   << ": cannot open the gene tree distribution file!"
                   << std::endl;
      continue;
    }
    if (family.mappingFile.size()) {
      // check that the mapping file exists if it is set
      std::ifstream isMap(family.mappingFile);
      if (!isMap || isMap.peek() == std::ifstream::traits_type::eof()) {
        Logger::info << "Excluding family " << family.name
                     << ": cannot open the mapping file!" << std::endl;
        continue;
      }
    }
    validFamilies.push_back(family);
  }
  families = validFamilies;
  checkEnoughFamilies(families);
  ParallelContext::barrier();
}

/**
 *  Callback for getSortedIndices
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
 *  the corresponding CCP objects and serialize them to a binary file.
 *  Each MPI rank will run this function on a subset of the families
 *
 *  @param families The gene families to process
 *  @param ccpRooting Mode to decide how to deal with rooted gene trees
 *  @param sampleFrequency Number used to subsample the gene trees (if set to 2,
 *  we only use half)
 *  @param ccpDimensionFile Output file with the per-family ccp sizes
 */
void generateCCPs(const AleArguments &args, const Families &families,
                  const std::string &ccpDimensionFile) {
  ParallelContext::barrier();
  Logger::timed << "Generating ccp files..." << std::endl;
  auto N = families.size();
  std::vector<unsigned int> localFamilyIndices;
  std::vector<unsigned int> localTreeSizes;
  std::vector<unsigned int> localCcpSizes;
  for (auto i = ParallelContext::getBegin(N); i < ParallelContext::getEnd(N);
       i++) {
    // read and seralize the ccps
    ConditionalClades ccp(families[i].startingGeneTree,
                          families[i].likelihoodFile, args.ccpRooting,
                          args.sampleFrequency);
    if (!ccp.isValid()) {
      Logger::error << "Error: invalid CCP for family " << families[i].name
                    << std::endl;
      ParallelContext::abort(1);
    }
    ccp.serialize(families[i].ccp);
    // record the dimensions for the subsequent sorting
    localFamilyIndices.push_back(i);
    localTreeSizes.push_back(ccp.getLeafNumber());
    localCcpSizes.push_back(ccp.getCladesNumber());
  }
  ParallelContext::barrier();
  // gather the family ccp dimensions from all MPI ranks
  std::vector<unsigned int> familyIndices;
  std::vector<unsigned int> treeSizes;
  std::vector<unsigned int> ccpSizes;
  ParallelContext::concatenateHetherogeneousUIntVectors(localFamilyIndices,
                                                        familyIndices);
  ParallelContext::concatenateHetherogeneousUIntVectors(localTreeSizes,
                                                        treeSizes);
  ParallelContext::concatenateHetherogeneousUIntVectors(localCcpSizes,
                                                        ccpSizes);
  assert(familyIndices.size() == families.size());
  // output the families and their cpp sizes, from the largest to the smallest
  auto sortedIndices = getSortedIndices(ccpSizes);
  ParallelOfstream os(ccpDimensionFile, true);
  os << "fam, leaves, ccps" << std::endl;
  for (unsigned int i = 0; i < familyIndices.size(); ++i) {
    auto j = sortedIndices[i];
    os << families[familyIndices[j]].name << ", " << treeSizes[j] << ", "
       << ccpSizes[j] << std::endl;
  }
  os.close();
  ParallelContext::barrier();
}

/**
 *  Remove the serialized CCP files to free some disk space
 */
void cleanupCCPs(const Families &families) {
  ParallelContext::barrier();
  Logger::timed << "Cleaning up ccp files..." << std::endl;
  for (const auto &family : families) {
    std::remove(family.ccp.c_str());
  }
  ParallelContext::barrier();
}

/**
 *  Exclude families that:
 *  - cover too few species
 *  - are too large (optional)
 *  - have too much uncertainty (optional)
 */
void trimFamilies(const AleArguments &args, Families &families) {
  Logger::timed << "Families: " << families.size() << std::endl;
  auto minSpecies = args.minCoveredSpecies;
  auto trimRatio = args.trimFamilyRatio;
  auto maxCladeSplitRatio = args.maxCladeSplitRatio;
  if (minSpecies > 0) {
    Logger::timed << "Trimming families covering less than " << minSpecies
                  << " species..." << std::endl;
    TrimFamilies::trimMinSpeciesCoverage(families, minSpecies);
    Logger::timed << "Remaining families: " << families.size() << std::endl;
  }
  if (trimRatio > 0.0) {
    Logger::timed << "Trimming (at least) " << trimRatio * 100.0
                  << "\% top families "
                  << "sorted by the number of observed clades..." << std::endl;
    TrimFamilies::trimHighCladesNumber(families, (1.0 - trimRatio));
    Logger::timed << "Remaining families: " << families.size() << std::endl;
  }
  if (maxCladeSplitRatio >= 0.0) {
    Logger::timed
        << "Trimming families with a ratio (observed clades)/(gene nodes) > "
        << maxCladeSplitRatio << "..." << std::endl;
    TrimFamilies::trimCladeSplitRatio(families, maxCladeSplitRatio);
    Logger::timed << "Remaining families: " << families.size() << std::endl;
  }
  checkEnoughFamilies(families);
  ParallelContext::barrier();
}

/**
 *  Compute an initial species tree from
 *  the given gene families and save it to a file
 */
void computeInitialSpeciesTree(const Families &families,
                               const SpeciesTreeAlgorithm &algo,
                               const std::string &outputFile) {
  switch (algo) {
  case SpeciesTreeAlgorithm::MiniNJ:
    MiniNJ::runMiniNJ(families)->save(outputFile);
    break;
  case SpeciesTreeAlgorithm::Random:
    std::make_unique<SpeciesTree>(families)->getTree().save(outputFile);
    break;
  case SpeciesTreeAlgorithm::User:
    assert(false);
    break;
  default:
    Logger::info << "\nError: The specified species tree building algorithm "
                 << "is not implemented in AleRax\n"
                 << std::endl;
    ParallelContext::abort(0);
  }
}

/**
 *  Generate or read the starting species tree and
 *  use it to initialize the current species tree
 *
 *  The tree can be either:
 *  - random
 *  - user-defined
 *  - computed with MiniNJ
 *
 *  @param args The program arguments
 *  @param families The gene family descriptions
 */
void initStartingSpeciesTree(AleArguments &args, const Families &families) {
  Logger::timed << "Initializing starting species tree..." << std::endl;
  // set the starting tree path
  auto startingSpeciesTree =
      Paths::getSpeciesTreeFile(args.output, "starting_species_tree.newick");
  if (args.speciesTreeAlgorithm == SpeciesTreeAlgorithm::User) {
    // ensure the initial tree file is readable
    if (ParallelContext::getRank() == 0) {
      try {
        SpeciesTree reader(args.speciesTree);
      } catch (const std::exception &e) {
        Logger::error << "Error while trying to parse the species tree:"
                      << std::endl;
        Logger::info << e.what() << std::endl;
        ParallelContext::abort(153);
      }
    }
    ParallelContext::barrier();
    // read the initial tree, provide it with internal node labels
    // and save it as the starting tree
    PLLRootedTree::labelRootedTree(args.speciesTree, startingSpeciesTree);
  } else {
    // compute the initial tree and save it as the starting tree
    computeInitialSpeciesTree(families, args.speciesTreeAlgorithm,
                              startingSpeciesTree);
  }
  ParallelContext::barrier();
  // set the current tree path for the further steps
  args.speciesTree =
      Paths::getSpeciesTreeFile(args.output, "inferred_species_tree.newick");
  if (ParallelContext::getRank() == 0) {
    // read the starting tree and save it as the current tree
    SpeciesTree currentTree(startingSpeciesTree);
    if (args.transferConstraint != TransferConstaint::RELDATED) {
      // the branch lengths of the starting tree are meaningless
      // unless it is dated, so we set all branch lengths to 1.0
      currentTree.getTree().equalizeBranchLengths();
    }
    currentTree.getTree().save(args.speciesTree);
  }
  ParallelContext::barrier();
}

/**
 *  Initialize the current species tree from the checkpoint
 */
void initSpeciesTreeFromCheckpoint(AleArguments &args,
                                   const std::string &ckpDir) {
  Logger::timed << "Initializing species tree from the checkpoint..."
                << std::endl;
  // set the current tree path for the further steps
  args.speciesTree =
      Paths::getSpeciesTreeFile(args.output, "inferred_species_tree.newick");
  // read the checkpointed tree and save it as the current tree
  auto newick = AleState::readCheckpointSpeciesTree(ckpDir);
  ParallelOfstream os(args.speciesTree, true);
  os << newick << std::endl;
  os.close();
  ParallelContext::barrier();
}

/**
 *  Check that a family conditional clade probability object
 *  is valid: check that the mapping gene <-> species is bijective
 */
bool checkCCP(const FamilyInfo &family, const ConditionalClades &ccp,
              const std::unordered_set<std::string> &allSpecies) {
  GeneSpeciesMapping mapping;
  mapping.fill(family.mappingFile, family.startingGeneTree);
  auto m = mapping.getMap();
  for (const auto &pair : ccp.getCidToLeaves()) {
    auto gene = pair.second;
    auto it = m.find(gene);
    if (it == m.end()) {
      Logger::error << "Error: gene " << gene << " from family " << family.name
                    << " is not mapped to any species" << std::endl;
      return false;
    }
    auto species = it->second;
    if (allSpecies.find(species) == allSpecies.end()) {
      Logger::error << "Error: gene " << gene << " from family " << family.name
                    << " is mapped to species " << species
                    << ", but this species is not in the species tree"
                    << std::endl;
      return false;
    }
  }
  return true;
}

/**
 *  Run checkCCP on all gene families
 */
void checkCCPAndSpeciesTree(const Families &families,
                            const std::string &speciesTreePath) {
  Logger::timed << "Checking that ccps and mappings are valid..." << std::endl;
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
  ParallelContext::parallelAnd(ok);
  if (!ok) {
    Logger::info << "\nProblems with CCPs, aborting." << std::endl;
    Logger::info << "Please check the files mentioned above\n" << std::endl;
    ParallelContext::abort(0);
  }
}

RecModelInfo buildRecModelInfo(const AleArguments &args) {
  return RecModelInfo(
      ArgumentsHelper::strToRecModel(args.reconciliationModelStr), args.recOpt,
      (args.modelParametrization ==
       ModelParametrization::PER_FAMILY), // per family rates
      args.gammaCategories, args.originationStrategy, args.pruneSpeciesTree,
      false, // rooted gene tree (option specific to GeneRax)
      false, // force gene tree root (option specific to GeneRax)
      false, // mad rooting (option specific to GeneRax)
      -1.0,  // branch length threshold
      args.transferConstraint, args.noDup, args.noDL, args.noTL,
      args.fractionMissingFile, args.memorySavings);
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
    // The initial value of the origination probabilities doesn't matter,
    // because they will be normalized during the likelihood computation.
    // We set it to the average of the other parameters such that all parameters
    // are in the same range (to facilitate gradient optimization)
    res[res.dimensions() - 1] = average;
  }
  return res;
}

void initCheckpoint(const AleArguments &args, const Families &families,
                    const std::string &ckpDir) {
  FileSystem::mkdir(ckpDir, true);
  ParallelContext::barrier();
  AleState::writeCheckpointCmd(args.getCommand(), ckpDir);
  AleState::writeCheckpointFamilies(families, ckpDir);
}

void runSpeciesTreeSearch(const AleArguments &args,
                          AleOptimizer &speciesTreeOptimizer) {
  Logger::info << std::endl;
  Logger::timed << "Start the species tree optimization..." << std::endl;
  if (args.randomSpeciesRoot) {
    Logger::info << std::endl;
    Logger::timed << "Random root position!" << std::endl;
    speciesTreeOptimizer.randomizeRoot();
  }
  // Species tree search or root search
  switch (args.speciesSearchStrategy) {
  case SpeciesSearchStrategy::HYBRID:
    speciesTreeOptimizer.optimize();
    break;
  case SpeciesSearchStrategy::REROOT:
    speciesTreeOptimizer.reroot();
    break;
  case SpeciesSearchStrategy::SKIP:
    Logger::info << std::endl;
    Logger::timed << "Optimization skipped!" << std::endl;
    break;
  default:
    assert(false); // not implemented yet
    break;
  }
  Logger::info << std::endl;
  Logger::timed << "End of the species tree optimization" << std::endl;
  speciesTreeOptimizer.saveSpeciesTree();
  Logger::timed << "Final species tree topology: " << args.speciesTree
                << std::endl;
}

void runFirstFinalRateOptimization(AleOptimizer &speciesTreeOptimizer) {
  Logger::info << std::endl;
  Logger::timed << "Final model rate optimization, non-thorough..."
                << std::endl;
  speciesTreeOptimizer.optimizeModelRates(false);
}

void runDateOptimization(const AleArguments &args,
                         AleOptimizer &speciesTreeOptimizer) {
  if (!args.inferSpeciationOrders) {
    return;
  }
  Logger::info << std::endl;
  Logger::timed << "Speciation order inference..." << std::endl;
  speciesTreeOptimizer.optimizeDates(true);
}

void runTransferHighwayInference(const AleArguments &args,
                                 AleOptimizer &speciesTreeOptimizer) {
  if (!args.highways) {
    return;
  }
  Logger::info << std::endl;
  Logger::timed << "Transfer highway inference..." << std::endl;
  speciesTreeOptimizer.inferHighways(args.highwayCandidateFile,
                                     args.highwayCandidatesStep1,
                                     args.highwayCandidatesStep2);
}

void runSecondFinalRateOptimization(const AleArguments &args,
                                    AleOptimizer &speciesTreeOptimizer) {
  if (!args.skipThoroughRates) {
    Logger::info << std::endl;
    Logger::timed << "Final model rate optimization, thorough..." << std::endl;
    speciesTreeOptimizer.optimizeModelRates(true);
  }
  speciesTreeOptimizer.saveRatesAndLL();
  auto finalLL = speciesTreeOptimizer.getEvaluator().computeLikelihood();
  Logger::timed << "Final species tree likelihood: ll=" << finalLL << std::endl;
}

void runReconciliationInference(const AleArguments &args,
                                AleOptimizer &speciesTreeOptimizer) {
  if (!args.geneTreeSamples) {
    return;
  }
  Logger::info << std::endl;
  Logger::timed << "Reconciling gene trees with the species tree..."
                << std::endl;
  speciesTreeOptimizer.reconcile(args.geneTreeSamples);
}

/**
 *  Main function of AleRax once the arguments have been parsed
 *
 *  @param args The program arguments
 */
void run(AleArguments &args) {
  Random::setSeed(static_cast<unsigned int>(args.seed));
  // initializing the output dir, ccp files, species tree, model info, etc.
  bool checkpointDetected = AleOptimizer::checkpointExists(args.output);
  auto ckpDir = AleOptimizer::getCheckpointDir(args.output);
  if (checkpointDetected) {
    Logger::info << "Checkpoint detected" << std::endl;
    AleState::checkCheckpointCmd(args.getCommand(), ckpDir);
  }
  auto ccpDir = FileSystem::joinPaths(args.output, "ccps");
  if (!checkpointDetected) {
    initAleRaxDirectories(args, ccpDir);
  }
  printInitialMessage(args);
  auto ccpDimensionFile = FileSystem::joinPaths(ccpDir, "ccpdim.txt");
  auto families = FamiliesFileParser::parseFamiliesFile(args.families);
  for (auto &family : families) {
    family.ccp = FileSystem::joinPaths(ccpDir, family.name + ".ccp");
  }
  if (!checkpointDetected) {
    filterInvalidFamilies(args, families);
    generateCCPs(args, families, ccpDimensionFile);
    trimFamilies(args, families);
    initStartingSpeciesTree(args, families);
  } else {
    AleState::filterCheckpointFamilies(families, ckpDir);
    initSpeciesTreeFromCheckpoint(args, ckpDir);
  }
  checkCCPAndSpeciesTree(families, args.speciesTree);
  auto coverageFile =
      FileSystem::joinPaths(args.output, "perSpeciesCoverage.txt");
  auto fractionMissingFile =
      FileSystem::joinPaths(args.output, "perSpeciesMissing.txt");
  Family::printStats(families, args.speciesTree, coverageFile,
                     fractionMissingFile);
  // initializing the optimizer
  auto info = buildRecModelInfo(args);
  auto startingRates = buildStartingRates(args, info);
  AleOptimizer speciesTreeOptimizer(
      args.speciesTree, families, info, args.modelParametrization,
      args.optimizationClassFile, startingRates, !args.fixRates,
      args.optVerbose, args.output);
  if (!checkpointDetected) {
    initCheckpoint(args, families, ckpDir);
    speciesTreeOptimizer.setCurrentStep(AleStep::SpeciesTreeOpt);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // species tree search
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::SpeciesTreeOpt) {
    runSpeciesTreeSearch(args, speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::ModelRateOpt1);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // model parameter optimization on the final species tree topology
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::ModelRateOpt1) {
    runFirstFinalRateOptimization(speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::RelDating);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // relative species dating on the final species tree topology
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::RelDating) {
    runDateOptimization(args, speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::Highways);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // highway search on the final species tree topology
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::Highways) {
    runTransferHighwayInference(args, speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::ModelRateOpt2);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // final model parameter optimization
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::ModelRateOpt2) {
    runSecondFinalRateOptimization(args, speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::Reconciliation);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // sampling reconciled gene trees
  if (speciesTreeOptimizer.getCurrentStep() <= AleStep::Reconciliation) {
    runReconciliationInference(args, speciesTreeOptimizer);
    speciesTreeOptimizer.setCurrentStep(AleStep::End);
    speciesTreeOptimizer.saveCheckpoint();
  }
  // end of the run
  Logger::info << std::endl;
  if (args.cleanupCCP) {
    cleanupCCPs(families);
  }
  args.printWarning();
  Logger::timed << "End of AleRax execution" << std::endl;
}

/**
 *  This main will be called by each process and is the entry
 *  point to AleRax. It inits the parallel context, the logger,
 *  reads the arguments, runs AleRax and closes everything.
 *
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

int main(int argc, char **argv) {
#ifdef WITH_MPI
  // the null communicator signals that we need to call MPI_Init
  return alerax_main(argc, argv, nullptr);
#else
  // this fake communicator signals that we're not using MPI
  int noMPIComm = -1;
  return alerax_main(argc, argv, &noMPIComm);
#endif
}
