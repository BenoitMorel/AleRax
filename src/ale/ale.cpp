#include "AleArguments.hpp"
#include <cstdio>
#include <ccp/ConditionalClades.hpp>
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Families.hpp>
#include "AleOptimizer.hpp"
#include "TrimFamilies.hpp"
#include <util/Paths.hpp>
#include <util/RecModelInfo.hpp>
#include <routines/Routines.hpp>
#include <routines/SlavesMain.hpp>
#include <IO/HighwayCandidateParser.hpp>

const char *version = "AleRax v1.0.0";


/**
 *  Check the validity of each gene family: 
 *  - check that the gene tree files exist
 *  - if set, check that the mapping files exist
 */
void filterInvalidFamilies(Families &families)
{
  Logger::timed << "Filtering families" << std::endl;
  Families validFamilies;
  for (const auto &family: families) {
    // check that the gene tree distribution file exists
    std::ifstream is(family.startingGeneTree);
    if (!is || is.peek() == std::ifstream::traits_type::eof()) {
      Logger::error << "Can't open input gene trees for family " << family.name << std::endl;
      continue;
    }
    if (family.mappingFile.size()) {
      // check that the mapping file exists if it is set
      std::ifstream isMap(family.mappingFile);
      if (!isMap || isMap.peek() == std::ifstream::traits_type::eof()) {
        Logger::error << "Can't open the mapping file for family " << family.name << std::endl;
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
bool compare(unsigned int a, unsigned int b, const std::vector<unsigned int> &data)
{
      return data[a]<data[b];
}

/**
 *  Return the indices of the input vector, sorted in descending order
 */
std::vector<unsigned int> getSortedIndices(const std::vector<unsigned int> &values)
{
  std::vector<unsigned int> indices(values.size());
  std::iota(std::begin(indices), std::end(indices), 0);
  std::sort(std::rbegin(indices), std::rend(indices), std::bind(compare,  std::placeholders::_1, std::placeholders::_2, values)); 
  return indices;
}

/**
 *  Read the gene tree distribution files (or the .ale files), compute
 *  the corresponding CCP objects, and serialize them to a binary file
 *  Each MPI rank will run this function on a subset of the families
 *
 *  @param ccpDir Output directory to store the seralized files
 *  @param ccpDimensionFile Output file with the per-family ccp sizes
 *  @param families The gene families to process
 *  @param ccpRooting Mode to decide how to deal with rooted gene trees
 *  @param sampleFrequency Number used to subsample the gene trees (if set to 2, we only use half)
 *
 *
 */
void generateCCPs(const std::string &ccpDir, 
    const std::string &ccpDimensionFile,
    Families &families, 
    CCPRooting ccpRooting,
    unsigned int sampleFrequency)
{
  ParallelContext::barrier();
  Logger::timed << "Generating ccp files..." << std::endl;
  for (auto &family: families) {
    family.ccp = FileSystem::joinPaths(ccpDir, family.name + ".ccp");
  } 
  auto N = families.size();
  
  std::vector<unsigned int> familyIndices;
  std::vector<unsigned int> treeSizes;
  std::vector<unsigned int> ccpSizes;
  for (auto i = ParallelContext::getBegin(N); i < ParallelContext::getEnd(N); i ++) {
    // Read and seralize the ccps
    ConditionalClades ccp(families[i].startingGeneTree, families[i].likelihoodFile, ccpRooting, sampleFrequency);
    ccp.serialize(families[i].ccp);
    // record the dimensions for subsequent sorting
    familyIndices.push_back(i);
    treeSizes.push_back(ccp.getLeafNumber());
    ccpSizes.push_back(ccp.getCladesNumber());
  }
  ParallelContext::barrier();
  // Gather the family ccp dimensions from all MPI ranks
  ParallelContext::concatenateHetherogeneousUIntVectors(familyIndices, familyIndices);
  ParallelContext::concatenateHetherogeneousUIntVectors(treeSizes, treeSizes);
  ParallelContext::concatenateHetherogeneousUIntVectors(ccpSizes, ccpSizes);
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
}

/**
 *  Remove the serialized CCP files to free some disk space
 */
void cleanupCCPs(Families &families) 
{
  ParallelContext::barrier();
  Logger::timed << "Cleaning up ccp files..." << std::endl;
  for (const auto &family: families) {
    std::remove(family.ccp.c_str());
  }
  ParallelContext::barrier();
}


/**
 *  Check that a family conditional clade probability object
 *  is valid: check that the mapping gene <-> species is bijective
 */
bool checkCCP(FamilyInfo &family,
    ConditionalClades &ccp, 
    std::unordered_set<std::string> &allSpecies)
{
  GeneSpeciesMapping mapping;
  mapping.fill(family.mappingFile, family.startingGeneTree);
  auto m = mapping.getMap();
  for (auto pair: ccp.getCidToLeaves()) {
    auto gene = pair.second;
    auto it = m.find(gene);
    if (it == m.end()) {
      Logger::error << "Gene " << gene << " not mapped to any species "
        << "in family "  << family.name << std::endl;
      return false;
    }
    auto species = it->second; 
    if (allSpecies.find(species) == allSpecies.end()) {
      Logger::error << "Gene " << gene << " from family " << family.name << " is mapped to species " << species << " but this species is not in the species tree" << std::endl;
      return false;
    }
  }
  return true;
}

/**
 *  Run checkCCP on all gene families
 */
void checkCCPAndSpeciesTree(Families &families,
    const std::string &speciesTreePath)
{
  PLLRootedTree speciesTree(speciesTreePath);
  auto labels = speciesTree.getLabels(true);
  auto N = families.size();
  bool ok = true;
  for (auto i = ParallelContext::getBegin(N); i < ParallelContext::getEnd(N); i ++) {
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
 *  @param minSpecies Species that cover less than minSpecies species will be trimmed
 *  @param trimRatio Proportion of families to trim, starting from the largest ones
 *  @param maxCladeSplitRatio Trim families with (ccp size) / (gene tree size) < maxCladeSplitRatio 
 */
void trimFamilies(Families &families, int minSpecies, double trimRatio,
    double maxCladeSplitRatio) 
{
  Logger::timed << "Families: " << families.size() << std::endl;
  if (minSpecies != -1) {
    Logger::timed << "Triming families covering less than " << minSpecies << " species " << std::endl;
    TrimFamilies::trimMinSpeciesCoverage(families, minSpecies);
    Logger::timed << "Remaining families: " << families.size() << std::endl;
  }
  if (trimRatio > 0.0) {
    Logger::timed << "Trimming families with too many clades (keeping " 
      << (1.0 - trimRatio) * 100.0 << "\% of the families) " << std::endl;
    TrimFamilies::trimHighCladesNumber(families, (1.0 - trimRatio));
    Logger::timed << "Remaining families: " << families.size() << std::endl;
  }
  if (maxCladeSplitRatio > 0) {
    Logger::timed << "Triming families with a ratio clades/nodes > " << maxCladeSplitRatio<< std::endl;
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
void initStartingSpeciesTree(AleArguments &args,
    Families &families)
{
  Logger::timed << "Initializing starting species tree..." << std::endl;
  auto startingSpeciesTree = Paths::getSpeciesTreeFile(
      args.output, 
      "starting_species_tree.newick");
  std::unique_ptr<PLLRootedTree> speciesTree(nullptr);
  if (args.speciesTreeAlgorithm == SpeciesTreeAlgorithm::User) {
    unsigned int canRead = 1;
    if (ParallelContext::getRank() == 0) {
      try {
        SpeciesTree reader(args.speciesTree);
      } catch (const std::exception &e) {
        Logger::info << "Error while trying to parse the species tree:" << std::endl;
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
    Routines::computeInitialSpeciesTree(families,
        args.output,
        args.speciesTreeAlgorithm)->save(startingSpeciesTree);

  }
  ParallelContext::barrier();
  args.speciesTree = Paths::getSpeciesTreeFile(args.output, "inferred_species_tree.newick");
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
    const std::string &speciesTreePath)
{
  PLLRootedTree speciesTree(speciesTreePath);
  ParallelOfstream os(perSpeciesRatesFile);
  for (auto label: speciesTree.getLabels(false)) {
    os << label << std::endl;
  }
  os.close();
  ParallelContext::barrier();
}


/**
 * Main function of AleRax once the arguments have been parsed
 *
 * @param args The program arguments
 */
void run( AleArguments &args)
{
  Random::setSeed(static_cast<unsigned int>(args.seed));
  FileSystem::mkdir(args.output, true);
  FileSystem::mkdir(args.output + "/species_trees", true);
  std::string ccpDir = FileSystem::joinPaths(args.output, "ccps");
  FileSystem::mkdir(ccpDir, true);
  Logger::initFileOutput(FileSystem::joinPaths(args.output, "alerax"));
  
  Logger::timed << version << std::endl; 
  args.printCommand();
  args.printSummary();
  
  auto families = FamiliesFileParser::parseFamiliesFile(args.families);
  if (!args.skipFamilyFiltering) {
    filterInvalidFamilies(families);
  }
  auto ccpDimensionFile = FileSystem::joinPaths(args.output, "ccpdim.txt");
  generateCCPs(ccpDir, ccpDimensionFile, families, args.ccpRooting, args.sampleFrequency);
  trimFamilies(families, args.minCoveredSpecies, args.trimFamilyRatio,
     args.maxCladeSplitRatio);
  if (families.size() == 0) {
    Logger::info << "No valid family, aborting" << std::endl;
    ParallelContext::abort(0);
  }
  initStartingSpeciesTree(args, families);
  Logger::timed << "Checking that ccp and mappings are valid..." << std::endl;
  checkCCPAndSpeciesTree(families, args.speciesTree); 
  RecModelInfo info(ArgumentsHelper::strToRecModel(args.reconciliationModelStr),
      args.perFamilyRates, // per family rates
      args.gammaCategories,
      args.originationStrategy,
      args.pruneSpeciesTree,
      false, // rooted gene tree
      false, // force gene tree root
      false, // mad rooting
      -1.0, // branch length threshold
      args.transferConstraint,
      false, // no dup
      args.noTL,
      args.fractionMissingFile,
      args.memorySavings);
  Parameters startingRates;
  switch (info.model) {
  case RecModel::UndatedDL:
    startingRates = Parameters(args.d, args.l);
    break;
  case RecModel::UndatedDTL:
    startingRates = Parameters(args.d, args.l, args.t);
    break;
  default:
    assert(false);
    break;
  }
  if (args.perSpeciesRates) {
    assert(args.speciesCategoryFile.size() == 0);
    args.speciesCategoryFile = FileSystem::joinPaths(args.output, "speciesRateCategories.txt");
    generatePerSpeciesRateFile(args.speciesCategoryFile, args.speciesTree);
  }
 
  std::string coverageFile(FileSystem::joinPaths(args.output, "fractionMissing.txt"));
  std::string fractionMissingFile(FileSystem::joinPaths(args.output, "perSpeciesCoverage.txt"));
  Family::printStats(families, args.speciesTree, coverageFile, fractionMissingFile);
  // init the optimizer
  AleOptimizer speciesTreeOptimizer(
      args.speciesTree,
      families,
      info,
      startingRates,
      !args.fixRates,
      args.verboseOptRates,
      args.speciesCategoryFile,
      args.output);
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
  speciesTreeOptimizer.optimizeModelRates(false);
  // Relative dating of the species tree
  if (args.inferSpeciationOrders) {
    speciesTreeOptimizer.optimizeDates(true);
    speciesTreeOptimizer.getEvaluator().computeLikelihood();
  }
  // Infer the highways
  if (args.highways) {
    auto highwaysOutputDir =  FileSystem::joinPaths(args.output, "highways");
  FileSystem::mkdir(highwaysOutputDir, true);
    // let's infer highways of transfers!
    auto highwayOutput = FileSystem::joinPaths(highwaysOutputDir,
      "highway_best_candidates.txt");
    std::vector<ScoredHighway> candidateHighways;
    // initial candidates
    if (args.highwayCandidateFile.size()) {
      // the user sets the candidates
      auto highways = HighwayCandidateParser::parse(args.highwayCandidateFile,
          speciesTreeOptimizer.getSpeciesTree().getTree());
      for (const auto &highway: highways) {
        candidateHighways.push_back(ScoredHighway(highway));
      }
    } else {
      // automatically search for candidates
      speciesTreeOptimizer.getCandidateHighways(candidateHighways, args.highwayCandidatesStep1);
    }
    // first filtering step: we add each highway candidate individually, set a small highway
    // probability, and keep the highway if the likelihood improves. We also sort the
    // highways per likelihood
    std::vector<ScoredHighway> filteredHighways;
    speciesTreeOptimizer.filterCandidateHighwaysFast(candidateHighways, filteredHighways);
    filteredHighways.resize(std::min(filteredHighways.size(), size_t(args.highwayCandidatesStep2)));
    std::vector<ScoredHighway> bestHighways;
    // optimize each highway probability individually. 
    // TODO: remove this step
    if (!args.highwaysSkipIndividualOptimization) {
      speciesTreeOptimizer.selectBestHighways(filteredHighways, bestHighways);
      speciesTreeOptimizer.saveBestHighways(bestHighways,
          highwayOutput);
    } else {
      bestHighways = filteredHighways;
    }
    // now optimize all highways together
    auto acceptedHighwayOutput = FileSystem::joinPaths(highwaysOutputDir,
      "highway_accepted_highways.txt");
    std::vector<ScoredHighway> acceptedHighways;
    speciesTreeOptimizer.optimizeAllHighways(bestHighways, acceptedHighways, true);
    speciesTreeOptimizer.saveBestHighways(acceptedHighways,
        acceptedHighwayOutput);
  }
  // one last round of DTL rate optimization
  if (!args.skipThoroughRates) {
    speciesTreeOptimizer.optimizeModelRates(true);
  }
  // sample reconciled gene trees
  Logger::timed <<"Sampling reconciled gene trees... (" << args.geneTreeSamples  << " samples)" << std::endl;
  speciesTreeOptimizer.reconcile(args.geneTreeSamples);
  speciesTreeOptimizer.saveSpeciesTree(); 
  speciesTreeOptimizer.saveRatesAndLL();
  if (args.cleanupCCP) {
    cleanupCCPs(families);
  }
  Logger::timed <<"End of the execution" << std::endl;
}

int alerax_main(int argc, char** argv, void* comm)
{
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

int internal_main(int argc, char** argv, void* comm)
{
  if (SlavesMain::isSlave(argc, argv)) {
    int slaveComm = -1; 
    return static_scheduled_main(argc, argv, &slaveComm);
  } else {
    return alerax_main(argc, argv, comm);
  }
}

int main(int argc, char** argv)
{
#ifdef WITH_MPI
  return internal_main(argc, argv, 0);
#else
  int noMPIComm = -1;
  return internal_main(argc, argv, &noMPIComm);
#endif
}


