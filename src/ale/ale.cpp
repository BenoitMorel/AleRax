#include "AleArguments.hpp"
#include <cstdio>
#include <ccp/ConditionalClades.hpp>
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/FamiliesFileParser.hpp>
#include "AleOptimizer.hpp"
#include "TrimFamilies.hpp"
#include <util/Paths.hpp>
#include <util/RecModelInfo.hpp>
#include <routines/Routines.hpp>
#include <routines/SlavesMain.hpp>
#include <IO/HighwayCandidateParser.hpp>

void filterInvalidFamilies(Families &families)
{
  Logger::timed << "Filtering families" << std::endl;
  Families validFamilies;
  for (const auto &family: families) {
    std::ifstream is(family.startingGeneTree);
    if (!is || is.peek() == std::ifstream::traits_type::eof()) {
      Logger::error << "Can't open input gene trees for family " << family.name << std::endl;
      continue;
    }
    validFamilies.push_back(family);  
  }
  families = validFamilies;
}


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


bool compare(unsigned int a, unsigned int b, const std::vector<unsigned int> &data)
{
      return data[a]<data[b];
}

std::vector<unsigned int> getSortedIndices(const std::vector<unsigned int> &values)
{
  std::vector<unsigned int> indices(values.size());
  std::iota(std::begin(indices), std::end(indices), 0);
  std::sort(std::rbegin(indices), std::rend(indices), std::bind(compare,  std::placeholders::_1, std::placeholders::_2, values)); 
  return indices;
}

void generateCCPs(const std::string &ccpDir, 
    const std::string &ccpDimensionFile,
    Families &families, 
    CCPRooting ccpRooting)
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
    ConditionalClades ccp(families[i].startingGeneTree, families[i].likelihoodFile, ccpRooting);
    ccp.serialize(families[i].ccp);
    familyIndices.push_back(i);
    treeSizes.push_back(ccp.getLeafNumber());
    ccpSizes.push_back(ccp.getCladesNumber());
  }
  ParallelContext::barrier();
  ParallelContext::concatenateHetherogeneousUIntVectors(familyIndices, familyIndices);
  ParallelContext::concatenateHetherogeneousUIntVectors(treeSizes, treeSizes);
  ParallelContext::concatenateHetherogeneousUIntVectors(ccpSizes, ccpSizes);
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

void cleanupCCPs(Families &families) 
{
  ParallelContext::barrier();
  Logger::timed << "Cleaning up ccp files..." << std::endl;
  for (const auto &family: families) {
    std::remove(family.ccp.c_str());
  }
  ParallelContext::barrier();
}

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

void run( AleArguments &args)
{
  Random::setSeed(static_cast<unsigned int>(args.seed));
  FileSystem::mkdir(args.output, true);
  FileSystem::mkdir(args.output + "/species_trees", true);
  std::string ccpDir = FileSystem::joinPaths(args.output, "ccps");
  FileSystem::mkdir(ccpDir, true);
  Logger::initFileOutput(FileSystem::joinPaths(args.output, "alerax"));
  auto families = FamiliesFileParser::parseFamiliesFile(args.families);
  filterInvalidFamilies(families);
  auto ccpDimensionFile = FileSystem::joinPaths(args.output, "ccpdim.txt");
  generateCCPs(ccpDir, ccpDimensionFile, families, args.ccpRooting);
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
      false, // per family rates
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
  AleOptimizer speciesTreeOptimizer(
      args.speciesTree,
      families,
      info,
      !args.fixRates,
      args.verboseOptRates,
      args.speciesCategoryFile,
      args.output);
  if (args.randomSpeciesRoot) {
    Logger::timed << "Random root position!" << std::endl;
    speciesTreeOptimizer.randomizeRoot();
  }
  switch (args.speciesSearchStrategy) {
  case SpeciesSearchStrategy::HYBRID:
    speciesTreeOptimizer.optimize();
    break;
  case SpeciesSearchStrategy::REROOT:
    speciesTreeOptimizer.optimizeModelRates(true);
    Logger::timed << "First root search, non thorough" << std::endl;
    speciesTreeOptimizer.rootSearch(2, false);
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
  if (args.inferSpeciationOrders) {
    speciesTreeOptimizer.optimizeDates(true);
    speciesTreeOptimizer.getEvaluator().computeLikelihood();
  }
  if (args.highways) {
    // let's infer highways of transfers!
    auto highwayOutput = FileSystem::joinPaths(args.output,
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
    // first filtering step
    std::vector<ScoredHighway> filteredHighways;
    speciesTreeOptimizer.filterCandidateHighwaysFast(candidateHighways, filteredHighways);
    filteredHighways.resize(std::min(filteredHighways.size(), size_t(args.highwayCandidatesStep2)));
    std::vector<ScoredHighway> bestHighways;
    speciesTreeOptimizer.selectBestHighways(filteredHighways, bestHighways);
    speciesTreeOptimizer.saveBestHighways(bestHighways,
        highwayOutput);
    auto acceptedHighwayOutput = FileSystem::joinPaths(args.output,
      "highway_accepted_highways.txt");
    std::vector<ScoredHighway> acceptedHighways;
    speciesTreeOptimizer.addHighways(bestHighways, acceptedHighways);
    speciesTreeOptimizer.saveBestHighways(acceptedHighways,
        acceptedHighwayOutput);
  }
  if (!args.skipThoroughRates) {
    speciesTreeOptimizer.optimizeModelRates(true);
  }
  Logger::timed <<"Sampling reconciled gene trees... (" << args.geneTreeSamples  << " samples)" << std::endl;
  speciesTreeOptimizer.reconcile(args.geneTreeSamples);
  speciesTreeOptimizer.saveSpeciesTree(); 
  if (args.cleanupCCP) {
    cleanupCCPs(families);
  }
  Logger::timed <<"End of the execution" << std::endl;
}

int genetegrator_main(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm); 
  Logger::init();
  Logger::timed << "AleRax v0.0.0" << std::endl; 
  AleArguments args(argc, argv); 
  args.printCommand();
  args.printSummary();
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
    return genetegrator_main(argc, argv, comm);
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


