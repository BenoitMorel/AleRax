#include "AleArguments.hpp"
#include <IO/Logger.hpp>
#include <climits>  


AleArguments::AleArguments(int iargc, char * iargv[]):
  argc(iargc),
  argv(iargv),
  reconciliationModelStr("UndatedDTL"),
  transferConstraint(TransferConstaint::PARENTS),
  originationStrategy(OriginationStrategy::UNIFORM),
  pruneSpeciesTree(false),
  noTL(false),
  gammaCategories(1),
  ccpRooting(CCPRooting::UNIFORM),
  speciesTreeAlgorithm(SpeciesTreeAlgorithm::User),
  speciesSearchStrategy(SpeciesSearchStrategy::EVAL),
  inferSpeciationOrders(false),
  fixRates(false),
  skipThoroughRates(false),
  highways(false),
  highwayCandidatesStep1(100),
  highwayCandidatesStep2(25),
  minCoveredSpecies(4),
  trimFamilyRatio(0.0),
  maxCladeSplitRatio(-1.0),
  geneTreeSamples(100),
  output("GeneTegrator"),
  cleanupCCP(true),
  seed(123),
  randomSpeciesRoot(false),
  verboseOptRates(false)
{
  if (argc == 1) {
    printHelp();
    ParallelContext::abort(0);
  }
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "-h" || arg == "--help") {
      printHelp();
      ParallelContext::abort(0);
    } else if (arg == "-f" || arg == "--families") {
      families = std::string(argv[++i]);
    } else if (arg == "-s" || arg == "--species-tree") {
      speciesTree = std::string(argv[++i]);
      speciesTreeAlgorithm = Enums::strToSpeciesTree(speciesTree);
    } else if (arg == "--infer-speciation-order") {
      inferSpeciationOrders = true;
    } else if (arg == "-r" || arg == "--rec-model") {
      reconciliationModelStr = std::string(argv[++i]);
    } else if (arg == "--skip-thorough-rates") {
      skipThoroughRates = true;
    } else if (arg == "--fix-rates") {
      fixRates = true;
    } else if (arg == "--highways") {
      highways = true;
    } else if (arg == "--highway-candidates-file") {
      highwayCandidateFile= std::string(argv[++i]);
    } else if (arg == "--highway-candidates-step1") {
      highwayCandidatesStep1 = atoi(argv[++i]);
    } else if (arg == "--highway-candidates-step2") {
      highwayCandidatesStep2 = atoi(argv[++i]);
    } else if (arg == "--transfer-constraint") {
      transferConstraint = ArgumentsHelper::strToTransferConstraint(std::string(argv[++i]));
    } else if (arg == "--origination") {
      originationStrategy = Enums::strToOrigination(std::string(argv[++i]));
    } else if (arg == "--prune-species-tree") {
      pruneSpeciesTree = true;
    } else if (arg == "--no-tl") {
      noTL = true;
    } else if (arg == "--trim-ratio") {
      trimFamilyRatio = atof(argv[++i]);
    } else if (arg == "--min-covered-species") {
      minCoveredSpecies = atof(argv[++i]);
    } else if (arg == "--max-clade-split-ratio") {
      maxCladeSplitRatio = atof(argv[++i]);
    } else if (arg == "--speciation-probability-categories") {
      gammaCategories = atoi(argv[++i]);
    } else if (arg == "--gene-tree-rooting") {
      ccpRooting = ArgumentsHelper::strToCCPRooting(std::string(argv[++i]));
    } else if (arg == "--fraction-missing") {
      fractionMissingFile = argv[++i];
    } else if (arg == "--species-categories") {
      speciesCategoryFile = argv[++i];
    } else if (arg == "--gene-tree-samples") {
      geneTreeSamples = atoi(argv[++i]);
    } else if (arg == "-p" || arg == "--prefix") {
      output = std::string(argv[++i]);
    } else if (arg == "--skip-cleanup-ccp") {
      cleanupCCP = false;
    } else if (arg == "--seed") {
      seed = atoi(argv[++i]);
    } else if (arg == "--random-species-root") {
      randomSpeciesRoot = true;
    } else if (arg == "--verbose-opt-rates") {
      verboseOptRates = true;
    } else {
      std::cerr << "Unknown argument " << arg << std::endl;
    }
  }
}

void AleArguments::printCommand() {
  Logger::info << "GeneRax was called as follow:" << std::endl;
  for (int i = 0; i < argc; ++i) {
    Logger::info << argv[i] << " ";
  }
  Logger::info << std::endl << std::endl;
}


void AleArguments::printHelp()
{
  Logger::info << "Printing help message:" << std::endl;
  Logger::info << "General options:" << std::endl;
  Logger::info << "\t-h, --help" << std::endl;
  Logger::info << "\t-f, --families <FAMILIES_INFORMATION>" << std::endl;
  Logger::info << "\t-s, --species-tree <SPECIES TREE>" << std::endl;
  Logger::info << "\t-p, --prefix <OUTPUTDIR>" << std::endl;
  Logger::info << "\t--gene-tree-samples <number of samples>" << std::endl;
  Logger::info << "\t--seed <seed>" << std::endl;
  
  Logger::info << "Reconciliation model options:" << std::endl; 
  Logger::info << "\t-r --rec-model <reconciliationModel>  {UndatedDL, UndatedDTL}" << std::endl;
  Logger::info << "\t--transfer-constraint {NONE, PARENTS, RELDATED}" << std::endl;
  Logger::info << "\t--prune-species-tree" << std::endl;
  Logger::info << "\t--speciation-probability-categories <VALUE>" << std::endl; 
  Logger::info << "\t--gene-tree-rooting {UNIFORM, ROOTED}" << std::endl;
  Logger::info << "\t--fraction-missing-file <filepath>" << std::endl;
  Logger::info << "\t--origination {UNIFORM, ROOT, LCA} "  << std::endl;

  Logger::info << "Search strategy options:" << std::endl; 
  Logger::info << "\t--species-tree-search {HYBRID, REROOT, EVAL, SKIP}" << std::endl;
  Logger::info << "\t--infer-speciation-order" << std::endl;
  
  
  Logger::info << "Trimming options" << std::endl;
  Logger::info << "\t--min-covered-species <value>" << std::endl;
  Logger::info << "\t--max-clade-split-ratio <value>" << std::endl;
  Logger::info << "\t--trim-family-ratio <proportion>" << std::endl;

  Logger::info << "Transfer highway options" << std::endl;
  Logger::info << "\t--highways" << std::endl;
  Logger::info << "\t--highway-candidates-file <filepath>" << std::endl;
  Logger::info << "\t--highway-candidates-step1 <value>" << std::endl;
  Logger::info << "\t--highway-candidates-step2 <value>" << std::endl;

  Logger::info << "For a more detailed description, please check the wiki on our github page" << std::endl;
}
