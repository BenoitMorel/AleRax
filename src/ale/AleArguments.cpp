#include "AleArguments.hpp"
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <climits>

const unsigned int DEFAULT_GENE_TREE_SAMPLES = 100;
const unsigned int DEFAULT_MIN_COVERED_SPECIES = 1;
const double DEFAULT_MAX_SPLIT_RATIO = -1.0;
const double DEFAULT_TRIM_FAMILY_RATIO = 0.0;
const double DEFAULT_DTL_RATES = 0.1;

AleArguments::AleArguments(int iargc, char *iargv[])
    : argc(iargc), argv(iargv), reconciliationModelStr("UndatedDTL"),
      transferConstraint(TransferConstaint::PARENTS),
      originationStrategy(OriginationStrategy::UNIFORM),
      pruneSpeciesTree(false), noTL(false), gammaCategories(1),
      ccpRooting(CCPRooting::UNIFORM),
      modelParametrization(ModelParametrization::GLOBAL), memorySavings(false),
      d(DEFAULT_DTL_RATES), l(DEFAULT_DTL_RATES), t(DEFAULT_DTL_RATES),
      speciesTreeAlgorithm(SpeciesTreeAlgorithm::User),
      speciesSearchStrategy(SpeciesSearchStrategy::SKIP),
      inferSpeciationOrders(false), fixRates(false), skipThoroughRates(false),
      recOpt(RecOpt::LBFGSB), highways(false), highwayCandidatesStep1(100),
      highwayCandidatesStep2(50), skipFamilyFiltering(false),
      minCoveredSpecies(DEFAULT_MIN_COVERED_SPECIES),
      trimFamilyRatio(DEFAULT_TRIM_FAMILY_RATIO),
      maxCladeSplitRatio(DEFAULT_MAX_SPLIT_RATIO), sampleFrequency(1),
      geneTreeSamples(DEFAULT_GENE_TREE_SAMPLES), output("alerax_output"),
      cleanupCCP(false), seed(123), randomSpeciesRoot(false),
      verboseOptRates(false) {
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
    } else if (arg == "--species-tree-search") {
      speciesSearchStrategy =
          ArgumentsHelper::strToSpeciesSearchStrategy(std::string(argv[++i]));
    } else if (arg == "--infer-speciation-order") {
      inferSpeciationOrders = true;
    } else if (arg == "-r" || arg == "--rec-model") {
      reconciliationModelStr = std::string(argv[++i]);
    } else if (arg == "--skip-thorough-rates") {
      skipThoroughRates = true;
    } else if (arg == "--rec-opt") {
      recOpt = ArgumentsHelper::strToRecOpt(std::string(argv[++i]));
    } else if (arg == "--fix-rates") {
      fixRates = true;
    } else if (arg == "--highways") {
      highways = true;
    } else if (arg == "--highway-candidate-file") {
      highwayCandidateFile = std::string(argv[++i]);
    } else if (arg == "--highway-candidate-step1") {
      highwayCandidatesStep1 = atoi(argv[++i]);
    } else if (arg == "--highway-candidate-step2") {
      highwayCandidatesStep2 = atoi(argv[++i]);
    } else if (arg == "--transfer-constraint") {
      transferConstraint =
          ArgumentsHelper::strToTransferConstraint(std::string(argv[++i]));
    } else if (arg == "--origination") {
      originationStrategy = Enums::strToOrigination(std::string(argv[++i]));
    } else if (arg == "--prune-species-tree") {
      pruneSpeciesTree = true;
    } else if (arg == "--no-tl") {
      noTL = true;
    } else if (arg == "--trim-ratio") {
      trimFamilyRatio = atof(argv[++i]);
    } else if (arg == "--skip-ratio") {
      trimFamilyRatio = atof(argv[++i]);
    } else if (arg == "--skip-family-filtering") {
      skipFamilyFiltering = true;
    } else if (arg == "--min-covered-species") {
      minCoveredSpecies = atof(argv[++i]);
    } else if (arg == "--max-clade-split-ratio") {
      maxCladeSplitRatio = atof(argv[++i]);
    } else if (arg == "--gene-tree-sample-frequency") {
      sampleFrequency = atoi(argv[++i]);
    } else if (arg == "--speciation-probability-categories") {
      gammaCategories = atoi(argv[++i]);
    } else if (arg == "--gene-tree-rooting") {
      ccpRooting = ArgumentsHelper::strToCCPRooting(std::string(argv[++i]));
    } else if (arg == "--fraction-missing-file") {
      fractionMissingFile = argv[++i];
    } else if (arg == "--gene-tree-samples") {
      geneTreeSamples = atoi(argv[++i]);
    } else if (arg == "--model-parametrization") {
      std::string temp(argv[++i]);
      modelParametrization = Enums::strToModelParametrization(temp);
      if (ModelParametrization::CUSTOM == modelParametrization) {
        optimizationClassFile = temp;
      }
    } else if (arg == "--memory-savings") {
      memorySavings = true;
    } else if (arg == "--d") {
      d = atof(argv[++i]);
    } else if (arg == "--l") {
      l = atof(argv[++i]);
    } else if (arg == "--t") {
      t = atof(argv[++i]);
    } else if (arg == "-p" || arg == "--prefix") {
      output = std::string(argv[++i]);
    } else if (arg == "--cleanup-ccp") {
      cleanupCCP = true;
    } else if (arg == "--seed") {
      seed = atoi(argv[++i]);
    } else if (arg == "--random-species-root") {
      randomSpeciesRoot = true;
    } else if (arg == "--verbose-opt-rates") {
      verboseOptRates = true;
    } else if (arg == "--per-family-rates" || arg == "--per-species-rates") {
      Logger::error
          << "Error: --per-family-rates and --per-species-rates are deprecated "
             "and have been replaced with --model-parametrization PER-FAMILY "
             "or --model-paramtetrization PER-SPECIES"
          << std::endl;
      ParallelContext::abort(10);
    } else {
      std::cerr << "Unknown argument " << arg << std::endl;
      ParallelContext::abort(10);
    }
  }
}

void AleArguments::printCommand() const {
  Logger::timed << "AleRax was called as follow:" << std::endl;
  for (int i = 0; i < argc; ++i) {
    Logger::info << argv[i] << " ";
  }
  Logger::info << std::endl << std::endl;
}

static std::string getOnOff(bool onOff) { return onOff ? "ON" : "OFF"; }

void AleArguments::printSummary() const {
  Logger::timed << "Run settings:" << std::endl;
  Logger::info << "\tFamily file: " << families << std::endl;
  switch (Enums::strToSpeciesTree(speciesTree)) {
  case SpeciesTreeAlgorithm::User:
    Logger::info << "\tStarting species tree: specified by user: "
                 << speciesTree << std::endl;
    break;
  case SpeciesTreeAlgorithm::Random:
    Logger::info << "\tStarting species tree: generated randomly " << std::endl;
    break;
  default:
    Logger::info
        << "\tStarting species tree: will be generated with the method "
        << speciesTree << std::endl;
    break;
  }
  Logger::info << "\tOutput directory: " << output << std::endl;
#ifdef WITH_MPI
  Logger::info << "\tMPI Ranks: " << ParallelContext::getSize() << std::endl;
#else
  Logger::info << "\tYou are running GeneRax without MPI (no parallelization)"
               << std::endl;
#endif
  Logger::info << "\tNumber of reconciled gene trees to sample: "
               << geneTreeSamples << std::endl;
  Logger::info << "\tRandom seed: " << seed << std::endl;
  Logger::info << "\tReconciliation model: " << reconciliationModelStr
               << std::endl;
  Logger::info << "\tModel parametrization: ";
  switch (modelParametrization) {
  case ModelParametrization::GLOBAL:
    Logger::info << " global to all species and families" << std::endl;
    break;
  case ModelParametrization::PER_SPECIES:
    Logger::info
        << " each species has a different set of rates, common to all families"
        << std::endl;
    break;
  case ModelParametrization::ORIGINATION_PER_SPECIES:
    Logger::info
        << " each species has a different set of origination probabilities. "
           "The other rates are global to all families and sepcies"
        << std::endl;
    break;
  case ModelParametrization::PER_FAMILY:
    Logger::info
        << " each family has a different set of rates, common to all species"
        << std::endl;
    break;
  case ModelParametrization::CUSTOM:
    Logger::info << " described in the user-given file: "
                 << optimizationClassFile << std::endl;
    break;
  };
  Logger::info << "Rate optimizer: " << ArgumentsHelper::recOptToStr(recOpt)
               << std::endl;
  Logger::info << "\tMemory savings: " << getOnOff(memorySavings) << std::endl;
  switch (transferConstraint) {
  case TransferConstaint::NONE:
    Logger::info << "\tTransfer constraints: no constraint" << std::endl;
    break;
  case TransferConstaint::PARENTS:
    Logger::info << "\tTransfer constraints: transfers to parents are forbidden"
                 << std::endl;
    break;
  case TransferConstaint::RELDATED:
    Logger::info
        << "\tTransfer constraints: transfers to the past are forbidden"
        << std::endl;
    break;
  }
  Logger::info << "\tPrune species mode is "
               << (pruneSpeciesTree ? "enabled" : "disabled") << std::endl;
  if (gammaCategories > 1) {
    Logger::info << "\tSpeciation probability categories: " << gammaCategories
                 << std::endl;
  }
  Logger::info << "\tGene tree rooting: ";
  switch (ccpRooting) {
  case CCPRooting::UNIFORM:
    Logger::info << "all gene tree root positions are considered with the same "
                    "probability"
                 << std::endl;
    break;
  case CCPRooting::ROOTED:
    Logger::info << "only the root of the input gene trees will be considered"
                 << std::endl;
    break;
  case CCPRooting::MAD:
    Logger::info << "all gene tree root positions are considered, with a "
                    "weight depending on MAD scores"
                 << std::endl;
    break;
  }
  if (fractionMissingFile.size()) {
    Logger::info << "\tFraction of missing gene file: " << fractionMissingFile
                 << std::endl;
  }
  Logger::info << "\tOrigination strategy: ";
  switch (originationStrategy) {
  case OriginationStrategy::UNIFORM:
    Logger::info << "gene families can originate from each species with the "
                    "same probability"
                 << std::endl;
    break;
  case OriginationStrategy::ROOT:
    Logger::info
        << "gene families only originate at the root of the species tree"
        << std::endl;
    break;
  case OriginationStrategy::LCA:
    Logger::info << "gene families only originate at the LCA of the species "
                    "that they cover"
                 << std::endl;
    break;
  case OriginationStrategy::OPTIMIZE:
    Logger::info << "gene families can originate from each species. The "
                    "origination probabilities will be estimated if "
                    "--per-species-rates or --param-opt-classes are set"
                 << std::endl;
    break;
  }
  Logger::info << "\tSpecies tree search: ";
  switch (speciesSearchStrategy) {
  case SpeciesSearchStrategy::SKIP:
  case SpeciesSearchStrategy::EVAL:
    Logger::info << "skipping species tree search" << std::endl;
    break;
  case SpeciesSearchStrategy::REROOT:
    Logger::info << "the starting species tree will only be rerooted"
                 << std::endl;
    break;
  case SpeciesSearchStrategy::HYBRID:
  case SpeciesSearchStrategy::SPR:
  case SpeciesSearchStrategy::TRANSFERS:
    Logger::info << "enabled (the species tree will be inferred from the gene "
                    "tree distributions)"
                 << std::endl;
    break;
  }
  if (inferSpeciationOrders) {
    Logger::info << "AleRax will estimate the relative order of speciation "
                    "events from the HGTs"
                 << std::endl;
  }
  Logger::info << "\tAleRax will exclude gene families covering less than "
               << minCoveredSpecies << " species" << std::endl;
  if (trimFamilyRatio != DEFAULT_TRIM_FAMILY_RATIO) {
    Logger::info << "\tAleRax will exclude a proportion of " << trimFamilyRatio
                 << " of the largest gene families" << std::endl;
  }
  if (maxCladeSplitRatio != DEFAULT_MAX_SPLIT_RATIO) {
    Logger::info << "\tAleRax will exclude the gene families for which the "
                    "ratio between the conditional clade probability size and "
                    "the number of nodes is greater than "
                 << maxCladeSplitRatio << std::endl;
  }
  Logger::info << std::endl;
}

void AleArguments::printHelp() {
  Logger::info << "Printing help message:" << std::endl;
  Logger::info << "General options:" << std::endl;
  Logger::info << "\t-h, --help" << std::endl;
  Logger::info << "\t-f, --families <FAMILIES_INFORMATION>" << std::endl;
  Logger::info << "\t-s, --species-tree <SPECIES TREE>" << std::endl;
  Logger::info << "\t-p, --prefix <OUTPUTDIR>" << std::endl;
  Logger::info << "\t--gene-tree-samples <number of samples>" << std::endl;
  Logger::info << "\t--seed <seed>" << std::endl;

  Logger::info << "Reconciliation model options:" << std::endl;
  Logger::info
      << "\t-r --rec-model <reconciliationModel>  {UndatedDL, UndatedDTL}"
      << std::endl;
  Logger::info << "\t---rec-opt <optimizer>  {Simplex, LBFGSB, Gradient}"
               << std::endl;
  Logger::info << "\t--transfer-constraint {NONE, PARENTS, RELDATED}"
               << std::endl;
  Logger::info << "\t--prune-species-tree" << std::endl;
  Logger::info << "\t--param-opt-classes filepath" << std::endl;
  Logger::info << "\t--speciation-probability-categories <VALUE>" << std::endl;
  Logger::info << "\t--gene-tree-rooting {UNIFORM, ROOTED}" << std::endl;
  Logger::info << "\t--fraction-missing-file <filepath>" << std::endl;
  Logger::info << "\t--origination {UNIFORM, ROOT, LCA} " << std::endl;
  Logger::info << "\t--fix-rates" << std::endl;
  Logger::info << "\t--per-family-rates" << std::endl;
  Logger::info << "\t--memory-savings" << std::endl;

  Logger::info << "Search strategy options:" << std::endl;
  Logger::info << "\t--species-tree-search {HYBRID, REROOT, SKIP}" << std::endl;
  Logger::info << "\t--infer-speciation-order" << std::endl;

  Logger::info << "Trimming options" << std::endl;
  Logger::info << "\t--min-covered-species <value>" << std::endl;
  Logger::info << "\t--max-clade-split-ratio <value>" << std::endl;
  Logger::info << "\t--skip-ratio <proportion>" << std::endl;
  Logger::info << "\t--gene-tree-sample-frequency <int>" << std::endl;

  Logger::info << "Transfer highway options" << std::endl;
  Logger::info << "\t--highways" << std::endl;
  Logger::info << "\t--highway-candidates-file <filepath>" << std::endl;
  Logger::info << "\t--highway-candidates-step1 <value>" << std::endl;
  Logger::info << "\t--highway-candidates-step2 <value>" << std::endl;

  Logger::info << "For a more detailed description, please check the wiki on "
                  "our github page"
               << std::endl;
}

void AleArguments::checkValid() {
  // this method is not really relevant anymore...
  bool ok = true;
  if (!ok) {
    Logger::error << "Error occured, aborting" << std::endl;
    ParallelContext::abort(1);
  }
}
