#include "AleArguments.hpp"

#include <climits>

#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>

const unsigned int DEFAULT_HIGHWAY_CANDIDATES_1 = 100;
const unsigned int DEFAULT_HIGHWAY_CANDIDATES_2 = 50;
const unsigned int DEFAULT_GENE_TREE_SAMPLES = 100;
const unsigned int DEFAULT_MIN_COVERED_SPECIES = 4;
const double DEFAULT_MAX_SPLIT_RATIO = -1.0;
const double DEFAULT_TRIM_FAMILY_RATIO = 0.0;
const double DEFAULT_DTL_RATES = 0.1;

AleArguments::AleArguments(int iargc, char *iargv[])
    : argc(iargc), argv(iargv), reconciliationModelStr("UndatedDTL"),
      transferConstraint(TransferConstaint::PARENTS), noTL(false),
      gammaCategories(1), pruneSpeciesTree(false),
      ccpRooting(CCPRooting::UNIFORM),
      originationStrategy(OriginationStrategy::UNIFORM), memorySavings(false),
      d(DEFAULT_DTL_RATES), l(DEFAULT_DTL_RATES), t(DEFAULT_DTL_RATES),
      speciesTreeAlgorithm(SpeciesTreeAlgorithm::User),
      speciesSearchStrategy(SpeciesSearchStrategy::SKIP),
      inferSpeciationOrders(false),
      modelParametrization(ModelParametrization::GLOBAL),
      recOpt(RecOpt::LBFGSB), fixRates(false), skipThoroughRates(false),
      highways(false), highwayCandidatesStep1(DEFAULT_HIGHWAY_CANDIDATES_1),
      highwayCandidatesStep2(DEFAULT_HIGHWAY_CANDIDATES_2),
      skipFamilyFiltering(false),
      minCoveredSpecies(DEFAULT_MIN_COVERED_SPECIES),
      trimFamilyRatio(DEFAULT_TRIM_FAMILY_RATIO),
      maxCladeSplitRatio(DEFAULT_MAX_SPLIT_RATIO), sampleFrequency(1),
      output("alerax_output"), geneTreeSamples(DEFAULT_GENE_TREE_SAMPLES),
      cleanupCCP(false), seed(123), randomSpeciesRoot(false),
      optVerbose(false) {
  if (argc == 1) {
    printHelp();
    ParallelContext::abort(0);
  }
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "-h" || arg == "--help") {
      printHelp();
      ParallelContext::abort(0);
    } else if (arg == "-p" || arg == "--prefix") {
      output = std::string(argv[++i]);
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
    } else if (arg == "--gene-tree-samples") {
      geneTreeSamples = atoi(argv[++i]);
    } else if (arg == "-r" || arg == "--rec-model") {
      reconciliationModelStr = std::string(argv[++i]);
    } else if (arg == "--transfer-constraint") {
      transferConstraint =
          ArgumentsHelper::strToTransferConstraint(std::string(argv[++i]));
    } else if (arg == "--no-tl") {
      noTL = true;
    } else if (arg == "--memory-savings") {
      memorySavings = true;
    } else if (arg == "--d") {
      d = atof(argv[++i]);
    } else if (arg == "--l") {
      l = atof(argv[++i]);
    } else if (arg == "--t") {
      t = atof(argv[++i]);
    } else if (arg == "--model-parametrization") {
      std::string temp(argv[++i]);
      modelParametrization = Enums::strToModelParametrization(temp);
      if (ModelParametrization::CUSTOM == modelParametrization) {
        optimizationClassFile = temp;
      }
    } else if (arg == "--rate-optimizer") {
      recOpt = ArgumentsHelper::strToRecOpt(std::string(argv[++i]));
    } else if (arg == "--fix-rates") {
      fixRates = true;
    } else if (arg == "--skip-thorough-rates") {
      skipThoroughRates = true;
    } else if (arg == "--speciation-gamma-categories") {
      gammaCategories = atoi(argv[++i]);
    } else if (arg == "--prune-species-tree") {
      pruneSpeciesTree = true;
    } else if (arg == "--fraction-missing-file") {
      fractionMissingFile = std::string(argv[++i]);
    } else if (arg == "--gene-tree-rooting") {
      ccpRooting = ArgumentsHelper::strToCCPRooting(std::string(argv[++i]));
    } else if (arg == "--origination") {
      originationStrategy = Enums::strToOrigination(std::string(argv[++i]));
    } else if (arg == "--highways") {
      highways = true;
    } else if (arg == "--highway-candidates-file") {
      highwayCandidateFile = std::string(argv[++i]);
    } else if (arg == "--highway-candidates-step1") {
      highwayCandidatesStep1 = atoi(argv[++i]);
    } else if (arg == "--highway-candidates-step2") {
      highwayCandidatesStep2 = atoi(argv[++i]);
    } else if (arg == "--skip-family-filtering") {
      skipFamilyFiltering = true;
    } else if (arg == "--min-covered-species") {
      minCoveredSpecies = atof(argv[++i]);
    } else if (arg == "--trim-ratio") {
      trimFamilyRatio = atof(argv[++i]);
    } else if (arg == "--max-clade-split-ratio") {
      maxCladeSplitRatio = atof(argv[++i]);
    } else if (arg == "--gene-tree-sample-frequency") {
      sampleFrequency = atoi(argv[++i]);
    } else if (arg == "--cleanup-ccp") {
      cleanupCCP = true;
    } else if (arg == "--seed") {
      seed = atoi(argv[++i]);
    } else if (arg == "--random-species-root") {
      randomSpeciesRoot = true;
    } else if (arg == "--verbose-opt") {
      optVerbose = true;
    } else if (arg == "--per-family-rates" || arg == "--per-species-rates") {
      Logger::info << "\nError: --per-family-rates and --per-species-rates are "
                      "deprecated and have been replaced "
                   << "with --model-parametrization PER-FAMILY or "
                      "--model-paramtetrization PER-SPECIES\n"
                   << std::endl;
      ParallelContext::abort(0);
    } else {
      Logger::info << "\nUnknown argument " << arg << "\n" << std::endl;
      ParallelContext::abort(0);
    }
  }
}

const std::string AleArguments::getCommand() const {
  auto cmd = std::string(argv[0]);
  for (int i = 1; i < argc; ++i) {
    cmd += " " + std::string(argv[i]);
  }
  return cmd;
}

void AleArguments::printCommand() const {
  Logger::timed << "AleRax was called as follows:" << std::endl;
  Logger::info << getCommand() << std::endl << std::endl;
}

static std::string getOnOff(bool onOff) { return onOff ? "ON" : "OFF"; }

void AleArguments::printSummary() const {
  Logger::timed << "Run settings:" << std::endl;
  Logger::info << "\tOutput directory: " << output << std::endl;
  Logger::info << "\tFamilies information: " << families << std::endl;
  Logger::info << "\tStarting species tree: ";
  switch (Enums::strToSpeciesTree(speciesTree)) {
  case SpeciesTreeAlgorithm::User:
    Logger::info << "will be imported from the user file: " << speciesTree
                 << std::endl;
    break;
  case SpeciesTreeAlgorithm::Random:
    Logger::info << "will be generated randomly" << std::endl;
    break;
  default:
    Logger::info << "will be generated with the method " << speciesTree
                 << std::endl;
    break;
  }
#ifdef WITH_MPI
  Logger::info << "\tMPI Ranks: " << ParallelContext::getSize() << std::endl;
#else
  Logger::info << "\tYou are running AleRax without MPI (no parallelization)"
               << std::endl;
#endif
  Logger::info << "\tRandom seed: " << seed << std::endl;
  Logger::info << "\tReconciliation model: " << reconciliationModelStr
               << std::endl;
  Logger::info << "\tTransfer constraint: ";
  switch (transferConstraint) {
  case TransferConstaint::NONE:
    Logger::info << "no constraint" << std::endl;
    break;
  case TransferConstaint::PARENTS:
    Logger::info << "transfers to parents are forbidden" << std::endl;
    break;
  case TransferConstaint::RELDATED:
    Logger::info << "transfers to the past are forbidden" << std::endl;
    break;
  }
  Logger::info << "\tMemory savings: " << getOnOff(memorySavings) << std::endl;
  Logger::info << "\tModel parametrization: ";
  switch (modelParametrization) {
  case ModelParametrization::GLOBAL:
    Logger::info << "rates are global to all species and families" << std::endl;
    break;
  case ModelParametrization::PER_FAMILY:
    Logger::info
        << "each family has a different set of rates, common to all species"
        << std::endl;
    break;
  case ModelParametrization::PER_SPECIES:
    Logger::info
        << "each species has a different set of rates, common to all families"
        << std::endl;
    break;
  case ModelParametrization::ORIGINATION_PER_SPECIES:
    Logger::info << "each species has a different set of origination "
                    "probabilities, other rates are global"
                 << std::endl;
    break;
  case ModelParametrization::CUSTOM:
    Logger::info << "clades of species sharing specific common rates are "
                    "described in the user file: "
                 << optimizationClassFile << std::endl;
    break;
  };
  Logger::info << "\tRate optimizer: " << ArgumentsHelper::recOptToStr(recOpt)
               << std::endl;
  if (gammaCategories > 1) {
    Logger::info << "\tSpeciation probability categories: " << gammaCategories
                 << std::endl;
  }
  Logger::info << "\tPrune species mode is "
               << (pruneSpeciesTree ? "enabled" : "disabled") << std::endl;
  if (fractionMissingFile.size()) {
    Logger::info << "\tFraction of missing gene file: " << fractionMissingFile
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
    Logger::info << "only the roots of the input gene trees are considered"
                 << std::endl;
    break;
  case CCPRooting::MAD:
    Logger::info << "all gene tree root positions are considered with weights "
                    "depending on their MAD scores"
                 << std::endl;
    break;
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
    Logger::info
        << "gene families can originate from each species, "
        << "the per-species origination probabilities will be estimated"
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
    Logger::info << "the starting species tree will be rerooted based on the "
                    "gene tree distributions"
                 << std::endl;
    break;
  case SpeciesSearchStrategy::HYBRID:
  case SpeciesSearchStrategy::SPR:
  case SpeciesSearchStrategy::TRANSFERS:
    Logger::info << "the optimal species tree topology will be inferred from "
                    "the gene tree distributions"
                 << std::endl;
    break;
  }
  if (inferSpeciationOrders) {
    Logger::info << "\tAleRax will estimate the relative order of speciation "
                    "events from the inferred HGTs"
                 << std::endl;
  }
  if (highways) {
    Logger::info << "\tAleRax will infer highways of horizontal gene transfer"
                 << std::endl;
    if (highwayCandidateFile.size()) {
      Logger::info
          << "\tCandidate highways will be obtained from the user file: "
          << highwayCandidateFile << std::endl;
    } else {
      Logger::info << "\tCandidate highways will be obtained from "
                      "reconciliations, max candidates: "
                   << highwayCandidatesStep1 << std::endl;
    }
    Logger::info << "\tMax number of the best highways to retain after "
                    "filtering the candidates: "
                 << highwayCandidatesStep2 << std::endl;
  }
  Logger::info << "\tNumber of reconciled gene trees to sample: "
               << geneTreeSamples << std::endl;
  Logger::info << "\tAleRax will exclude gene families covering less than "
               << minCoveredSpecies << " species" << std::endl;
  if (trimFamilyRatio != DEFAULT_TRIM_FAMILY_RATIO) {
    Logger::info
        << "\tAleRax will exclude top " << trimFamilyRatio * 100
        << "\% of the gene families with the highest numbers of observed clades"
        << std::endl;
  }
  if (maxCladeSplitRatio != DEFAULT_MAX_SPLIT_RATIO) {
    Logger::info << "\tAleRax will exclude gene families with the observed "
                    "clades to gene nodes ratio greater than "
                 << maxCladeSplitRatio << std::endl;
  }
  if (sampleFrequency != 1) {
    Logger::info << "\tAleRax will use only every i-th tree from the given "
                    "gene tree distributions to obtain "
                 << "observed clades, i = " << sampleFrequency << std::endl;
  }
  printWarning();
  Logger::info << std::endl;
}

void AleArguments::printHelp() const {
  Logger::info << "Printing help message:" << std::endl;
  Logger::info << "General options:" << std::endl;
  Logger::info << "\t-h, --help" << std::endl;
  Logger::info << "\t-f, --families <FAMILIES FILE>" << std::endl;
  Logger::info
      << "\t-s, --species-tree <SPECIES TREE> {Random, MiniNJ, <filepath>}"
      << std::endl;
  Logger::info << "\t-p, --prefix <OUTPUT DIR>" << std::endl;
  Logger::info << "\t--gene-tree-samples <number of samples>" << std::endl;
  Logger::info << "\t--seed <seed>" << std::endl;

  Logger::info << "Reconciliation model options:" << std::endl;
  Logger::info
      << "\t-r --rec-model <reconciliationModel> {UndatedDL, UndatedDTL}"
      << std::endl;
  Logger::info << "\t--transfer-constraint {NONE, PARENTS, RELDATED}"
               << std::endl;
  Logger::info << "\t--model-parametrization {GLOBAL, PER-FAMILY, PER-SPECIES, "
                  "ORIGINATION-PER-SPECIES, <filepath>}"
               << std::endl;
  Logger::info << "\t--speciation-gamma-categories <number of categories>"
               << std::endl;
  Logger::info << "\t--origination {UNIFORM, ROOT, LCA, OPTIMIZE}" << std::endl;
  Logger::info << "\t--rate-optimizer {LBFGSB, GRADIENT, SIMPLEX, GSL_SIMPLEX}"
               << std::endl;
  Logger::info << "\t--fix-rates" << std::endl;
  Logger::info << "\t--verbose-opt" << std::endl;
  Logger::info << "\t--gene-tree-rooting {UNIFORM, MAD, ROOTED}" << std::endl;
  Logger::info << "\t--prune-species-tree" << std::endl;
  Logger::info << "\t--fraction-missing-file <filepath>" << std::endl;
  Logger::info << "\t--memory-savings" << std::endl;

  Logger::info << "Search strategy options:" << std::endl;
  Logger::info << "\t--species-tree-search {HYBRID, REROOT, SKIP}" << std::endl;
  Logger::info << "\t--infer-speciation-order" << std::endl;

  Logger::info << "Trimming options:" << std::endl;
  Logger::info << "\t--min-covered-species <int>" << std::endl;
  Logger::info << "\t--trim-ratio <proportion>" << std::endl;
  Logger::info << "\t--max-clade-split-ratio <value>" << std::endl;
  Logger::info << "\t--gene-tree-sample-frequency <int>" << std::endl;

  Logger::info << "Transfer highway options:" << std::endl;
  Logger::info << "\t--highways" << std::endl;
  Logger::info << "\t--highway-candidates-file <filepath>" << std::endl;
  Logger::info << "\t--highway-candidates-step1 <int>" << std::endl;
  Logger::info << "\t--highway-candidates-step2 <int>" << std::endl;

  Logger::info << "For a more detailed description please check the wiki on "
                  "our github page:"
               << std::endl;
  Logger::info << "https://github.com/BenoitMorel/AleRax/wiki" << std::endl;
  Logger::info << std::endl;
}

void AleArguments::printWarning() const {
  if (transferConstraint == TransferConstaint::NONE) {
    Logger::info << std::string(112, '!') << std::endl;
    Logger::info << "! WARNING: --transfer-constraint NONE is biologically "
                    "nonsensical and "
                 << "should be used for testing purposes only !" << std::endl;
    Logger::info << std::string(112, '!') << std::endl;
  }
}

void AleArguments::checkValid() const {
  bool ok = true;
  // species tree block
  if (speciesSearchStrategy == SpeciesSearchStrategy::HYBRID) {
    if (transferConstraint == TransferConstaint::RELDATED) {
      ok = false;
      Logger::info << "\nError: species tree search is not compatible "
                   << "with --transfer-constraint RELDATED\n"
                   << std::endl;
    }
    if (modelParametrization == ModelParametrization::CUSTOM) {
      ok = false;
      Logger::info << "\nError: species tree search is not compatible "
                   << "with --model-parametrization <filename>\n"
                   << std::endl;
    }
  }
  if (speciesSearchStrategy == SpeciesSearchStrategy::REROOT) {
    if (transferConstraint == TransferConstaint::RELDATED &&
        !inferSpeciationOrders) {
      ok = false;
      Logger::info << "\nError: species tree root inference is not compatible "
                   << "with --transfer-constraint RELDATED unless "
                   << "using --infer-speciation-order\n"
                   << std::endl;
    }
    if (modelParametrization == ModelParametrization::CUSTOM) {
      ok = false;
      Logger::info << "\nError: species tree root inference is not compatible "
                   << "with --model-parametrization <filename>\n"
                   << std::endl;
    }
  }
  // recmodel block
  if (reconciliationModelStr != "UndatedDTL") {
    if (transferConstraint != TransferConstaint::PARENTS) {
      ok = false;
      Logger::info << "\nError: --transfer-constraint can only be used "
                   << "with --rec-model UndatedDTL\n"
                   << std::endl;
    }
    if (highways) {
      ok = false;
      Logger::info << "\nError: --highways can only be used "
                   << "with --rec-model UndatedDTL\n"
                   << std::endl;
    }
  }
  // missing species block
  if (pruneSpeciesTree && fractionMissingFile.size()) {
    ok = false;
    Logger::info << "\nError: --fraction-missing-file cannot be used "
                 << "with --prune-species-tree\n"
                 << std::endl;
  }
  // model parametrization block
  if ((modelParametrization == ModelParametrization::GLOBAL ||
       modelParametrization == ModelParametrization::PER_FAMILY) &&
      originationStrategy == OriginationStrategy::OPTIMIZE) {
    ok = false;
    Logger::info << "\nError: --origination OPTIMIZE cannot be used "
                 << "with --model-parametrization GLOBAL or "
                 << "--model-parametrization PER-FAMILY\n"
                 << std::endl;
  }
  if (modelParametrization == ModelParametrization::ORIGINATION_PER_SPECIES &&
      originationStrategy != OriginationStrategy::OPTIMIZE) {
    ok = false;
    Logger::info << "\nError: --model-parametrization ORIGINATION-PER-SPECIES "
                    "should only be used "
                 << "with --origination OPTIMIZE\n"
                 << std::endl;
  }
  // reldating block
  if (transferConstraint != TransferConstaint::RELDATED &&
      inferSpeciationOrders) {
    ok = false;
    Logger::info << "\nError: --infer-speciation-order can only be used "
                 << "with --transfer-constraint RELDATED\n"
                 << std::endl;
  }
  // highway block
  if (!highways) {
    if (highwayCandidateFile.size()) {
      ok = false;
      Logger::info << "\nError: --highway-candidates-file can only be used "
                   << "with --highways\n"
                   << std::endl;
    }
    if (highwayCandidatesStep1 != DEFAULT_HIGHWAY_CANDIDATES_1) {
      ok = false;
      Logger::info << "\nError: --highway-candidate-step1 can only be used "
                   << "with --highways\n"
                   << std::endl;
    }
    if (highwayCandidatesStep2 != DEFAULT_HIGHWAY_CANDIDATES_2) {
      ok = false;
      Logger::info << "\nError: --highway-candidate-step2 can only be used "
                   << "with --highways\n"
                   << std::endl;
    }
  }
  if (highwayCandidateFile.size() &&
      highwayCandidatesStep1 != DEFAULT_HIGHWAY_CANDIDATES_1) {
    ok = false;
    Logger::info << "\nError: --highway-candidate-step1 cannot be used "
                 << "with --highway-candidates-file\n"
                 << std::endl;
  }
  // final
  if (!ok) {
    ParallelContext::abort(0);
  }
}
