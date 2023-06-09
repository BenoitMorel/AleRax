#pragma once

#include "util/enums.hpp"
#include <string>
#include "IO/ArgumentsHelper.hpp"

  

class AleArguments {
public:
  AleArguments(int argc, char * argv[]);
  void printHelp();
  void printCommand();
  void printSummary();
public:
  
  int argc;
  char ** argv;

  // input data
  std::string families;
  std::string speciesTree;
  std::string reconciliationModelStr;
  
  // model
  TransferConstaint transferConstraint;
  OriginationStrategy originationStrategy;
  bool pruneSpeciesTree;
  bool noTL;
  unsigned int gammaCategories;
  CCPRooting ccpRooting;
  std::string fractionMissingFile;
  std::string speciesCategoryFile;
  bool memorySavings;

  // search
  SpeciesTreeAlgorithm speciesTreeAlgorithm;
  SpeciesSearchStrategy speciesSearchStrategy;
  bool inferSpeciationOrders;
  bool fixRates;
  bool skipThoroughRates;

  // highways
  bool highways;
  std::string highwayCandidateFile;
  unsigned int highwayCandidatesStep1;
  unsigned int highwayCandidatesStep2;

  // trimming
  int minCoveredSpecies;
  double trimFamilyRatio;
  double maxCladeSplitRatio;

  // output
  unsigned int geneTreeSamples;
  std::string output;
  bool cleanupCCP;

  // random seed
  int seed;
  
  // experimental
  bool randomSpeciesRoot;
  bool verboseOptRates;
};

