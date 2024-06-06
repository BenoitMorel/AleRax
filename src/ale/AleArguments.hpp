#pragma once

#include "IO/ArgumentsHelper.hpp"
#include "util/enums.hpp"
#include <string>

/**
 *  Parse and store the program arguments
 */
class AleArguments {
public:
  /**
   *  Parse the arguments from main()
   */
  AleArguments(int argc, char *argv[]);

  /**
   *  Print AleRax' help message
   */
  void printHelp();

  /**
   *  Print the command line used to call the AleRax' executable
   */
  void printCommand() const;

  /**
   *  Write a user-friendly summary of the most important
   *  parameters set by the user
   */
  void printSummary() const;

  /**
   *  Check that the arguments are compatible with each other
   *  If not, terminates the program with an explicit error message
   */
  void checkValid();

public:
  // the parameters of the main() function
  int argc;
  char **argv;

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
  std::string optimizationClassFile;
  ModelParametrization modelParametrization;
  bool memorySavings;
  double d;
  double l;
  double t;

  // search
  SpeciesTreeAlgorithm speciesTreeAlgorithm;
  SpeciesSearchStrategy speciesSearchStrategy;
  bool inferSpeciationOrders;
  bool fixRates;
  bool skipThoroughRates;
  RecOpt recOpt;

  // highways
  bool highways;
  std::string highwayCandidateFile;
  unsigned int highwayCandidatesStep1;
  unsigned int highwayCandidatesStep2;

  // trimming
  bool skipFamilyFiltering;
  int minCoveredSpecies;
  double trimFamilyRatio;
  double maxCladeSplitRatio;
  unsigned int sampleFrequency;

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
