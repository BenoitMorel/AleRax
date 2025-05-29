#pragma once

#include <string>

#include "IO/ArgumentsHelper.hpp"
#include "util/enums.hpp"

/**
 *  Parses and stores the program arguments
 */
class AleArguments {
public:
  /**
   *  Parse the arguments from main()
   */
  AleArguments(int iargc, char *iargv[]);

  /**
   *  Print AleRax' help message
   */
  void printHelp() const;

  /**
   *  Return the command line used to call the AleRax' executable
   */
  const std::string getCommand() const;

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
   *  Print warnings if some of the arguments used are valid, but
   *  still might result in undesired behaviour
   */
  void printWarning() const;

  /**
   *  Check that the arguments are compatible with each other.
   *  If not, terminates the program with an explicit error message
   */
  void checkValid() const;

public:
  // the parameters of the main() function
  int argc;
  char **argv;

  // input data
  std::string families;
  std::string speciesTree;

  // model
  std::string reconciliationModelStr;
  TransferConstaint transferConstraint;
  bool noDL;
  bool noTL;
  unsigned int gammaCategories;
  bool pruneSpeciesTree;
  std::string fractionMissingFile;
  CCPRooting ccpRooting;
  OriginationStrategy originationStrategy;
  bool memorySavings;
  double d;
  double l;
  double t;

  // search
  SpeciesTreeAlgorithm speciesTreeAlgorithm;
  SpeciesSearchStrategy speciesSearchStrategy;
  bool inferSpeciationOrders;
  ModelParametrization modelParametrization;
  std::string optimizationClassFile;
  RecOpt recOpt;
  bool fixRates;
  bool skipThoroughRates;

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
  std::string output;
  unsigned int geneTreeSamples;
  bool cleanupCCP;

  // random seed
  int seed;

  // experimental
  bool randomSpeciesRoot;
  bool optVerbose;
};
