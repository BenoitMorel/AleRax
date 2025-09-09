#pragma once

#include <string>
#include <vector>

#include <IO/Families.hpp>
#include <IO/HighwayCandidateParser.hpp>
#include <trees/SpeciesTree.hpp>

#include "AleModelParameters.hpp"

/**
 *  Represents the current step in the AleRax pipeline,
 *  to restart from the right place after a checkpoint.
 *  The values should be ordered and contiguous
 */
enum class AleStep {
  Init = 0,
  SpeciesTreeOpt = 1,
  ModelRateOpt1 = 2,
  RelDating = 3,
  Highways = 4,
  ModelRateOpt2 = 5,
  Reconciliation = 6,
  End = 7
};

/**
 *  Stores all information about the current state of an AleRax run:
 *  - the current step in the pipeline
 *  - the current species tree
 *  - the current model parameters
 *
 *  Implements functions to handle checkpoints
 */
struct AleState {
  /**
   *  Constructor
   */
  AleState(const std::string &speciesTreePath)
      : currentStep(AleStep::Init),
        speciesTree(std::make_unique<SpeciesTree>(speciesTreePath)),
        mixtureAlpha(1.0) {}

  /**
   *  Dump the current run arguments to the checkpoint directory
   */
  static void writeCheckpointCmd(const std::string &currentCmd,
                                 const std::string &checkpointDir);

  /**
   *  Make sure the checkpoint used the same arguments as
   *  the current run
   */
  static void checkCheckpointCmd(const std::string &currentCmd,
                                 const std::string &checkpointDir);

  /**
   *  Dump the list of the accepted family names to the checkpoint
   *  directory
   */
  static void writeCheckpointFamilies(const Families &families,
                                      const std::string &checkpointDir);

  /**
   *  Retain only the families listed in the checkpoint
   */
  static void filterCheckpointFamilies(Families &families,
                                       const std::string &checkpointDir);

  /**
   *  Dump the current state to the checkpoint directory
   */
  void serialize(const std::string &checkpointDir) const;

  /**
   *  Load the current state from the checkpoint directory
   */
  void unserialize(const std::string &checkpointDir,
                   const std::vector<std::string> &perCoreFamilyNames);

  /**
   *  Read the species tree newick from the dumped current state
   *  (from the checkpoint directory)
   */
  static std::string
  readCheckpointSpeciesTree(const std::string &checkpointDir);

  // the running step
  AleStep currentStep;
  // the current species tree
  std::unique_ptr<SpeciesTree> speciesTree;
  // the names of the families, to map the model parameters to their
  // respective families
  std::vector<std::string> localFamilyNames;
  // the current model parameters
  double mixtureAlpha;
  std::vector<AleModelParameters> perLocalFamilyModelParams;
  std::vector<Highway> transferHighways;
};
