#pragma once

#include <string>
#include <vector>

#include "AleModelParameters.hpp"
#include <trees/SpeciesTree.hpp>

/**
 *  Represents the current step in the AleRax pipeline,
 *  to restart from the right place after a checkpoint
 *  The values should be ordered and contiguous
 */
enum class AleStep {
  Init = 0,
  SpeciesTreeOpt = 1,
  ModelRateOpt = 2,
  RelDating = 3,
  Highways = 4,
  ModelRateOpt2 = 5,
  Reconciliation = 6
};

/**
 *  Store all information about the current state of an AleRax run:
 *  - the current step in the pipeline
 *  - the current species tree
 *  - the current model parameters
 */
struct AleState {

  /**
   *  Constructor to build the AleState when not restarting
   *  from a checkpoint
   */
  AleState(const std::string &speciesTreePath)
      : currentStep(AleStep::Init),
        speciesTree(std::make_unique<SpeciesTree>(speciesTreePath)) {}

  /**
   *  Dump the current state to an output directory (checkpoint)
   */
  void serialize(const std::string &outputDir) const;

  /**
   *  Load the current state from an input directory (checkpoint
   */
  void unserialize(const std::string &inputDir,
                   const std::vector<std::string> &perCoreFamilyNames);

  /**
   *  Read the species tree newick from the dumped current state
   *  (from the checkpoint directory)
   */
  static std::string readSpeciesTreeNewick(const std::string inputDir);

  // running step
  AleStep currentStep;
  // current species tree
  std::unique_ptr<SpeciesTree> speciesTree;
  // the name of the families, to map the model parameters to their
  // respective families
  std::vector<std::string> familyNames;
  // current model parameters
  std::vector<AleModelParameters> perFamilyModelParameters;
};
