#pragma once

#include <vector>
#include <string>

#include <trees/SpeciesTree.hpp>
#include "AleModelParameters.hpp"


enum class AleStep {
  Init = 0, SpeciesTreeOpt, ModelRateOpt, RelDating, Highways, Reconciliation
};

/**
 *  Store all information about the current state of an AleRax run:
 *  - the current species tree
 *  - the current model parameters
 */ 
struct AleState {
 
  /**
   *  Constructor
   */
  AleState(const std::string &speciesTreePath):
    currentStep(AleStep::Init),
    speciesTree(std::make_unique<SpeciesTree>(speciesTreePath))
  {}

  /**
   *  Dump the current state to an output file
   */
  void serialize(const std::string &outputFilePath);
  
  /**
   *  Load the current state from an input file
   */
  void unserialize(const std::string &inputFilePath);

  // running step
  AleStep currentStep;
  // current species tree
  std::unique_ptr<SpeciesTree> speciesTree;
  // current model parameters
  std::vector<AleModelParameters> perFamilyModelParameters;
};
