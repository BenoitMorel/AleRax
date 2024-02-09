#pragma once

#include <trees/SpeciesTree.hpp>
#include "AleModelParameters.hpp"

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
    speciesTree(std::make_unique<SpeciesTree>(speciesTreePath))
  {}

  // current species tree
  std::unique_ptr<SpeciesTree> speciesTree;
  // current model parameters
  AleModelParameters modelParameters;

};
