#pragma once

#include <vector>

#include "AleOptimizer.hpp"

class Highways {
public:
  /**
   *  Infer potential highway candidates from the undated reconciliations
   */
  static void getCandidateHighways(AleOptimizer &optimizer,
                                   std::vector<ScoredHighway> &highways,
                                   unsigned int maxCandidates);

  /**
   *  Filter the highway candidates by testing them with a small hardcoded
   * probability
   */
  static void
  filterCandidateHighwaysFast(AleOptimizer &optimizer,
                              const std::vector<ScoredHighway> &highways,
                              std::vector<ScoredHighway> &filteredHighways);

  static void
  optimizeAllHighways(AleOptimizer &optimizer,
                      const std::vector<ScoredHighway> &candidateHighways,
                      std::vector<ScoredHighway> &acceptedHighways,
                      bool thorough);
};
