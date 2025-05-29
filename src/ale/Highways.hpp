#pragma once

#include <vector>

#include <IO/HighwayCandidateParser.hpp>

class AleOptimizer;

struct ScoredHighway {
  ScoredHighway() {}
  ScoredHighway(const Highway &highway, double score = 0.0)
      : highway(highway), score(score) {}
  Highway highway;
  double score;
  bool operator<(const ScoredHighway &other) const {
    return score < other.score;
  }
  bool operator==(const Highway &other) const {
    return (highway.src == other.src) && (highway.dest == other.dest);
  }
};

/**
 *  Implements the main functions to infer transfer highways, to compute their
 *  highway probabilities and to add them to the reconciliation model
 */
class Highways {
public:
  /**
   *  Include the provided highways to the potential candidates in the order
   *  defined by the frequencies of the transfer directions they represent
   *  in the undated reconciliations
   */
  static void
  getSortedCandidatesFromList(AleOptimizer &optimizer,
                              const std::vector<Highway> &highways,
                              std::vector<ScoredHighway> &candidateHighways);

  /**
   *  Infer potential highway candidates from the undated reconciliations
   *  by sampling the most frequent transfer directions
   */
  static void
  getCandidateHighways(AleOptimizer &optimizer,
                       std::vector<ScoredHighway> &candidateHighways,
                       unsigned int maxCandidates);

  /**
   *  Filter the highway candidates by separately testing each of them with
   *  a small hardcoded highway proba. Highways resulting in higher likelihood
   *  increase are retained
   */
  static void
  filterCandidateHighways(AleOptimizer &optimizer,
                          const std::vector<ScoredHighway> &candidateHighways,
                          std::vector<ScoredHighway> &filteredHighways,
                          unsigned int maxCandidates);

  /**
   *  Optimize jointly the highway probas of the filtered candidates. Highways
   *  with significant resulting highway probas are added to the recmodel
   */
  static void
  optimizeAllHighways(AleOptimizer &optimizer,
                      const std::vector<ScoredHighway> &filteredHighways,
                      std::vector<ScoredHighway> &acceptedHighways,
                      bool thorough);
};
