#pragma once

#include <search/SpeciesRootSearch.hpp>
#include <trees/SpeciesTree.hpp>
#include <IO/FamiliesFileParser.hpp>
#include "UndatedDLMultiModel.hpp"
#include "UndatedDTLMultiModel.hpp"
#include <trees/PLLRootedTree.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <memory>
#include <vector>
#include <maths/ModelParameters.hpp>
#include "AleEvaluator.hpp"

struct ScoredHighway {
  ScoredHighway() 
  {}

  
  ScoredHighway(const Highway &highway, double score = 0.0, double scoreDiff = 0.0):
    highway(highway),
    score(score),
    scoreDiff(scoreDiff)
  {}

  Highway highway;
  double score;
  double scoreDiff;
  bool operator < (const ScoredHighway &other) {
    return score < other.score;
  }
};


class AleOptimizer: public SpeciesTree::Listener {
public:
  AleOptimizer(const std::string speciesTreeFile, 
      const Families &families, 
      const RecModelInfo &info,
      bool optimizeRates,
      bool optimizeVerbose,
      const std::string &speciesCategoryFile,
      const std::string &outputDir);

  void optimize();
  double sprSearch(unsigned int radius);
  double rootSearch(unsigned int maxDepth, bool thorough = false);
  double transferSearch();
  void onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);
  void reconcile(unsigned int samples);
  double optimizeModelRates(bool thorough = false);
  void optimizeDates(bool thorough = true);
  SpeciesTree &getSpeciesTree() {return *_speciesTree;}
  void randomizeRoot();
  void saveSpeciesTree();
  void saveSupportTree();
  SpeciesTreeLikelihoodEvaluatorInterface &getEvaluator() {return *_evaluator;}
  void getCandidateHighways(std::vector<ScoredHighway> &highways, unsigned int toTest);
  void filterCandidateHighwaysFast(const std::vector<ScoredHighway> &highways, std::vector<ScoredHighway> &filteredHighways);
  void selectBestHighways(const std::vector<ScoredHighway> &highways, std::vector<ScoredHighway> &bestHighways);
  void saveBestHighways(const std::vector<ScoredHighway> &highways,
      const std::string &output);
  void addHighways(const std::vector<ScoredHighway> &candidateHighways,
      std::vector<ScoredHighway> &acceptedHighways);
  void saveRELLSupports();
private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  PerCoreGeneTrees _geneTrees;
  RecModelInfo _info;
  std::unique_ptr<AleEvaluator> _evaluator;
  AleModelParameters _modelRates;
  std::string _outputDir;
  SpeciesSearchState _searchState;
  RootLikelihoods _rootLikelihoods;
  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
};



