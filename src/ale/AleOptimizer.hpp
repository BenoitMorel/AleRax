#pragma once

#include "AleEvaluator.hpp"
#include "AleState.hpp"
#include "UndatedDLMultiModel.hpp"
#include "UndatedDTLMultiModel.hpp"
#include <IO/FamiliesFileParser.hpp>
#include <IO/FileSystem.hpp>
#include <maths/ModelParameters.hpp>
#include <memory>
#include <optimizers/DTLOptimizer.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <search/SpeciesRootSearch.hpp>
#include <trees/PLLRootedTree.hpp>
#include <trees/SpeciesTree.hpp>
#include <vector>

struct ScoredHighway {
  ScoredHighway() {}

  ScoredHighway(const Highway &highway, double score = 0.0,
                double scoreDiff = 0.0)
      : highway(highway), score(score), scoreDiff(scoreDiff) {}

  Highway highway;
  double score;
  double scoreDiff;
  bool operator<(const ScoredHighway &other) const {
    return score < other.score;
  }
  bool operator==(const Highway &other) const {
    return highway.src == other.src && highway.dest == other.dest;
  }
};

bool cmpHighwayByProbability(const ScoredHighway &a, const ScoredHighway &b);

class AleOptimizer : public SpeciesTree::Listener,
                     public SpeciesSearchState::Listener,
                     public DTLOptimizerListener {
public:
  AleOptimizer(const std::string speciesTreeFile, const Families &families,
               const RecModelInfo &info,
               ModelParametrization modelParametrization,
               const Parameters &startingRates, bool optimizeRates,
               bool optimizeVerbose, const std::string &optimizationClassFile,
               const std::string &outputDir);

  /**
   *  Optimize the species tree topology
   */
  void optimize();

  void enableCheckpoints(bool enable) { _enableCheckpoints = enable; }

  /**
   *  Optize the species tree root
   */
  double rootSearch(unsigned int maxDepth, bool thorough = false);

  /**
   *  Callback called when the species tree topology changes
   */
  void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);

  /**
   *  Callback called when the species tree improves (the likelihood)
   */
  virtual void betterTreeCallback();

  /**
   *  Callback called when better model parameters have been found
   */
  virtual void onBetterParametersFoundCallback();

  /**
   *  Sample reconciliations and generate many output files
   */
  void reconcile(unsigned int samples);

  /**
   *  Optimize the model parameters
   */
  double optimizeModelRates(bool thorough = false);

  /**
   *  Optimize the relative order of speciation events
   */
  void optimizeDates(bool thorough = true);

  /**
   *  Randomly pick a branch in the species tree and reroot it
   *  at this branch
   */
  void randomizeRoot();

  /**
   *  Save the species tree
   */
  void saveSpeciesTree();

  /**
   *  Save the species tree with its support values
   */
  void saveSupportTree();

  /**
   *  Save the DTL rates and the per-family likelihoods
   */
  void saveRatesAndLL();

  /**
   *  Accessor
   */
  AleEvaluator &getEvaluator() { return *_evaluator; }

  /**
   *  Accessor
   */
  SpeciesTree &getSpeciesTree() { return *_state.speciesTree; }

  /**
   *  Accessor
   */
  std::vector<AleModelParameters> &getModelParameters() {
    return _state.perFamilyModelParameters;
  }
  const std::vector<AleModelParameters> &getModelParameters() const {
    return _state.perFamilyModelParameters;
  }

  const RecModelInfo &getRecModelInfo() const { return _info; }

  void saveBestHighways(const std::vector<ScoredHighway> &highways,
                        const std::string &output);
  void saveRELLSupports();
  std::string getHighwaysOutputDir() const;

  void saveCheckpoint() const;
  void loadCheckpoint();
  bool checkpointExists() const { return checkpointExists(_outputDir); }
  static bool checkpointExists(const std::string &outputDir);
  static std::string getCheckpointDir(const std::string &outputDir) {
    return FileSystem::joinPaths(outputDir, "checkpoint");
  }
  AleStep getCurrentStep() const { return _state.currentStep; }
  void setCurrentStep(AleStep step) { _state.currentStep = step; }

private:
  AleState _state;
  const Families &_families;
  PerCoreGeneTrees _geneTrees;
  RecModelInfo _info;
  std::unique_ptr<AleEvaluator> _evaluator;
  std::string _outputDir;
  std::string _checkpointDir;
  std::unique_ptr<SpeciesSearchState> _speciesTreeSearchState;
  RootLikelihoods _rootLikelihoods;
  bool _enableCheckpoints;
  double sprSearch(unsigned int radius);
  double transferSearch();
  std::string
  saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick",
                           bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str,
                                  bool masterRankOnly = true);
};
