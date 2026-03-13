#pragma once

#include <memory>

#include <IO/FileSystem.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <util/RecModelInfo.hpp>

#include "AleEvaluator.hpp"
#include "AleState.hpp"
#include "Highways.hpp"

class AleOptimizer : public SpeciesTree::Listener,
                     public SpeciesSearchState::Listener,
                     public DTLOptimizerListener {
public:
  /**
   *  Constructor
   */
  AleOptimizer(const std::string &speciesTreeFile, const Families &families,
               const RecModelInfo &info,
               const ModelParametrization &modelParametrization,
               const std::string &optimizationClassFile,
               const Parameters &startingRates, bool optimizeRates,
               bool optimizeVerbose, const std::string &outputDir);

  /**
   *  Species tree search functions
   */
  // Randomly pick a branch in the species tree and reroot it
  // at this branch
  void randomizeRoot();
  // Optimize the species tree topology
  void optimize();
  // Optimize the species tree rooting
  void reroot();

  /**
   *  Parameter optimization functions
   */
  // Optimize the model parameters
  double optimizeModelRates(bool thorough = false);
  // Optimize the relative order of speciation events
  void optimizeDates(bool thorough = true);

  /**
   *  Callback functions
   */
  // Callback called when the species tree node dates change
  virtual void onSpeciesDatesChange();
  // Callback called when the species tree topology changes
  virtual void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);
  // Callback called when the species tree improves (the likelihood)
  virtual void betterTreeCallback();
  // Callback called when better model parameters have been found
  virtual void onBetterParametersFoundCallback();

  /**
   *  Saving functions
   */
  // Save the species tree
  void saveSpeciesTree();
  // Save the species tree with its branch support values
  void saveSpeciesBranchSupports();
  // Save the species tree with its root support values
  void saveSpeciesRootSupports();
  // Save the DTL rates and the per-family likelihoods
  void saveRatesAndLL();

  /**
   *  Sample reconciliations and generate many output files
   */
  void reconcile(unsigned int samples);

  /**
   *  Functions to handle transfer highways
   */
  void inferHighways(const std::string &highwayCandidateFile,
                     unsigned int highwayCandidatesStep1,
                     unsigned int highwayCandidatesStep2);
  std::string getHighwaysOutputDir() const {
    return FileSystem::joinPaths(_outputDir, "highways");
  }

  /**
   *  Functions to handle checkpoints
   */
  void saveCheckpoint() { _state.serialize(getCheckpointDir(_outputDir)); };
  bool checkpointExists() const { return checkpointExists(_outputDir); }
  static bool checkpointExists(const std::string &outputDir) {
    return FileSystem::dirExists(getCheckpointDir(outputDir));
  }
  static std::string getCheckpointDir(const std::string &outputDir) {
    return FileSystem::joinPaths(outputDir, "checkpoint");
  }
  void enableCheckpoints(bool enable) { _enableCheckpoints = enable; }

  /**
   *  Functions to interact with the AleRax workflow
   */
  AleStep getCurrentStep() const { return _state.currentStep; }
  void setCurrentStep(AleStep step) { _state.currentStep = step; }

  /**
   *  Accessors
   */
  // general info
  unsigned int getLocalFamilyNumber() const {
    return _geneTrees.getTrees().size();
  }
  const RecModelInfo &getRecModelInfo() const { return _info; }
  // species tree and model parameters
  SpeciesTree &getSpeciesTree() { return *_state.speciesTree; }
  double &getMixtureAlpha() { return _state.mixtureAlpha; }
  std::vector<Highway> &getTransferHighways() {
    return _state.transferHighways;
  }
  std::vector<AleModelParameters> &getModelParameters() {
    return _state.perLocalFamilyModelParams;
  }
  // likelihood evaluator
  AleEvaluator &getEvaluator() { return *_evaluator; }

private:
  double rootSearch(unsigned int maxDepth, bool thorough = false);
  double transferSearch();
  double sprSearch(unsigned int radius);
  void saveGeneConsensusTree(const std::string &geneSampleFile,
                             const std::string &outputFile);
  void saveFamiliesTakingHighway(const Highway &highway,
                                 const VectorDouble &perFamilyTransfers,
                                 const std::string &directory);
  void saveBestHighways(const std::vector<ScoredHighway> &highways,
                        const std::string &outputFile);

private:
  // all families
  const Families _families;
  // families of the current MPI rank
  const PerCoreGeneTrees _geneTrees;
  const RecModelInfo _info;
  RootLikelihoods _rootLikelihoods;
  const std::string _outputDir;
  bool _enableCheckpoints;
  // helper classes
  AleState _state;
  std::unique_ptr<AleEvaluator> _evaluator;
  std::unique_ptr<SpeciesSearchState> _speciesTreeSearchState;
};
