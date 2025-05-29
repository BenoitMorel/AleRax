#pragma once

#include <memory>
#include <vector>

#include <IO/Families.hpp>
#include <maths/ModelParameters.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <search/SpeciesSearchCommon.hpp>
#include <trees/PLLRootedTree.hpp>
#include <trees/SpeciesTree.hpp>

#include "AleModelParameters.hpp"
#include "OptimizationClasses.hpp"
#include "UndatedDLMultiModel.hpp"
#include "UndatedDTLMultiModel.hpp"

class RecModelInfo;
using MultiEvaluation = MultiModel;
using MultiEvaluationPtr = std::shared_ptr<MultiEvaluation>;
using PerCoreMultiEvaluations = std::vector<MultiEvaluationPtr>;
class AleOptimizer;

struct ScoredFamily {
  ScoredFamily() {}
  ScoredFamily(const std::string &familyName, double score = 0.0)
      : familyName(familyName), score(score) {}
  std::string familyName;
  double score;
  bool operator<(const ScoredFamily &other) const {
    if (score == other.score) {
      return familyName < other.familyName;
    }
    return score < other.score;
  }
};

/**
 *  Implements the main functions to evaluate the species tree likelihood
 *  and to sample reconciliations
 */
class AleEvaluator : public SpeciesTreeLikelihoodEvaluatorInterface {
public:
  /**
   *  Constructor
   *  @param optimizer
   *  @param speciesTree Reference to the current AleState species tree
   *  @param info Description of the reconciliation model
   *  @param modelParametrization Describes how parameters should be optimized
   *  @param optimizationClassFile Path to the file defining parameter
   * optimization classes
   *  @param mixtureAlpha Reference to the current AleState model alpha
   *  @param perLocalFamilyModelParams Reference to the current AleState model
   * parameters
   *  @param transferHighways Reference to the current AleState transfer
   * highways
   *  @param optimizeRates If set to false, model parameter optimization will be
   * skipped
   *  @param optimizeVerbose If set to true, the optimization routines will
   * print more logs
   *  @param families List of all gene families
   *  @param geneTrees List of the gene families assigned to the MPI rank
   *  @param outputDir AleRax output directory
   */
  AleEvaluator(AleOptimizer &optimizer, SpeciesTree &speciesTree,
               const RecModelInfo &info,
               const ModelParametrization &modelParametrization,
               const std::string &optimizationClassFile, double &mixtureAlpha,
               std::vector<AleModelParameters> &perLocalFamilyModelParams,
               std::vector<Highway> &transferHighways, bool optimizeRates,
               bool optimizeVerbose, const Families &families,
               const PerCoreGeneTrees &geneTrees, const std::string &outputDir);

  /**
   *  Destructor
   */
  virtual ~AleEvaluator() {}

  /**
   *  See description in the parent class.
   *  We don't need to implement those functions in AleRax
   */
  virtual void pushRollback() {}
  virtual void popAndApplyRollback() {}

  /**
   *  Compute the species tree likelihood
   *
   *  @param perFamLL Structure to store the per-family likelihoods (if set)
   *  @return The log likelihood
   */
  virtual double computeLikelihood(PerFamLL *perFamLL = nullptr);

  /**
   *  Fast approximated version of computeLikelihood, used by some
   *  of the optimization heuristics to speedup the search, unless
   *  providesFastLikelihoodImpl returns false (and so far it returns false)
   */
  virtual double computeLikelihoodFast();
  virtual bool providesFastLikelihoodImpl() const { return false; }

  /**
   *  Compute the likelihood of family i (using the _geneTrees indexing)
   */
  double computeFamilyLikelihood(unsigned int i);

  /**
   *  Functions to set the recmodel parameters
   */
  void setAlpha(double alpha);
  void setFamilyParameters(unsigned int family, const Parameters &parameters);

  /**
   *  Optimize the model rates
   */
  virtual double optimizeModelRates(bool thorough);

  /**
   *  Update node information after modifying the species tree
   */
  virtual void onSpeciesDatesChange();
  virtual void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);

  /**
   *  Sample reconciliations for family i (using the _geneTrees indexing)
   */
  void sampleFamilyScenarios(unsigned int i, unsigned int samples,
                             std::vector<std::shared_ptr<Scenario>> &scenarios);

  /**
   *  Sample reconciliations and return the HGT frequencies and
   *  per-species events
   *
   *  Obtain possible transfers to be used in transfer-guided
   *  species tree search, relative species dating and
   *  transfer highway inference
   */
  virtual void
  getTransferInformation(SpeciesTree &speciesTree,
                         TransferFrequencies &frequencies,
                         PerSpeciesEvents &perSpeciesEvents,
                         PerCorePotentialTransfers &potentialTransfers);

  /**
   *  Functions to handle highways
   */
  void addHighway(const Highway &highway);
  void removeHighway();
  void saveSnapshotPerFamilyLL();
  void savePerFamilyLikelihoodDiff(const std::string &outputFile);
  unsigned int getInputTreesNumber() const;

  /**
   *  Is the recmodel dated?
   */
  virtual bool isDated() const { return _info.isDated(); }

  /**
   *  Are we in the prune species tree mode?
   */
  virtual bool pruneSpeciesTree() const { return _info.pruneSpeciesTree; }

  /**
   *  Should the optimization routines print verbose logs?
   */
  virtual bool isVerbose() const { return _optimizeVerbose; }

  /**
   *  Accessors
   */
  const RecModelInfo &getRecModelInfo() const { return _info; }
  std::vector<AleModelParameters> &getModelParameters() {
    return _modelParameters;
  }
  const std::vector<Highway> &getHighways() const { return _highways; }
  MultiModel &getEvaluation(unsigned int i) { return *_evaluations[i]; }
  unsigned int getLocalFamilyNumber() const {
    return _geneTrees.getTrees().size();
  }
  std::string getOutputDir() const { return _outputDir; }
  const OptimizationClasses &getOptimizationClasses() const {
    return _optimizationClasses;
  }

private:
  // optimize the alpha parameter for speciation probability gamma categories
  double optimizeGammaRates();
  // set precision of the evaluator of the i-th family
  void resetEvaluation(unsigned int i, bool highPrecision);
  // try to set all evaluators to low precision
  void resetAllPrecisions();
  // print the number of evaluators having high precision
  void printHighPrecisionCount();

private:
  AleOptimizer &_optimizer;
  SpeciesTree &_speciesTree;
  const RecModelInfo &_info;
  // _optimizationClasses["D"][5] is the index of the optimization class
  // of the duplication probability of species branch with index 5.
  // All parameters that belong to the same optimization class share
  // the same value
  const OptimizationClasses _optimizationClasses;
  double &_mixtureAlpha;
  std::vector<AleModelParameters> &_modelParameters;
  std::vector<Highway> &_highways;
  bool _optimizeRates;
  bool _optimizeVerbose;
  // global families
  const Families &_families;
  // local families
  const PerCoreGeneTrees &_geneTrees;
  // evaluations for local families
  PerCoreMultiEvaluations _evaluations;
  PerCoreMultiEvaluations _approxEvaluations;
  std::vector<int> _highPrecisions;
  std::vector<double> _snapshotPerFamilyLL;
  const std::string _outputDir;
};
