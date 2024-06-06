#pragma once

#include "UndatedDLMultiModel.hpp"
#include "UndatedDTLMultiModel.hpp"
#include <IO/FamiliesFileParser.hpp>
#include <maths/ModelParameters.hpp>
#include <memory>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <search/SpeciesRootSearch.hpp>
#include <trees/PLLRootedTree.hpp>
#include <trees/SpeciesTree.hpp>
#include <vector>

#include "AleModelParameters.hpp"
#include "OptimizationClasses.hpp"
#include <search/UFBoot.hpp>

class RecModelInfo;
using MultiEvaluation = MultiModel;
using MultiEvaluationPtr = std::shared_ptr<MultiEvaluation>;
using PerCoreMultiEvaluation = std::vector<MultiEvaluationPtr>;
using ParamOptGlasses = std::unordered_map<char, std::vector<unsigned int>>;
class AleOptimizer;

/**
 *  Implement the main functions to evaluate the species tree likelihood
 *  and to sample reconciliations
 */
class AleEvaluator : public SpeciesTreeLikelihoodEvaluatorInterface {
public:
  /**
   *  Constructor
   *  @param optimizer
   *  @param speciesTree The instance of the species tree
   *  @param info
   *  @param modelParametrization Describes how parameters should be optimized
   *  @param modelRates The instance of the model parameters
   *  @param optimizeRates If set to false, model parameter optimization will be
   * skipped
   *  @param optimizeVerbose If set to true, the optimization routines will
   * print more logs
   *  @param families The gene family descriptions
   *  @param geneTrees List of the gene families assigned to the MPI rank
   *  @param outputDir AleRax' output directory
   */
  AleEvaluator(AleOptimizer &optimizer, SpeciesTree &speciesTree,
               const RecModelInfo &info,
               ModelParametrization modelParametrization,
               std::vector<AleModelParameters> &modelParameters,
               bool optimizeRates, bool optimizeVerbose,
               const Families &families, PerCoreGeneTrees &geneTrees,
               const std::string &speciesCategoryFile,
               const std::string &outputDir);
  /**
   *  Destructor
   */
  virtual ~AleEvaluator() {}

  /**
   *  Compute the species tree likelihood
   *
   *  @brief perFamLL Structure to store the per-family likelihoods (if set)
   *  @return The log likelihood
   */
  virtual double computeLikelihood(PerFamLL *perFamLL = nullptr);

  /**
   *  Fast approximated version of computeLikelihood, used by some
   *  of the optimization heuristics to speedup the search, unless
   *  providesFastLikelihoodImpl return false (and so far it returns false)
   */
  virtual double computeLikelihoodFast();
  virtual bool providesFastLikelihoodImpl() const { return false; }

  /**
   * Compute the likelihood of family i (using the _geneTrees indexing)
   */
  double computeFamilyLikelihood(unsigned int i);

  /**
   *  Optimize the model rates
   */
  virtual double optimizeModelRates(bool thorough);

  /**
   *  Return true if the reconciliation model is dated
   */
  virtual bool isDated() const { return _info.isDated(); }

  /**
   *  Sample reconciliations and returns the HGT frequencies and
   *  per-species events (see
   *
   *  See description in the parent class
   */
  virtual void
  getTransferInformation(SpeciesTree &speciesTree,
                         TransferFrequencies &frequencies,
                         PerSpeciesEvents &perSpeciesEvents,
                         PerCorePotentialTransfers &potentialTransfers);

  /**
   *  Are we in prune species tree mode?
   */
  virtual bool pruneSpeciesTree() const {
    return getRecModelInfo().pruneSpeciesTree;
  }

  const RecModelInfo &getRecModelInfo() const { return _info; }

  /**
   * Set the alpha parameter of the gamma function for the family rate
   * categories
   */
  virtual void setAlpha(double alpha);
  virtual void setParameters(const std::vector<Parameters> &parameters);
  virtual void setFamilyParameters(unsigned int family,
                                   const Parameters &parameters);
  virtual void onSpeciesDatesChange();
  virtual void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);
  /**
   *  See description in the parent class
   *  We don't need to implement those functions in AleRax
   */
  virtual void pushRollback() {}
  virtual void popAndApplyRollback() {}
  void printHightPrecisionCount();
  MultiModel &getEvaluation(unsigned int i) { return *_evaluations[i]; }

  void addHighway(const Highway &highway);
  const std::vector<Highway> &getHighways() const { return _highways; }
  void removeHighway();
  void sampleScenarios(unsigned int family, unsigned int samples,
                       std::vector<std::shared_ptr<Scenario>> &scenarios);
  std::string getOutputDir() const { return _outputDir; }
  void savePerFamilyLikelihoodDiff(const std::string &output);
  void saveSnapshotPerFamilyLL();
  std::vector<AleModelParameters> &getModelParameters() {
    return _modelParameters;
  }
  unsigned int getLocalFamilyNumber() const {
    return _geneTrees.getTrees().size();
  }
  const std::vector<unsigned int> &getSpeciesToCat() const {
    return _speciesToCat;
  }
  const OptimizationClasses &getOptimizationClasses() const {
    return _optimizationClasses;
  }

protected:
  virtual double optimizeGammaRates();
  void resetEvaluation(unsigned int i, bool highPrecision);
  /**
   *  Tries to set all precisions to low precision
   */
  void resetAllPrecisions();

private:
  AleOptimizer &_optimizer;
  SpeciesTree &_speciesTree;
  const RecModelInfo &_info;
  std::vector<AleModelParameters> &_modelParameters;
  bool _optimizeRates;
  std::vector<Highway> _highways;
  const Families &_families;
  PerCoreMultiEvaluation _evaluations;
  PerCoreMultiEvaluation _approxEvaluations;
  PerCoreGeneTrees &_geneTrees;
  std::vector<int> _highPrecisions;
  std::string _outputDir;
  std::vector<double> _snapshotPerFamilyLL;
  bool _optimizeVerbose;

  // _optimizationClasses["D"][5] is the index of the optimization class
  // of the duplication probability of species branch with index 5.
  // All parameters that belong to the same optimization class share
  // the same value
  OptimizationClasses _optimizationClasses;

  std::vector<unsigned int> _speciesToCat;
};
