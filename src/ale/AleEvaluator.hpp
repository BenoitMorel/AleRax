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

#include "AleModelParameters.hpp"
#include <search/UFBoot.hpp>

class RecModelInfo;
using MultiEvaluation = MultiModel;
using MultiEvaluationPtr = 
  std::shared_ptr<MultiEvaluation>;
using PerCoreMultiEvaluation = std::vector<MultiEvaluationPtr>;



class AleEvaluator: public SpeciesTreeLikelihoodEvaluatorInterface {
public:
  AleEvaluator(SpeciesTree &speciesTree,
      AleModelParameters &modelRates, 
      bool optimizeRates,
      bool optimizeVerbose,
      const Families &families,
      PerCoreGeneTrees &geneTrees,
      const std::string &outputDir);
  virtual ~AleEvaluator() {}
  virtual double computeLikelihood(PerFamLL *perFamLL = nullptr); 
  virtual double computeLikelihoodFast();
  virtual bool providesFastLikelihoodImpl() const {return false;}
  virtual bool isDated() const {return _modelRates.getInfo().isDated();}
  virtual double optimizeModelRates(bool thorough);
  virtual void pushRollback() {}
  virtual void popAndApplyRollback() {}
  virtual void getTransferInformation(SpeciesTree &speciesTree,
    TransferFrequencies &frequencies,
    PerSpeciesEvents &perSpeciesEvents,
    PerCorePotentialTransfers &potentialTransfers);
  virtual bool pruneSpeciesTree() const {return _modelRates.getInfo().pruneSpeciesTree;}
  virtual void setAlpha(double alpha);
  virtual void setParameters(Parameters &parameters);
  virtual void onSpeciesDatesChange();  
  virtual void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);
  void printHightPrecisionCount();
  MultiModel &getEvaluation(unsigned int i) {return *_evaluations[i];}

  void addHighway(const Highway &highway);
  void removeHighway();
  void sampleScenarios(unsigned int family, unsigned int samples,
      std::vector< std::shared_ptr<Scenario> > &scenarios);
  std::string getOutputDir() const {return _outputDir;}
  void savePerFamilyLikelihoodDiff(const std::string &output); 
  void saveSnapshotPerFamilyLL();
protected:
  virtual double optimizeGammaRates();
  void resetEvaluation(unsigned int i, bool highPrecision);
private:
  SpeciesTree &_speciesTree;
  AleModelParameters &_modelRates;
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
};


