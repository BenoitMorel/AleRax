#include "AleOptimizer.hpp"

#include <algorithm>
#include <fstream>
#include <memory>

#include <IO/Logger.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <parallelization/ParallelContext.hpp>
#include <search/SpeciesTransferSearch.hpp>

static MultiEvaluationPtr
createModel(SpeciesTree &speciesTree, const FamilyInfo &family,
            const RecModelInfo &info, const double alpha,
            const AleModelParameters &modelParameters,
            const std::vector<Highway> &highways, bool highPrecision) {
  std::shared_ptr<MultiModelInterface> model;
  GeneSpeciesMapping mapping;
  mapping.fill(family.mappingFile, family.startingGeneTree);
  switch (info.model) {
  case RecModel::UndatedDL:
    if (highPrecision) {
      model = std::make_shared<UndatedDLMultiModel<ScaledValue>>(
          speciesTree.getTree(), mapping, info, family.ccp);
    } else {
      model = std::make_shared<UndatedDLMultiModel<double>>(
          speciesTree.getTree(), mapping, info, family.ccp);
    }
    break;
  case RecModel::UndatedDTL:
    if (highPrecision) {
      model = std::make_shared<UndatedDTLMultiModel<ScaledValue>>(
          speciesTree.getDatedTree(), mapping, info, family.ccp);
    } else {
      model = std::make_shared<UndatedDTLMultiModel<double>>(
          speciesTree.getDatedTree(), mapping, info, family.ccp);
    }
    break;
  default:
    assert(false);
  }
  model->setAlpha(alpha);
  RatesVector rates;
  modelParameters.getRateVector(rates);
  model->setRates(rates);
  model->setHighways(highways);
  return model;
}

AleEvaluator::AleEvaluator(
    AleOptimizer &optimizer, SpeciesTree &speciesTree, const RecModelInfo &info,
    const ModelParametrization &modelParametrization,
    const std::string &optimizationClassFile, double &mixtureAlpha,
    std::vector<AleModelParameters> &perLocalFamilyModelParams,
    std::vector<Highway> &transferHighways, bool optimizeRates,
    bool optimizeVerbose, const Families &families,
    const PerCoreGeneTrees &geneTrees, const std::string &outputDir)
    : _optimizer(optimizer), _speciesTree(speciesTree), _info(info),
      _optimizationClasses(_speciesTree.getTree(), modelParametrization,
                           optimizationClassFile, _info),
      _mixtureAlpha(mixtureAlpha), _modelParameters(perLocalFamilyModelParams),
      _highways(transferHighways), _optimizeRates(optimizeRates),
      _optimizeVerbose(optimizeVerbose), _families(families),
      _geneTrees(geneTrees), _highPrecisions(getLocalFamilyNumber(), -1),
      _outputDir(outputDir) {
  Logger::timed << "Initializing ccps and evaluators..." << std::endl;
  _evaluations.resize(getLocalFamilyNumber());
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    resetEvaluation(i, false);
  }
  ParallelContext::barrier();
  unsigned int totalCladesNumber = 0;
  unsigned int worstFamilyCladesNumber = 0;
  for (auto &evaluation : _evaluations) {
    auto ccpSize = evaluation->getCCP().getCladesNumber();
    totalCladesNumber += ccpSize;
    worstFamilyCladesNumber = std::max(worstFamilyCladesNumber, ccpSize);
  }
  ParallelContext::barrier();
  ParallelContext::sumUInt(totalCladesNumber);
  ParallelContext::maxUInt(worstFamilyCladesNumber);
  double perCoreCladesNumber =
      double(totalCladesNumber) / double(ParallelContext::getSize());
  double loadBalancing = std::min(1.0, double(perCoreCladesNumber) /
                                           double(worstFamilyCladesNumber));
  unsigned int effectiveFamiliesNumber =
      totalCladesNumber / worstFamilyCladesNumber;
  Logger::timed << "Initializing ccps finished" << std::endl;
  Logger::timed << "Total number of clades: " << totalCladesNumber << std::endl;
  Logger::timed << "Load balancing: " << loadBalancing << std::endl;
  Logger::timed << "Recommended maximum number of cores: "
                << effectiveFamiliesNumber << std::endl;
}

void AleEvaluator::resetEvaluation(unsigned int i, bool highPrecision) {
  auto famIndex = _geneTrees.getTrees()[i].familyIndex;
  _evaluations[i] =
      createModel(_speciesTree, _families[famIndex], _info, _mixtureAlpha,
                  _modelParameters[i], _highways, highPrecision);
  auto ll = _evaluations[i]->computeLogLikelihood();
  if (highPrecision) {
    _highPrecisions[i] = 1;
  } else {
    _highPrecisions[i] = -1;
    if (!std::isnormal(ll)) {
      resetEvaluation(i, true);
    }
  }
}

void AleEvaluator::resetAllPrecisions() {
  auto llBefore = computeLikelihoodFast();
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    if (_highPrecisions[i] != -1) {
      resetEvaluation(i, false);
    }
  }
  auto llAfter = computeLikelihoodFast();
  if (fabs(llBefore - llAfter) > 0.1) {
    Logger::info << "Likelihood changed after lowering the precision: "
                 << std::endl;
    Logger::info << "Before: ll=" << llBefore << std::endl;
    Logger::info << "After:  ll=" << llAfter << std::endl;
  }
}

void AleEvaluator::printHighPrecisionCount() {
  unsigned int high = 0;
  unsigned int low = 0;
  for (auto v : _highPrecisions) {
    if (v >= 0) {
      high++;
    } else {
      low++;
    }
  }
  ParallelContext::barrier();
  ParallelContext::sumUInt(high);
  ParallelContext::sumUInt(low);
  Logger::info << "Families in low-precision mode: " << low
               << ", in high-precision mode: " << high << std::endl;
}

double AleEvaluator::computeLikelihoodFast() {
  // printHighPrecisionCount();
  return computeLikelihood();
}

double AleEvaluator::computeLikelihood(PerFamLL *perFamLL) {
  if (perFamLL) {
    perFamLL->clear();
  }
  double sumLL = 0.0;
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    auto ll = computeFamilyLikelihood(i);
    if (perFamLL) {
      perFamLL->push_back(ll);
    }
    sumLL += ll;
  }
  ParallelContext::barrier();
  ParallelContext::sumDouble(sumLL);
  return sumLL;
}

double AleEvaluator::computeFamilyLikelihood(unsigned int i) {
  auto ll = _evaluations[i]->computeLogLikelihood();
  if (_highPrecisions[i] == -1 && !std::isnormal(ll)) {
    // We are in the low precision mode (we use double)
    // and it's not accurate enough, switch to the high
    // precision mode and recompute the ll
    resetEvaluation(i, true);
    ll = _evaluations[i]->computeLogLikelihood();
  }
  if (!std::isnormal(ll)) {
    // Bad ll even in the high precision mode!
    Logger::error << "Error: ll=" << ll << " for family "
                  << _geneTrees.getTrees()[i].name << std::endl;
    assert(false);
  }
  /*
  if (_highPrecisions[i] >= 0 && _highPrecisions[i] % 20 == 0) {
    // We are in the high precision mode, we now check if we can
    // switch to the low precision mode to make computations faster
    resetEvaluation(i, false);
  }
  */
  if (_highPrecisions[i] >= 0) {
    _highPrecisions[i]++;
  }
  return ll;
}

void AleEvaluator::setAlpha(double alpha) {
  _mixtureAlpha = alpha;
  for (auto &evaluation : _evaluations) {
    evaluation->setAlpha(_mixtureAlpha);
  }
}

void AleEvaluator::setFamilyParameters(unsigned int family,
                                       const Parameters &parameters) {
  RatesVector rateVector;
  _modelParameters[family].setParameters(parameters);
  _modelParameters[family].getRateVector(rateVector);
  _evaluations[family]->setRates(rateVector);
}

/**
 *  Optimizes a set of DTL parameters that are shared among gene families
 */
class DTLGlobalParametersOptimizer : public FunctionToOptimize {
public:
  DTLGlobalParametersOptimizer(AleEvaluator &evaluator)
      : _evaluator(evaluator) {}
  void setParameters(Parameters &parameters) {
    parameters.ensurePositivity();
    auto fullParameters =
        _evaluator.getOptimizationClasses().getFullParameters(parameters);
    for (unsigned int i = 0; i < _evaluator.getLocalFamilyNumber(); ++i) {
      _evaluator.setFamilyParameters(i, fullParameters);
    }
  }
  virtual double evaluate(Parameters &parameters) {
    setParameters(parameters);
    auto res = _evaluator.computeLikelihood();
    parameters.setScore(res);
    return res;
  }

private:
  AleEvaluator &_evaluator;
};

/**
 *  Optimizes DTL parameters for a given family
 */
class DTLFamilyParametersOptimizer : public FunctionToOptimize {
public:
  DTLFamilyParametersOptimizer(AleEvaluator &evaluator, unsigned int family)
      : _evaluator(evaluator), _family(family) {}
  void setParameters(Parameters &parameters) {
    parameters.ensurePositivity();
    auto fullParameters =
        _evaluator.getOptimizationClasses().getFullParameters(parameters);
    _evaluator.setFamilyParameters(_family, fullParameters);
  }
  virtual double evaluate(Parameters &parameters) {
    setParameters(parameters);
    auto res = _evaluator.computeFamilyLikelihood(_family);
    parameters.setScore(res);
    return res;
  }

private:
  AleEvaluator &_evaluator;
  unsigned int _family;
};

double AleEvaluator::optimizeModelRates(bool thorough) {
  auto ll = computeLikelihood();
  if (_optimizeRates) {
    Logger::timed << "[Species search] Optimizing model rates ";
    OptimizationSettings settings;
    settings.listeners.push_back(&_optimizer);
    settings.verbose = _optimizeVerbose;
    settings.strategy = _info.recOpt;
    settings.lineSearchMinImprovement = std::max(0.1, -ll / 10000.0);
    if (!thorough) {
      Logger::info << "(light), ll=" << ll << std::endl;
      settings.minAlpha = 0.01;
      settings.startingAlpha = 0.5;
    } else {
      Logger::info << "(thorough), ll=" << ll << std::endl;
      if (-ll < 100.0) {
        settings.lineSearchMinImprovement = 0.01;
      }
      settings.minAlpha = 0.005;
      settings.startingAlpha = 0.01;
    }
    settings.optimizationMinImprovement = settings.lineSearchMinImprovement;
    settings.factr = LBFGSBPrecision::MEDIUM;
    if (_info.perFamilyRates) {
      Logger::timed << "[Species search]   Free parameters: "
                    << _optimizationClasses.getFreeParameters() << " per family"
                    << std::endl;
      _optimizer.enableCheckpoints(false); // to avoid MPI issues
      for (unsigned int family = 0; family < getLocalFamilyNumber(); ++family) {
        DTLFamilyParametersOptimizer function(*this, family);
        auto categorizedParameters =
            getOptimizationClasses().getCompressedParameters(
                _modelParameters[family].getParameters());
        auto bestParameters = DTLOptimizer::optimizeParameters(
            function, categorizedParameters, settings);
        // set the found best parameters to the family
        function.setParameters(bestParameters);
      }
      _optimizer.enableCheckpoints(true); // to avoid MPI issues
      ll = computeLikelihood();
    } else {
      Logger::timed << "[Species search]   Free parameters: "
                    << _optimizationClasses.getFreeParameters() << std::endl;
      DTLGlobalParametersOptimizer function(*this);
      auto categorizedParameters =
          getOptimizationClasses().getCompressedParameters(
              _modelParameters[0].getParameters());
      ParallelContext::barrier();
      auto bestParameters = DTLOptimizer::optimizeParameters(
          function, categorizedParameters, settings);
      // set the found best parameters to all families
      function.setParameters(bestParameters);
      ll = computeLikelihood();
    }
    Logger::timed << "[Species search]   After model rate opt, ll=" << ll
                  << std::endl;
    //            << ", rates=" << _modelParameters << std::endl;
  }
  ll = optimizeGammaRates();
  resetAllPrecisions();
  ll = computeLikelihood();
  return ll;
}

static double callback(void *p, double x) {
  auto *evaluator = (AleEvaluator *)p;
  evaluator->setAlpha(x);
  auto ll = evaluator->computeLikelihood();
  return -ll;
}

double AleEvaluator::optimizeGammaRates() {
  auto ll = computeLikelihood();
  if (_info.gammaCategories == 1) {
    return ll;
  }
  Logger::timed << "[Species search] Optimizing gamma categories" << std::endl;
  double minAlpha = CORAX_OPT_MIN_ALPHA;
  double maxAlpha = CORAX_OPT_MAX_ALPHA;
  double startingAlpha = _mixtureAlpha;
  double tolerance = 0.1;
  double fx = -ll;
  double f2x = 1.0;
  ParallelContext::barrier();
  double alpha =
      corax_opt_minimize_brent(minAlpha, startingAlpha, maxAlpha, tolerance,
                               &fx, &f2x, (void *)this, &callback);
  ll = -fx;
  setAlpha(alpha);
  std::vector<double> categories(_info.gammaCategories);
  corax_compute_gamma_cats(alpha, categories.size(), &categories[0],
                           CORAX_GAMMA_RATES_MEAN);
  Logger::timed << "[Species search]   After gamma cat  opt, ll=" << ll
                << std::endl;
  if (_optimizeVerbose) {
    Logger::info << "alpha=" << alpha << std::endl;
    Logger::info << "speciation rate categories: ";
    for (auto c : categories) {
      Logger::info << c << " ";
    }
    Logger::info << std::endl;
  }
  return ll;
}

void AleEvaluator::onSpeciesDatesChange() {
  for (auto &evaluation : _evaluations) {
    evaluation->onSpeciesDatesChange();
  }
}

void AleEvaluator::onSpeciesTreeChange(
    const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) {
  for (auto &evaluation : _evaluations) {
    evaluation->onSpeciesTreeChange(nodesToInvalidate);
  }
}

void AleEvaluator::sampleFamilyScenarios(
    unsigned int i, unsigned int samples,
    std::vector<std::shared_ptr<Scenario>> &scenarios) {
  assert(i < getLocalFamilyNumber());
  scenarios.clear();
  bool ok = getEvaluation(i).sampleReconciliations(samples, scenarios);
  if (_highPrecisions[i] == -1 && !ok) {
    // We are in the low precision mode (we use double)
    // and it's not accurate enough, switch to the high
    // precision mode and resample
    scenarios.clear();
    resetEvaluation(i, true);
    ok = getEvaluation(i).sampleReconciliations(samples, scenarios);
  }
  if (!ok) {
    // Couldn't sample even in the high precision mode!
    Logger::error << "Error: cannot sample reconciliations for family "
                  << _geneTrees.getTrees()[i].name << std::endl;
    assert(false);
  }
}

void AleEvaluator::getTransferInformation(
    SpeciesTree &speciesTree, TransferFrequencies &transferFrequencies,
    PerSpeciesEvents &perSpeciesEvents,
    PerCorePotentialTransfers &potentialTransfers) {
  // this is duplicated code from Routines...
  const auto labelToId = speciesTree.getTree().getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getTree().getDeterministicIdToLabel();
  const unsigned int labelsNumber = idToLabel.size();
  transferFrequencies.count =
      MatrixUint(labelsNumber, VectorUint(labelsNumber, 0));
  transferFrequencies.idToLabel = idToLabel;
  perSpeciesEvents = PerSpeciesEvents(speciesTree.getTree().getNodeNumber());
  auto infoCopy = _info;
  infoCopy.model = RecModel::UndatedDTL;
  infoCopy.originationStrategy = OriginationStrategy::UNIFORM;
  infoCopy.transferConstraint = TransferConstaint::PARENTS;
  for (const auto &geneTree : _geneTrees.getTrees()) {
    const auto &family = _families[geneTree.familyIndex];
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    UndatedDTLMultiModel<ScaledValue> evaluation(speciesTree.getDatedTree(),
                                                 mapping, infoCopy, family.ccp);
    std::vector<std::shared_ptr<Scenario>> scenarios;
    // Warning:
    // Using Random::getProba() in the sampling function makes
    // the random state inconsistent between the MPI ranks.
    // Call ParallelContext::makeRandConsistent() right after
    // all MPI ranks passed the loop
    bool ok = evaluation.sampleReconciliations(1, scenarios);
    assert(ok);
    assert(scenarios.size() == 1);
    auto &scenario = *scenarios[0];
    scenario.countTransfers(labelToId, transferFrequencies.count);
    scenario.gatherReconciliationStatistics(perSpeciesEvents);
    potentialTransfers.addScenario(scenario);
  }
  ParallelContext::barrier();
  ParallelContext::makeRandConsistent();
  for (unsigned int i = 0; i < labelsNumber; ++i) {
    ParallelContext::sumVectorUInt(transferFrequencies.count[i]);
  }
  perSpeciesEvents.parallelSum();
  assert(ParallelContext::isRandConsistent());
}

void AleEvaluator::addHighway(const Highway &highway) {
  _highways.push_back(highway);
  for (auto &evaluation : _evaluations) {
    evaluation->setHighways(_highways);
  }
}

void AleEvaluator::removeHighway() {
  _highways.pop_back();
  for (auto &evaluation : _evaluations) {
    evaluation->setHighways(_highways);
  }
}

void AleEvaluator::saveSnapshotPerFamilyLL() {
  std::vector<double> localLikelihoods;
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    auto ll = computeFamilyLikelihood(i);
    localLikelihoods.push_back(ll);
  }
  ParallelContext::barrier();
  ParallelContext::concatenateHetherogeneousDoubleVectors(localLikelihoods,
                                                          _snapshotPerFamilyLL);
  assert(_snapshotPerFamilyLL.size() == _families.size());
}

void AleEvaluator::savePerFamilyLikelihoodDiff(const std::string &outputFile) {
  std::vector<unsigned int> localIndices;
  std::vector<double> localLikelihoods;
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    auto famIndex = _geneTrees.getTrees()[i].familyIndex;
    auto ll = computeFamilyLikelihood(i);
    localIndices.push_back(famIndex);
    localLikelihoods.push_back(ll);
  }
  ParallelContext::barrier();
  std::vector<unsigned int> indices;
  std::vector<double> likelihoods;
  ParallelContext::concatenateHetherogeneousUIntVectors(localIndices, indices);
  ParallelContext::concatenateHetherogeneousDoubleVectors(localLikelihoods,
                                                          likelihoods);
  assert(indices.size() == _snapshotPerFamilyLL.size());
  if (ParallelContext::getRank() == 0) {
    std::vector<ScoredFamily> scoredFamilies;
    for (unsigned int i = 0; i < indices.size(); ++i) {
      const auto &family = _families[indices[i]];
      auto ll = likelihoods[i];
      scoredFamilies.push_back(
          ScoredFamily(family.name, ll - _snapshotPerFamilyLL[i]));
    }
    std::sort(scoredFamilies.begin(), scoredFamilies.end());
    std::ofstream os(outputFile);
    os << "fam, llDiff" << std::endl;
    for (const auto &sf : scoredFamilies) {
      os << sf.familyName << ", " << sf.score << std::endl;
    }
    os.close();
  }
  ParallelContext::barrier();
}

unsigned int AleEvaluator::getInputTreesNumber() const {
  unsigned int totalInputTrees = 0;
  for (auto &evaluation : _evaluations) {
    totalInputTrees += evaluation->getCCP().getInputTreesNumber();
  }
  ParallelContext::barrier();
  ParallelContext::sumUInt(totalInputTrees);
  return totalInputTrees;
}
