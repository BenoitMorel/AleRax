#include "AleOptimizer.hpp"
#include <IO/FileSystem.hpp>
#include <IO/IO.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <maths/Random.hpp>
#include <memory>
#include <optimizers/DTLOptimizer.hpp>
#include <search/DatedSpeciesTreeSearch.hpp>
#include <search/SpeciesSPRSearch.hpp>
#include <search/SpeciesTransferSearch.hpp>
#include <util/Paths.hpp>

static std::shared_ptr<MultiModel>
createModel(SpeciesTree &speciesTree, const FamilyInfo &family,
            const RecModelInfo &info, const AleModelParameters &modelParameters,
            const std::vector<Highway> &highways, bool highPrecision) {
  std::shared_ptr<MultiModel> model;
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
    break;
  }
  RatesVector rates;
  modelParameters.getRateVector(rates);
  model->setRates(rates);
  model->setHighways(highways);
  return model;
}

AleEvaluator::AleEvaluator(
    AleOptimizer &optimizer, SpeciesTree &speciesTree, const RecModelInfo &info,
    ModelParametrization modelParametrization,
    std::vector<AleModelParameters> &modelParameters, bool optimizeRates,
    bool optimizeVerbose, const Families &families, PerCoreGeneTrees &geneTrees,
    const std::string &speciesCategoryFile, const std::string &outputDir)
    : _optimizer(optimizer), _speciesTree(speciesTree), _info(info),
      _modelParameters(modelParameters), _optimizeRates(optimizeRates),
      _families(families), _geneTrees(geneTrees),
      _highPrecisions(_geneTrees.getTrees().size(), -1), _outputDir(outputDir),
      _optimizeVerbose(optimizeVerbose),
      _optimizationClasses(_speciesTree.getTree(), modelParametrization,
                           speciesCategoryFile, _info) {
  Logger::timed << "Initializing ccps and evaluators..." << std::endl;
  _evaluations.resize(_geneTrees.getTrees().size());
  for (unsigned int i = 0; i < _geneTrees.getTrees().size(); ++i) {
    resetEvaluation(i, false);
  }
  ParallelContext::barrier();
  unsigned int cladeNumber = 0;
  unsigned int worstFamily = 0;
  for (auto &evaluation : _evaluations) {
    cladeNumber += evaluation->getCCP().getCladesNumber();
    worstFamily = std::max(worstFamily, evaluation->getCCP().getCladesNumber());
  }
  unsigned int totalCladesNumber = cladeNumber;
  ParallelContext::maxUInt(worstFamily);
  ParallelContext::sumUInt(totalCladesNumber);
  double averageCladesNumber =
      double(totalCladesNumber) / double(ParallelContext::getSize());
  Logger::timed << "Initializing ccps finished" << std::endl;
  Logger::timed << "Total number of clades: " << totalCladesNumber << std::endl;
  Logger::timed << "Load balancing: "
                << std::min(1.0,
                            double(averageCladesNumber) / double(worstFamily))
                << std::endl;
  Logger::timed << "Recommended maximum number of cores: "
                << totalCladesNumber / worstFamily << std::endl;
}

void AleEvaluator::resetEvaluation(unsigned int i, bool highPrecision) {
  auto famIndex = _geneTrees.getTrees()[i].familyIndex;
  auto &family = _families[famIndex];
  _evaluations[i] = createModel(_speciesTree, family, _info,
                                _modelParameters[i], _highways, highPrecision);
  _highPrecisions[i] = highPrecision;
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
  for (unsigned int i = 0; i < _geneTrees.getTrees().size(); ++i) {
    if (_highPrecisions[i] != -1) {
      resetEvaluation(i, false);
    }
  }
  auto llAfter = computeLikelihoodFast();
  if (fabs(llBefore - llAfter) > 0.1) {
    Logger::info << "Likelihood changed after a reset of the precision: "
                 << std::endl;
    Logger::info << "Before: ll=" << llBefore << std::endl;
    Logger::info << "After:  ll=" << llAfter << std::endl;
  }
}

double AleEvaluator::computeLikelihoodFast() { return computeLikelihood(); }

double AleEvaluator::computeLikelihood(PerFamLL *perFamLL) {
  std::vector<double> localLL;
  if (perFamLL) {
    perFamLL->clear();
  }
  double sumLL = 0.0;
  for (unsigned int i = 0; i < _evaluations.size(); ++i) {
    auto ll = computeFamilyLikelihood(i);
    sumLL += ll;
    if (perFamLL) {
      perFamLL->push_back(ll);
    }
  }
  // printHightPrecisionCount();
  ParallelContext::sumDouble(sumLL);
  return sumLL;
}

double AleEvaluator::computeFamilyLikelihood(unsigned int i) {
  auto famIndex = _geneTrees.getTrees()[i].familyIndex;
  auto ll = _evaluations[i]->computeLogLikelihood();
  auto &family = _families[famIndex];
  if (_highPrecisions[i] == -1 && !std::isnormal(ll)) {
    // we are in low precision mode (we use double)
    // and it's not accurate enough, switch to
    // high precision mode

    resetEvaluation(i, true);
    ll = _evaluations[i]->computeLogLikelihood();
  }
  if (!std::isnormal(ll)) {
    std::cerr << "Error: ll=" << ll << " for family " << family.name
              << std::endl;
  }
  assert(std::isnormal(ll));
  /*
  if (_highPrecisions[i] >= 0 && _highPrecisions[i] % 20 == 0) {
    // we are in high precision mode, we now check if we can
    // switch to low precision mode to make computations faster
      resetEvaluation(i, false);
  }
  */
  if (_highPrecisions[i] >= 0) {
    _highPrecisions[i]++;
  }
  return ll;
}

void AleEvaluator::setAlpha(double alpha) {
  for (auto evaluation : _evaluations) {
    evaluation->setAlpha(alpha);
  }
}

void AleEvaluator::printHightPrecisionCount() {
  unsigned int high = 0;
  unsigned int low = 0;
  for (auto v : _highPrecisions) {
    if (v >= 0) {
      high++;
    } else {
      low++;
    }
  }
  ParallelContext::sumUInt(high);
  ParallelContext::sumUInt(low);
  // Logger::info << " Double: " << low << " scaledvalue: " << high <<
  // std::endl;
}

void AleEvaluator::onSpeciesTreeChange(
    const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) {
  for (auto &evaluation : _evaluations) {
    evaluation->onSpeciesTreeChange(nodesToInvalidate);
  }
}

void AleEvaluator::onSpeciesDatesChange() {
  for (auto &evaluation : _evaluations) {
    evaluation->onSpeciesDatesChange();
  }
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

/**
 *  Optimize a set of DTL parameters that are shared among gene families
 */
class DTLParametersOptimizerGlobal : public FunctionToOptimize {
public:
  DTLParametersOptimizerGlobal(AleEvaluator &evaluator)
      : _evaluator(evaluator) {}

  virtual double evaluate(Parameters &parameters) {
    if (0 !=
        parameters
            .dimensions()) { // happens if no family is assigned to this core
      parameters.ensurePositivity();
      auto fullParameters =
          _evaluator.getOptimizationClasses().getFullParameters(parameters);
      for (unsigned int i = 0; i < _evaluator.getLocalFamilyNumber(); ++i) {
        _evaluator.setFamilyParameters(i, fullParameters);
      }
    }
    auto res = _evaluator.computeLikelihood();
    parameters.setScore(res);
    return res;
  }

private:
  AleEvaluator &_evaluator;
};

class DTLFamilyParametersOptimizer : public FunctionToOptimize {
public:
  DTLFamilyParametersOptimizer(AleEvaluator &evaluator, unsigned int family)
      : _evaluator(evaluator), _family(family) {}

  void setParameters(Parameters &parameters) {
    parameters.ensurePositivity();
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

void AleEvaluator::setParameters(const std::vector<Parameters> &parameters) {
  for (unsigned int i = 0; i < getLocalFamilyNumber(); ++i) {
    setFamilyParameters(i, parameters[i]);
  }
}

void AleEvaluator::setFamilyParameters(unsigned int family,
                                       const Parameters &parameters) {
  RatesVector rateVector;
  _modelParameters[family].setParameters(parameters);
  _modelParameters[family].getRateVector(rateVector);
  _evaluations[family]->setRates(rateVector);
}

double AleEvaluator::optimizeModelRates(bool thorough) {
  double ll = 0.0;
  OptimizationSettings settings;
  settings.listeners.push_back(&_optimizer);
  if (_optimizeRates) {
    settings.strategy = _info.recOpt;
    settings.verbose = _optimizeVerbose;
    settings.factr = LBFGSBPrecision::MEDIUM;
    ll = computeLikelihood();
    Logger::timed << "[Species search] Optimizing model rates ";
    if (thorough) {
      Logger::info << "(thorough)" << std::endl;
    } else {
      Logger::info << "(light)" << std::endl;
    }
    if (!thorough) {
      settings.lineSearchMinImprovement = std::max(0.1, ll / 10000.0);
      settings.minAlpha = 0.01;
      settings.startingAlpha = 0.5;
      settings.optimizationMinImprovement = settings.lineSearchMinImprovement;
    } else {
      settings.lineSearchMinImprovement = std::max(0.1, ll / 10000.0);
      if (ll < 100.0) {
        settings.lineSearchMinImprovement = 0.01;
        settings.optimizationMinImprovement = settings.lineSearchMinImprovement;
      }
      settings.startingAlpha = 0.01;
      settings.minAlpha = 0.005;
      settings.optimizationMinImprovement = settings.lineSearchMinImprovement;
    }
    if (_info.perFamilyRates) {
      _optimizer.enableCheckpoints(false); // to avoid MPI issues
      for (unsigned int family = 0; family < _evaluations.size(); ++family) {
        DTLFamilyParametersOptimizer function(*this, family);
        auto categorizedParameters =
            getOptimizationClasses().getCompressedParameters(
                _modelParameters[family].getParameters());
        auto bestParameters = DTLOptimizer::optimizeParameters(
            function, categorizedParameters, settings);
        function.setParameters(bestParameters);
      }
      _optimizer.enableCheckpoints(true); // to avoid MPI issues
      ll = computeLikelihood();
      Logger::timed << "[Species search]   After model rate opt, ll=" << ll
                    << std::endl;
    } else {
      Logger::timed << "Free parameters: "
                    << _optimizationClasses.getFreeParameters() << std::endl;
      DTLParametersOptimizerGlobal function(*this);

      // fake parameters, in case there are no family assigned to this core
      // we can't use empty parameter vector, otherwise the optimizer returns
      // and this creates inconsistancy between MPI nodes
      Parameters categorizedParameters(3);
      if (!_modelParameters.empty()) {
        categorizedParameters =
            getOptimizationClasses().getCompressedParameters(
                _modelParameters[0].getParameters());
      }
      ParallelContext::barrier();
      auto bestParameters = DTLOptimizer::optimizeParameters(
          function, categorizedParameters, settings);
      function.evaluate(bestParameters); // set the parameters
      ll = computeLikelihood();
      Logger::timed
          << "[Species search]   After model rate opt, ll=" << ll
          << std::endl; // " rates: " << _modelParameters << std::endl;
    }
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
  auto gammaCategories = _info.gammaCategories;
  auto ll = computeLikelihood();
  if (gammaCategories == 1) {
    return ll;
  }
  double minAlpha = CORAX_OPT_MIN_ALPHA;
  double maxAlpha = CORAX_OPT_MAX_ALPHA;
  double startingAlpha = 1.0;
  double tolerance = 0.1;
  double f2x = 1.0;
  double alpha =
      corax_opt_minimize_brent(minAlpha, startingAlpha, maxAlpha, tolerance,
                               &ll, &f2x, (void *)this, &callback);
  setAlpha(alpha);
  std::vector<double> categories(_info.gammaCategories);
  corax_compute_gamma_cats(alpha, categories.size(), &categories[0],
                           CORAX_GAMMA_RATES_MEAN);
  Logger::timed << "[Species search]   After gamma cat  opt, ll=" << ll
                << std::endl;
  Logger::info << "alpha = " << alpha << std::endl;
  Logger::info << "rate categories: ";
  for (auto c : categories) {
    Logger::info << c << " ";
  }
  Logger::info << std::endl;
  return ll;
}

void AleEvaluator::getTransferInformation(
    SpeciesTree &speciesTree, TransferFrequencies &transferFrequencies,
    PerSpeciesEvents &perSpeciesEvents,
    PerCorePotentialTransfers &potentialTransfers) {
  // this is duplicated code from Routines...
  const auto labelToId = speciesTree.getTree().getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getTree().getDeterministicIdToLabel();
  const unsigned int labelsNumber = idToLabel.size();
  const VectorUint zeros(labelsNumber, 0);
  transferFrequencies.count = MatrixUint(labelsNumber, zeros);
  transferFrequencies.idToLabel = idToLabel;
  perSpeciesEvents = PerSpeciesEvents(speciesTree.getTree().getNodeNumber());
  auto infoCopy = _info;
  infoCopy.model = RecModel::UndatedDTL;
  infoCopy.originationStrategy = OriginationStrategy::UNIFORM;
  infoCopy.transferConstraint = TransferConstaint::PARENTS;
  for (const auto &geneTree : _geneTrees.getTrees()) {
    auto &family = (_families)[geneTree.familyIndex];
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    UndatedDTLMultiModel<ScaledValue> evaluation(speciesTree.getDatedTree(),
                                                 mapping, infoCopy, family.ccp);

    evaluation.computeLogLikelihood();
    // warning, this might make the random state
    // inconsistent between the MPI ranks
    // ParallelContext::makeRandConsistent() needs to be called
    // right after the loop
    std::vector<std::shared_ptr<Scenario>> scenarios;
    bool ok = evaluation.sampleReconciliations(1, scenarios);
    assert(ok);
    assert(scenarios.size() == 1);
    auto &scenario = *scenarios[0];
    scenario.countTransfers(labelToId, transferFrequencies.count);
    scenario.gatherReconciliationStatistics(perSpeciesEvents);
    potentialTransfers.addScenario(scenario);
  }
  ParallelContext::makeRandConsistent();
  for (unsigned int i = 0; i < labelsNumber; ++i) {
    ParallelContext::sumVectorUInt(transferFrequencies.count[i]);
  }
  perSpeciesEvents.parallelSum();
  assert(ParallelContext::isRandConsistent());
}

void AleEvaluator::sampleScenarios(
    unsigned int family, unsigned int samples,
    std::vector<std::shared_ptr<Scenario>> &scenarios) {
  assert(family < _evaluations.size());
  scenarios.clear();
  getEvaluation(family).computeLogLikelihood();
  bool ok = getEvaluation(family).sampleReconciliations(samples, scenarios);
  if (!ok) {
    scenarios.clear();
    resetEvaluation(family, true);
    ok = getEvaluation(family).sampleReconciliations(samples, scenarios);
    if (!ok) {
      std::cerr << "Error: cannot sample reconciliations for family "
                << _families[_geneTrees.getTrees()[family].familyIndex].name
                << std::endl;
      assert(ok);
    }
  }
}

struct ScoredString {
  ScoredString(const std::string str, double score) : str(str), score(score) {}

  bool operator<(const ScoredString &other) const {
    if (score == other.score) {
      return str < other.str;
    }
    return score < other.score;
  }
  std::string str;
  double score;
};

void AleEvaluator::savePerFamilyLikelihoodDiff(const std::string &output) {
  std::vector<unsigned int> indices;
  std::vector<double> likelihoods;
  for (unsigned int i = 0; i < _evaluations.size(); ++i) {
    auto famIndex = _geneTrees.getTrees()[i].familyIndex;
    auto ll = _evaluations[i]->computeLogLikelihood();
    indices.push_back(famIndex);
    likelihoods.push_back(ll);
  }
  std::vector<unsigned int> allIndices;
  std::vector<double> allLikelihoods;
  ParallelContext::concatenateHetherogeneousDoubleVectors(likelihoods,
                                                          allLikelihoods);
  ParallelContext::concatenateHetherogeneousUIntVectors(indices, allIndices);
  assert(allLikelihoods.size() == _snapshotPerFamilyLL.size());
  ParallelOfstream os(output);
  std::vector<ScoredString> scoredFamilies;
  for (unsigned int i = 0; i < allLikelihoods.size(); ++i) {
    auto &family = _families[allIndices[i]];
    auto ll = allLikelihoods[i];
    scoredFamilies.push_back(
        ScoredString(family.name, ll - _snapshotPerFamilyLL[i]));
  }
  std::sort(scoredFamilies.begin(), scoredFamilies.end());
  for (const auto &scoredFamily : scoredFamilies) {
    os << scoredFamily.score << " " << scoredFamily.str << std::endl;
  }
}

void AleEvaluator::saveSnapshotPerFamilyLL() {
  std::vector<double> likelihoods;
  for (unsigned int i = 0; i < _evaluations.size(); ++i) {
    auto ll = _evaluations[i]->computeLogLikelihood();
    likelihoods.push_back(ll);
  }
  ParallelContext::concatenateHetherogeneousDoubleVectors(likelihoods,
                                                          _snapshotPerFamilyLL);
}
