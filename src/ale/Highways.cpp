#include "Highways.hpp"

#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <cmath>
#include <cstdio>
#include <optimizers/DTLOptimizer.hpp>
#include <search/SpeciesTransferSearch.hpp>

// const double MIN_PH = 0.00000001;
// const double MAX_PH = 0.8;

class HighwayFunction : public FunctionToOptimize {
public:
  HighwayFunction(AleEvaluator &evaluator,
                  const std::vector<Highway *> &highways)
      : _highways(highways), _evaluator(evaluator) {}

  virtual double evaluate(Parameters &parameters) {
    double v = evaluatePrint(parameters, false);
    // Logger::timed << "Evaluate transfer " << std::setprecision(17) <<
    // parameters << std::endl;
    return v;
  }

  virtual double evaluatePrint(Parameters &parameters, bool print,
                               const std::string outputDir = "") {
    assert(parameters.dimensions() == _highways.size());
    // parameters.constrain(MIN_PH, MAX_PH);
    for (unsigned int i = 0; i < _highways.size(); ++i) {
      Highway highwayCopy = *_highways[i];
      highwayCopy.proba = parameters[i];
      _evaluator.addHighway(highwayCopy);
    }
    auto ll = _evaluator.computeLikelihood();
    if (print) {
      assert(outputDir.size());
      std::string out = FileSystem::joinPaths(
          outputDir,
          std::string("transferll_") + std::string(_highways[0]->src->label) +
              std::string("_") + std::string(_highways[0]->dest->label));
      _evaluator.savePerFamilyLikelihoodDiff(out);
    }
    for (auto highway : _highways) {
      (void)(highway);
      _evaluator.removeHighway();
    }
    parameters.setScore(ll);
    return ll;
  }

private:
  const std::vector<Highway *> &_highways;
  AleEvaluator &_evaluator;
};

static Parameters testHighwayFast(AleEvaluator &evaluator,
                                  const Highway &highway,
                                  const std::string &highwaysOutputDir,
                                  double startingProbability = 0.01) {
  std::vector<Highway *> highways;
  auto copy = highway;
  highways.push_back(&copy);
  HighwayFunction f(evaluator, highways);
  Parameters parameters(1);
  parameters[0] = startingProbability;
  f.evaluatePrint(parameters, true, highwaysOutputDir);
  return parameters;
}

static Parameters testHighways(AleEvaluator &evaluator,
                               const std::vector<Highway *> &highways,
                               const Parameters &startingProbabilities,
                               bool optimize, bool thorough) {
  assert(highways.size() == startingProbabilities.dimensions());
  HighwayFunction f(evaluator, highways);
  if (optimize) {
    OptimizationSettings settings;
    settings.strategy = evaluator.getRecModelInfo().recOpt;
    settings.minAlpha = 0.001;
    settings.epsilon = 0.000001;
    settings.factr = LBFGSBPrecision::MEDIUM;
    // settings.verbose = true;
    if (thorough) {
      settings.individualParamOpt = true;
      settings.individualParamOptMinImprovement = 10000.0;
    }
    auto res =
        DTLOptimizer::optimizeParameters(f, startingProbabilities, settings);
    // res.constrain(MIN_PH, MAX_PH);
    return res;
  } else {
    auto parameters = startingProbabilities;
    f.evaluate(parameters);
    return parameters;
  }
}

static Parameters optimizeSingleHighway(AleEvaluator &evaluator,
                                        Highway &highway,
                                        double startingProbability) {
  std::vector<Highway *> highways;
  auto copy = highway;
  highways.push_back(&copy);
  HighwayFunction f(evaluator, highways);
  Parameters startingProbabilities(1);
  startingProbabilities[0] = startingProbability;
  OptimizationSettings settings;
  settings.strategy = RecOpt::LBFGSB;
  settings.minAlpha = 0.001;
  settings.epsilon = 0.000001;
  // settings.verbose = true;
  settings.factr = LBFGSBPrecision::LOW;
  auto res =
      DTLOptimizer::optimizeParameters(f, startingProbabilities, settings);
  // res.constrain(MIN_PH, MAX_PH);
  return res;
}

static bool isHighwayCompatible(Highway &highway, const RecModelInfo &info,
                                const DatedTree &tree) {
  auto from = highway.src;
  auto to = highway.dest;
  switch (info.transferConstraint) {
  case TransferConstaint::NONE:
    return true;
  case TransferConstaint::PARENTS:
    while (from) {
      if (to == from) {
        return false;
      }
      from = from->parent;
    }
    return true;
  case TransferConstaint::RELDATED:
    return tree.canTransferUnderRelDated(from->node_index, to->node_index);
  }
  assert(false);
  return false;
}

void Highways::getCandidateHighways(AleOptimizer &optimizer,
                                    std::vector<ScoredHighway> &scoredHighways,
                                    unsigned int maxCandidates) {
  auto &speciesTree = optimizer.getSpeciesTree();
  unsigned int minTransfers = 2;
  MovesBlackList blacklist;
  std::vector<TransferMove> transferMoves;
  SpeciesTransferSearch::getSortedTransferList(
      speciesTree, optimizer.getEvaluator(), minTransfers, blacklist,
      transferMoves);

  for (const auto &transferMove : transferMoves) {
    auto prune = speciesTree.getNode(transferMove.prune);
    auto regraft = speciesTree.getNode(transferMove.regraft);
    Highway highway(regraft, prune);
    auto lca = speciesTree.getTree().getLCA(prune, regraft);
    unsigned int distance = 0;
    while (prune != lca) {
      distance += 1;
      prune = prune->parent;
    }
    while (regraft != lca) {
      distance += 1;
      regraft = regraft->parent;
    }
    if (distance >= 5) {
      scoredHighways.push_back(ScoredHighway(highway, 0.0));
    } else {
      Logger::timed << "Rejecting (dist) candidate: " << highway.src->label
                    << "->" << highway.dest->label << " d = " << distance
                    << std::endl;
    }
    if (scoredHighways.size() >= maxCandidates) {
      break;
    }
  }
}

std::vector<ScoredHighway>
Highways::getSortedCandidatesFromList(AleOptimizer &optimizer,
                                      std::vector<Highway> &candidateHighways) {
  auto &speciesTree = optimizer.getSpeciesTree();
  unsigned int minTransfers = 1;
  MovesBlackList blacklist;
  std::vector<ScoredHighway> scoredHighways;
  std::vector<TransferMove> transferMoves;
  SpeciesTransferSearch::getSortedTransferList(
      speciesTree, optimizer.getEvaluator(), minTransfers, blacklist,
      transferMoves);

  for (const auto &transferMove : transferMoves) {
    auto prune = speciesTree.getNode(transferMove.prune);
    auto regraft = speciesTree.getNode(transferMove.regraft);
    Highway highway(regraft, prune);
    if (std::find(candidateHighways.begin(), candidateHighways.end(),
                  highway) != candidateHighways.end()) {
      scoredHighways.push_back(ScoredHighway(highway, 0.0));
    }
  }
  for (const auto &highway : candidateHighways) {
    if (std::find(scoredHighways.begin(), scoredHighways.end(), highway) ==
        scoredHighways.end()) {
      scoredHighways.push_back(ScoredHighway(highway, 0.0));
    }
  }
  return scoredHighways;
}

void Highways::filterCandidateHighwaysFast(
    AleOptimizer &optimizer, const std::vector<ScoredHighway> &highways,
    std::vector<ScoredHighway> &filteredHighways, size_t sample_size) {
  auto &evaluator = optimizer.getEvaluator();
  auto &speciesTree = optimizer.getSpeciesTree();
  Logger::timed << "Filering " << highways.size() << " candidate highways"
                << std::endl;
  double initialLL = evaluator.computeLikelihood();
  Logger::timed << "initial ll=" << initialLL << std::endl;
  evaluator.saveSnapshotPerFamilyLL();
  for (const auto &scoredHighway : highways) {
    double proba = 0.01;
    auto highway = scoredHighway.highway;
    if (!isHighwayCompatible(highway, optimizer.getRecModelInfo(),
                             speciesTree.getDatedTree())) {
      Logger::info << "Incompatible highway " << highway.src->label << "->"
                   << highway.dest->label << std::endl;
      continue;
    }
    Logger::timed << "Testing candidate: " << highway.src->label << "->"
                  << highway.dest->label << " with p = " << proba << std::endl;
    auto plausibility_params = testHighwayFast(
        evaluator, highway, optimizer.getHighwaysOutputDir(), proba);
    auto llDiff = plausibility_params.getScore() - initialLL;
    if (llDiff < 0.01) {
      proba = 0.1;
      Logger::timed
          << "No improvement with small probability! Trying again with p = "
          << proba << std::endl;
      plausibility_params = testHighwayFast(
          evaluator, highway, optimizer.getHighwaysOutputDir(), proba);
      llDiff = plausibility_params.getScore() - initialLL;
    }

    if (llDiff > 0.01) {
      auto parameters = optimizeSingleHighway(evaluator, highway, proba);
      llDiff = parameters.getScore() - initialLL;
      if (2 * llDiff > log(sample_size)) {
        evaluator.addHighway(highway);
        initialLL = parameters.getScore();
        Logger::timed << "Accepting candidate: ";
        highway.proba = parameters[0];
        filteredHighways.push_back(ScoredHighway(highway, -llDiff));
        evaluator.saveSnapshotPerFamilyLL();
      } else {
        Logger::timed << "Rejecting (BIC) candidate: ";
      }
    } else {
      Logger::timed << "Rejecting (noImprov) candidate: ";
    }
    Logger::info << highway.src->label << "->" << highway.dest->label
                 << " ll diff = " << llDiff << " best proba = " << highway.proba
                 << std::endl;
  }
  for (auto &highway : filteredHighways) {
    (void)(highway);
    evaluator.removeHighway();
  }
  std::sort(filteredHighways.begin(), filteredHighways.end());
}

void Highways::optimizeAllHighways(
    AleOptimizer &optimizer,
    const std::vector<ScoredHighway> &candidateHighways,
    std::vector<ScoredHighway> &acceptedHighways, bool thorough) {
  auto &evaluator = optimizer.getEvaluator();
  Logger::timed << "Trying to add all candidate highways simultaneously"
                << std::endl;
  std::vector<Highway> highways;
  std::vector<Highway *> highwaysPtr;
  Parameters startingProbabilities;
  for (const auto candidate : candidateHighways) {
    highways.push_back(candidate.highway);
    startingProbabilities.addValue(candidate.highway.proba);
  }
  for (auto &highway : highways) {
    highwaysPtr.push_back(&highway);
  }
  auto parameters = testHighways(evaluator, highwaysPtr, startingProbabilities,
                                 true, thorough);
  Logger::info << parameters << std::endl;
  for (unsigned int i = 0; i < candidateHighways.size(); ++i) {
    ScoredHighway sh(candidateHighways[i]);
    sh.highway.proba = parameters[i];
    acceptedHighways.push_back(sh);
    evaluator.addHighway(sh.highway);
  }
  std::sort(acceptedHighways.rbegin(), acceptedHighways.rend(),
            cmpHighwayByProbability);
}
