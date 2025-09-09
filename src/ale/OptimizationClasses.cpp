#include "OptimizationClasses.hpp"

#include <fstream>
#include <map>
#include <sstream>
#include <unordered_set>

#include <parallelization/ParallelContext.hpp>
#include <trees/PLLRootedTree.hpp>
#include <util/RecModelInfo.hpp>

static void extractSpeciesToCatRecursive(
    corax_rnode_t *node, std::vector<unsigned int> &speciesToCat,
    unsigned int cat, const std::map<std::string, unsigned int> &labelToCat) {
  std::string label = node->label;
  auto it = labelToCat.find(label);
  if (it != labelToCat.end()) {
    cat = it->second;
  }
  speciesToCat[node->node_index] = cat;
  if (node->left) {
    extractSpeciesToCatRecursive(node->left, speciesToCat, cat, labelToCat);
    extractSpeciesToCatRecursive(node->right, speciesToCat, cat, labelToCat);
  }
}

static void fillFromOptimizationClassFile(
    const std::string &optimizationClassFile,
    const std::unordered_set<std::string> &allSpeciesLabels,
    const std::vector<char> &allTypes, const std::string &rootLabel,
    std::map<char, std::map<std::string, unsigned int>> &labelToCat,
    unsigned int &classNumber) {
  std::ifstream is(optimizationClassFile);
  if (!is) {
    Logger::error << "Error: cannot open parametrization file "
                  << optimizationClassFile << std::endl;
    ParallelContext::abort(30);
  }
  std::string line;
  while (std::getline(is, line)) {
    if (line.size() == 0) {
      continue;
    }
    std::stringstream iss(line);
    std::string speciesLabel;
    std::string paramTypes;
    iss >> speciesLabel;
    iss >> paramTypes;
    if (speciesLabel == rootLabel) {
      // the root has already been added
      continue;
    }
    if (allSpeciesLabels.find(speciesLabel) == allSpeciesLabels.end()) {
      Logger::error << "Error: label " << speciesLabel << " in "
                    << optimizationClassFile << " is not in the species tree"
                    << std::endl;
      ParallelContext::abort(30);
    }
    if (paramTypes.size() == 0) {
      Logger::error
          << "Error: no parameter type (e.g. D, T, L...) specified for species "
          << speciesLabel << " in " << optimizationClassFile << std::endl;
      ParallelContext::abort(30);
    }
    for (unsigned int i = 0; i < paramTypes.size(); ++i) {
      auto paramType = paramTypes[i];
      if (std::find(allTypes.begin(), allTypes.end(), paramType) ==
          allTypes.end()) {
        Logger::error << "Error: the parameter type encoded with the letter "
                      << paramType << " in " << optimizationClassFile
                      << " does not correspond to any parameter "
                      << "supported by the current reconciliation model"
                      << std::endl;
        ParallelContext::abort(30);
      }
      labelToCat[paramType].insert({speciesLabel, classNumber});
      classNumber++;
    }
  }
  is.close();
}

static void fillFromAllSpecies(
    const std::unordered_set<std::string> &allSpeciesLabels,
    const std::vector<char> &allTypes, const std::string &rootLabel,
    std::map<char, std::map<std::string, unsigned int>> &labelToCat,
    unsigned int &classNumber) {
  for (auto &speciesLabel : allSpeciesLabels) {
    if (speciesLabel == rootLabel) {
      // the root has already been added
      continue;
    }
    for (auto paramType : allTypes) {
      labelToCat[paramType].insert({speciesLabel, classNumber});
      classNumber++;
    }
  }
}

static void fillFromOriginations(
    const std::unordered_set<std::string> &allSpeciesLabels,
    const std::vector<char> &allTypes, const std::string &rootLabel,
    std::map<char, std::map<std::string, unsigned int>> &labelToCat,
    unsigned int &classNumber) {
  const auto ORIGINATION_TYPE = 'O';
  if (std::find(allTypes.begin(), allTypes.end(), ORIGINATION_TYPE) ==
      allTypes.end()) {
    // origination probabilities do not belong to the types to optimize,
    // we shouldn't get here
    assert(false);
  }
  for (auto &speciesLabel : allSpeciesLabels) {
    if (speciesLabel == rootLabel) {
      // the root has already been added
      continue;
    }
    labelToCat[ORIGINATION_TYPE].insert({speciesLabel, classNumber});
    classNumber++;
  }
}

OptimizationClasses::OptimizationClasses(
    const PLLRootedTree &speciesTree,
    const ModelParametrization &parametrization,
    const std::string &optimizationClassFile, const RecModelInfo &info)
    : _allTypes(info.getParamTypes()), _classNumber(0) {
  unsigned int N = speciesTree.getNodeNumber();
  const std::string rootLabel(speciesTree.getRoot()->label);
  const auto allSpeciesLabels = speciesTree.getLabels(false);
  std::map<char, std::map<std::string, unsigned int>> labelToCat;
  for (unsigned int i = 0; i < _allTypes.size(); ++i) {
    auto paramType = _allTypes[i];
    _allTypeIndices.insert({paramType, i});
    labelToCat[paramType] = std::map<std::string, unsigned int>();
    labelToCat[paramType].insert({rootLabel, _classNumber});
    _classNumber++;
    _classes[paramType] = std::vector<unsigned int>(N, i);
  }
  switch (parametrization) {
  case ModelParametrization::GLOBAL:
  case ModelParametrization::PER_FAMILY:
    // nothing to be done
    break;
  case ModelParametrization::PER_SPECIES:
    fillFromAllSpecies(allSpeciesLabels, _allTypes, rootLabel, labelToCat,
                       _classNumber);
    break;
  case ModelParametrization::ORIGINATION_PER_SPECIES:
    fillFromOriginations(allSpeciesLabels, _allTypes, rootLabel, labelToCat,
                         _classNumber);
    break;
  case ModelParametrization::CUSTOM:
    fillFromOptimizationClassFile(optimizationClassFile, allSpeciesLabels,
                                  _allTypes, rootLabel, labelToCat,
                                  _classNumber);
    break;
  default:
    assert(false);
  }
  for (auto paramType : _allTypes) {
    extractSpeciesToCatRecursive(speciesTree.getRoot(), _classes[paramType],
                                 labelToCat[paramType][rootLabel],
                                 labelToCat[paramType]);
  }
}

static bool normalizeParamType(const char paramType) {
  return paramType == 'O';
}

Parameters OptimizationClasses::getCompressedParameters(
    const Parameters &fullParameters) const {
  Parameters res(_classNumber);
  for (auto type : _allTypes) {
    const auto typeIndex = _allTypeIndices.at(type);
    const auto &typeClasses = _classes.at(type);
    double norm = 1.0;
    if (normalizeParamType(type)) {
      norm = 0.0;
      for (unsigned int species = 0; species < typeClasses.size(); ++species) {
        norm += fullParameters[species * _allTypes.size() + typeIndex];
      }
    }
    for (unsigned int species = 0; species < typeClasses.size(); ++species) {
      auto cat = typeClasses[species];
      res[cat] = fullParameters[species * _allTypes.size() + typeIndex] / norm;
    }
  }
  return res;
}

Parameters OptimizationClasses::getFullParameters(
    const Parameters &compressedParameters) const {
  Parameters res(_allTypes.size() * _classes.at(_allTypes[0]).size());
  for (auto type : _allTypes) {
    const auto typeIndex = _allTypeIndices.at(type);
    const auto &typeClasses = _classes.at(type);
    double norm = 1.0;
    if (normalizeParamType(type)) {
      norm = 0.0;
      for (unsigned int species = 0; species < typeClasses.size(); ++species) {
        auto cat = typeClasses[species];
        norm += compressedParameters[cat];
      }
    }
    for (unsigned int species = 0; species < typeClasses.size(); ++species) {
      auto cat = typeClasses[species];
      res[species * _allTypes.size() + typeIndex] =
          compressedParameters[cat] / norm;
    }
  }
  return res;
}
