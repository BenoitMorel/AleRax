#include "OptimizationClasses.hpp"

#include <sstream>

#include <IO/IO.hpp>
#include <trees/PLLRootedTree.hpp>
#include <util/RecModelInfo.hpp>

// helper function for getSpeciesToCat
static void
extractSpeciesToCatRec(corax_rnode_t *node,
                       std::vector<unsigned int> &speciesToCat,
                       unsigned int _classNumber,
                       const std::map<std::string, unsigned int> &labelToCat) {
  std::string label = node->label;
  auto it = labelToCat.find(label);
  if (it != labelToCat.end()) {
    _classNumber = it->second;
  }
  speciesToCat[node->node_index] = _classNumber;
  if (node->left) {
    extractSpeciesToCatRec(node->left, speciesToCat, _classNumber, labelToCat);
    extractSpeciesToCatRec(node->right, speciesToCat, _classNumber, labelToCat);
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
      // already added by default
      continue;
    }
    if (allSpeciesLabels.find(speciesLabel) == allSpeciesLabels.end()) {
      Logger::error << "Error, label " << line << " from "
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
        Logger::error
            << "Error: the parameter type encoded with the letter " << paramType
            << " does not correspond to any parameter supported by the current "
               "reconciliation model. Please edit the file "
            << optimizationClassFile << std::endl;
        ParallelContext::abort(30);
      }
      labelToCat[paramType].insert({speciesLabel, classNumber});
      classNumber++;
    }
  }
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
    // origination probabilities do not belong to the types to optimize, return
    return;
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
    const PLLRootedTree &speciesTree, ModelParametrization parametrization,
    const std::string &optimizationClassFile, const RecModelInfo &info)
    : _classNumber(0) {
  unsigned int N = speciesTree.getNodeNumber();
  std::map<char, std::map<std::string, unsigned int>> labelToCat;
  const std::string rootLabel(speciesTree.getRoot()->label);
  auto allSpeciesLabels = speciesTree.getLabels(false);
  _allTypes = info.getParamTypes();
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
  };

  for (auto paramType : _allTypes) {
    extractSpeciesToCatRec(speciesTree.getRoot(), _classes[paramType],
                           labelToCat[paramType][rootLabel],
                           labelToCat[paramType]);
  }
}

static bool normalizeParamType(char paramType) { return paramType == 'O'; }

Parameters OptimizationClasses::getCompressedParameters(
    const Parameters &fullParameters) const {
  Parameters res(_classNumber);
  for (auto type : _allTypes) {
    const auto &typeClasses = _classes.at(type);
    double norm = 1.0;
    if (normalizeParamType(type)) {
      norm = 0;
      for (unsigned int species = 0; species < typeClasses.size(); ++species) {
        auto typeIndex = _allTypeIndices.at(type);
        norm += fullParameters[species * _allTypes.size() + typeIndex];
      }
    }
    for (unsigned int species = 0; species < typeClasses.size(); ++species) {
      auto cat = typeClasses[species];
      auto typeIndex = _allTypeIndices.at(type);
      res[cat] = fullParameters[species * _allTypes.size() + typeIndex] / norm;
    }
  }
  return res;
}

Parameters OptimizationClasses::getFullParameters(
    const Parameters &compressedParameters) const {
  Parameters res(_allTypes.size() * _classes.at(_allTypes[0]).size());
  for (auto type : _allTypes) {
    const auto &typeClasses = _classes.at(type);
    double norm = 1.0;
    if (normalizeParamType(type)) {
      norm = 0.0;
      for (unsigned int species = 0; species < typeClasses.size(); ++species) {
        auto cat = typeClasses[species];
        auto typeIndex = _allTypeIndices.at(type);
        norm += res[species * _allTypes.size() + typeIndex] =
            compressedParameters[cat];
      }
    }
    for (unsigned int species = 0; species < typeClasses.size(); ++species) {
      auto cat = typeClasses[species];
      auto typeIndex = _allTypeIndices.at(type);
      res[species * _allTypes.size() + typeIndex] =
          compressedParameters[cat] / norm;
    }
  }
  return res;
}
