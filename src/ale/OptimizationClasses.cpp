#include "OptimizationClasses.hpp"

#include <sstream>

#include <util/RecModelInfo.hpp>
#include <trees/PLLRootedTree.hpp>
#include <IO/IO.hpp>

// helper function for getSpeciesToCat
static void extractSpeciesToCatRec(corax_rnode_t *node,
    std::vector<unsigned int> &speciesToCat,
    unsigned int _classNumber,
    const std::map<std::string, unsigned int> &labelToCat)
{
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

OptimizationClasses::OptimizationClasses(const PLLRootedTree &speciesTree,
  const std::string &optimizationClassFile,
  const RecModelInfo &info):
  _classNumber(0)
{
  unsigned int N = speciesTree.getNodeNumber();
  std::map<char, std::map<std::string, unsigned int> >labelToCat;
  const std::string rootLabel(speciesTree.getRoot()->label);
  std::string line;
  auto allLabels = speciesTree.getLabels(false);
  _allTypes = info.getParamTypes();
  for (unsigned int i = 0; i < _allTypes.size(); ++i) {
    auto paramType = _allTypes[i];
    _allTypeIndices.insert({paramType, i});
    labelToCat[paramType] = std::map<std::string, unsigned int>();
    labelToCat[paramType].insert({rootLabel, _classNumber});
    _classNumber++;
    _classes[paramType] = std::vector<unsigned int>(N, i);
  }
  std::ifstream is(optimizationClassFile);
  if (is) {
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
      if (allLabels.find(speciesLabel) == allLabels.end()) {
        Logger::error << "Warning, label " << line << " from " << 
          optimizationClassFile << " is not in the species tree" << std::endl;
        ParallelContext::abort(30);
      }
      if (paramTypes.size() == 0) {
        Logger::error << "Error: no parameter type (e.g. D, T, L...) specified for species " << speciesLabel << " in " << optimizationClassFile << std::endl;
        ParallelContext::abort(30);
      }
      for (unsigned int i = 0; i < paramTypes.size(); ++i) {
        auto paramType = paramTypes[i];
        if (std::find(_allTypes.begin(), _allTypes.end(), paramType) == _allTypes.end()) {
          Logger::error << "Error: the parameter type encoded with the letter " << paramType << " does not correspond to any parameter supported by the current reconciliation model. Please edit the file " << optimizationClassFile << std::endl;
          ParallelContext::abort(30);
        }
        labelToCat[paramType].insert({speciesLabel, _classNumber});
        _classNumber++;
      }
    }
  }
  for (auto paramType: _allTypes) {
    extractSpeciesToCatRec(speciesTree.getRoot(),
        _classes[paramType],
        labelToCat[paramType][rootLabel],
        labelToCat[paramType]);
  }
}


Parameters OptimizationClasses::getCompressedParameters(const Parameters &fullParameters) const
{
  Parameters res(_classNumber);
  for (auto type: _allTypes) {
    const auto &typeClasses = _classes.at(type);
    for (unsigned int species = 0; species < typeClasses.size(); ++species) {
      auto cat = typeClasses[species];
      auto typeIndex = _allTypeIndices.at(type);
      res[cat] = fullParameters[species * _allTypes.size() + typeIndex]; 
    }
  }
  return res;
}

Parameters OptimizationClasses::getFullParameters(const Parameters &compressedParameters) const
{
  Parameters res(_allTypes.size() *_classes.at(_allTypes[0]).size());
  for (auto type: _allTypes) {
    const auto &typeClasses = _classes.at(type);
    for (unsigned int species = 0; species < typeClasses.size(); ++species) {
      auto cat = typeClasses[species];
      auto typeIndex = _allTypeIndices.at(type);
      res[species * _allTypes.size() + typeIndex] = compressedParameters[cat];
    }
  }
  return res;
}

