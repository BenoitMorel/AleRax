#pragma once

#include <unordered_map>
#include <string>
#include <vector>
#include <maths/Parameters.hpp>

class PLLRootedTree;
struct RecModelInfo;

class OptimizationClasses {
public:
  /**
   *  Constructor
   *  @param speciesTree 
   *  @param optimizationClassFile
   *  @param info 
   */
  OptimizationClasses(const PLLRootedTree &speciesTree,
    const std::string &optimizationClassFile,
    const RecModelInfo &info);

  /**
   *  Accessor
   *  @param paramType e.g. 'D', 'T', 'L' etc.
   *  @param speciesNodeIndex The index of the corresponding species branch
   */
  unsigned int getClass(char paramType, unsigned int speciesNodeIndex) const {
    return _classes.at(paramType)[speciesNodeIndex];
  }

  Parameters getCompressedParameters(const Parameters &fullParameters) const;
  Parameters getFullParameters(const Parameters &compressedParameters) const;

private:
  std::unordered_map<char, std::vector<unsigned int> > _classes;
  std::vector<char> _allTypes; // e.g. [D, T, L]
  std::unordered_map<char, unsigned int> _allTypeIndices; // indices of the elements of _allTypes
  unsigned int _classNumber;
};

