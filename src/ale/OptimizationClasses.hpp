#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include <maths/Parameters.hpp>
#include <util/enums.hpp>

class PLLRootedTree;
struct RecModelInfo;

/**
 *  Describes which DTLO parameters should be shared among which species
 * branches. An optimization class is a set of pairs of (paramType,
 * speciesNodeIndex) that share the same value when optimizing the model
 * parameters. Each optimization class adds one free parameter in the whole
 * model.
 */
class OptimizationClasses {
public:
  /**
   *  Constructor
   *  @param speciesTree
   *  @param parametrization
   *  @param optimizationClassFile
   *  @param info
   */
  OptimizationClasses(const PLLRootedTree &speciesTree,
                      const ModelParametrization &parametrization,
                      const std::string &optimizationClassFile,
                      const RecModelInfo &info);

  /**
   *  The two following methods convert the compressed parameters to the full
   * parameters and vice versa. Applying both methods results in the same
   * initial parameter vector.
   *
   *  -Full parameters: raw parameter vector, with one value for each species
   * and each parameter type, the order of the elements is the same as in the
   * AleModelParameters class. -Compressed parameters: one value per
   * optimization class
   */
  Parameters getCompressedParameters(const Parameters &fullParameters) const;
  Parameters getFullParameters(const Parameters &compressedParameters) const;

  /**
   *  Accessor
   *  @param paramType e.g. 'D', 'T', 'L' etc.
   *  @param speciesNodeIndex The index of the corresponding species branch
   */
  unsigned int getClass(char paramType, unsigned int speciesNodeIndex) const {
    return _classes.at(paramType)[speciesNodeIndex];
  }

  /**
   *  Accessor
   */
  unsigned int getFreeParameters() const { return _classNumber; }

private:
  // _classes[paramType][speciesNodeIndex] = integer representing the
  // corresponding class
  std::unordered_map<char, std::vector<unsigned int>> _classes;
  // list of the represented parameter types (e.g. [D, T, L])
  std::vector<char> _allTypes;
  // indices of the elements of _allTypes
  std::unordered_map<char, unsigned int> _allTypeIndices;
  // number of optimization classes (number of free parameters in the model)
  unsigned int _classNumber;
};
