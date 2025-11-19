#pragma once

#include <algorithm>
#include <vector>

#include <maths/Parameters.hpp>
#include <util/types.hpp>

/**
 * Stores the model parameters for a gene family
 */
class AleModelParameters {
public:
  /**
   *  @brief Default constructor
   */
  AleModelParameters() : _speciesBranchNumber(0) {}

  /**
   *  @brief Constructor
   *  @param startingRates Set of rates to assign to each category
   *  @param speciesBranchNumber Number of species branches in the species tree
   */
  AleModelParameters(const Parameters &startingRates,
                     unsigned int speciesBranchNumber)
      : _paramTypeNumber(startingRates.dimensions()),
        _speciesBranchNumber(speciesBranchNumber),
        _parameters(_paramTypeNumber * _speciesBranchNumber) {
    for (unsigned int i = 0; i < getFreeParameters(); ++i) {
      _parameters[i] = startingRates[i % getParamTypeNumber()];
    }
  }

  /**
   *  @brief Constructor
   *  @param paramTypeNumber Number of different rates per species branch
   *  @param speciesBranchNumber Number of species branches in the species tree
   *  @param startingValue Starting value for all rates
   */
  AleModelParameters(unsigned int paramTypeNumber,
                     unsigned int speciesBranchNumber, double startingValue)
      : _paramTypeNumber(paramTypeNumber),
        _speciesBranchNumber(speciesBranchNumber),
        _parameters(_paramTypeNumber * _speciesBranchNumber) {
    for (unsigned int i = 0; i < getFreeParameters(); ++i) {
      _parameters[i] = startingValue;
    }
  }

  /**
   *  Number of free parameters for a single species category
   */
  unsigned int getParamTypeNumber() const { return _paramTypeNumber; }

  /**
   *  Number of species branches
   */
  unsigned int getSpeciesBranchNumber() const { return _speciesBranchNumber; }

  /**
   *  Set rates from a parameter vector
   */
  void setParameters(const Parameters &parameters) {
    assert(parameters.dimensions() == getFreeParameters());
    _parameters = parameters;
  }

  /**
   *  Get rates as a parameter vector
   */
  const Parameters &getParameters() const { return _parameters; }
  Parameters &getParameters() { return _parameters; }

  /**
   *  Set a specific rate for a specific species
   */
  void setParameter(unsigned int species, unsigned int rate, double val) {
    _parameters[species * getParamTypeNumber() + rate] = val;
  }

  /**
   *  Get a specific rate for a specific species
   */
  double getParameter(unsigned int species, unsigned int rate) const {
    return _parameters[species * getParamTypeNumber() + rate];
  }
  double &getParameter(unsigned int species, unsigned int rate) {
    return _parameters[species * getParamTypeNumber() + rate];
  }

  /**
   *  Fill rates for a family. RatesVector is indexed as:
   *  rates[paramType][speciesNodeIndex]
   */
  void getRateVector(RatesVector &rates) const {
    rates.resize(getParamTypeNumber());
    for (unsigned int t = 0; t < getParamTypeNumber(); ++t) {
      rates[t].resize(getSpeciesBranchNumber());
      for (unsigned int s = 0; s < getSpeciesBranchNumber(); ++s) {
        rates[t][s] = getParameter(s, t);
      }
    }
  }

private:
  /**
   *  Number of tree parameters for all categories
   */
  unsigned int getFreeParameters() const {
    return getSpeciesBranchNumber() * getParamTypeNumber();
  }

private:
  unsigned int _paramTypeNumber;
  unsigned int _speciesBranchNumber;
  Parameters _parameters;
};
