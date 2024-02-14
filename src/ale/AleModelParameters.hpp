#pragma once

#include <util/RecModelInfo.hpp>
#include <util/types.hpp>
#include <maths/ModelParameters.hpp>

#include <vector>
#include <algorithm>
#include <fstream>

/**
 * Stores the model parameters for a gene family
 */
class AleModelParameters {
public:
  /**
   *  @brief Default constructor
   */
  AleModelParameters():
    _speciesBranchNumber(0)
  {
    
  }
 
  /**
   *  @brief Constructor
   *  @param startingRates the set of rates to assign to each category
   *  @param speciesBranchNumber
   *  @param info Description of the model
   */
  AleModelParameters(const Parameters &startingRates,
      unsigned int speciesBranchNumber):
    _paramTypeNumber(startingRates.dimensions()),
    _speciesBranchNumber(speciesBranchNumber),
    _parameters(_speciesBranchNumber, startingRates)
  {
    for (unsigned int i = 0; i < getFreeParameters(); ++i) {
      _parameters[i] = startingRates[i % getParamTypeNumber()];
    }
  }

  /**
   *  Number of free parameters for a single species category
   */
  unsigned int getParamTypeNumber() const {return _paramTypeNumber;}
  

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
  const Parameters &getParameters() {return _parameters;}
  
  /**
   *  Fill rates for a family. Rates is indexed as:
   *  rates[event][species * C + category]
   */
  void getRateVector(RatesVector &rates) const
  {
    rates.resize(getParamTypeNumber());
    for (unsigned int d = 0; d < getParamTypeNumber(); ++d) { 
      rates[d].resize(getSpeciesBranchNumber());
      for (unsigned int s = 0; s < getSpeciesBranchNumber(); ++s) {
        rates[d][s] = getParameter(s, d);
      }
    }
  }

  Parameters getCategorizedParameters(std::vector<unsigned int> &speciesToCategories)
  {
    auto catNumber = *std::max_element(speciesToCategories.begin(), 
        speciesToCategories.end()) + 1; 
    Parameters res(catNumber * getParamTypeNumber());
    for (unsigned int type = 0; type < getParamTypeNumber(); ++type) {
      for (unsigned int species = 0; species < speciesToCategories.size(); ++species) {
        auto cat = speciesToCategories[species];
        res[cat * getParamTypeNumber() + type] = getParameter(species, type);
      }
    }
    return res;
  }

  static Parameters getParametersFromCategorized(const Parameters &categorizedParameters,
      const std::vector<unsigned int> &speciesToCategories,
      unsigned int paramTypeNumber)
  {
    Parameters res(paramTypeNumber * speciesToCategories.size());
    for (unsigned int type = 0; type < paramTypeNumber; ++type) {
      for (unsigned int species = 0; species < speciesToCategories.size(); ++species) {
        auto cat = speciesToCategories[species];
        res[species * paramTypeNumber +  type] = 
          categorizedParameters[cat * paramTypeNumber + type];
      }
    }
    return res;
  }
  
  /**
   *  Get a specific rate for a specific species
   */
  double getParameter(unsigned int species, unsigned int rate) const {
    return _parameters[species * getParamTypeNumber() + rate];
  }
  
  /**
   *  Set a specific rate for a specific species
   */
  void setParameter(unsigned int species, unsigned int rate, double val) {
    _parameters[species * getParamTypeNumber() + rate] = val;
  }
private:
  unsigned int getSpeciesBranchNumber() const {return _speciesBranchNumber;}

  /**
   * Number of tree parameters for all categories
   */
  unsigned int getFreeParameters() const {return getSpeciesBranchNumber() * getParamTypeNumber();} 
  
  
  
  unsigned int _paramTypeNumber;
  unsigned int _speciesBranchNumber;
  Parameters _parameters;
};



