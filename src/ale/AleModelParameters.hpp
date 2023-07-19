#pragma once

#include <util/RecModelInfo.hpp>

#include <vector>
#include <algorithm>
#include <fstream>

/**
 * Stores the model parameters: DTL rates, global or per species category
 * A category is a set of species that shares the same DTL rates
 */
class AleModelParameters {
public:
  /**
   *  @brief Default constructor
   */
  AleModelParameters():
    _catNumber(0)
  {
    
  }
 
  /**
   *  @brief Constructor
   *  @param startingRates the set of rates to assign to each category
   *  @param speciesToCategory mapping species node index to category
   *  @param cat index to cat label
   *  @param info Description of the model
   */
  AleModelParameters(const Parameters &startingRates,
      const std::vector<unsigned int> &speciesToCategory,
      const std::vector<std::string> &catToLabel,
      unsigned int families,
      const RecModelInfo &info):
    _info(info),
    _speciesToCat(speciesToCategory),
    _catToLabel(catToLabel),
    _catNumber(*std::max_element(std::begin(speciesToCategory), std::end(speciesToCategory)) + 1),
    _familyCatNumber(info.perFamilyRates ? families : 1),
    _parameters(_catNumber * _familyCatNumber, startingRates)
  {
    for (unsigned int i = 0; i < getTotalFreeParameters(); ++i) {
      _parameters[i] = startingRates[i % perCategoryFreeParameters()];
    }
  }
  
  /**
   *  Number of free parameters for a single species category
   */
  unsigned int perCategoryFreeParameters() const {return _info.modelFreeParameters();}
  
  

  /**
   *  Get model description
   */
  const RecModelInfo &getInfo() const {return _info;}

  /**
   *  Set rates from a parameter vector
   */
  void setParameters(const Parameters &parameters) {
    assert(parameters.dimensions() == getTotalFreeParameters());
    _parameters = parameters;
  }

  void setParametersForFamily(unsigned int family, const Parameters &parameters) {
    assert(parameters.dimensions() == getPerFamilyFreeParameters());
    for (unsigned int i = 0; i < getPerFamilyFreeParameters(); ++i) {
      _parameters[i + family * getPerFamilyFreeParameters()] = parameters[i];
    }
  }

  /**
   *  Fill rates for a family. Rates is indexed as:
   *  rates[event][species * C + category]
   */
  void getRatesForFamily(unsigned int family, RatesVector &rates) const
  {
    if (_familyCatNumber == 1) {
      family = 0;
    }
    unsigned int N = getSpecesNumber();
    unsigned int F = perCategoryFreeParameters();
    rates.resize(F);
    for (unsigned int d = 0; d < F; ++d) { 
      rates[d].resize(N);
      for (unsigned int s = 0; s < N; ++s) {
        rates[d][s] = getRateFromSpecies(family, s, d);
      }
    }
  }

  /**
   *  Get rates as a parameter vector
   */
  const Parameters &getParameters() {return _parameters;}
  
  Parameters getParametersForFamily(unsigned int family) {
    Parameters parameters(getPerFamilyFreeParameters());
    for (unsigned int i = 0; i < parameters.dimensions(); ++i) {
      parameters[i] = _parameters[i + family * getPerFamilyFreeParameters()];   
    }
    return parameters;
  }

  friend std::ostream &operator<<(std::ostream &os, const AleModelParameters &mp)  {
    os << "[";
    for (unsigned int c = 0; c < mp.categoryNumber(); ++c) {
      os << mp._catToLabel[c] << " ";
      os << "(";
      for (unsigned int r = 0; r < mp.perCategoryFreeParameters(); ++r) {
        os << mp.getRateFromCat(0, c, r);
        if (r !=  mp.perCategoryFreeParameters() -1) {
          os << ",";
        }
      }
      os << ")";
      if (c != mp.categoryNumber() - 1) {
        os << ",";
      }
    }
    os << "]";
    return os;
  }

private:
  /**
   *  Number of species categories
   */
  unsigned int categoryNumber() const {return _catNumber;}
  unsigned int getSpecesNumber() const {return _speciesToCat.size();}

  /**
   * Number of tree parameters for all categories
   */
  unsigned int getTotalFreeParameters() const {return getPerFamilyFreeParameters() * _familyCatNumber;} 
  unsigned int getPerFamilyFreeParameters() const {return categoryNumber() * perCategoryFreeParameters();} 
  
  /**
   *  Get a specific rate for a specific species
   */
  double getRateFromSpecies(unsigned int family, unsigned int species, unsigned int rate) const {
    return getRateFromCat(family, _speciesToCat[species], rate);
  }
  
  double getRateFromCat(unsigned int family, unsigned int cat, unsigned int rate) const {
    return _parameters[family * getPerFamilyFreeParameters() + cat * perCategoryFreeParameters() + rate];
  }
  RecModelInfo _info;
  std::vector<unsigned int> _speciesToCat;
  std::vector<std::string> _catToLabel;
  unsigned int _catNumber;
  unsigned int _familyCatNumber;
  Parameters _parameters;
};



