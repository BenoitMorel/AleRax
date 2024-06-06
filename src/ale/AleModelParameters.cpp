#include "AleModelParameters.hpp"

void AleModelParameters::getRateVector(RatesVector &rates) const {
  rates.resize(getParamTypeNumber());
  for (unsigned int d = 0; d < getParamTypeNumber(); ++d) {
    rates[d].resize(getSpeciesBranchNumber());
    for (unsigned int s = 0; s < getSpeciesBranchNumber(); ++s) {
      rates[d][s] = getParameter(s, d);
    }
  }
}
