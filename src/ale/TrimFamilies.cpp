#include "TrimFamilies.hpp"

#include <algorithm>
#include <vector>

#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <ccp/ConditionalClades.hpp>
#include <parallelization/ParallelContext.hpp>
#include <util/types.hpp>

void TrimFamilies::trimMinSpeciesCoverage(Families &families,
                                          unsigned int minCoverage) {
  auto N = families.size();
  // filter families in parallel
  std::vector<unsigned int> localToKeep;
  for (auto i = ParallelContext::getBegin(N); i < ParallelContext::getEnd(N);
       ++i) {
    const auto &family = families[i];
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    if (mapping.getCoveredSpecies().size() >= minCoverage) {
      localToKeep.push_back(i);
    }
  }
  ParallelContext::barrier();
  // gather family information from all MPI ranks
  std::vector<unsigned int> toKeep;
  ParallelContext::concatenateHetherogeneousUIntVectors(localToKeep, toKeep);
  // trim families
  Families familiesCopy = families;
  families.clear();
  for (auto i : toKeep) {
    assert(ParallelContext::isIntEqual(i));
    families.push_back(familiesCopy[i]);
  }
}

void TrimFamilies::trimHighCladesNumber(Families &families, double keepRatio) {
  assert(keepRatio >= 0.0);
  auto N = families.size();
  // get family observed clade numbers in parallel
  std::vector<unsigned int> localCcpSizes;
  std::vector<unsigned int> localCheck;
  for (auto i = ParallelContext::getBegin(N); i < ParallelContext::getEnd(N);
       ++i) {
    ConditionalClades ccp;
    ccp.unserialize(families[i].ccp);
    localCcpSizes.push_back(ccp.getCladesNumber());
    localCheck.push_back(i);
  }
  ParallelContext::barrier();
  // gather family information from all MPI ranks
  std::vector<unsigned int> ccpSizes;
  std::vector<unsigned int> check;
  ParallelContext::concatenateHetherogeneousUIntVectors(localCcpSizes,
                                                        ccpSizes);
  ParallelContext::concatenateHetherogeneousUIntVectors(localCheck, check);
  // sort families by the number of observed clades
  std::vector<PairUInt> sizeToIndex(N);
  for (unsigned int i = 0; i < N; ++i) {
    sizeToIndex[i].first = ccpSizes[i];
    sizeToIndex[i].second = i;
  }
  std::sort(sizeToIndex.begin(), sizeToIndex.end());
  // trim families
  Families familiesCopy = families;
  families.clear();
  unsigned int cutFrom = (unsigned int)(double(N) * keepRatio);
  assert(cutFrom <= N);
  unsigned int maxCcpSizeKept = 0;
  for (unsigned int i = 0; i < cutFrom; ++i) {
    assert(ParallelContext::isIntEqual(sizeToIndex[i].second));
    families.push_back(familiesCopy[sizeToIndex[i].second]);
    maxCcpSizeKept = sizeToIndex[i].first;
  }
  Logger::info << "Trimming families with more than " << maxCcpSizeKept
               << " clades" << std::endl;
}

void TrimFamilies::trimCladeSplitRatio(Families &families, double maxRatio) {
  assert(maxRatio >= 0.0);
  auto N = families.size();
  // filter families in parallel
  std::vector<unsigned int> localToKeep;
  unsigned int allClades = 0;
  unsigned int keptClades = 0;
  for (auto i = ParallelContext::getBegin(N); i < ParallelContext::getEnd(N);
       ++i) {
    ConditionalClades ccp;
    ccp.unserialize(families[i].ccp);
    auto c = ccp.getCladesNumber();
    auto n = ccp.getLeafNumber() * 2 - 3;
    allClades += c;
    if (maxRatio * double(n) >= double(c)) {
      localToKeep.push_back(i);
      keptClades += c;
    }
  }
  ParallelContext::barrier();
  // gather family information from all MPI ranks
  std::vector<unsigned int> toKeep;
  ParallelContext::concatenateHetherogeneousUIntVectors(localToKeep, toKeep);
  ParallelContext::sumUInt(allClades);
  ParallelContext::sumUInt(keptClades);
  Logger::info << "Clades before trimming: " << allClades << std::endl;
  Logger::info << "Clades after trimming: " << keptClades << std::endl;
  // trim families
  Families familiesCopy = families;
  families.clear();
  for (auto i : toKeep) {
    assert(ParallelContext::isIntEqual(i));
    families.push_back(familiesCopy[i]);
  }
}
