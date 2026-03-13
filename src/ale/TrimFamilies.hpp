#pragma once

#include <IO/Families.hpp>

/**
 *  A collection of routines to trim "bad" families (that are too large,
 *  that have too much uncertainty, or that cover not enough species)
 */
class TrimFamilies {
public:
  TrimFamilies() = delete;

  /**
   *  Remove the families that cover less than minCoverage species
   *
   *  @param families Input and output families
   *  @param minCoverage The minimum number of species that must be covered
   */
  static void trimMinSpeciesCoverage(Families &families,
                                     unsigned int minCoverage);

  /**
   *  Sort the families from the largest to the smallest (in terms of
   *  CCP size) and keep the smallest ones
   *
   *  @param families Input and output families
   *  @param keepRatio The proportion of families to keep
   */
  static void trimHighCladesNumber(Families &families, double keepRatio);

  /**
   *  @brief Trim the families with high gene tree uncertainty
   *
   *  For each family, compute the ratio between the CCP size and
   *  the gene tree size, and remove the families with a ratio greater
   *  than maxRatio
   *
   *  @param families Input and output families
   *  @param maxRatio Threshold ratio
   *
   */
  static void trimCladeSplitRatio(Families &families, double maxRatio);
};
