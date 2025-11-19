#pragma once

#include <IO/HighwayCandidateParser.hpp>
#include <ccp/ConditionalClades.hpp>
#include <likelihoods/reconciliation_models/BaseReconciliationModel.hpp>

/**
 *  Temp storage for a clade event sampled in
 *  the reconciliation space
 */
template <class REAL> struct ReconciliationCell {
  Scenario::Event event;
  REAL maxProba;
  double blLeft;
  double blRight;
};

/**
 *  Base class for reconciliation models that take as input
 *  gene tree distributions (represented by conditional clade
 *  probabilities)
 */
class MultiModelInterface : public BaseReconciliationModel {
public:
  MultiModelInterface(PLLRootedTree &speciesTree,
                      const GeneSpeciesMapping &geneSpeciesMapping,
                      const RecModelInfo &info, const std::string &ccpFile)
      : BaseReconciliationModel(speciesTree, geneSpeciesMapping, info),
        _memorySavings(info.memorySavings) {
    _ccp.unserialize(ccpFile);
    mapGenesToSpecies();
  }

  virtual ~MultiModelInterface() {}

  const ConditionalClades &getCCP() const { return _ccp; }

  virtual void setAlpha(double alpha) = 0;
  virtual void setRates(const RatesVector &rates) = 0;
  virtual void setHighways(const std::vector<Highway> &highways) {
    (void)(highways);
  }

  virtual double computeLogLikelihood() = 0;

  virtual bool inferMLScenario(Scenario &scenario) = 0;

  virtual bool
  sampleReconciliations(unsigned int samples,
                        std::vector<std::shared_ptr<Scenario>> &scenarios) = 0;

protected:
  // should be called on changing the model rates (alpha/DTLO/
  // highways), as the LLs computed earlier become irrelevant
  void resetCache() { _llCache.clear(); }

private:
  virtual void mapGenesToSpecies();

protected:
  const bool _memorySavings;
  ConditionalClades _ccp;
  std::unordered_map<size_t, double> _llCache;
};

inline void MultiModelInterface::mapGenesToSpecies() {
  // fill _geneToSpecies, _speciesCoverage and _numberOfCoveredSpecies
  this->_geneToSpecies.clear();
  this->_speciesCoverage.resize(this->_speciesTree.getLeafNumber(), 0);
  this->_numberOfCoveredSpecies = 0;
  const auto &cidToLeaves = _ccp.getCidToLeaves();
  for (const auto &p : cidToLeaves) {
    const auto cid = p.first;
    const auto &geneName = p.second;
    const auto &speciesName = this->_geneNameToSpeciesName[geneName];
    const auto spid = this->_speciesNameToId[speciesName];
    this->_geneToSpecies[cid] = spid;
    if (!this->_speciesCoverage[spid]) {
      this->_numberOfCoveredSpecies++;
    }
    this->_speciesCoverage[spid]++;
  }
  // build the pruned species tree representation based
  // on the species coverage
  onSpeciesTreeChange(nullptr);
}

/**
 *  Implements all the methods that require the REAL template
 *  and are not dependent on modelling HGT
 */
template <class REAL> class MultiModel : public MultiModelInterface {
public:
  MultiModel(PLLRootedTree &speciesTree,
             const GeneSpeciesMapping &geneSpeciesMapping,
             const RecModelInfo &info, const std::string &ccpFile)
      : MultiModelInterface(speciesTree, geneSpeciesMapping, info, ccpFile) {}

  virtual ~MultiModel() {}

  /**
   *  Compute the LL of the current species tree
   */
  virtual double computeLogLikelihood();

  /**
   *  Find the ML reconciled gene tree and store it in
   *  the scenario variable
   */
  virtual bool inferMLScenario(Scenario &scenario);

  /**
   *  Sample reconciled gene trees and add them to
   *  the scenarios vector
   */
  virtual bool
  sampleReconciliations(unsigned int samples,
                        std::vector<std::shared_ptr<Scenario>> &scenarios);

protected:
  // functions to work with CLVs
  virtual void allocateMemory() = 0;
  virtual void deallocateMemory() = 0;
  virtual void updateCLV(CID cid) = 0;

  /**
   *  Probability of an arbitrary gene family to originate
   *  somewhere in the species tree and survive (at least as
   *  a single copy on a single leaf).
   *  All possible reconciliations of the current family are
   *  included by this scenario, thus it is the conditioning
   *  probability
   */
  virtual double getLikelihoodFactor(unsigned int category) = 0;

  /**
   *  Probability of the current family to evolve starting
   *  from a given species branch.
   *  Sum over all scenarios stemming from the species branch
   */
  virtual REAL getRootCladeLikelihood(corax_rnode_t *speciesNode,
                                      unsigned int category) = 0;

  /**
   *  Compute the CLV value for a given cid (clade id) and a given
   *  species node. Return true in case of success.
   *  If recCell is set, sample the next event in the reconciliation
   *  space
   */
  virtual bool
  computeProbability(CID cid, corax_rnode_t *speciesNode, unsigned int category,
                     REAL &proba,
                     ReconciliationCell<REAL> *recCell = nullptr) = 0;

  /**
   *  Return a uniform random REAL from the [0,max) interval
   */
  REAL getRandom(REAL max) {
    auto rand = max * Random::getProba();
    scale(rand);
    if (rand == max) { // may occur due to rounding issues
      rand = REAL();
    }
    return rand;
  }

private:
  /**
   *  Return an undated or dated species tree hash
   */
  virtual size_t getHash() = 0;

  /**
   *  Number of mixture categories for the per-branch
   *  probabilities of elementary events
   */
  unsigned int getGammaCatNumber() { return this->_info.gammaCategories; }

  /**
   *  Fill CLV for each observed gene clade
   */
  void updateCLVs();

  /**
   *  Jointly sample a category and a species branch based on
   *  the probability of the gene tree to be rooted in it
   */
  corax_rnode_t *sampleOriginationSpecies(unsigned int &category);

  /**
   *  Wrapper for the backtrace function
   */
  bool computeScenario(Scenario &scenario, bool stochastic);

  /**
   *  Recursively sample a reconciled gene tree
   */
  bool backtrace(CID cid, corax_rnode_t *speciesNode, corax_unode_t *geneNode,
                 unsigned int category, Scenario &scenario, bool stochastic);
};

template <class REAL> double MultiModel<REAL>::computeLogLikelihood() {
  if (this->_ccp.skip()) {
    return 0.0;
  }
  auto hash = getHash();
  // check if LL is already computed for this species tree
  auto cacheIt = this->_llCache.find(hash);
  if (cacheIt != this->_llCache.end()) {
    // too many collisions for dated getHash(),
    // delete the condition line upon fixing
    if (!this->_info.isDated())
      return cacheIt->second;
  }
  if (this->_memorySavings) {
    // initialize CLVs, as they have not been created with the model
    allocateMemory();
  }
  // compute CLVs based on recomputed per-branch probas
  updateCLVs();
  // compute the species tree LL
  REAL res = REAL();
  std::vector<REAL> categoryLikelihoods(getGammaCatNumber(), REAL());
  for (unsigned int c = 0; c < getGammaCatNumber(); ++c) {
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      categoryLikelihoods[c] += getRootCladeLikelihood(speciesNode, c);
    }
    // condition on survival
    categoryLikelihoods[c] /=
        getLikelihoodFactor(c); // no need to scale here as factor <= 1.0
  }
  // sum over the categories
  res = std::accumulate(categoryLikelihoods.begin(), categoryLikelihoods.end(),
                        REAL());
  // normalize by the number of categories
  res /= static_cast<double>(getGammaCatNumber());
  scale(res);
  auto ll = getLog(res);
  // remember the computed LL for this species tree
  this->_llCache[hash] = ll;
  if (this->_memorySavings) {
    // delete CLVs to use less RAM
    deallocateMemory();
  }
  return ll;
}

template <class REAL>
bool MultiModel<REAL>::inferMLScenario(Scenario &scenario) {
  assert(false); // not implemented for MultiModel
  return computeScenario(scenario, false);
}

template <class REAL>
bool MultiModel<REAL>::sampleReconciliations(
    unsigned int samples, std::vector<std::shared_ptr<Scenario>> &scenarios) {
  if (this->_memorySavings) {
    // initialize CLVs, as they have not been created with the model
    allocateMemory();
  }
  // compute CLVs based on recomputed per-branch probas
  updateCLVs();
  // sample and save scenarios
  bool ok = true;
  for (unsigned int i = 0; i < samples; ++i) {
    scenarios.push_back(std::make_shared<Scenario>());
    ok &= computeScenario(*scenarios.back(), true);
    if (!ok)
      break;
  }
  if (this->_memorySavings) {
    // delete CLVs to use less RAM
    deallocateMemory();
  }
  return ok;
}

template <class REAL> void MultiModel<REAL>::updateCLVs() {
  // recompute the per-branch probas of elementary events
  this->beforeComputeCLVs();
  // perform updateCLV() for each observed gene clade:
  // postorder gene clade traversal is granted
  for (CID cid = 0; cid < this->_ccp.getCladesNumber(); ++cid) {
    updateCLV(cid);
  }
  this->_allSpeciesNodesInvalid = false;
  this->_invalidatedSpeciesNodes.clear();
}

template <class REAL>
corax_rnode_t *
MultiModel<REAL>::sampleOriginationSpecies(unsigned int &category) {
  auto totalLikelihood = REAL();
  for (unsigned int c = 0; c < getGammaCatNumber(); ++c) {
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      totalLikelihood += getRootCladeLikelihood(speciesNode, c);
    }
  }
  auto toSample = getRandom(totalLikelihood);
  auto sumLikelihood = REAL();
  for (unsigned int c = 0; c < getGammaCatNumber(); ++c) {
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      sumLikelihood += getRootCladeLikelihood(speciesNode, c);
      if (sumLikelihood > toSample) {
        category = c;
        return speciesNode;
      }
    }
  }
  assert(false);
  return nullptr;
}

template <class REAL>
bool MultiModel<REAL>::computeScenario(Scenario &scenario, bool stochastic) {
  assert(stochastic);
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  // sample rate category and origination species
  unsigned int category = 0;
  auto originationSpecies = sampleOriginationSpecies(category);
  // init scenario
  scenario.setSpeciesTree(&this->_speciesTree);
  auto geneRoot = scenario.generateGeneRoot();
  scenario.setGeneRoot(geneRoot);
  auto virtualRootIndex = 2 * this->_ccp.getLeafNumber();
  scenario.setVirtualRootIndex(virtualRootIndex);
  // run sampling recursion
  auto ok = backtrace(rootCID, originationSpecies, geneRoot, category, scenario,
                      stochastic);
  return ok;
}

template <class REAL>
bool MultiModel<REAL>::backtrace(CID cid, corax_rnode_t *speciesNode,
                                 corax_unode_t *geneNode, unsigned int category,
                                 Scenario &scenario, bool stochastic) {
  auto c = category;
  REAL proba;
  // compute the probability of the clade on the species branch
  // to obtain the sampling probability from it further
  if (!computeProbability(cid, speciesNode, c, proba)) {
    return false;
  }
  // create a recCell, set its sampling probability and sample
  // an event of the clade on the species branch
  ReconciliationCell<REAL> recCell;
  recCell.maxProba = getRandom(proba);
  if (!computeProbability(cid, speciesNode, c, proba, &recCell)) {
    return false;
  }
  // fill the recCell's event fields:
  // - the current node index of the reconciled gene tree
  // - the current species node index
  // - the type of the previous event in the scenario
  if (scenario.getGeneNodeBuffer().size() == 1) {
    recCell.event.geneNode = scenario.getVirtualRootIndex();
  } else {
    recCell.event.geneNode = geneNode->node_index;
  }
  recCell.event.speciesNode = speciesNode->node_index;
  recCell.event.previousType = scenario.getLastEventType();
  // overwrite the recCell's event child node indices:
  // replace the cids of the sampled child clades with
  // the accordingly built child nodes of the reconciled gene tree
  CID leftCid = 0;
  CID rightCid = 0;
  corax_unode_t *leftGeneNode = nullptr;
  corax_unode_t *rightGeneNode = nullptr;
  if (recCell.event.type == ReconciliationEventType::EVENT_S ||
      recCell.event.type == ReconciliationEventType::EVENT_D ||
      recCell.event.type == ReconciliationEventType::EVENT_T) {
    // save the original clade indices
    leftCid = recCell.event.leftGeneIndex;
    rightCid = recCell.event.rightGeneIndex;
    // make child nodes for the scenario
    scenario.generateGeneChildren(geneNode, leftGeneNode, rightGeneNode);
    leftGeneNode->length = recCell.blLeft;
    rightGeneNode->length = recCell.blRight;
    recCell.event.leftGeneIndex = leftGeneNode->node_index;
    recCell.event.rightGeneIndex = rightGeneNode->node_index;
  }
  // add the overwritten recCell's event to the scenario
  bool addEvent = true;
  if (recCell.event.type == ReconciliationEventType::EVENT_DL ||
      (recCell.event.type == ReconciliationEventType::EVENT_TL &&
       recCell.event.pllDestSpeciesNode == nullptr)) {
    addEvent = false;
  }
  if (addEvent) {
    scenario.addEvent(recCell.event);
  }
  // recursion to sample for the child clades or to resample
  bool ok = true;
  std::string geneLabel;
  switch (recCell.event.type) {
  case ReconciliationEventType::EVENT_S:
    scenario.setLastEventType(recCell.event.type);
    ok &= backtrace(leftCid, this->getSpeciesLeft(speciesNode), leftGeneNode, c,
                    scenario, stochastic);
    scenario.setLastEventType(recCell.event.type);
    ok &= backtrace(rightCid, this->getSpeciesRight(speciesNode), rightGeneNode,
                    c, scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_D:
    scenario.setLastEventType(recCell.event.type);
    ok &=
        backtrace(leftCid, speciesNode, leftGeneNode, c, scenario, stochastic);
    scenario.setLastEventType(recCell.event.type);
    ok &= backtrace(rightCid, speciesNode, rightGeneNode, c, scenario,
                    stochastic);
    break;
  case ReconciliationEventType::EVENT_T:
    // source species
    ok &=
        backtrace(leftCid, speciesNode, leftGeneNode, c, scenario, stochastic);
    // receiving species
    scenario.setLastEventType(ReconciliationEventType::EVENT_T);
    ok &= backtrace(rightCid, recCell.event.pllDestSpeciesNode, rightGeneNode,
                    c, scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_SL:
    scenario.setLastEventType(recCell.event.type);
    ok &= backtrace(cid, recCell.event.pllDestSpeciesNode, geneNode, c,
                    scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_DL:
    // the gene got duplicated, but one copy was lost, we resample again
    ok &= backtrace(cid, speciesNode, geneNode, c, scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_TL:
    if (recCell.event.pllDestSpeciesNode == nullptr) {
      // the gene was lost in the receiving species, we resample again
      ok &= backtrace(cid, speciesNode, geneNode, c, scenario, stochastic);
    } else {
      scenario.setLastEventType(ReconciliationEventType::EVENT_TL);
      ok &= backtrace(cid, recCell.event.pllDestSpeciesNode, geneNode, c,
                      scenario, stochastic);
    }
    break;
  case ReconciliationEventType::EVENT_None:
    // on leaf speciation
    geneLabel = recCell.event.label;
    geneNode->label = new char[geneLabel.size() + 1];
    memcpy(geneNode->label, geneLabel.c_str(), geneLabel.size() + 1);
    break;
  default:
    ok = false;
  }
  return ok;
}
