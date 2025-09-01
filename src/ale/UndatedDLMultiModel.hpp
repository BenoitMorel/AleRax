#pragma once

#include "MultiModel.hpp"

/**
 *  Implements all HGT modelling-dependent functions of
 *  the MultiModel class in the context of the UndatedDL
 *  model
 */
template <class REAL> class UndatedDLMultiModel : public MultiModel<REAL> {
public:
  UndatedDLMultiModel(PLLRootedTree &speciesTree,
                      const GeneSpeciesMapping &geneSpeciesMapping,
                      const RecModelInfo &info, const std::string &ccpFile);

  virtual ~UndatedDLMultiModel() {}

  virtual void setAlpha(double alpha);
  virtual void setRates(const RatesVector &rates);

private:
  unsigned int _gammaCatNumber;
  std::vector<double> _gammaScalers;
  RatesVector _dlRates;
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _OP; // Origination probability, per species branch
  std::vector<REAL> _uE;   // Extinction probability, per species branch
  OriginationStrategy _originationStrategy;

  // Element e of the gene clade's DLCLV stores the probability of the clade,
  // given the clade is mapped to the species branch e.
  // In the paper: Pi_{e,gamma} of a clade gamma for each branch e
  using DLCLV = std::vector<REAL>;
  // vector of DLCLVs for all observed clades
  std::vector<DLCLV> _dlclvs;

  // functions to work with CLVs
  virtual void allocateMemory();
  virtual void deallocateMemory();
  virtual void updateCLV(CID cid);

  // functions to work with probabilities
  virtual void recomputeSpeciesProbabilities();
  virtual double getLikelihoodFactor(unsigned int category);
  virtual REAL getRootCladeLikelihood(corax_rnode_t *speciesNode,
                                      unsigned int category);
  virtual bool computeProbability(CID cid, corax_rnode_t *speciesNode,
                                  unsigned int category, REAL &proba,
                                  ReconciliationCell<REAL> *recCell = nullptr);

  // functions to work with _llCache
  virtual size_t getHash() { return this->getSpeciesTreeHash(); }
};

/**
 *  Constructor
 */
template <class REAL>
UndatedDLMultiModel<REAL>::UndatedDLMultiModel(
    PLLRootedTree &speciesTree, const GeneSpeciesMapping &geneSpeciesMapping,
    const RecModelInfo &info, const std::string &ccpFile)
    : MultiModel<REAL>(speciesTree, geneSpeciesMapping, info, ccpFile),
      _gammaCatNumber(info.gammaCategories),
      _gammaScalers(_gammaCatNumber, 1.0),
      _PD(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 0.2),
      _PL(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 0.2),
      _PS(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 1.0),
      _OP(this->getAllSpeciesNodeNumber(),
          1.0 / static_cast<double>(this->getAllSpeciesNodeNumber())),
      _uE(this->getAllSpeciesNodeNumber() * _gammaCatNumber, REAL()),
      _originationStrategy(info.originationStrategy) {
  auto N = this->getAllSpeciesNodeNumber();
  // set gamma scalers with the default alpha
  setAlpha(1.0);
  // set all DLO rates to the default value
  _dlRates.resize(this->_info.modelFreeParameters(),
                  std::vector<double>(N, 0.2));
  // initialize DLCLVs if needed
  if (!this->_memorySavings) {
    allocateMemory();
  }
}

template <class REAL> void UndatedDLMultiModel<REAL>::setAlpha(double alpha) {
  corax_compute_gamma_cats(alpha, _gammaScalers.size(), &_gammaScalers[0],
                           CORAX_GAMMA_RATES_MEAN);
  this->invalidateAllSpeciesNodes();
  this->resetCache();
}

template <class REAL>
void UndatedDLMultiModel<REAL>::setRates(const RatesVector &rates) {
  assert(rates.size() == this->_info.modelFreeParameters());
  _dlRates = rates;
  this->invalidateAllSpeciesNodes();
  this->resetCache();
}

/**
 *  Allocate memory to the CLVs
 */
template <class REAL> void UndatedDLMultiModel<REAL>::allocateMemory() {
  DLCLV nullCLV(this->getAllSpeciesNodeNumber() * _gammaCatNumber, REAL());
  _dlclvs = std::vector<DLCLV>(this->_ccp.getCladesNumber(), nullCLV);
}

/**
 *  Free memory allocated to the CLVs
 */
template <class REAL> void UndatedDLMultiModel<REAL>::deallocateMemory() {
  _dlclvs = std::vector<DLCLV>();
}

/**
 *  Compute the CLV for a given clade
 */
template <class REAL> void UndatedDLMultiModel<REAL>::updateCLV(CID cid) {
  auto &uq = _dlclvs[cid];
  std::fill(uq.begin(), uq.end(), REAL());
  bool ok;
  for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
    // postorder species tree traversal is granted
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      auto e = speciesNode->node_index;
      auto ec = e * _gammaCatNumber + c;
      REAL p = REAL();
      ok = computeProbability(cid, speciesNode, c, p);
      assert(ok);
      uq[ec] = p;
    }
  }
}

/**
 *  Compute the per species branch probabilities of
 *  the elementary events of clade evolution
 */
template <class REAL>
void UndatedDLMultiModel<REAL>::recomputeSpeciesProbabilities() {
  auto allSpeciesNumber = this->getAllSpeciesNodeNumber();
  // recompute _PD, _PL, _PS
  auto &dupRates = _dlRates[0];
  auto &lossRates = _dlRates[1];
  assert(allSpeciesNumber == dupRates.size());
  assert(allSpeciesNumber == lossRates.size());
  std::fill(_PD.begin(), _PD.end(), 0.0);
  std::fill(_PL.begin(), _PL.end(), 0.0);
  std::fill(_PS.begin(), _PS.end(), 0.0);
  for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      auto e = speciesNode->node_index;
      auto ec = e * _gammaCatNumber + c;
      _PD[ec] = dupRates[e];
      _PL[ec] = lossRates[e];
      _PS[ec] = _gammaScalers[c];
      if (this->_info.noDup) {
        _PD[ec] = 0.0;
      }
      auto sum = _PD[ec] + _PL[ec] + _PS[ec];
      _PD[ec] /= sum;
      _PL[ec] /= sum;
      _PS[ec] /= sum;
    }
  }
  // recompute _OP
  std::vector<corax_rnode_t *> speciesNodesBuffer;
  std::vector<corax_rnode_t *> *possibleOriginationSpeciesNodes = nullptr;
  switch (_originationStrategy) {
  case OriginationStrategy::UNIFORM:
  case OriginationStrategy::OPTIMIZE:
    possibleOriginationSpeciesNodes = &(this->getPrunedSpeciesNodes());
    break;
  case OriginationStrategy::ROOT:
    speciesNodesBuffer.push_back(this->_speciesTree.getRoot());
    possibleOriginationSpeciesNodes = &speciesNodesBuffer;
    break;
  case OriginationStrategy::LCA:
    speciesNodesBuffer.push_back(this->getCoveredSpeciesLCA());
    possibleOriginationSpeciesNodes = &speciesNodesBuffer;
    break;
  }
  std::fill(_OP.begin(), _OP.end(), 0.0);
  if (_originationStrategy == OriginationStrategy::OPTIMIZE) {
    auto &oriRates = _dlRates[2];
    assert(allSpeciesNumber == oriRates.size());
    double sum = 0.0;
    for (auto speciesNode : *possibleOriginationSpeciesNodes) {
      auto e = speciesNode->node_index;
      sum += oriRates[e];
    }
    for (auto speciesNode : *possibleOriginationSpeciesNodes) {
      auto e = speciesNode->node_index;
      _OP[e] = oriRates[e] / sum;
    }
  } else {
    double sum = static_cast<double>(possibleOriginationSpeciesNodes->size());
    for (auto speciesNode : *possibleOriginationSpeciesNodes) {
      auto e = speciesNode->node_index;
      _OP[e] = 1.0 / sum;
    }
  }
  // recompute _uE
  std::fill(_uE.begin(), _uE.end(), REAL());
  // iterate several times to resolve _uE probas with
  // fixed point optimization
  unsigned int maxIt = 4;
  for (unsigned int it = 0; it < maxIt; ++it) {
    for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
      // postorder species tree traversal is granted
      for (auto speciesNode : this->getPrunedSpeciesNodes()) {
        auto e = speciesNode->node_index;
        auto ec = e * _gammaCatNumber + c;
        REAL temp;
        REAL proba = REAL();
        // L scenario
        temp = REAL(_PL[ec]);
        scale(temp);
        proba += temp;
        // S scenario
        if (this->getSpeciesLeft(speciesNode)) {
          // internal branch
          auto f = this->getSpeciesLeft(speciesNode)->node_index;
          auto g = this->getSpeciesRight(speciesNode)->node_index;
          auto fc = f * _gammaCatNumber + c;
          auto gc = g * _gammaCatNumber + c;
          temp = _uE[fc] * (_uE[gc] * _PS[ec]); // SEE scenario
        } else {
          // terminal branch
          temp = REAL(_PS[ec] * this->_fm[e]); // S but not observed scenario
        }
        scale(temp);
        proba += temp;
        // DEE scenario
        temp = _uE[ec] * (_uE[ec] * _PD[ec]);
        scale(temp);
        proba += temp;
        assert(proba < REAL(1.000001));
        _uE[ec] = proba;
      }
    }
  } // end of iteration
}

/**
 *  Correction factor to the species tree likelihood,
 *  because we condition on survival
 */
template <class REAL>
double UndatedDLMultiModel<REAL>::getLikelihoodFactor(unsigned int category) {
  double factor(0.0);
  auto c = category;
  for (auto speciesNode : this->getPrunedSpeciesNodes()) {
    auto e = speciesNode->node_index;
    auto ec = e * _gammaCatNumber + c;
    factor += (1.0 - _uE[ec]) * _OP[e];
  }
  return factor;
}

/**
 *  Probability of the current family to evolve starting
 *  from a given species branch
 */
template <class REAL>
REAL UndatedDLMultiModel<REAL>::getRootCladeLikelihood(
    corax_rnode_t *speciesNode, unsigned int category) {
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto c = category;
  auto e = speciesNode->node_index;
  auto ec = e * _gammaCatNumber + c;
  REAL likelihood = _dlclvs[rootCID][ec] * _OP[e];
  scale(likelihood);
  return likelihood;
}

/**
 *  Compute the CLV value for a given cid (clade id) and a given
 *  species node and write it to the proba variable
 */
template <class REAL>
bool UndatedDLMultiModel<REAL>::computeProbability(
    CID cid, corax_rnode_t *speciesNode, unsigned int category, REAL &proba,
    ReconciliationCell<REAL> *recCell) {
  proba = REAL();
  REAL temp;
  bool isGeneLeaf = this->_ccp.isLeaf(cid);
  bool isSpeciesLeaf = !this->getSpeciesLeft(speciesNode);
  auto c = category;
  auto e = speciesNode->node_index;
  auto ec = e * _gammaCatNumber + c;
  unsigned int f = 0;
  unsigned int g = 0;
  unsigned int fc = 0;
  unsigned int gc = 0;
  if (!isSpeciesLeaf) {
    f = this->getSpeciesLeft(speciesNode)->node_index;
    g = this->getSpeciesRight(speciesNode)->node_index;
    fc = f * _gammaCatNumber + c;
    gc = g * _gammaCatNumber + c;
  }
  REAL maxProba = REAL();
  if (recCell) {
    recCell->event.geneNode = cid;
    recCell->event.speciesNode = e;
    maxProba = recCell->maxProba;
  }
  // S events on terminal species branches can happen
  // for terminal gene nodes only:
  if (isGeneLeaf) {
    // - S event on a terminal species branch (only for compatible genes and
    // species)
    if (isSpeciesLeaf && this->_geneToSpecies[cid] == e) {
      temp = REAL(_PS[ec]);
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_None;
        recCell->event.label = this->_ccp.getLeafLabel(cid);
        return true;
      }
    }
  }
  // S events on internal species branches and D events can happen
  // for ancestral gene nodes only:
  for (const auto &cladeSplit : this->_ccp.getCladeSplits(cid)) {
    auto cidLeft = cladeSplit.left;
    auto cidRight = cladeSplit.right;
    auto freq = cladeSplit.frequency;
    // - S event on an internal species branch
    if (!isSpeciesLeaf) {
      temp = _dlclvs[cidLeft][fc] * _dlclvs[cidRight][gc] * (_PS[ec] * freq);
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_S;
        recCell->event.leftGeneIndex = cidLeft;
        recCell->event.rightGeneIndex = cidRight;
        recCell->blLeft = cladeSplit.blLeft;
        recCell->blRight = cladeSplit.blRight;
        return true;
      }
      temp = _dlclvs[cidRight][fc] * _dlclvs[cidLeft][gc] * (_PS[ec] * freq);
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_S;
        recCell->event.leftGeneIndex = cidRight;
        recCell->event.rightGeneIndex = cidLeft;
        recCell->blLeft = cladeSplit.blRight;
        recCell->blRight = cladeSplit.blLeft;
        return true;
      }
    }
    // - D event
    temp = _dlclvs[cidLeft][ec] * _dlclvs[cidRight][ec] * (_PD[ec] * freq);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_D;
      recCell->event.leftGeneIndex = cidLeft;
      recCell->event.rightGeneIndex = cidRight;
      recCell->blLeft = cladeSplit.blLeft;
      recCell->blRight = cladeSplit.blRight;
      return true;
    }
  }
  // SL events and DL events can happen
  // for any of gene nodes:
  // - SL event (only on an internal species branch)
  if (!isSpeciesLeaf) {
    temp = _dlclvs[cid][fc] * (_uE[gc] * _PS[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.lostSpeciesNode = g;
      recCell->event.pllDestSpeciesNode = this->getSpeciesLeft(speciesNode);
      recCell->event.pllLostSpeciesNode = this->getSpeciesRight(speciesNode);
      return true;
    }
    temp = _dlclvs[cid][gc] * (_uE[fc] * _PS[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.lostSpeciesNode = f;
      recCell->event.pllDestSpeciesNode = this->getSpeciesRight(speciesNode);
      recCell->event.pllLostSpeciesNode = this->getSpeciesLeft(speciesNode);
      return true;
    }
  }
  if (!this->_info.noDL) {
    // - DL event
    temp = proba / (1.0 - (_uE[ec] * _PD[ec] * 2.0));
    scale(temp);
    proba = temp;
    if (recCell && proba > maxProba) {
      // in fact, nothing happens, we'll have to resample
      recCell->event.type = ReconciliationEventType::EVENT_DL;
      return true;
    }
  }
  if (recCell) {
    Logger::error << "error: proba=" << proba << ", maxProba=" << maxProba
                  << " (proba < maxProba)" << std::endl;
    return false; // we haven't sampled any event, this should not happen
  }
  if (proba > REAL(1.0)) {
    Logger::error << "error: proba=" << proba << " (proba > 1)" << std::endl;
    return false;
  }
  return true;
}
