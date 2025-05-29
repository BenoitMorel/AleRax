#pragma once

#include <trees/DatedTree.hpp>

#include "MultiModel.hpp"

/**
 *  Highway with the proba scaled per rate category
 */
struct WeightedHighway {
  Highway highway;
  std::vector<double> proba;
};

/**
 *  Implements all HGT modelling-dependent functions of
 *  the MultiModel class in the context of the UndatedDTL
 *  model
 */
template <class REAL> class UndatedDTLMultiModel : public MultiModel<REAL> {
public:
  UndatedDTLMultiModel(DatedTree &speciesTree,
                       const GeneSpeciesMapping &geneSpeciesMapping,
                       const RecModelInfo &info, const std::string &ccpFile);

  virtual ~UndatedDTLMultiModel() {}

  virtual void setAlpha(double alpha);
  virtual void setRates(const RatesVector &rates);
  virtual void setHighways(const std::vector<Highway> &highways);

private:
  DatedTree &_datedTree;
  unsigned int _gammaCatNumber;
  std::vector<double> _gammaScalers;
  RatesVector _dtlRates;
  std::vector<std::vector<WeightedHighway>>
      _highways;           // Highways, per species branch
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PT; // Transfer probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _OP; // Origination probability, per species branch
  std::vector<REAL> _uE;   // Extinction probability, per species branch
  std::vector<REAL> _tE; // Transfer-extinction probability, per species branch
  OriginationStrategy _originationStrategy;
  TransferConstaint _transferConstraint;
  // Allowed transfers, per species branch
  std::vector<std::vector<corax_rnode_t *>> _transferCandidateSpeciesNodes;

  struct DTLCLV {
    // Element e of the gene clade's _uq stores the probability of the clade,
    // given the clade is mapped to the species branch e.
    // In the paper: Pi_{e,gamma} of a clade gamma for each branch e
    std::vector<REAL> _uq;
    // Element e of the gene clade's _ut stores the probability of the clade
    // to be transferred from the species node e to some other species node.
    // In the paper: \bar{Pi}_{e,gamma} of a clade gamma for each branch e
    std::vector<REAL> _tq;
    DTLCLV() : _uq(0), _tq(0) {}
    DTLCLV(unsigned int speciesNumber, unsigned int gammaCategories)
        : _uq(speciesNumber * gammaCategories, REAL()),
          _tq(speciesNumber * gammaCategories, REAL()) {}
  };
  // vector of DTLCLVs for all observed clades
  std::vector<DTLCLV> _dtlclvs;

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

  // functions to deal with transfers
  void updateTransferCandidates();
  bool sampleTransferEvent(CID cid, corax_rnode_t *srcSpeciesNode,
                           unsigned int category, Scenario::Event &event);

  double getTransferWeightNorm(unsigned int speciesNodeIndex) const {
    auto e = speciesNodeIndex;
    return static_cast<double>(_transferCandidateSpeciesNodes[e].size());
  }

  // functions to work with _llCache
  virtual size_t getHash();
};

/**
 *  Constructor
 */
template <class REAL>
UndatedDTLMultiModel<REAL>::UndatedDTLMultiModel(
    DatedTree &speciesTree, const GeneSpeciesMapping &geneSpeciesMapping,
    const RecModelInfo &info, const std::string &ccpFile)
    : MultiModel<REAL>(speciesTree.getRootedTree(), geneSpeciesMapping, info,
                       ccpFile),
      _datedTree(speciesTree), _gammaCatNumber(info.gammaCategories),
      _gammaScalers(_gammaCatNumber, 1.0),
      _PD(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 0.2),
      _PL(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 0.2),
      _PT(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 0.1),
      _PS(this->getAllSpeciesNodeNumber() * _gammaCatNumber, 1.0),
      _OP(this->getAllSpeciesNodeNumber(),
          1.0 / static_cast<double>(this->getAllSpeciesNodeNumber())),
      _uE(this->getAllSpeciesNodeNumber() * _gammaCatNumber, REAL()),
      _tE(this->getAllSpeciesNodeNumber() * _gammaCatNumber, REAL()),
      _originationStrategy(info.originationStrategy),
      _transferConstraint(info.transferConstraint) {
  auto N = this->getAllSpeciesNodeNumber();
  // set gamma scalers with the default alpha
  setAlpha(1.0);
  // set all DTLO rates to the default value
  _dtlRates.resize(this->_info.modelFreeParameters(),
                   std::vector<double>(N, 0.2));
  // initialize highways and per-species transfer candidates
  _highways.resize(N);
  _transferCandidateSpeciesNodes.resize(N);
  // initialize DTLCLVs if needed
  if (!this->_memorySavings) {
    allocateMemory();
  }
}

template <class REAL> void UndatedDTLMultiModel<REAL>::setAlpha(double alpha) {
  corax_compute_gamma_cats(alpha, _gammaScalers.size(), &_gammaScalers[0],
                           CORAX_GAMMA_RATES_MEAN);
  this->invalidateAllSpeciesNodes();
  this->resetCache();
}

template <class REAL>
void UndatedDTLMultiModel<REAL>::setRates(const RatesVector &rates) {
  assert(rates.size() == this->_info.modelFreeParameters());
  _dtlRates = rates;
  this->invalidateAllSpeciesNodes();
  this->resetCache();
}

template <class REAL>
void UndatedDTLMultiModel<REAL>::setHighways(
    const std::vector<Highway> &highways) {
  for (auto &speciesWeightedHighways : _highways) {
    speciesWeightedHighways.clear();
  }
  for (auto highway : highways) {
    WeightedHighway hp;
    hp.highway = highway;
    // map the highway to the pruned species tree
    if (this->prunedMode()) {
      hp.highway.src = this->_speciesToPrunedNode[highway.src->node_index];
      hp.highway.dest = this->_speciesToPrunedNode[highway.dest->node_index];
    } else {
      hp.highway.src = highway.src;
      hp.highway.dest = highway.dest;
    }
    // in the pruned mode a highway from/to a species not covered by the gene
    // family may have the src or dest species absent or may look like a
    // self-transfer
    if (!hp.highway.src || !hp.highway.dest ||
        (hp.highway.src == hp.highway.dest)) {
      // this highway should not affect this gene family
      continue;
    }
    // this value will be normalized later on, separately per category
    hp.proba.resize(_gammaCatNumber, highway.proba);
    // update _highways
    _highways[hp.highway.src->node_index].push_back(hp);
  }
  this->invalidateAllSpeciesNodes();
  this->resetCache();
}

/**
 *  Allocate memory to the CLVs
 */
template <class REAL> void UndatedDTLMultiModel<REAL>::allocateMemory() {
  DTLCLV nullCLV(this->getAllSpeciesNodeNumber(), _gammaCatNumber);
  _dtlclvs = std::vector<DTLCLV>(this->_ccp.getCladesNumber(), nullCLV);
}

/**
 *  Free memory allocated to the CLVs
 */
template <class REAL> void UndatedDTLMultiModel<REAL>::deallocateMemory() {
  _dtlclvs = std::vector<DTLCLV>();
}

/**
 *  Compute the CLV for a given clade
 */
template <class REAL> void UndatedDTLMultiModel<REAL>::updateCLV(CID cid) {
  auto &uq = _dtlclvs[cid]._uq;
  auto &tq = _dtlclvs[cid]._tq;
  std::fill(uq.begin(), uq.end(), REAL());
  std::fill(tq.begin(), tq.end(), REAL());
  // iterate several times to resolve the DL and TL terms with
  // fixed point optimization: not needed if we don't model TL
  unsigned int maxIt = this->_info.noTL ? 1 : 4;
  for (unsigned int it = 0; it < maxIt; ++it) {
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
      // now that we've got the clade proba on every species branch,
      // we can compute the clade transfer probas
      for (auto speciesNode : this->getPrunedSpeciesNodes()) {
        auto e = speciesNode->node_index;
        auto ec = e * _gammaCatNumber + c;
        REAL p = REAL();
        for (auto destSpeciesNode : _transferCandidateSpeciesNodes[e]) {
          auto d = destSpeciesNode->node_index;
          auto dc = d * _gammaCatNumber + c;
          p += uq[dc];
        }
        p /= getTransferWeightNorm(e);
        scale(p);
        tq[ec] = p;
      }
    }
  }
}

/**
 *  Update the list of the allowed transfer-receiving species
 *  for each species branch
 */
template <class REAL>
void UndatedDTLMultiModel<REAL>::updateTransferCandidates() {
  for (auto &speciesTransferCandidates : _transferCandidateSpeciesNodes) {
    speciesTransferCandidates.clear();
  }
  for (auto speciesNode : this->getPrunedSpeciesNodes()) {
    auto e = speciesNode->node_index;
    if (_transferConstraint == TransferConstaint::NONE) {
      // include all the species nodes except for the self
      for (auto destSpeciesNode : this->getPrunedSpeciesNodes()) {
        auto d = destSpeciesNode->node_index;
        if (d == e) {
          continue;
        }
        _transferCandidateSpeciesNodes[e].push_back(destSpeciesNode);
      }
    }
    if (_transferConstraint == TransferConstaint::PARENTS) {
      // identify all parents of the self
      std::unordered_set<unsigned int> parents;
      auto parent = speciesNode;
      while (parent) {
        parents.insert(parent->node_index);
        parent = parent->parent;
      }
      // include all the species nodes except for the parents
      for (auto destSpeciesNode : this->getPrunedSpeciesNodes()) {
        auto d = destSpeciesNode->node_index;
        if (parents.end() != parents.find(d)) {
          continue;
        }
        _transferCandidateSpeciesNodes[e].push_back(destSpeciesNode);
      }
    }
    if (_transferConstraint == TransferConstaint::RELDATED) {
      // include all the species nodes younger than the self's parent
      for (auto destSpeciesNode : this->getPrunedSpeciesNodes()) {
        auto d = destSpeciesNode->node_index;
        if (!_datedTree.canTransferUnderRelDated(e, d)) {
          continue;
        }
        _transferCandidateSpeciesNodes[e].push_back(destSpeciesNode);
      }
    }
  }
}

/**
 *  Compute the per species branch probabilities of
 *  the elementary events of clade evolution
 */
template <class REAL>
void UndatedDTLMultiModel<REAL>::recomputeSpeciesProbabilities() {
  auto allSpeciesNumber = this->getAllSpeciesNodeNumber();
  // recompute _PD, _PL, _PT, _PS and highway.proba
  auto &dupRates = _dtlRates[0];
  auto &lossRates = _dtlRates[1];
  auto &transferRates = _dtlRates[2];
  assert(allSpeciesNumber == dupRates.size());
  assert(allSpeciesNumber == lossRates.size());
  assert(allSpeciesNumber == transferRates.size());
  std::fill(_PD.begin(), _PD.end(), 0.0);
  std::fill(_PL.begin(), _PL.end(), 0.0);
  std::fill(_PT.begin(), _PT.end(), 0.0);
  std::fill(_PS.begin(), _PS.end(), 0.0);
  for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      auto e = speciesNode->node_index;
      auto ec = e * _gammaCatNumber + c;
      _PD[ec] = dupRates[e];
      _PL[ec] = lossRates[e];
      _PT[ec] = transferRates[e];
      _PS[ec] = _gammaScalers[c];
      if (this->_info.noDup) {
        _PD[ec] = 0.0;
      }
      auto sum = _PD[ec] + _PL[ec] + _PT[ec] + _PS[ec];
      for (const auto &highway : _highways[e]) {
        sum += highway.highway.proba;
      }
      _PD[ec] /= sum;
      _PL[ec] /= sum;
      _PT[ec] /= sum;
      _PS[ec] /= sum;
      for (auto &highway : _highways[e]) {
        assert(highway.highway.proba >= 0.0);
        // highway proba, normalized per category
        highway.proba[c] = highway.highway.proba / sum;
        assert(highway.proba[c] < 1.0);
      }
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
    auto &oriRates = _dtlRates[3];
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
  // recompute _uE and _tE
  updateTransferCandidates();
  std::fill(_uE.begin(), _uE.end(), REAL());
  std::fill(_tE.begin(), _tE.end(), REAL());
  // iterate several times to resolve _uE and _tE probas with
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
        // TEE scenario
        temp = _uE[ec] * (_tE[ec] * _PT[ec]);
        scale(temp);
        proba += temp;
        // highway TEE scenario
        for (const auto &highway : _highways[e]) {
          auto d = highway.highway.dest->node_index;
          auto dc = d * _gammaCatNumber + c;
          temp = _uE[ec] * (_uE[dc] * highway.proba[c]);
          scale(temp);
          proba += temp;
        }
        assert(proba < REAL(1.000001));
        _uE[ec] = proba;
      }
      // now that we've got extinction probas for every species branch,
      // we can compute transfer-extinction probas
      for (auto speciesNode : this->getPrunedSpeciesNodes()) {
        auto e = speciesNode->node_index;
        auto ec = e * _gammaCatNumber + c;
        REAL proba = REAL();
        for (auto destSpeciesNode : _transferCandidateSpeciesNodes[e]) {
          auto d = destSpeciesNode->node_index;
          auto dc = d * _gammaCatNumber + c;
          proba += _uE[dc];
        }
        proba /= getTransferWeightNorm(e);
        scale(proba);
        assert(proba < REAL(1.000001));
        _tE[ec] = proba;
      }
    }
  } // end of iteration
}

/**
 *  Correction factor to the species tree likelihood,
 *  because we condition on survival
 */
template <class REAL>
double UndatedDTLMultiModel<REAL>::getLikelihoodFactor(unsigned int category) {
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
REAL UndatedDTLMultiModel<REAL>::getRootCladeLikelihood(
    corax_rnode_t *speciesNode, unsigned int category) {
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto c = category;
  auto e = speciesNode->node_index;
  auto ec = e * _gammaCatNumber + c;
  REAL likelihood = _dtlclvs[rootCID]._uq[ec] * _OP[e];
  scale(likelihood);
  return likelihood;
}

/**
 *  Sample a transfer destination species branch and
 *  write it into the given recCell event
 */
template <class REAL>
bool UndatedDTLMultiModel<REAL>::sampleTransferEvent(
    unsigned int cid, corax_rnode_t *srcSpeciesNode, unsigned int category,
    Scenario::Event &event) {
  auto c = category;
  auto e = srcSpeciesNode->node_index;
  auto ec = e * _gammaCatNumber + c;
  auto survivingTransferSum = _dtlclvs[cid]._tq[ec] * getTransferWeightNorm(e);
  auto toSample = this->getRandom(survivingTransferSum);
  auto sumProba = REAL();
  for (auto destSpeciesNode : _transferCandidateSpeciesNodes[e]) {
    auto d = destSpeciesNode->node_index;
    auto dc = d * _gammaCatNumber + c;
    sumProba += _dtlclvs[cid]._uq[dc];
    if (sumProba > toSample) {
      event.destSpeciesNode = d;
      event.pllDestSpeciesNode = destSpeciesNode;
      return true;
    }
  }
  return false;
}

/**
 *  Compute the CLV value for a given cid (clade id) and a given
 *  species node and write it to the proba variable
 */
template <class REAL>
bool UndatedDTLMultiModel<REAL>::computeProbability(
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
  // S events on internal species branches, D events and T events can happen
  // for ancestral gene nodes only:
  for (const auto &cladeSplit : this->_ccp.getCladeSplits(cid)) {
    auto cidLeft = cladeSplit.left;
    auto cidRight = cladeSplit.right;
    auto freq = cladeSplit.frequency;
    // - S event on an internal species branch
    if (!isSpeciesLeaf) {
      temp = _dtlclvs[cidLeft]._uq[fc] * _dtlclvs[cidRight]._uq[gc] *
             (_PS[ec] * freq);
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
      temp = _dtlclvs[cidRight]._uq[fc] * _dtlclvs[cidLeft]._uq[gc] *
             (_PS[ec] * freq);
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
    temp = _dtlclvs[cidLeft]._uq[ec] * _dtlclvs[cidRight]._uq[ec] *
           (_PD[ec] * freq);
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
    // - T event
    temp = _dtlclvs[cidLeft]._uq[ec] * _dtlclvs[cidRight]._tq[ec] *
           (_PT[ec] * freq);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_T;
      if (!sampleTransferEvent(cidRight, speciesNode, c, recCell->event)) {
        return false;
      }
      recCell->event.leftGeneIndex = cidLeft;
      recCell->event.rightGeneIndex = cidRight;
      recCell->blLeft = cladeSplit.blLeft;
      recCell->blRight = cladeSplit.blRight;
      return true;
    }
    temp = _dtlclvs[cidRight]._uq[ec] * _dtlclvs[cidLeft]._tq[ec] *
           (_PT[ec] * freq);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_T;
      if (!sampleTransferEvent(cidLeft, speciesNode, c, recCell->event)) {
        return false;
      }
      recCell->event.leftGeneIndex = cidRight;
      recCell->event.rightGeneIndex = cidLeft;
      recCell->blLeft = cladeSplit.blRight;
      recCell->blRight = cladeSplit.blLeft;
      return true;
    }
    // - highway T event
    for (const auto &highway : _highways[e]) {
      auto d = highway.highway.dest->node_index;
      auto dc = d * _gammaCatNumber + c;
      temp = _dtlclvs[cidLeft]._uq[ec] * _dtlclvs[cidRight]._uq[dc] *
             (highway.proba[c] * freq);
      scale(temp);
      proba += temp;
      if (proba > REAL(1.0)) {
        std::cerr << "error " << _dtlclvs[cidLeft]._uq[ec] << " "
                  << _dtlclvs[cidRight]._uq[dc] << " " << highway.proba[c]
                  << " " << freq << std::endl;
      }
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_T;
        recCell->event.destSpeciesNode = d;
        recCell->event.pllDestSpeciesNode = highway.highway.dest;
        recCell->event.leftGeneIndex = cidLeft;
        recCell->event.rightGeneIndex = cidRight;
        recCell->blLeft = cladeSplit.blLeft;
        recCell->blRight = cladeSplit.blRight;
        return true;
      }
      temp = _dtlclvs[cidRight]._uq[ec] * _dtlclvs[cidLeft]._uq[dc] *
             (highway.proba[c] * freq);
      scale(temp);
      proba += temp;
      if (proba > REAL(1.0)) {
        std::cerr << "error " << _dtlclvs[cidRight]._uq[ec] << " "
                  << _dtlclvs[cidLeft]._uq[dc] << " " << highway.proba[c] << " "
                  << freq << std::endl;
      }
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_T;
        recCell->event.destSpeciesNode = d;
        recCell->event.pllDestSpeciesNode = highway.highway.dest;
        recCell->event.leftGeneIndex = cidRight;
        recCell->event.rightGeneIndex = cidLeft;
        recCell->blLeft = cladeSplit.blRight;
        recCell->blRight = cladeSplit.blLeft;
        return true;
      }
    }
  }
  // SL events, DL events and TL events can happen
  // for any of gene nodes:
  // - SL event (only on an internal species branch)
  if (!isSpeciesLeaf) {
    temp = _dtlclvs[cid]._uq[fc] * (_uE[gc] * _PS[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.lostSpeciesNode = g;
      recCell->event.pllDestSpeciesNode = this->getSpeciesLeft(speciesNode);
      recCell->event.pllLostSpeciesNode = this->getSpeciesRight(speciesNode);
      return true;
    }
    temp = _dtlclvs[cid]._uq[gc] * (_uE[fc] * _PS[ec]);
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
  if (this->_info.noTL) { // no TL events allowed
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
  } else { // TL events allowed
    if (!this->_info.noDL) {
      // - DL event
      temp = _dtlclvs[cid]._uq[ec] * (_uE[ec] * _PD[ec] * 2.0);
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        // in fact, nothing happens, we'll have to resample
        recCell->event.type = ReconciliationEventType::EVENT_DL;
        return true;
      }
    }
    // - TL event
    // we transfer, but the gene gets extinct in the receiving species
    temp = _dtlclvs[cid]._uq[ec] * (_tE[ec] * _PT[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      // in fact, nothing happens, we'll have to resample
      recCell->event.type = ReconciliationEventType::EVENT_TL;
      recCell->event.pllDestSpeciesNode = nullptr;
      return true;
    }
    // we transfer, and the gene gets extinct in the sending species
    temp = _dtlclvs[cid]._tq[ec] * (_uE[ec] * _PT[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_TL;
      if (!sampleTransferEvent(cid, speciesNode, c, recCell->event)) {
        return false;
      }
      return true;
    }
    // - highway TL event
    for (const auto &highway : _highways[e]) {
      auto d = highway.highway.dest->node_index;
      auto dc = d * _gammaCatNumber + c;
      // we transfer, but the gene gets extinct in the receiving species
      temp = _dtlclvs[cid]._uq[ec] * (_uE[dc] * highway.proba[c]);
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        // in fact, nothing happens, we'll have to resample
        recCell->event.type = ReconciliationEventType::EVENT_TL;
        recCell->event.pllDestSpeciesNode = nullptr;
        return true;
      }
      // we transfer, and the gene gets extinct in the sending species
      temp = _dtlclvs[cid]._uq[dc] * (_uE[ec] * highway.proba[c]);
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_TL;
        recCell->event.destSpeciesNode = d;
        recCell->event.pllDestSpeciesNode = highway.highway.dest;
        return true;
      }
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

template <class REAL> size_t UndatedDTLMultiModel<REAL>::getHash() {
  auto hash = this->getSpeciesTreeHash();
  switch (_transferConstraint) {
  case TransferConstaint::NONE:
  case TransferConstaint::PARENTS:
    return hash;
  case TransferConstaint::RELDATED:
    return this->_datedTree.getOrderingHash(hash);
  default:
    assert(false);
  }
  return hash;
}
