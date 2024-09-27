#pragma once

#include "MultiModel.hpp"
#include <IO/GeneSpeciesMapping.hpp>
#include <ccp/ConditionalClades.hpp>
#include <maths/ScaledValue.hpp>
#include <trees/DatedTree.hpp>
#include <trees/PLLRootedTree.hpp>

class RecModelInfo;
double log(ScaledValue v);

struct WeightedHighway {
  Highway highway;
  double proba;
};

template <class REAL>
class UndatedDTLMultiModel : public MultiModelTemplate<REAL> {

public:
  UndatedDTLMultiModel(DatedTree &speciesTree,
                       const GeneSpeciesMapping &geneSpeciesMapping,
                       const RecModelInfo &info, const std::string &ccpFile);

  virtual ~UndatedDTLMultiModel() {}

  virtual void setRates(const RatesVector &);
  virtual void setAlpha(double alpha);
  virtual double computeLogLikelihood();
  virtual corax_rnode_t *sampleSpeciesNode(unsigned int &category);
  virtual void setHighways(const std::vector<Highway> &highways) {
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
      if (/*!_speciesToPrunedNode[highway.dest->node_index] ||*/ !hp.highway
              .src ||
          !hp.highway.dest) {
        // this highway should not affect this family
        continue;
      }
      hp.proba = highway.proba; // this value will be normalized later on
      if (hp.highway.src != hp.highway.dest) {
        // do not count transfers to self
        _highways[hp.highway.src->node_index].push_back(hp);
      }
    }
    resetCache();
    recomputeSpeciesProbabilities();
  }

  virtual void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) {
    MultiModel::onSpeciesTreeChange(nodesToInvalidate);
    resetCache();
  }

protected:
  virtual void updateCLVs();
  virtual void allocateMemory();
  virtual void deallocateMemory();

private:
  REAL getTransferSum(unsigned int cid, unsigned int e, unsigned int c);

  DatedTree &_datedTree;
  size_t _gammaCatNumber;
  std::vector<double> _gammaScalers;
  RatesVector _dtlRates;
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PT; // Transfer probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _OP; // Origination probability, per species branch
  std::vector<REAL> _uE;   // Extinction probability, per species branch
  std::vector<REAL> _uEBar;
  std::unordered_map<size_t, double> _llCache;
  std::vector<std::vector<WeightedHighway>> _highways;
  TransferConstaint _transferConstraint;
  OriginationStrategy _originationStrategy;
  /**
   *  All intermediate results needed to compute the reconciliation likelihood
   *  each gene node has one DTLCLV object
   *  Each DTLCLV gene  object is a function of the DTLCLVs of the direct
   * children genes
   */
  struct DTLCLV {
    DTLCLV() : _survivingTransferSum(REAL()) {}

    DTLCLV(size_t speciesNumber, size_t gammaCategories)
        : _uq(speciesNumber * gammaCategories, REAL()),
          _correctionSum(speciesNumber * gammaCategories, REAL()),
          _correctionNorm(speciesNumber * gammaCategories, 0.0),
          _survivingTransferSum(gammaCategories, REAL()) {}
    // probability of a gene node rooted at a species node
    std::vector<REAL> _uq;
    std::vector<REAL> _correctionSum;
    std::vector<double> _correctionNorm;

    // sum of transfer probabilities. Can be computed only once
    // for all species, to reduce computation complexity
    // (each entry corresponds to one category)
    std::vector<REAL> _survivingTransferSum;
  };

  std::vector<DTLCLV> _dtlclvs;

  void updateCLV(CID cid);
  double getLikelihoodFactor(unsigned int category);
  virtual void recomputeSpeciesProbabilities();
  virtual bool computeProbability(CID cid, corax_rnode_t *speciesNode,
                                  size_t category, REAL &proba,
                                  ReconciliationCell<REAL> *recCell = nullptr);
  bool sampleTransferEvent(unsigned int cid, corax_rnode_t *originSpeciesNode,
                           size_t category, Scenario::Event &event);
  size_t getHash();
  void resetCache() { _llCache.clear(); }

  double getTransferWeightNorm() const {
    return double(this->getPrunedSpeciesNodeNumber());
  }
};

template <class REAL>
UndatedDTLMultiModel<REAL>::UndatedDTLMultiModel(
    DatedTree &speciesTree, const GeneSpeciesMapping &geneSpeciesMapping,
    const RecModelInfo &info, const std::string &ccpFile)
    : MultiModelTemplate<REAL>(speciesTree.getRootedTree(), geneSpeciesMapping,
                               info, ccpFile),
      _datedTree(speciesTree), _gammaCatNumber(info.gammaCategories),
      _gammaScalers(_gammaCatNumber, 1.0),
      _PD(this->_speciesTree.getNodeNumber() * _gammaCatNumber, 0.2),
      _PL(this->_speciesTree.getNodeNumber() * _gammaCatNumber, 0.2),
      _PT(this->_speciesTree.getNodeNumber() * _gammaCatNumber, 0.1),
      _PS(this->_speciesTree.getNodeNumber() * _gammaCatNumber, 1.0),
      _OP(this->_speciesTree.getNodeNumber(),
          1.0 / static_cast<double>(this->_speciesTree.getNodeNumber())),
      _uE(this->_speciesTree.getNodeNumber() * _gammaCatNumber, REAL()),
      _uEBar(this->_speciesTree.getNodeNumber() * _gammaCatNumber, REAL()),
      _transferConstraint(info.transferConstraint),
      _originationStrategy(info.originationStrategy) {
  auto N = this->_speciesTree.getNodeNumber();
  _highways.resize(N);
  _dtlRates.resize(this->_info.modelFreeParameters(),
                   std::vector<double>(N, 0.2));
  this->onSpeciesTreeChange(nullptr);
  setAlpha(1.0);
  if (!this->_memorySavings) {
    allocateMemory();
  }
}

template <class REAL> void UndatedDTLMultiModel<REAL>::allocateMemory() {
  DTLCLV nullCLV(this->_speciesTree.getNodeNumber(), _gammaCatNumber);
  _dtlclvs = std::vector<DTLCLV>(2 * (this->_ccp.getCladesNumber()), nullCLV);
}

template <class REAL> void UndatedDTLMultiModel<REAL>::deallocateMemory() {
  _dtlclvs.clear();
}

template <class REAL>
void UndatedDTLMultiModel<REAL>::setRates(const RatesVector &rates) {
  resetCache();
  assert(rates.size() == this->_info.modelFreeParameters());
  _dtlRates = rates;
  recomputeSpeciesProbabilities();
}
template <class REAL> void UndatedDTLMultiModel<REAL>::setAlpha(double alpha) {
  resetCache();
  corax_compute_gamma_cats(alpha, _gammaScalers.size(), &_gammaScalers[0],
                           CORAX_GAMMA_RATES_MEAN);
  recomputeSpeciesProbabilities();
}

/*
 *  We fill the intermediate likelihood table (P_{e,u} in the paper) for a given
 *  gene CCP
 */
template <class REAL> void UndatedDTLMultiModel<REAL>::updateCLV(CID cid) {
  auto &clv = _dtlclvs[cid];
  auto &uq = clv._uq;
  auto &correctionSum = clv._correctionSum;
  auto &correctionNorm = clv._correctionNorm;

  std::fill(uq.begin(), uq.end(), REAL());
  auto tempUq = uq;
  auto N = this->getPrunedSpeciesNodeNumber();
  unsigned int maxIt = this->_info.noTL ? 1 : 4;
  std::fill(correctionSum.begin(), correctionSum.end(), REAL());
  std::fill(correctionNorm.begin(), correctionNorm.end(), N);

  if (_transferConstraint == TransferConstaint::PARENTS) {
    auto postOrder = this->_speciesTree.getPostOrderNodes();
    for (auto it = postOrder.rbegin(); it != postOrder.rend(); ++it) {
      auto speciesNode = *it;
      auto e = speciesNode->node_index;
      for (size_t c = 0; c < _gammaCatNumber; ++c) {
        auto parent = speciesNode;
        auto ec = e * _gammaCatNumber + c;
        while (parent) {
          correctionNorm[ec] -= 1.0;
          auto p = parent->node_index;
          parent = parent->parent;
        }
      }
    }
  }
  std::fill(clv._survivingTransferSum.begin(), clv._survivingTransferSum.end(),
            REAL());
  // iterate several times to resolve the TL term with
  // fixed point optimizaiton
  for (unsigned int it = 0; it < maxIt; ++it) {
    std::vector<REAL> sums(_gammaCatNumber, REAL());
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      auto e = speciesNode->node_index;
      for (size_t c = 0; c < _gammaCatNumber; ++c) {
        auto ec = e * _gammaCatNumber + c;
        REAL v = REAL();
        computeProbability(cid, speciesNode, c, v);
        tempUq[ec] = v;
        scale(v);
        sums[c] += v;
      }
    }
    std::swap(tempUq, uq);
    std::fill(correctionSum.begin(), correctionSum.end(), REAL());
    std::fill(clv._survivingTransferSum.begin(),
              clv._survivingTransferSum.end(), REAL());
    // precompute ancestral correction sum
    // (to forbid transfers to parents)
    if (_transferConstraint == TransferConstaint::PARENTS) {
      auto postOrder = this->_speciesTree.getPostOrderNodes();
      for (auto it = postOrder.rbegin(); it != postOrder.rend(); ++it) {
        auto speciesNode = *it;
        auto e = speciesNode->node_index;
        for (size_t c = 0; c < _gammaCatNumber; ++c) {
          auto parent = speciesNode;
          auto ec = e * _gammaCatNumber + c;
          while (parent) {
            auto p = parent->node_index;
            auto pc = p * _gammaCatNumber + c;
            auto temp = uq[pc];
            scale(temp);
            correctionSum[ec] += temp;
            parent = parent->parent;
          }
        }
      }
    }
    // precompute the related correction sum
    // (to forbit transfers to older lineages)
    if (_transferConstraint == TransferConstaint::RELDATED) {
      std::vector<REAL> softDatedSums(N * _gammaCatNumber, REAL());
      std::vector<REAL> softDatedSum(_gammaCatNumber, REAL());
      for (auto leaf : this->_speciesTree.getLeaves()) {
        auto e = leaf->node_index;
        for (size_t c = 0; c < _gammaCatNumber; ++c) {
          auto ec = e * _gammaCatNumber + c;
          softDatedSum[c] += uq[ec] / getTransferWeightNorm();
        }
      }
      for (auto it = _datedTree.getOrderedSpeciations().rbegin();
           it != _datedTree.getOrderedSpeciations().rend(); ++it) {
        auto node = (*it);
        auto e = node->node_index;
        for (size_t c = 0; c < _gammaCatNumber; ++c) {
          auto ec = e * _gammaCatNumber + c;
          softDatedSums[ec] = softDatedSum[c];
          auto temp = uq[ec] / getTransferWeightNorm();
          scale(temp);
          softDatedSum[c] += temp;
        }
      }
      for (auto node : this->getPrunedSpeciesNodes()) {
        auto e = node->node_index;
        auto p = node->parent ? node->parent->node_index : e;
        if (e != p) {
          for (size_t c = 0; c < _gammaCatNumber; ++c) {
            auto pc = p * _gammaCatNumber + c;
            auto ec = e * _gammaCatNumber + c;
            correctionSum[ec] = softDatedSums[pc];
          }
        }
      }
    }
    for (size_t c = 0; c < _gammaCatNumber; ++c) {
      clv._survivingTransferSum[c] = sums[c];
    }
  }
}

template <class REAL> void UndatedDTLMultiModel<REAL>::updateCLVs() {
  for (CID cid = 0; cid < this->_ccp.getCladesNumber(); ++cid) {
    updateCLV(cid);
  }
}

template <class REAL>
double UndatedDTLMultiModel<REAL>::computeLogLikelihood() {
  if (this->_ccp.skip()) {
    return 0.0;
  }
  auto hash = this->getHash();
  auto cacheIt = _llCache.find(hash);
  if (cacheIt != _llCache.end()) {
    return cacheIt->second;
  }

  this->beforeComputeLogLikelihood();
  if (this->_memorySavings) {
    allocateMemory();
  }
  updateCLVs();
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  std::vector<REAL> categoryLikelihoods(_gammaCatNumber, REAL());

  for (auto speciesNode : this->getPrunedSpeciesNodes()) {
    auto e = speciesNode->node_index;
    for (size_t c = 0; c < _gammaCatNumber; ++c) {
      categoryLikelihoods[c] +=
          _dtlclvs[rootCID]._uq[e * _gammaCatNumber + c] * _OP[e];
    }
  }
  // condition on survival
  for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
    categoryLikelihoods[c] /= getLikelihoodFactor(c);
  }
  // sum over the categories
  REAL res = std::accumulate(categoryLikelihoods.begin(),
                             categoryLikelihoods.end(), REAL());
  // normalize by the number of categories
  res /= double(_gammaCatNumber);
  // ths root correction makes sure that UndatedDTLMultiModel and
  // UndatedDTL model are equivalent when there is one tree per
  // family: the UndatedDTLMultiModel integrates over all possible
  // roots and adds a 1/numberOfGeneRoots weight that is not
  // present un the UndatedDTL, so we multiply back here
  auto rootCorrection = double(this->getPrunedSpeciesNodeNumber());
  res *= rootCorrection;
  auto ret = log(res);
  _llCache[hash] = ret;
  if (this->_memorySavings) {
    deallocateMemory();
  }
  return ret;
}

template <class REAL>
bool UndatedDTLMultiModel<REAL>::sampleTransferEvent(
    unsigned int cid, corax_rnode_t *originSpeciesNode, size_t category,
    Scenario::Event &event) {
  auto e = originSpeciesNode->node_index;
  auto N = this->getPrunedSpeciesNodeNumber();
  auto c = category;
  auto ec = e * _gammaCatNumber + c;
  auto &clv = _dtlclvs[cid];
  auto &survivingTransferSum = clv._survivingTransferSum;
  auto &correctionSum = clv._correctionSum;
  REAL max = REAL();
  switch (_transferConstraint) {
  case TransferConstaint::NONE:
    max = survivingTransferSum[c];
    break;
  case TransferConstaint::PARENTS:
    max = survivingTransferSum[c] - correctionSum[ec];
    break;
  case TransferConstaint::RELDATED:
    max = correctionSum[ec] * static_cast<REAL>(N);
    break;
  default:
    assert(false);
  }
  auto samplingProba = Random::getProba();
  max *= samplingProba;
  REAL sum = REAL();

  std::unordered_set<unsigned int> parents;
  if (_transferConstraint == TransferConstaint::PARENTS) {
    auto parent = originSpeciesNode;
    while (parent) {
      parents.insert(parent->node_index);
      parent = parent->parent;
    }
  }

  for (auto speciesNode : this->getPrunedSpeciesNodes()) {
    auto h = speciesNode->node_index;
    if (_transferConstraint == TransferConstaint::NONE) {
    }
    // parent mode: do not continue if the receiving species
    // is a parent of the source species
    if (_transferConstraint == TransferConstaint::PARENTS) {
      if (parents.end() != parents.find(h)) {
        continue;
      }
    }
    if (_transferConstraint == TransferConstaint::RELDATED) {
      if (originSpeciesNode->parent) {
        auto p = originSpeciesNode->parent->node_index;
        if (_datedTree.getRank(p) >= _datedTree.getRank(h)) {
          continue;
        }
      }
    }
    auto hc = h * _gammaCatNumber + c;
    sum += _dtlclvs[cid]._uq[hc];
    if (sum > max) {
      event.pllDestSpeciesNode = speciesNode;
      ;
      event.destSpeciesNode = h;
      return true;
    }
  }
  return false;
}

template <class REAL>
void UndatedDTLMultiModel<REAL>::recomputeSpeciesProbabilities() {
  auto &dupRates = _dtlRates[0];
  auto &lossRates = _dtlRates[1];
  auto &transferRates = _dtlRates[2];
  auto maxSpeciesId = this->_speciesTree.getNodeNumber();
  assert(maxSpeciesId == dupRates.size());
  assert(maxSpeciesId == lossRates.size());
  assert(maxSpeciesId == transferRates.size());
  for (auto speciesNode : this->getPrunedSpeciesNodes()) {
    auto e = speciesNode->node_index;
    for (size_t c = 0; c < _gammaCatNumber; ++c) {
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
        highway.proba = highway.highway.proba / sum;
        assert(highway.proba < 1.0);
      }
    }
  }

  std::vector<corax_rnode_t *> speciesNodesBuffer;
  std::vector<corax_rnode_t *> *possibleSpeciesRootNodes = nullptr;
  switch (_originationStrategy) {
  case OriginationStrategy::UNIFORM:
  case OriginationStrategy::OPTIMIZE:
    possibleSpeciesRootNodes = &(this->getPrunedSpeciesNodes());
    break;
  case OriginationStrategy::ROOT:
    speciesNodesBuffer.push_back(this->_speciesTree.getRoot());
    possibleSpeciesRootNodes = &speciesNodesBuffer;
    break;
  case OriginationStrategy::LCA:
    speciesNodesBuffer.push_back(this->getSpeciesLCA());
    possibleSpeciesRootNodes = &speciesNodesBuffer;
    break;
  }
  if (_originationStrategy == OriginationStrategy::OPTIMIZE) {
    double sum = 0.0;
    for (auto speciesNode : *possibleSpeciesRootNodes) {
      sum += _dtlRates[3][speciesNode->node_index];
    }
    sum /= static_cast<double>(this->getPrunedSpeciesNodes().size());
    for (auto speciesNode : *possibleSpeciesRootNodes) {
      _OP[speciesNode->node_index] =
          _dtlRates[3][speciesNode->node_index] / sum;
    }
  } else {
    std::fill(_OP.begin(), _OP.end(), 0.0);
    for (auto speciesNode : *possibleSpeciesRootNodes) {
      _OP[speciesNode->node_index] = 1.0;
    }
  }
  std::fill(_uE.begin(), _uE.end(), REAL());
  auto transferSum = std::vector<REAL>(_gammaCatNumber, REAL());
  unsigned int maxIt = 4;
  for (unsigned int it = 0; it < maxIt; ++it) {
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      auto e = speciesNode->node_index;
      for (size_t c = 0; c < _gammaCatNumber; ++c) {
        auto ec = e * _gammaCatNumber + c;
        if (it == maxIt - 1 && !speciesNode->left) {
          _uE[ec] = _uE[ec] * (1.0 - this->_fm[e]) + REAL(this->_fm[e]);
          scale(_uE[ec]);
        } else {
          // L
          REAL proba(_PL[ec]);
          REAL temp;
          // D
          temp = _uE[ec] * _uE[ec] * _PD[ec];
          scale(temp);
          proba += temp;
          // T
          temp = _uEBar[ec] * _uE[ec] * _PT[ec];
          scale(temp);
          proba += temp;
          // S
          if (this->getSpeciesLeft(speciesNode)) {
            auto g = this->getSpeciesLeft(speciesNode)->node_index;
            auto h = this->getSpeciesRight(speciesNode)->node_index;
            auto gc = g * _gammaCatNumber + c;
            auto hc = h * _gammaCatNumber + c;
            temp = _uE[gc] * _uE[hc] * _PS[ec];
            scale(temp);
            proba += temp;
          }
          // transfer highway
          for (const auto &highway : _highways[e]) {
            auto d = highway.highway.dest->node_index;
            auto dc = d * _gammaCatNumber + c;
            temp = _uE[ec] * _uE[dc] * highway.proba;
            scale(temp);
            proba += temp;
          }
          _uE[ec] = proba;
          if (!(proba < REAL(1.0))) {
            std::cerr << "hey " << proba << " < " << REAL(1.0) << " is wrong"
                      << std::endl;
          }
          assert(proba < REAL(1.000001));
        }
      }
    }
    // now compute transfer sum for the next iteration
    if (it < maxIt - 1) {
      // precompute ancestral correction sum
      // (to forbid transfers to parents)
      std::fill(transferSum.begin(), transferSum.end(), REAL());
      for (auto speciesNode : this->getPrunedSpeciesNodes()) {
        auto e = speciesNode->node_index;
        for (size_t c = 0; c < _gammaCatNumber; ++c) {
          auto ec = e * _gammaCatNumber + c;
          transferSum[c] += _uE[ec];
        }
      }
      std::fill(_uEBar.begin(), _uEBar.end(), REAL());
      std::vector<double> correctionNorm(this->_speciesTree.getNodeNumber() *
                                             _gammaCatNumber,
                                         getTransferWeightNorm());
      if (_transferConstraint == TransferConstaint::PARENTS) {
        auto postOrder = this->_speciesTree.getPostOrderNodes();
        for (auto it = postOrder.rbegin(); it != postOrder.rend(); ++it) {
          auto speciesNode = *it;
          auto e = speciesNode->node_index;
          for (size_t c = 0; c < _gammaCatNumber; ++c) {
            auto parent = speciesNode;
            auto ec = e * _gammaCatNumber + c;
            while (parent) {
              correctionNorm[ec] -= 1.0;
              auto p = parent->node_index;
              auto pc = p * _gammaCatNumber + c;
              auto temp = _uE[pc];
              scale(temp);
              _uEBar[ec] += temp;
              parent = parent->parent;
            }
          }
        }
      } else {
        assert(false);
      }
      for (auto speciesNode : this->getPrunedSpeciesNodes()) {
        auto e = speciesNode->node_index;
        for (size_t c = 0; c < _gammaCatNumber; ++c) {
          auto ec = e * _gammaCatNumber + c;
          _uEBar[ec] = (transferSum[c] - _uEBar[ec]) / correctionNorm[ec];
        }
      }
    }
  } // iterations to account for TL
}

template <class REAL>
REAL UndatedDTLMultiModel<REAL>::getTransferSum(unsigned int cid,
                                                unsigned int e,
                                                unsigned int c) {
  auto ec = e * _gammaCatNumber + c;
  switch (_transferConstraint) {
  case TransferConstaint::NONE:
    return _dtlclvs[cid]._survivingTransferSum[c] / getTransferWeightNorm();
  case TransferConstaint::PARENTS:
    return (_dtlclvs[cid]._survivingTransferSum[c] -
            _dtlclvs[cid]._correctionSum[ec]) /
           _dtlclvs[cid]._correctionNorm[ec];
  case TransferConstaint::RELDATED:
    return _dtlclvs[cid]._correctionSum[ec];
  default:
    assert(false);
  }
}

template <class REAL>
bool UndatedDTLMultiModel<REAL>::computeProbability(
    CID cid, corax_rnode_t *speciesNode, size_t category, REAL &proba,
    ReconciliationCell<REAL> *recCell) {
  proba = REAL();
  bool isSpeciesLeaf = !this->getSpeciesLeft(speciesNode);
  auto e = speciesNode->node_index;
  auto c = category;
  auto ec = e * _gammaCatNumber + c;
  REAL maxProba = REAL();
  if (recCell) {
    recCell->event.geneNode = cid;
    recCell->event.speciesNode = e;
    recCell->event.type = ReconciliationEventType::EVENT_None;
    maxProba = recCell->maxProba;
  }
  if (this->_ccp.isLeaf(cid) && isSpeciesLeaf) {
    if (this->_geneToSpecies[cid] == e) {
      proba = REAL(_PS[ec]);
      if (recCell) {
        recCell->event.label = this->_ccp.getLeafLabel(cid);
      }
    }
    return true;
  }
  REAL temp;
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

  // iterate over all gene CCPs
  for (const auto &cladeSplit : this->_ccp.getCladeSplits(cid)) {
    auto cidLeft = cladeSplit.left;
    auto cidRight = cladeSplit.right;
    auto freq = cladeSplit.frequency;
    if (not isSpeciesLeaf) {
      // S event;
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
        recCell->blLeft = cladeSplit.blLeft;
        recCell->blRight = cladeSplit.blRight;
        return true;
      }
    }
    // D event
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
    // T event
    temp = getTransferSum(cidLeft, e, c) * (_PT[ec] * freq);
    temp *= _dtlclvs[cidRight]._uq[ec];
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_T;

      if (!sampleTransferEvent(cidLeft, speciesNode, c, recCell->event)) {
        return false;
      }
      recCell->event.destSpeciesNode =
          recCell->event.pllDestSpeciesNode->node_index;
      recCell->event.leftGeneIndex = cidRight;
      recCell->event.rightGeneIndex = cidLeft;
      recCell->blLeft = cladeSplit.blLeft;
      recCell->blRight = cladeSplit.blRight;
      return true;
    }
    temp = getTransferSum(cidRight, e, c) * (_PT[ec] * freq);
    temp *= _dtlclvs[cidLeft]._uq[ec];
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_T;
      if (!sampleTransferEvent(cidRight, speciesNode, c, recCell->event)) {
        return false;
      }
      recCell->event.destSpeciesNode =
          recCell->event.pllDestSpeciesNode->node_index;
      recCell->event.leftGeneIndex = cidLeft;
      recCell->event.rightGeneIndex = cidRight;
      recCell->blLeft = cladeSplit.blLeft;
      recCell->blRight = cladeSplit.blRight;
      return true;
    }

    // highway transfers
    for (const auto &highway : _highways[e]) {
      auto d = highway.highway.dest->node_index;
      auto dc = d * _gammaCatNumber + c;
      temp = (_dtlclvs[cidLeft]._uq[ec] * _dtlclvs[cidRight]._uq[dc]) *
             (highway.proba * freq);
      scale(temp);
      proba += temp;
      if (proba > REAL(1.0)) {
        std::cerr << "error " << _dtlclvs[cidLeft]._uq[ec] << " "
                  << _dtlclvs[cidRight]._uq[dc] << " " << highway.proba << " "
                  << freq << std::endl;
      }
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_T;
        recCell->event.pllDestSpeciesNode = this->_speciesTree.getNode(d);
        recCell->event.destSpeciesNode = d;
        recCell->event.leftGeneIndex = cidLeft;
        recCell->event.rightGeneIndex = cidRight;
        recCell->blLeft = cladeSplit.blLeft;
        recCell->blRight = cladeSplit.blRight;
        return true;
      }
      temp = (_dtlclvs[cidRight]._uq[ec] * _dtlclvs[cidLeft]._uq[dc]) *
             (highway.proba * freq);
      scale(temp);
      proba += temp;
      if (proba > REAL(1.0)) {
        std::cerr << "error " << _dtlclvs[cidRight]._uq[ec] << " "
                  << _dtlclvs[cidLeft]._uq[dc] << " " << highway.proba << " "
                  << freq << std::endl;
      }
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_T;
        recCell->event.pllDestSpeciesNode = this->_speciesTree.getNode(d);
        recCell->event.destSpeciesNode = d;
        recCell->event.leftGeneIndex = cidRight;
        recCell->event.rightGeneIndex = cidLeft;
        recCell->blLeft = cladeSplit.blLeft;
        recCell->blRight = cladeSplit.blRight;
        return true;
      }
    }
  } // end of the iteraiton over the gene CCPs

  if (not isSpeciesLeaf) {
    // SL event
    temp = _dtlclvs[cid]._uq[fc] * (_uE[gc] * _PS[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.destSpeciesNode = f;
      recCell->event.pllDestSpeciesNode = this->getSpeciesLeft(speciesNode);
      recCell->event.pllLostSpeciesNode = this->getSpeciesRight(speciesNode);
      return true;
    }
    temp = _dtlclvs[cid]._uq[gc] * (_uE[fc] * _PS[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.destSpeciesNode = g;
      recCell->event.pllDestSpeciesNode = this->getSpeciesRight(speciesNode);
      recCell->event.pllLostSpeciesNode = this->getSpeciesLeft(speciesNode);
      return true;
    }
  }
  // DL event
  temp = _dtlclvs[cid]._uq[ec] * (2.0 * _uE[ec]) * _PD[ec];
  scale(temp);
  proba += temp;
  if (recCell && proba > maxProba) {
    recCell->event.type = ReconciliationEventType::EVENT_DL;
    return true;
  }
  // TL event
  if (!this->_info.noTL) {
    // the gene is transfered to the dest species and goes extinct in
    // the src species
    temp = getTransferSum(cid, e, c) * (_PT[ec]);
    temp *= _uE[ec];
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_TL;
      if (!sampleTransferEvent(cid, speciesNode, c, recCell->event)) {
        return false;
      }
      recCell->event.destSpeciesNode =
          recCell->event.pllDestSpeciesNode->node_index;
      return true;
    }
    // the gene is transfered to the dest species and goes extinct in
    // the dest species
    temp = _dtlclvs[cid]._uq[ec] * (_PT[ec]);
    temp *= _uEBar[ec];
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      // in fact, nothing happens, we'll have to resample
      recCell->event.type = ReconciliationEventType::EVENT_TL;
      recCell->event.pllDestSpeciesNode = nullptr;
      return true;
    }
    // TL from a highway
    for (const auto &highway : _highways[e]) {
      auto d = highway.highway.dest->node_index;
      auto dc = d * _gammaCatNumber + c;
      // we transfer but the gene gets extinct in the receiving species
      temp = (_dtlclvs[cid]._uq[ec] * _uE[dc]) * highway.proba;
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        // in fact, nothing happens, we'll have to resample
        recCell->event.type = ReconciliationEventType::EVENT_TL;
        recCell->event.pllDestSpeciesNode = nullptr;
        return true;
      }
      // we transfer and the gene gets extinct in the sending species
      temp = (_dtlclvs[cid]._uq[dc] * _uE[ec]) * highway.proba;
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_TL;
        recCell->event.pllDestSpeciesNode = highway.highway.dest;
        recCell->event.destSpeciesNode = d;
        return true;
      }
    }
  }

  if (recCell) {
    std::cerr << "boum " << proba << " " << maxProba << std::endl;
    return false;
  }
  if (proba > REAL(1.0)) {
    return false;
  }
  return true;
}

/**
 *  Correction factor because we condition on survival
 */
template <class REAL>
double UndatedDTLMultiModel<REAL>::getLikelihoodFactor(unsigned int category) {
  double factor(0.0);
  for (auto speciesNode : this->getPrunedSpeciesNodes()) {
    auto e = speciesNode->node_index;
    factor += _OP[e * _gammaCatNumber + category] *
              (1.0 - _uE[e * _gammaCatNumber + category]);
  }
  return factor;
}

template <class REAL>
corax_rnode_t *
UndatedDTLMultiModel<REAL>::sampleSpeciesNode(unsigned int &category) {
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto &uq = _dtlclvs[rootCID]._uq;
  REAL totalLL = REAL();
  for (auto node : this->getPrunedSpeciesNodes()) {
    auto e = node->node_index;
    for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
      totalLL += uq[e * _gammaCatNumber + c] * _OP[e];
    }
  }
  auto toSample = totalLL * Random::getProba();
  auto sumLL = REAL();
  for (auto node : this->getPrunedSpeciesNodes()) {
    auto e = node->node_index;
    for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
      sumLL += uq[e * _gammaCatNumber + c] * _OP[e];
      if (sumLL >= toSample) {
        category = c;
        return node;
      }
    }
  }
  assert(false);
  return nullptr;
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
