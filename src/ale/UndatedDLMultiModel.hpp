#pragma once

#include "MultiModel.hpp"
#include <IO/GeneSpeciesMapping.hpp>
#include <ccp/ConditionalClades.hpp>
#include <maths/ScaledValue.hpp>
#include <trees/PLLRootedTree.hpp>

class RecModelInfo;
double log(ScaledValue v);

template <class REAL>
class UndatedDLMultiModel : public MultiModelTemplate<REAL> {
public:
  UndatedDLMultiModel(PLLRootedTree &speciesTree,
                      const GeneSpeciesMapping &geneSpeciesMapping,
                      const RecModelInfo &info, const std::string &ccpFile);

  virtual ~UndatedDLMultiModel() {}

  virtual void setRates(const RatesVector &);
  virtual double computeLogLikelihood();
  virtual corax_rnode_t *sampleSpeciesNode(unsigned int &category);

private:
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _OP; // Origination probability, per species branch
  std::vector<double> _uE; // Extinction probability, per species branch
  OriginationStrategy _originationStrategy;
  using DLCLV = std::vector<REAL>;
  std::vector<DLCLV> _dlclvs;
  RatesVector _dlRates;

  REAL getLikelihoodFactor();
  virtual void recomputeSpeciesProbabilities();
  virtual bool computeProbability(CID cid, corax_rnode_t *speciesNode,
                                  size_t category, REAL &proba,
                                  ReconciliationCell<REAL> *recCell = nullptr);
};

template <class REAL>
UndatedDLMultiModel<REAL>::UndatedDLMultiModel(
    PLLRootedTree &speciesTree, const GeneSpeciesMapping &geneSpeciesMapping,
    const RecModelInfo &info, const std::string &ccpFile)
    : MultiModelTemplate<REAL>(speciesTree, geneSpeciesMapping, info, ccpFile),
      _PD(this->getAllSpeciesNodeNumber(), 0.2),
      _PL(this->getAllSpeciesNodeNumber(), 0.2),
      _PS(this->getAllSpeciesNodeNumber(), 1.0),
      _OP(this->getAllSpeciesNodeNumber(),
          1.0 / static_cast<double>(this->getAllSpeciesNodeNumber())),
      _uE(this->getAllSpeciesNodeNumber(), 0.0),
      _originationStrategy(info.originationStrategy) {
  std::vector<REAL> zeros(this->getAllSpeciesNodeNumber(), REAL());
  _dlclvs = std::vector<std::vector<REAL>>(this->_ccp.getCladesNumber(), zeros);
  for (unsigned int e = 0; e < this->_speciesTree.getNodeNumber(); ++e) {
    double sum = _PD[e] + _PL[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PS[e] /= sum;
  }
  this->onSpeciesTreeChange(nullptr);
}

template <class REAL>
void UndatedDLMultiModel<REAL>::setRates(const RatesVector &rates) {
  assert(rates.size() == this->_info.modelFreeParameters());
  _dlRates = rates;
  recomputeSpeciesProbabilities();
}

template <class REAL> double UndatedDLMultiModel<REAL>::computeLogLikelihood() {
  if (this->_ccp.skip()) {
    return 0.0;
  }
  this->beforeComputeLogLikelihood();
  std::vector<REAL> zeros(this->_speciesTree.getNodeNumber(), REAL());
  _dlclvs = std::vector<std::vector<REAL>>(this->_ccp.getCladesNumber(), zeros);
  for (CID cid = 0; cid < this->_ccp.getCladesNumber(); ++cid) {
    for (auto speciesNode : this->getPrunedSpeciesNodes()) {
      auto category = 0;
      computeProbability(cid, speciesNode, category,
                         _dlclvs[cid][speciesNode->node_index]);
    }
  }
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  REAL res = REAL();
  for (auto speciesNode : this->getPrunedSpeciesNodes()) {
    auto e = speciesNode->node_index;
    res += _dlclvs[rootCID][e] * _OP[e];
  }
  // the root correction makes sure that UndatedDLMultiModel and
  // UndatedDL model are equivalent when there is one tree per
  // family: the UndatedDLMultiModel integrates over all possible
  // roots and adds a 1/numberOfGeneRoots weight that is not
  // present un the UndatedDL, so we multiply back here
  auto rootCorrection = double(this->_ccp.getLeafNumber() * 2 - 3);
  return log(res) - log(getLikelihoodFactor()) + log(rootCorrection);
}

template <class REAL>
void UndatedDLMultiModel<REAL>::recomputeSpeciesProbabilities() {
  auto &dupRates = _dlRates[0];
  auto &lossRates = _dlRates[1];
  assert(this->getAllSpeciesNodeNumber() == dupRates.size());
  assert(this->getAllSpeciesNodeNumber() == lossRates.size());
  _PD = dupRates;
  _PL = lossRates;
  _PS.resize(this->getAllSpeciesNodeNumber());
  for (auto node : this->getPrunedSpeciesNodes()) {
    auto e = node->node_index;
    if (this->_info.noDup) {
      _PD[e] = 0.0;
    }
    auto sum = _PD[e] + _PL[e] + 1.0;
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PS[e] = 1.0 / sum;
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
      sum += _dlRates[2][speciesNode->node_index];
    }
    sum /= static_cast<double>(this->getPrunedSpeciesNodes().size());
    for (auto speciesNode : *possibleSpeciesRootNodes) {
      _OP[speciesNode->node_index] = _dlRates[2][speciesNode->node_index] / sum;
    }
  } else {
    std::fill(_OP.begin(), _OP.end(), 0.0);
    for (auto speciesNode : *possibleSpeciesRootNodes) {
      _OP[speciesNode->node_index] = 1.0;
    }
  }
  for (auto speciesNode : this->getPrunedSpeciesNodes()) {
    auto e = speciesNode->node_index;
    double a = _PD[e];
    double b = -1.0;
    double c = _PL[e];
    if (this->getSpeciesLeft(speciesNode)) {
      c += _PS[e] * _uE[this->getSpeciesLeft(speciesNode)->node_index] *
           _uE[this->getSpeciesRight(speciesNode)->node_index];
    }
    double proba = solveSecondDegreePolynome(a, b, c);
    assert(proba >= 0.0 && proba <= 1.0);
    _uE[speciesNode->node_index] = proba;
  }
}

template <class REAL>
bool UndatedDLMultiModel<REAL>::computeProbability(
    CID cid, corax_rnode_t *speciesNode, size_t, REAL &proba,
    ReconciliationCell<REAL> *recCell) {
  proba = REAL();
  bool isSpeciesLeaf = !this->getSpeciesLeft(speciesNode);
  auto e = speciesNode->node_index;
  REAL maxProba = REAL();
  if (recCell) {
    recCell->event.geneNode = cid;
    recCell->event.speciesNode = e;
    recCell->event.type = ReconciliationEventType::EVENT_None;
    maxProba = recCell->maxProba;
  }
  // terminal gene and species nodes
  if (this->_ccp.isLeaf(cid) && isSpeciesLeaf) {
    if (this->_geneToSpecies[cid] == e) {
      proba = REAL(_PS[e]);
      if (recCell) {
        recCell->event.label = this->_ccp.getLeafLabel(cid);
      }
    }
    return true;
  }
  REAL temp;
  unsigned int f = 0;
  unsigned int g = 0;
  if (!isSpeciesLeaf) {
    f = this->getSpeciesLeft(speciesNode)->node_index;
    g = this->getSpeciesRight(speciesNode)->node_index;
  }

  for (const auto &cladeSplit : this->_ccp.getCladeSplits(cid)) {
    auto cidLeft = cladeSplit.left;
    auto cidRight = cladeSplit.right;
    auto freq = cladeSplit.frequency;
    if (not isSpeciesLeaf) {
      // S event;
      temp = _dlclvs[cidLeft][f] * _dlclvs[cidRight][g] * (_PS[e] * freq);
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
      temp = _dlclvs[cidRight][f] * _dlclvs[cidLeft][g] * (_PS[e] * freq);
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
    // D events
    temp = _dlclvs[cidLeft][e] * _dlclvs[cidRight][e] * (_PD[e] * freq);
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
  if (not isSpeciesLeaf) {
    // SL event
    temp = _dlclvs[cid][f] * (_uE[g] * _PS[e]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.destSpeciesNode = f;
      recCell->event.pllDestSpeciesNode = this->getSpeciesLeft(speciesNode);
      recCell->event.pllLostSpeciesNode = this->getSpeciesRight(speciesNode);
      return true;
    }

    temp = _dlclvs[cid][g] * (_uE[f] * _PS[e]);
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
  // proba /= (1.0 - 2.0 * _PD[e] * _uE[e]);
  if (recCell) {
    // std::cerr << "cerr " << proba << " " << maxProba << " " << (proba >
    // maxProba) << std::endl;
    return false; // we haven't sampled any event
  }
  return true;
}

template <class REAL> REAL UndatedDLMultiModel<REAL>::getLikelihoodFactor() {
  REAL factor(0.0);
  for (auto speciesNode : this->getPrunedSpeciesNodes()) {
    auto e = speciesNode->node_index;
    factor += (REAL(1.0) - REAL(_uE[e]));
  }
  return factor;
}

template <class REAL>
corax_rnode_t *
UndatedDLMultiModel<REAL>::sampleSpeciesNode(unsigned int &category) {
  category = 0;
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto &uq = _dlclvs[rootCID];
  auto totalLL = std::accumulate(uq.begin(), uq.end(), REAL());
  REAL toSample = totalLL * Random::getProba();
  REAL sumLL = REAL();
  for (auto speciesNode : this->getPrunedSpeciesNodes()) {
    sumLL += uq[speciesNode->node_index];
    if (sumLL > toSample) {
      return speciesNode;
    }
  }
  assert(false);
  return nullptr;
}
