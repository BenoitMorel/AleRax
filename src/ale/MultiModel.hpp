#pragma once

#include "util/enums.hpp"
#include <IO/HighwayCandidateParser.hpp>
#include <IO/LibpllParsers.hpp>
#include <ccp/ConditionalClades.hpp>
#include <likelihoods/reconciliation_models/BaseReconciliationModel.hpp>
#include <memory>

class GeneSpeciesMapping;
class PLLRootedTree;

template <class REAL> struct ReconciliationCell {
  Scenario::Event event;
  REAL maxProba;
  double blLeft;
  double blRight;
};

#define EXCLUDE_ABOVE_PRUNED
#define EXCLUDE_DEAD_NODES
#define COMPUTE_TL

// helper function for updateSpeciesToPrunedNode
static void
auxUntilPrunedRoot(corax_rnode_t *speciesNode, corax_rnode_t *prunedRoot,
                   std::vector<corax_rnode_t *> &speciesToPrunedNode) {
  assert(speciesNode);
  assert(prunedRoot);
  speciesToPrunedNode[speciesNode->node_index] = prunedRoot;
  if (speciesNode == prunedRoot) {
    return;
  }
  if (speciesNode->left) {
    assert(speciesNode->right);
    auxUntilPrunedRoot(speciesNode->left, prunedRoot, speciesToPrunedNode);
    auxUntilPrunedRoot(speciesNode->right, prunedRoot, speciesToPrunedNode);
  }
}

// helper function for updateSpeciesToPrunedNode
static void
fillUnsampledSpeciesRec(corax_rnode_t *unsampledNode,
                        corax_rnode_t *prunedNodeToAssign,
                        std::vector<corax_rnode_t *> &speciesToPrunedNode) {
  speciesToPrunedNode[unsampledNode->node_index] = prunedNodeToAssign;
  if (unsampledNode->left) {
    fillUnsampledSpeciesRec(unsampledNode->left, prunedNodeToAssign,
                            speciesToPrunedNode);
    fillUnsampledSpeciesRec(unsampledNode->right, prunedNodeToAssign,
                            speciesToPrunedNode);
  }
}

/**
 * Base class for reconciliation models that take as input
 * gene tree distributions (represented by conditional clade
 * probabilities)
 */
class MultiModel : public BaseReconciliationModel {
public:
  MultiModel(PLLRootedTree &speciesTree,
             const GeneSpeciesMapping &geneSpeciesMapping,
             const RecModelInfo &info, const std::string &ccpFile)
      : BaseReconciliationModel(speciesTree, geneSpeciesMapping, info),
        _memorySavings(info.memorySavings) {
    _ccp.unserialize(ccpFile);
    mapGenesToSpecies();
  }
  virtual ~MultiModel() {}

  const ConditionalClades &getCCP() const { return _ccp; }

  virtual void setHighways(const std::vector<Highway> &highways) {
    (void)(highways);
  }

  virtual void onSpeciesDatesChange() {}

  virtual void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) {
    BaseReconciliationModel::onSpeciesTreeChange(nodesToInvalidate);
    updateSpeciesToPrunedNode();
  }

  virtual bool
  sampleReconciliations(unsigned int samples,
                        std::vector<std::shared_ptr<Scenario>> &scenarios) = 0;

protected:
  // methods used to handle memory in memory saving mode
  virtual void allocateMemory() {}
  virtual void deallocateMemory() {}
  virtual void updateCLVs() {}
  corax_rnode_t *getSpeciesLCA() {
    auto tree = &this->_speciesTree;
    corax_rnode_t *lca = nullptr;
    for (auto it : this->_speciesNameToId) {
      auto spid = it.second;
      auto node = tree->getNode(spid);
      lca = tree->getLCA(lca, node);
    }
    return _speciesToPrunedNode[lca->node_index];
  }

  virtual void updateSpeciesToPrunedNode() {
    if (!_speciesToPrunedNode.size()) {
      _speciesToPrunedNode.resize(this->getAllSpeciesNodeNumber());
    }
    std::fill(_speciesToPrunedNode.begin(), _speciesToPrunedNode.end(),
              nullptr);
    for (auto speciesNode : this->getAllSpeciesNodes()) {
      auto e = speciesNode->node_index;
      if (speciesNode->left) {
        auto left = speciesNode->left->node_index;
        auto right = speciesNode->right->node_index;
        if (_speciesToPrunedNode[left] && _speciesToPrunedNode[right]) {
          // this node belongs to the pruned nodes
          _speciesToPrunedNode[e] = speciesNode;
        } else if (_speciesToPrunedNode[left]) {
          _speciesToPrunedNode[e] = _speciesToPrunedNode[left];
#ifndef EXCLUDE_DEAD_NODES
          fillUnsampledSpeciesRec(speciesNode->right, _speciesToPrunedNode[e],
                                  _speciesToPrunedNode);
#endif
        } else if (_speciesToPrunedNode[right]) {
          _speciesToPrunedNode[e] = _speciesToPrunedNode[right];
#ifndef EXCLUDE_DEAD_NODES
          fillUnsampledSpeciesRec(speciesNode->left, _speciesToPrunedNode[e],
                                  _speciesToPrunedNode);
#endif
        } // else do nothing
      } else {
        if (this->_speciesCoverage[e]) {
          _speciesToPrunedNode[e] = speciesNode;
        }
      }
    }
    // if the  root of the pruned species tree is not the root, we need
    // to map all parents and siblings of the pruned root to the pruned root
#ifndef EXCLUDE_ABOVE_PRUNED
    auxUntilPrunedRoot(this->getSpeciesTree().getRoot(), this->getPrunedRoot(),
                       _speciesToPrunedNode);
#endif
  }

  void mapGenesToSpecies() {
    const auto &cidToLeaves = _ccp.getCidToLeaves();
    this->_speciesNameToId.clear();
    this->_geneToSpecies.resize(_ccp.getCladesNumber());
    for (auto node : this->_allSpeciesNodes) {
      if (!node->left) {
        this->_speciesNameToId[node->label] = node->node_index;
      }
    }
    this->_speciesCoverage =
        std::vector<unsigned int>(this->getAllSpeciesNodeNumber(), 0);
    for (auto p : cidToLeaves) {
      auto cid = p.first;
      const auto &geneName = cidToLeaves.at(cid);
      const auto &speciesName = this->_geneNameToSpeciesName[geneName];
      this->_geneToSpecies[cid] = this->_speciesNameToId[speciesName];
      this->_speciesCoverage[this->_geneToSpecies[cid]]++;
    }
  }
  ConditionalClades _ccp;
  bool _memorySavings;
  std::vector<corax_rnode_t *> _speciesToPrunedNode;
};

/**
 * Implements all basic methods required by the child classes
 * and that require the REAL template
 */
template <class REAL> class MultiModelTemplate : public MultiModel {
public:
  MultiModelTemplate(PLLRootedTree &speciesTree,
                     const GeneSpeciesMapping &geneSpeciesMapping,
                     const RecModelInfo &info, const std::string &ccpFile)
      : MultiModel(speciesTree, geneSpeciesMapping, info, ccpFile) {}
  virtual ~MultiModelTemplate() {}
  /**
   * Samples a reconciled gene tree and stores it into
   * scenario
   */
  virtual bool inferMLScenario(Scenario &scenario);

  /**
   *  Sample scenarios and add them to the scenarios vector
   */
  virtual bool
  sampleReconciliations(unsigned int samples,
                        std::vector<std::shared_ptr<Scenario>> &scenarios);

protected:
  /**
   * Sample the species node from with the sampled gene tree
   * originates. The implementation depends on the origination
   * mode (unifor, from the root etc.)
   */
  virtual corax_rnode_t *sampleSpeciesNode(unsigned int &category) = 0;

private:
  /**
   * Compute the CLV for a given cid (clade id) and a
   * given species node. Returns true in case of sucess
   * If recCell is set, the function samples the next
   * event in the reconciliation space
   */
  virtual bool
  computeProbability(CID cid, corax_rnode_t *speciesNode, size_t category,
                     REAL &proba,
                     ReconciliationCell<REAL> *recCell = nullptr) = 0;

  /**
   * Recursively sample a reconciled gene tree
   */
  bool backtrace(unsigned int cid, corax_rnode_t *speciesRoot,
                 corax_unode_t *geneNode, unsigned int category,
                 Scenario &scenario, bool stochastic);
  bool _computeScenario(Scenario &scenario, bool stochastic);
};

template <class REAL>
bool MultiModelTemplate<REAL>::sampleReconciliations(
    unsigned int samples, std::vector<std::shared_ptr<Scenario>> &scenarios) {
  allocateMemory();
  updateCLVs();
  for (unsigned int i = 0; i < samples; ++i) {
    scenarios.push_back(std::make_shared<Scenario>());
    if (!_computeScenario(*scenarios.back(), true)) {
      deallocateMemory();
      return false;
    }
  }
  deallocateMemory();
  return true;
}

template <class REAL>
bool MultiModelTemplate<REAL>::inferMLScenario(Scenario &scenario) {
  assert(false); // not implemented for MultiModel
  return _computeScenario(scenario, false);
}

template <class REAL>
bool MultiModelTemplate<REAL>::_computeScenario(Scenario &scenario,
                                                bool stochastic) {
  assert(stochastic);
  unsigned int category = 0;
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto speciesRoot = sampleSpeciesNode(category);
  auto rootIndex = 2 * _ccp.getLeafNumber();
  scenario.setVirtualRootIndex(rootIndex);
  scenario.setSpeciesTree(&this->_speciesTree);
  auto virtualGeneRoot = scenario.generateVirtualGeneRoot();
  scenario.setGeneRoot(virtualGeneRoot);
  auto res = backtrace(rootCID, speciesRoot, virtualGeneRoot, category,
                       scenario, stochastic);
  return res;
}

template <class REAL>
bool MultiModelTemplate<REAL>::backtrace(unsigned int cid,
                                         corax_rnode_t *speciesNode,
                                         corax_unode_t *geneNode,
                                         unsigned int category,
                                         Scenario &scenario, bool stochastic) {
  REAL proba;
  if (!computeProbability(cid, speciesNode, category, proba)) {
    return false;
  }
  ReconciliationCell<REAL> recCell;
  recCell.maxProba = proba * Random::getProba();
  if (!computeProbability(cid, speciesNode, category, proba, &recCell)) {
    return false;
  }
  if (scenario.getGeneNodeBuffer().size() == 1) {
    recCell.event.geneNode = scenario.getVirtualRootIndex();
  } else {
    recCell.event.geneNode = geneNode->node_index;
  }

  corax_unode_t *leftGeneNode;
  corax_unode_t *rightGeneNode;
  auto leftCid = recCell.event.leftGeneIndex;
  auto rightCid = recCell.event.rightGeneIndex;
  if (recCell.event.type == ReconciliationEventType::EVENT_S ||
      recCell.event.type == ReconciliationEventType::EVENT_D ||
      recCell.event.type == ReconciliationEventType::EVENT_T) {
    scenario.generateGeneChildren(geneNode, leftGeneNode, rightGeneNode);
    recCell.event.leftGeneIndex = leftGeneNode->node_index;
    recCell.event.rightGeneIndex = rightGeneNode->node_index;
    leftGeneNode->length = recCell.blLeft;
    rightGeneNode->length = recCell.blRight;
    recCell.event.previous_event_type = scenario.getLastEventType();
  }
  bool addEvent = true;
  if ((recCell.event.type == ReconciliationEventType::EVENT_TL &&
       recCell.event.pllDestSpeciesNode == nullptr) ||
      recCell.event.type == ReconciliationEventType::EVENT_DL) {
    addEvent = false;
  }

  if (addEvent) {
    scenario.addEvent(recCell.event);
  }

  bool ok = true;
  std::string label;
  auto c = category;

  scenario.setLastEventType(recCell.event.type);
  switch (recCell.event.type) {
  case ReconciliationEventType::EVENT_S:
    ok &= backtrace(leftCid, this->getSpeciesLeft(speciesNode), leftGeneNode, c,
                    scenario, stochastic);
    scenario.setLastEventType(recCell.event.type);
    ok &= backtrace(rightCid, this->getSpeciesRight(speciesNode), rightGeneNode,
                    c, scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_D:
    ok &=
        backtrace(leftCid, speciesNode, leftGeneNode, c, scenario, stochastic);
    scenario.setLastEventType(recCell.event.type);
    ok &= backtrace(rightCid, speciesNode, rightGeneNode, c, scenario,
                    stochastic);
    break;
  case ReconciliationEventType::EVENT_DL:
    scenario.setLastEventType(ReconciliationEventType::EVENT_S);
    ok &= backtrace(cid, speciesNode, geneNode, c, scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_SL:
    ok &= backtrace(cid, recCell.event.pllDestSpeciesNode, geneNode, c,
                    scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_T:
    // source species
    scenario.setLastEventType(ReconciliationEventType::EVENT_S);
    ok &=
        backtrace(leftCid, speciesNode, leftGeneNode, c, scenario, stochastic);
    // dest species
    scenario.setLastEventType(ReconciliationEventType::EVENT_T);
    ok &= backtrace(rightCid, recCell.event.pllDestSpeciesNode, rightGeneNode,
                    c, scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_TL:
    if (recCell.event.pllDestSpeciesNode == nullptr) {
      // the gene was lost in the recieving species, we resample again
      scenario.setLastEventType(ReconciliationEventType::EVENT_S);
      ok &= backtrace(cid, speciesNode, geneNode, c, scenario, stochastic);
    } else {
      scenario.setLastEventType(ReconciliationEventType::EVENT_TL);
      ok &= backtrace(cid, recCell.event.pllDestSpeciesNode, geneNode, c,
                      scenario, stochastic);
    }
    break;
  case ReconciliationEventType::EVENT_None:
    label = _ccp.getCidToLeaves().at(cid);
    geneNode->label = new char[label.size() + 1];
    memcpy(geneNode->label, label.c_str(), label.size() + 1);
    break;
  default:
    ok = false;
  }
  return ok;
}
