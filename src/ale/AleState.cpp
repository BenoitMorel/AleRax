#include "AleState.hpp"

#include <cassert>
#include <iostream>
#include <sstream>

#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <IO/ParallelOfstream.hpp>
#include <parallelization/ParallelContext.hpp>

const std::string CKP_INCOMPLETE = std::string(10, '*');
const std::string CKP_COMPLETE = std::string(10, '#');

static void finishMessageAndExit() {
  Logger::info << "To run anew, please rename or remove or change "
               << "the output directory\n"
               << std::endl;
  ParallelContext::abort(0);
}

bool AleState::checkpointExists(const std::string &checkpointDir) {
  if (FileSystem::dirExists(checkpointDir)) {
    Logger::info << "Checkpoint detected" << std::endl;
    // The checkpoint directory exists. Now check that the header of
    // the mainCheckpoint.txt file is complete, which would imply
    // that all checkpoint files exist and are updated, since
    // - args.txt and fams.txt are written before mainCheckpoint.txt
    // - the header is complete only after all family files have been updated
    // Note that the current algorithm:
    // - does control the order of writing files to the OS buffer (i.e.
    //   safe against software interruptions: process kill)
    // - does NOT control the order of writing files to the disk (i.e.
    //   NOT safe against hardware interruptions: power loss or kernel crash)
    auto checkpointPath =
        FileSystem::joinPaths(checkpointDir, "mainCheckpoint.txt");
    std::string headerStr;
    std::ifstream is(checkpointPath);
    if (!is || !std::getline(is, headerStr) || headerStr != CKP_COMPLETE) {
      Logger::info << "\nThe previous run failed to properly save "
                   << "the checkpoint." << std::endl;
      finishMessageAndExit();
    }
    is.close();
    return true;
  }
  return false;
}

void AleState::writeCheckpointCmd(const std::string &currentCmd,
                                  const std::string &checkpointDir) {
  auto cmdPath = FileSystem::joinPaths(checkpointDir, "args.txt");
  ParallelOfstream os(cmdPath, true);
  os << "mpiRanks:" << std::endl;
  os << ParallelContext::getSize() << std::endl;
  os << "command:" << std::endl;
  os << currentCmd << std::endl;
  os.close();
  ParallelContext::barrier();
}

void AleState::checkCheckpointCmd(const std::string &currentCmd,
                                  const std::string &checkpointDir) {
  auto cmdPath = FileSystem::joinPaths(checkpointDir, "args.txt");
  std::ifstream is(cmdPath);
  if (!is.good()) {
    Logger::error << "Error: cannot read cmd checkpoint file " << cmdPath
                  << std::endl;
    ParallelContext::abort(32);
  }
  unsigned int currentRanks = ParallelContext::getSize();
  unsigned int checkpointRanks = 0;
  std::string checkpointCmd;
  std::string bufferStr;
  is >> bufferStr;
  is >> checkpointRanks;
  is >> bufferStr;
  is >> std::ws;
  std::getline(is, checkpointCmd);
  is.close();
  std::string errorStr;
  if (checkpointRanks != currentRanks) {
    errorStr = "number of MPI ranks";
  } else if (checkpointCmd != currentCmd) {
    errorStr = "command";
  }
  if (errorStr.size()) {
    Logger::info << "\nThe " << errorStr << " used to run the program is "
                 << "different from the checkpoint." << std::endl;
    Logger::info << "To restart from the checkpoint please run the program "
                 << "exactly as given in the " << cmdPath << " file."
                 << std::endl;
    finishMessageAndExit();
  }
}

void AleState::writeCheckpointFamilies(const Families &families,
                                       const std::string &checkpointDir) {
  auto famPath = FileSystem::joinPaths(checkpointDir, "fams.txt");
  ParallelOfstream os(famPath, true);
  for (const auto &family : families) {
    os << family.name << std::endl;
  }
  os.close();
  ParallelContext::barrier();
}

void AleState::filterCheckpointFamilies(Families &families,
                                        const std::string &checkpointDir) {
  Families validFamilies;
  std::unordered_set<std::string> checkpointNames;
  auto famPath = FileSystem::joinPaths(checkpointDir, "fams.txt");
  std::ifstream is(famPath);
  if (!is.good()) {
    Logger::error << "Error: cannot read family checkpoint file " << famPath
                  << std::endl;
    ParallelContext::abort(32);
  }
  std::string name;
  while (std::getline(is, name)) {
    if (name.size()) {
      checkpointNames.insert(name);
    }
  }
  is.close();
  for (const auto &family : families) {
    if (checkpointNames.find(family.name) != checkpointNames.end()) {
      validFamilies.push_back(family);
    }
  }
  families = validFamilies;
  ParallelContext::barrier();
}

void AleState::serialize(const std::string &checkpointDir) const {
  assert(currentStep > AleStep::Init);
  assert(!localFamilyNames.empty());
  ParallelContext::barrier();
  auto checkpointPath =
      FileSystem::joinPaths(checkpointDir, "mainCheckpoint.txt");
  ParallelOfstream os(checkpointPath, true);
  os << CKP_INCOMPLETE << std::endl;
  // current step
  os << static_cast<unsigned int>(currentStep) << std::endl;
  // dated species tree
  os << speciesTree->toString() << std::endl;
  // mixture alpha
  os << mixtureAlpha << std::endl;
  // transfer highways
  for (const auto &highway : transferHighways) {
    os << highway.src->label << " " << highway.dest->label << " "
       << highway.proba << std::endl;
  }
  os.close();
  ParallelContext::barrier();
  // get the new indexing scheme of the unserialized species tree
  PLLRootedTree referenceTree(speciesTree->toString(), false);
  auto spidNewToOld = referenceTree.getNodeIndexMapping(speciesTree->getTree());
  // param vectors
  assert(perLocalFamilyModelParams.size() == localFamilyNames.size());
  for (unsigned int f = 0; f < localFamilyNames.size(); ++f) {
    const auto &family = localFamilyNames[f];
    const auto &mp = perLocalFamilyModelParams[f];
    auto paramTypeNumber = mp.getParamTypeNumber();
    auto speciesBranchNumber = mp.getSpeciesBranchNumber();
    auto paramPath = FileSystem::joinPaths(checkpointDir, family + ".txt");
    std::ofstream os(paramPath);
    os << paramTypeNumber << " " << speciesBranchNumber << std::endl;
    for (unsigned int e = 0; e < speciesBranchNumber; ++e) {
      for (unsigned int r = 0; r < paramTypeNumber; ++r) {
        os << mp.getParameter(spidNewToOld[e], r) << " ";
      }
    }
    os << std::endl;
    os.close();
  }
  ParallelContext::barrier();
  // reset the header (to ensure all checkpoint files have been updated)
  if (ParallelContext::getRank() == 0) {
    std::string headerStr;
    std::fstream fs(checkpointPath, std::ios::in | std::ios::out);
    assert(fs && std::getline(fs, headerStr) && headerStr == CKP_INCOMPLETE);
    fs.seekp(0, std::ios::beg);
    fs << CKP_COMPLETE << std::endl;
    fs.close();
  }
  ParallelContext::barrier();
}

void AleState::unserialize(const std::string &checkpointDir) {
  assert(currentStep == AleStep::Init);
  assert(!localFamilyNames.empty());
  ParallelContext::barrier();
  auto checkpointPath =
      FileSystem::joinPaths(checkpointDir, "mainCheckpoint.txt");
  std::ifstream is(checkpointPath);
  if (!is.good()) {
    Logger::error << "Error: cannot read checkpoint file " << checkpointPath
                  << std::endl;
    ParallelContext::abort(32);
  }
  std::string headerStr;
  std::getline(is, headerStr);
  // current step
  unsigned int bufferUint = 0;
  is >> bufferUint;
  currentStep = static_cast<AleStep>(bufferUint);
  assert(currentStep > AleStep::Init);
  if (currentStep == AleStep::End) {
    Logger::info << "\nThe previous run finished successfully according to "
                 << "the checkpoint." << std::endl;
    finishMessageAndExit();
  }
  // dated species tree (but we already know it)
  std::string bufferStr;
  is >> bufferStr;
  // mixture alpha
  is >> mixtureAlpha;
  // transfer highways
  auto labelToNode = speciesTree->getTree().getLabelToNode(false);
  is >> std::ws;
  std::string line;
  while (std::getline(is, line)) {
    std::string src;
    std::string dest;
    double proba = -1.0;
    std::istringstream iss(line);
    iss >> src;
    iss >> dest;
    iss >> proba;
    assert(proba >= 0.0);
    Highway highway(labelToNode.find(src)->second,
                    labelToNode.find(dest)->second);
    highway.proba = proba;
    transferHighways.push_back(highway);
  }
  is.close();
  // param vectors
  for (const auto &family : localFamilyNames) {
    auto paramPath = FileSystem::joinPaths(checkpointDir, family + ".txt");
    std::ifstream is(paramPath);
    if (!is.good()) {
      Logger::error << "Error: cannot read parameter checkpoint file "
                    << paramPath << std::endl;
      ParallelContext::abort(32);
    }
    unsigned int paramTypeNumber = 0;
    unsigned int speciesBranchNumber = 0;
    is >> paramTypeNumber;
    is >> speciesBranchNumber;
    auto mp = AleModelParameters(paramTypeNumber, speciesBranchNumber, 0.0);
    for (unsigned int i = 0; i < mp.getParameters().dimensions(); ++i) {
      is >> mp.getParameters()[i];
    }
    is.close();
    perLocalFamilyModelParams.push_back(mp);
  }
  ParallelContext::barrier();
}

std::string
AleState::readCheckpointSpeciesTree(const std::string &checkpointDir) {
  auto checkpointPath =
      FileSystem::joinPaths(checkpointDir, "mainCheckpoint.txt");
  std::ifstream is(checkpointPath);
  if (!is.good()) {
    Logger::error << "Error: cannot read checkpoint file " << checkpointPath
                  << std::endl;
    ParallelContext::abort(32);
  }
  std::string headerStr;
  std::getline(is, headerStr);
  unsigned int bufferUint = 0;
  is >> bufferUint;
  std::string speciesTreeNewick;
  is >> speciesTreeNewick;
  is.close();
  return speciesTreeNewick;
}
