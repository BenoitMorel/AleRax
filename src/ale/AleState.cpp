#include "AleState.hpp"

#include <unordered_set>

#include <IO/FileSystem.hpp>
#include <IO/ParallelOfstream.hpp>
#include <parallelization/ParallelContext.hpp>

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
  std::string bufferStr;
  unsigned int checkpointRanks = 0;
  std::string checkpointCmd;
  is >> bufferStr;
  is >> checkpointRanks;
  is >> bufferStr;
  is >> std::ws;
  std::getline(is, checkpointCmd);
  is.close();
  if (checkpointRanks != currentRanks) {
    Logger::info << "\nThe number of MPI ranks used to run the program is "
                    "different from the checkpoint."
                 << std::endl;
    Logger::info
        << "To run anew please rename or remove or change the output directory."
        << std::endl;
    Logger::info << "To restart from the checkpoint please run the program "
                    "exactly as given in the "
                 << cmdPath << " file\n"
                 << std::endl;
    ParallelContext::abort(0);
  }
  if (checkpointCmd != currentCmd) {
    Logger::info << "\nThe command used to run the program is different from "
                    "the checkpoint."
                 << std::endl;
    Logger::info
        << "To run anew please rename or remove or change the output directory."
        << std::endl;
    Logger::info << "To restart from the checkpoint please run the program "
                    "exactly as given in the "
                 << cmdPath << " file\n"
                 << std::endl;
    ParallelContext::abort(0);
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
    if (checkpointNames.find(family.name) == checkpointNames.end()) {
      continue;
    }
    validFamilies.push_back(family);
  }
  families = validFamilies;
  ParallelContext::barrier();
}

void AleState::serialize(const std::string &checkpointDir) const {
  ParallelContext::barrier();
  auto checkpointPath =
      FileSystem::joinPaths(checkpointDir, "mainCheckpoint.txt");
  ParallelOfstream os(checkpointPath, true);
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
}

void AleState::unserialize(const std::string &checkpointDir,
                           const std::vector<std::string> &perCoreFamilyNames) {
  ParallelContext::barrier();
  localFamilyNames = perCoreFamilyNames;
  auto checkpointPath =
      FileSystem::joinPaths(checkpointDir, "mainCheckpoint.txt");
  std::ifstream is(checkpointPath);
  if (!is.good()) {
    Logger::error << "Error: cannot read checkpoint file " << checkpointPath
                  << std::endl;
    ParallelContext::abort(32);
  }
  // current step
  unsigned int bufferUint = 0;
  is >> bufferUint;
  currentStep = static_cast<AleStep>(bufferUint);
  if (currentStep == AleStep::End) {
    Logger::info << "\nThe previous run finished successfully according to the "
                    "checkpoint."
                 << std::endl;
    Logger::info << "To run anew please rename or remove or change the output "
                    "directory\n"
                 << std::endl;
    ParallelContext::abort(0);
  }
  // dated species tree
  std::string bufferStr;
  is >> bufferStr;
  speciesTree = std::make_unique<SpeciesTree>(bufferStr, false);
  // mixture alpha
  is >> mixtureAlpha;
  // transfer highways
  auto labelToNode = speciesTree->getTree().getLabelToNode(false);
  is >> std::ws;
  std::string line;
  while (std::getline(is, line)) {
    if (line.size()) {
      std::string src;
      std::string dest;
      double proba;
      std::stringstream iss(line);
      iss >> src;
      iss >> dest;
      iss >> proba;
      Highway highway(labelToNode.find(src)->second,
                      labelToNode.find(dest)->second);
      highway.proba = proba;
      transferHighways.push_back(highway);
    }
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
  assert(perLocalFamilyModelParams.size() == perCoreFamilyNames.size());
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
  unsigned int bufferUint = 0;
  is >> bufferUint;
  std::string speciesTreeNewick;
  is >> speciesTreeNewick;
  is.close();
  return speciesTreeNewick;
}
