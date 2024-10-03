#include "AleState.hpp"

#include <IO/FileSystem.hpp>
#include <fstream>

void AleState::serialize(const std::string &outputDir) const {
  ParallelContext::barrier();
  if (ParallelContext::getRank() == 0) {
    std::string checkpointPath =
        FileSystem::joinPaths(outputDir, "mainCheckpoint.txt");
    std::ofstream os(checkpointPath);
    // current step
    os << static_cast<unsigned int>(currentStep) << std::endl;
    // dated species tree
    // we rescale the branch lengths because the relative dates will be
    // unserialized from the branch lengths
    speciesTree->getDatedTree().rescaleBranchLengths();
    os << speciesTree->getTree().getNewickString() << std::endl;
  }
  assert(perFamilyModelParameters.size() == familyNames.size());
  for (unsigned int f = 0; f < perFamilyModelParameters.size(); ++f) {
    // param vectors
    auto family = familyNames[f];
    auto &mp = perFamilyModelParameters[f];
    std::string outputFilePath =
        FileSystem::joinPaths(outputDir, family + ".txt");
    std::ofstream os(outputFilePath);
    os << mp.getParamTypeNumber() << " " << mp.getSpeciesBranchNumber()
       << std::endl;
    for (unsigned int i = 0; i < mp.getParameters().dimensions(); ++i) {
      os << mp.getParameters()[i] << " ";
    }
    os << std::endl;
  }
  ParallelContext::barrier();
}

std::string AleState::readSpeciesTreeNewick(const std::string inputDir) {
  std::string checkpointPath =
      FileSystem::joinPaths(inputDir, "mainCheckpoint.txt");
  std::ifstream is(checkpointPath);
  if (!is.good()) {
    Logger::error << "Error, can't read checkpoint file " << checkpointPath
                  << std::endl;
    ParallelContext::abort(32);
  }
  unsigned int bufferUint = 0;
  is >> bufferUint;
  std::string speciesTreeNewick;
  is >> speciesTreeNewick;
  return speciesTreeNewick;
}

void AleState::unserialize(const std::string &inputDir,
                           const std::vector<std::string> &perCoreFamilyNames) {
  ParallelContext::barrier();
  familyNames = perCoreFamilyNames;
  std::string checkpointPath =
      FileSystem::joinPaths(inputDir, "mainCheckpoint.txt");
  std::ifstream is(checkpointPath);
  if (!is.good()) {
    Logger::error << "Error, can't read checkpoint file " << checkpointPath
                  << std::endl;
    ParallelContext::abort(32);
  }
  unsigned int bufferUint = 0;
  std::string bufferStr;
  // current step
  is >> bufferUint;
  currentStep = static_cast<AleStep>(bufferUint);
  // dated species tree
  is >> bufferStr;
  speciesTree = std::make_unique<SpeciesTree>(bufferStr, false);
  // param vectors
  for (auto &family : perCoreFamilyNames) {
    std::string paramsPath = FileSystem::joinPaths(inputDir, family + ".txt");
    std::ifstream is(paramsPath);
    if (!is.good()) {
      Logger::error << "Can't open parameter checkpoint file " << paramsPath
                    << std::endl;
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
    perFamilyModelParameters.push_back(mp);
  }
  assert(perCoreFamilyNames.size() == perFamilyModelParameters.size());
  ParallelContext::barrier();
}
