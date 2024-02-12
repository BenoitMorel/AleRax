#include "AleState.hpp"

#include <fstream>

void AleState::serialize(const std::string &outputFilePath)
{
  std::ofstream os(outputFilePath);
  // current step
  os << static_cast<unsigned int>(currentStep) << std::endl;
  // dated species tree
  // we rescale the branch lengths because the relative dates will be 
  // unserialized from the branch lengths
  speciesTree->getDatedTree().rescaleBranchLengths();
  os << speciesTree->getTree().getNewickString() << std::endl;
  // param vector

}
  

void AleState::unserialize(const std::string &inputFilePath)
{
  std::ifstream is(inputFilePath);
  unsigned int bufferUint = 0;
  std::string bufferStr;
  // current step
  is >> bufferUint;
  currentStep = static_cast<AleStep>(bufferUint);
  // dated species tree
  is >> bufferStr;
  speciesTree = std::make_unique<SpeciesTree>(bufferStr, true);
}


