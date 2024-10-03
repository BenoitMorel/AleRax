#include <fstream>
#include <iostream>
#include <trees/PLLUnrootedTree.hpp>

std::vector<std::string> getLines(const std::string path) {
  std::ifstream is(path);
  std::string line;
  std::vector<std::string> res;
  while (std::getline(is, line)) {
    if (line.size() > 2) {
      res.push_back(line);
    }
  }
  return res;
}

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Error, syntax is:" << std::endl;
    std::cerr << "consensus input_trees threshold output_consensus"
              << std::endl;
    return 1;
  }
  std::string input(argv[1]);
  double threshold = atof(argv[2]);
  std::string output(argv[3]);

  auto newicks = getLines(input);
  auto outNewick = PLLUnrootedTree::buildConsensusTree(newicks, threshold);
  std::ofstream os(output);
  os << outNewick << std::endl;
  std::cout << outNewick << std::endl;

  return 0;
}
