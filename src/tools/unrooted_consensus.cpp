#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <trees/PLLUnrootedTree.hpp>

std::vector<std::string> getLines(const std::string &path) {
  std::ifstream is(path);
  std::string line;
  std::vector<std::string> lines;
  while (std::getline(is, line)) {
    if (line.size() > 2) {
      lines.push_back(line);
    }
  }
  is.close();
  return lines;
}

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Error, syntax is:" << std::endl;
    std::cerr << "consensus input_trees threshold output_consensus"
              << std::endl;
    return 1;
  }
  std::string inputFile(argv[1]);
  double threshold = atof(argv[2]);
  std::string outputFile(argv[3]);
  auto utreeStrs = getLines(inputFile);
  auto consUtreeStr = PLLUnrootedTree::buildConsensusTree(utreeStrs, threshold);
  std::ofstream os(outputFile);
  os << consUtreeStr << std::endl;
  std::cout << consUtreeStr << std::endl;
  os.close();
  return 0;
}
