#include <fstream>
#include <iostream>
#include <string>

#include <trees/DatedTree.hpp>

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Error, syntax is:" << std::endl;
    std::cerr << "order_speciations input_tree output_tree" << std::endl;
    return 1;
  }
  std::string inputFile(argv[1]);
  std::string outputFile(argv[2]);
  PLLRootedTree rootedTree(inputFile, true);
  DatedTree datedTree(rootedTree, true);
  datedTree.rescaleBranchLengths();
  std::ofstream os(outputFile);
  os << rootedTree.getNewickString() << std::endl;
  os.close();
  return 0;
}
