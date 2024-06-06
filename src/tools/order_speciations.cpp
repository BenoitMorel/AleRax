#include <fstream>
#include <iostream>
#include <trees/DatedTree.hpp>

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Error, syntax is:" << std::endl;
    std::cerr << "order_speciation input_tree output_tree" << std::endl;
    return 1;
  }
  std::string input(argv[1]);
  std::string output(argv[2]);
  PLLRootedTree rootedTree(input, true);
  DatedTree datedTree(&rootedTree, true);
  datedTree.rescaleBranchLengths();
  std::ofstream os(output);
  os << rootedTree.getNewickString() << std::endl;
  return 0;
}
