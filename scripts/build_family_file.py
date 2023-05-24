#!/usr/bin/env python
import os
import sys

"""
Script that generates a family file for GeneRax.
"""

def get_family(filename):
  return filename.split(".")[0]

def join_abs(directory, filename):
  return os.path.abspath(os.path.join(directory, filename))

def build_families(trees_dir, mappings_dir, output_file):
  families = {}
  family_names = {}
  trees = {}
  mappings = {}

  if (trees_dir != "NONE"):
    for tree in os.listdir(trees_dir):
      trees[get_family(tree)] = join_abs(trees_dir, tree)
      family_names[get_family(tree)] = True
  if (mappings_dir != "NONE"):
    for mapping in os.listdir(mappings_dir):
      mappings[get_family(mapping)] = join_abs(mappings_dir, mapping)
      family_names[get_family(mapping)] = True
  for family in family_names:
    cell = {}
    if (family in trees):
      cell["starting_gene_tree"] = trees[family]
    if (family in mappings):
      cell["mapping"] = mappings[family]
    if (family in families):
      print("[Error]: family " + family + " was found twice in the alignments directory.")
      exit(1)
    families[family] = cell
  
  with open(output_file, "w") as writer:
    writer.write("[FAMILIES]\n")
    for family in families:
      cell = families[family]
      writer.write("- " + family + "\n")
      for param in cell:
        writer.write(param + " = " + cell[param] + "\n")
  print("Output families file" + output_file)

if (__name__ == "__main__"): 
  if (len(sys.argv) != 4): 
     print(" [Error] Invalid syntax")
     print(" Usage: python " + os.path.basename(__file__) + " trees_dir mappings_dir output_file")
     print(" Example: python "  + os.path.basename(__file__) + " /home/mytrees/ /home/mymappings/ /home/families.txt")
     exit(1)
  trees_dir = sys.argv[1]
  mappings_dir = sys.argv[2]
  output_file = sys.argv[3]
  build_families(trees_dir, mappings_dir, output_file)




