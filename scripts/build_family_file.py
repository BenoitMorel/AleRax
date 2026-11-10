#!/usr/bin/env python3
import os
import sys

"""
Script that generates a family file for AleRax
"""

def get_family_name(filename):
  return filename.split(".")[0]

def join_abs(directory, filename):
  return os.path.abspath(os.path.join(directory, filename))

def build_families(trees_dir, mappings_dir, output_file):
  names = []
  trees = {}
  mappings = {}
  # read the directories
  if (trees_dir not in ["-", "NONE"]):
    for tree in os.listdir(trees_dir):
      name = get_family_name(tree)
      names.append(name)
      if (name in trees):
        print("[Error] Family " + name + " found twice in the trees directory")
        exit(1)
      trees[name] = join_abs(trees_dir, tree)
  if (mappings_dir not in ["-", "NONE"]):
    for mapping in os.listdir(mappings_dir):
      name = get_family_name(mapping)
      if (name not in names):
        print("[Warning] Mapping " + mapping + " does not correspond to any tree and will be ignored")
      if (name in mappings):
        print("[Error] Family " + name + " found twice in the mappings directory")
        exit(1)
      mappings[name] = join_abs(mappings_dir, mapping)
  # fill the families
  families = {}
  for name in names:
    family = {}
    if (name in trees):
      family["gene_tree"] = trees[name]
    if (name in mappings):
      family["mapping"] = mappings[name]
    families[name] = family
  # write the families file
  with open(output_file, "w") as writer:
    writer.write("[FAMILIES]\n")
    for name in families:
      family = families[name]
      writer.write("- " + name + "\n")
      for param in family:
        writer.write(param + " = " + family[param] + "\n")
  print(f"Output families file {output_file} with {len(families)} families")

# main
if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print(" [Error] Invalid syntax")
    print(" Usage: python3 " + os.path.basename(__file__) + " trees_dir mappings_dir output_file")
    print(" Example: python3 " + os.path.basename(__file__) + " /home/mytrees/ /home/mymappings/ /home/families.txt")
    exit(1)
  trees_dir = sys.argv[1]
  mappings_dir = sys.argv[2]
  output_file = sys.argv[3]
  build_families(trees_dir, mappings_dir, output_file)
