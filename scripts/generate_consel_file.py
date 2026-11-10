#!/usr/bin/env python3
import os
import sys

def generate_consel_file(output_file, rundirs):
  families = []
  trees = []
  tree_lls = {}
  tree_rundir = {}
  index = 0
  for rundir in rundirs:
    tree = str(index + 1)
    trees.append(tree)
    tree_rundir[tree] = rundir
    lls = {}
    # read the likelihoods
    ll_file = os.path.join(rundir, "per_fam_likelihoods.tsv")
    print("Reading " + ll_file)
    with open(ll_file, "r") as reader:
      # skip the header
      reader.readline()
      for line in reader:
        items = line.split()
        family = items[0]
        ll = float(items[1])
        lls[family] = ll
        if (index == 0):
          families.append(family)
    if (index != 0):
      if (len(lls) != len(families)):
        print(f"[Error] {rundir} has {len(lls)} per-family likelihoods while {tree_rundir[trees[0]]} has {len(families)} per-family likelihoods")
        exit(1)
    tree_lls[tree] = lls
    index += 1
  assert(len(families) > 0)
  print("Number of families: " + str(len(families)))
  print("Number of trees: " + str(len(trees)))
  print("Mapping tree items to rundirs:")
  print("  item\trundir")
  for tree in trees:
    print("  " + tree + "\t" + tree_rundir[tree])
  # write the consel file
  print("Writing the per-family likelihoods into " + output_file)
  with open(output_file, "w") as writer:
    writer.write("  " + str(len(trees)) + " " + str(len(families)) + "\n")
    for tree in trees:
      writer.write(tree)
      lls = tree_lls[tree]
      for family in families:
        writer.write(" " + str(lls[family]))
      writer.write("\n")
  print("You can now analyze the output file with the consel suite by running:")
  print("  makermt --puzzle " + output_file)
  print("  consel " + output_file.split(".")[0] + ".rmt #WARNING: the exact filename might be slightly different")
  print("  catpv " + output_file.split(".")[0] + ".pv")

# main
if (__name__ == "__main__"):
  if (len(sys.argv) < 4):
    print(" [Error] Invalid syntax")
    print(" Usage: python3 " + os.path.basename(__file__) + " output_file aleraxrundir1 aleraxrundir2 [aleraxrundir3 ...]")
    exit(1)
  output_file = sys.argv[1]
  rundirs = sys.argv[2:]
  generate_consel_file(output_file, rundirs)
