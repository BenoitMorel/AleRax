#!/usr/bin/env python
import os
import sys


def generate_consel_file(output, rundirs):
  tree_to_lls = {}
  tree_file = {}
  index = 0
  columns = 0
  first_tree = ""
  families = []
  trees = []
  for rundir in rundirs:
    tree = "item " + str(index + 1)
    trees.append(tree)
    if (index == 0):
      first_tree = tree
    lls = {}
    tree_file[tree] = rundir
    
    # read the likelihoods
    ll_file = os.path.join(rundir, "per_fam_likelihoods.txt")
    print("reading " + ll_file)
    for line in open(ll_file).readlines():
      sp = line.split()
      fam = sp[0]
      ll = float(sp[1])
      lls[fam] = ll
      if (index == 0):
        families.append(fam)
    if (index == 0):
      columns = len(lls)
    else:
      if (columns != len(lls)):
        print("Error: " + rundir + " has " + str(len(lls)) + " per-family likelihoods while " + tree_file[tree] + " has " + str(columns) + " per-family likelihoods")
        print("Aborting")
        sys.exit(1)
 
    tree_to_lls[tree] = lls
    index += 1
        
  assert(columns == len(families))
  print("Number of per-family likelihoods: " + str(columns))
  print("Number of trees:" + str(len(trees)))
  print("Mapping tree-name to rundir:")
  for tree in trees:
    print("    " + tree + " " + tree_file[tree])
  print("Writing the per-family likelihoods into " + output)
  
  with open(output, "w") as writer:
    writer.write("  " + str(len(trees)) + " " + str(columns) + "\n")
    for tree in trees:
      writer.write(tree + " ")
      lls = tree_to_lls[tree]
      for fam in families:
        writer.write(str(str(lls[fam]) + " "))
      writer.write("\n")
  print("You can now analyze the output file with the consel suite by running:")
  print("  makermt --puzzle " + output) 
  print("  consel " + output.split(".")[0] + ".rmt #WARNING: the exact filename might be slightly different")
  print("  catpv " + output.split(".")[0] + ".pv")

if (__name__ == "__main__"): 
  if (len(sys.argv) < 4): 
    print(" Usage: python " + os.path.basename(__file__) + " output_consel_file aleraxrundir1 aleraleraxrundir2 [aleraxrundir3 ...]")
    sys.exit(1)
  output = sys.argv[1]
  rundirs = sys.argv[2:]
  generate_consel_file(output, rundirs)

