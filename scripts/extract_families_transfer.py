#!/usr/bin/env python3
import os
import sys

def get_frequency(transfer_file, src, dest):
  with open(transfer_file, "r") as reader:
    # skip the header
    reader.readline()
    for line in reader:
      items = line.split()
      assert(len(items) == 3)
      if (src == items[0] and dest == items[1]):
        return float(items[2])
  return 0.0

def extract(alerax_dir, src, dest, min_freq, output_file):
  res = []
  summaries = os.path.join(alerax_dir, "reconciliations", "summaries")
  for file in os.listdir(summaries):
    if (file.endswith("_meanTransfers.tsv")):
      transfer_file = os.path.join(summaries, file)
      freq = get_frequency(transfer_file, src, dest)
      if (freq > min_freq):
        family = file.replace("_meanTransfers.tsv", "")
        res.append((freq, family))
  res.sort(reverse = True)
  total = 0.0
  for p in res:
    total += p[0]
  print(f"Writing families with more than {min_freq} T events from {src} to {dest} into {output_file}")
  with open(output_file, "w") as writer:
    writer.write(str(total) + "\n")
    for p in res:
      writer.write(str(p[0]) + " " + p[1] + "\n")
  print(f"{len(res)} families written")

# main
if (__name__ == "__main__"):
  if (len(sys.argv) != 6):
    print("This script takes as input an AleRax output directory, two species names, and a threshold value")
    print("It will store in the output file the list of families for which the frequency of transfers between those two species is above the threshold, as well as the frequencies")
    print(" [Error] Invalid syntax")
    print(" Usage: python3 " + os.path.basename(__file__) + " alerax_dir src dest threshold output_file")
    exit(1)
  alerax_dir = sys.argv[1]
  src = sys.argv[2]
  dest = sys.argv[3]
  min_freq = float(sys.argv[4])
  output_file = sys.argv[5]
  extract(alerax_dir, src, dest, min_freq, output_file)
