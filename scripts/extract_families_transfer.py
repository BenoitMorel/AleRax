#!/usr/bin/env python
import os
import sys




def get_frequency(transfer_file, species_from, species_to):
  reader = open(transfer_file)
  while(True):
    line = reader.readline()
    if (not line):
      break
    sp = line.split(" ")
    if (len(sp) < 2):
      continue
    if (species_from != sp[0]):
      continue
    if (species_to != sp[1]):
      continue
    print(sp[2])
    return float(sp[2])
  return 0.0

def extract(alerax_dir, from_species, to_species, min_value, output_file):
  res = []
  summaries = os.path.join(alerax_dir, "reconciliations", "summaries")
  for f in os.listdir(summaries):
    if (f.endswith("_transfers.txt")):
      family = f.replace("_transfers.txt", "")
      event_file = os.path.join(summaries, f)
      freq = get_frequency(event_file, from_species, to_species)
      if (freq > min_value):
        res.append((freq, family))
  res.sort(reverse = True)
  total = 0.0
  for item in res:
    total += item[0]
  with open(output_file, "w") as writer:
    writer.write(str(total) + "\n")
    for item in res:
      writer.write(str(item[0]) + " " + item[1] + "\n")

if (__name__ == "__main__"): 
  if (len(sys.argv) != 6): 
     print("This script takes as input an AleRax output directory, two species,and a threshold value")
     print("It will store in a file the list of families for which the frequency of the transfer between those two species  is above this threshold, as well as the frequencies")
     print(" [Error] Invalid syntax")
     print(" Usage: python " + os.path.basename(__file__) + " alerax_dir from to min_value output_file")
     exit(1)
  alerax_dir = sys.argv[1]
  from_species = sys.argv[2]
  to_species = sys.argv[3]
  min_value = float(sys.argv[4])
  output_file = sys.argv[5]
  
  extract(alerax_dir, from_species, to_species, min_value, output_file)







