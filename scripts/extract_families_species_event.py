#!/usr/bin/env python
import os
import sys



def get_event_index(event):
  if (event == "S"):
    return 1
  if (event == "D"):
    return 2
  if (event == "L"):
    return 3
  if (event == "T"):
    return 4
  assert(False)

def get_frequency(event_file, species, event):
  reader = open(event_file)
  # skip the first line
  reader.readline() 
  event_index = get_event_index(event)
  while(True):
    line = reader.readline()
    if (not line):
      break
    sp = line.split(",")
    if (species == sp[0]):
      return float(sp[event_index])
  return 0.0

def extract(alerax_dir, species, event, min_value, output_file):
  res = []
  summaries = os.path.join(alerax_dir, "reconciliations", "summaries")
  for f in os.listdir(summaries):
    if (f.endswith("_perspecies_eventcount.txt")):
      event_file = os.path.join(summaries, f)
      freq = get_frequency(event_file, species, event)
      if (freq > min_value):
        family = f.replace("_perspecies_eventcount.txt", "")
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
     print("This script takes as input an AleRax output directory, the name of a species, an event type (D, T, L, S), and a threshold value")
     print("It will store in a file the list of families for which the frequency of this event is above this threshold, as well as the frequencies")
     print(" [Error] Invalid syntax")
     print(" Usage: python " + os.path.basename(__file__) + " alerax_dir species event_type min_value output_file")
     print(" Event can be S, D, L, or T")
     exit(1)
  alerax_dir = sys.argv[1]
  species = sys.argv[2]
  event = sys.argv[3]
  min_value = float(sys.argv[4])
  output_file = sys.argv[5]
  
  if (not event in ["S", "D", "L", "T"]):
     print(" Event can be S, D, L, or T")
     exit(1)
  extract(alerax_dir, species, event, min_value, output_file)





