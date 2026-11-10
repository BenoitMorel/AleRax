#!/usr/bin/env python3
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
  if (event == "O"):
    return 6
  assert(False)
  return 0

def get_frequency(event_file, species, event):
  event_index = get_event_index(event)
  with open(event_file, "r") as reader:
    # skip the header
    reader.readline()
    for line in reader:
      items = line.split()
      assert(len(items) == 10)
      if (species == items[0]):
        return float(items[event_index])
  return 0.0

def extract(alerax_dir, species, event, min_freq, output_file):
  res = []
  summaries = os.path.join(alerax_dir, "reconciliations", "summaries")
  for file in os.listdir(summaries):
    if (file.endswith("_meanSpeciesEventCounts.tsv")):
      event_file = os.path.join(summaries, file)
      freq = get_frequency(event_file, species, event)
      if (freq > min_freq):
        family = file.replace("_meanSpeciesEventCounts.tsv", "")
        res.append((freq, family))
  res.sort(reverse = True)
  total = 0.0
  for p in res:
    total += p[0]
  print(f"Writing families with more than {min_freq} {event} events in {species} into {output_file}")
  with open(output_file, "w") as writer:
    writer.write(str(total) + "\n")
    for p in res:
      writer.write(str(p[0]) + " " + p[1] + "\n")
  print(f"{len(res)} families written")

# main
if (__name__ == "__main__"):
  if (len(sys.argv) != 6):
    print("This script takes as input an AleRax output directory, a species name, an event type, and a threshold value")
    print("It will store in the output file the list of families for which the frequency of this event in that species is above the threshold, as well as the frequencies")
    print(" [Error] Invalid syntax")
    print(" Usage: python3 " + os.path.basename(__file__) + " alerax_dir species event threshold output_file")
    print(" Event can be one of S, D, L, T, O")
    exit(1)
  alerax_dir = sys.argv[1]
  species = sys.argv[2]
  event = sys.argv[3]
  min_freq = float(sys.argv[4])
  output_file = sys.argv[5]
  if (event not in ["S", "D", "L", "T", "O"]):
    print("[Error] Invalid event type: " + event)
    print("Event can be one of S, D, L, T, O")
    exit(1)
  extract(alerax_dir, species, event, min_freq, output_file)
