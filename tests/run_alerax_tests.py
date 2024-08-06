#!/usr/bin/env python3
import os
import sys
import subprocess
import shutil
from  distutils.spawn import find_executable

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
REPO_DIR = os.path.realpath(os.path.join(SCRIPT_DIR, os.pardir))
ALERAX = os.path.join(REPO_DIR, "build", "bin", "alerax")
FAMILIES_SCRIPT = os.path.join(REPO_DIR, "scripts", "build_family_file.py")
DATA_DIR = os.path.join(REPO_DIR, "data", "test_data")
OUTPUT = os.path.join(REPO_DIR, "tests", "outputs")

try:
  os.makedirs(OUTPUT)
except:
  pass

def is_mpi_installed():
  return find_executable("mpiexec") is not None

def get_test_name(dataset, model, heterogeneity, cores):
  test_name = dataset
  test_name += "_" + model
  test_name += "_" + heterogeneity
  test_name += "_" + str(cores)
  return test_name

def reset_dir(directory):
  shutil.rmtree(directory, ignore_errors=True)
  os.makedirs(directory)

def generate_families_file_data(test_data, test_output):
  families_file = os.path.join(test_output, "families.txt")
  command = []
  command.append("python3")
  command.append(FAMILIES_SCRIPT)
  command.append(os.path.join(test_data, "mrbayes_trees"))
  mappings_dir = os.path.join(test_data, "mappings")
  if (os.path.isdir(mappings_dir)):
    command.append(mappings_dir)
  else:
    command.append("NONE")
  command.append(families_file)
  logs_file_path = os.path.join(test_data, "families_script_logs.txt")
  with open(logs_file_path, "w") as writer:
    subprocess.check_call(command, stdout = writer, stderr = writer)
  return families_file

def generate_families_file(test_output, starting_trees = "NONE", mappings = "NONE"):
  families_file = os.path.join(test_output, "families.txt")
  command = []
  command.append("python3")
  command.append(FAMILIES_SCRIPT)
  command.append(starting_trees)
  command.append(mappings)
  command.append(families_file)
  logs_file_path = os.path.join(test_output, "families_script_logs.txt")
  with open(logs_file_path, "w") as writer:
    subprocess.check_call(command, stdout = writer, stderr = writer)
  return families_file


def is_string_in_file(string, file_name):
  return string in open(file_name).read()

def count_string_in_file(string, file_name):
  return open(file_name).read().count(string)

def check_reconciliation(test_output, model):
  reconciliations_path = os.path.join(test_output, "alerax", "reconciliations")
  nhx_dup_A = os.path.join(reconciliations_path, "gene_dup_A_reconciliated.nhx")
  nhx_dup_AB = os.path.join(reconciliations_path, "gene_dup_AB_reconciliated.nhx")
  nhx_transfer_A_D = os.path.join(reconciliations_path, "gene_transfer_A_D_reconciliated.nhx")
  if (not is_string_in_file("[&&NHX:S=A:D=Y:H=N:B=0]", nhx_dup_A)):
    print("Failed to infer a duplication in species A (" + nhx_dup_A + ")")
    return False
  if (count_string_in_file("=Y", nhx_dup_A) != 1):
    print("Inferred to many events in " + nhx_dup_A)
    return False;
  if (not is_string_in_file("[&&NHX:S=AB:D=Y:H=N:B=0]", nhx_dup_AB)):
    print("Failed to infer a duplication in species AB (" + nhx_dup_AB + ")")
    return False
  if (count_string_in_file("=Y", nhx_dup_AB) != 1):
    print("Inferred to many events in " + nhx_dup_AB)
    return False;
  if (model == "UndatedDTL"):
    if (not is_string_in_file("[&&NHX:S=A:D=N:H=Y@A@D:B=0]", nhx_transfer_A_D)):
      print("Failed to infer a transfer from A to D (" + nhx_transfer_A_D + ")")
      return False
    if (count_string_in_file("=Y", nhx_transfer_A_D) != 1):
      print("Inferred to many events in " + nhx_transfer_A_D)
      return False;
  return True

def run_alerax(test_data, test_output, families_file, model, heterogeneity, cores):
  command = []
  if (cores > 1):
    command.append("mpiexec")
    command.append("-np")
    command.append(str(cores))
  command.append(ALERAX)
  command.append("-f")
  command.append(families_file)
  command.append("-s")
  command.append(os.path.join(test_data, "speciesTree.newick"))
  command.append("--rec-model")
  command.append(model)
  if (heterogeneity == "species"):
      command.append("--model-parametrization")
      command.append("PER-SPECIES")
  if (heterogeneity == "families"):
      command.append("--model-parametrization")
      command.append("PER-FAMILY")
  command.append("-p")
  command.append(os.path.join(test_output, "alerax"))
  logs_file_path = os.path.join(test_output, "tests_logs.txt")
  with open(logs_file_path, "w") as writer:
    try:
      subprocess.check_call(command, stdout = writer, stderr = writer)
    except Exception as inst:
      print("The following command failed: ")
      print(" ".join(command))
      raise inst


def run_test(dataset, model, heterogeneity, cores):
  test_name = get_test_name(dataset, model, heterogeneity, cores)
  test_output = os.path.join(OUTPUT, test_name)
  reset_dir(test_output)
  test_data = os.path.join(DATA_DIR, dataset)
  #ok = check_reconciliation(test_output, model)
  try:
    families_file = generate_families_file_data(test_data, test_output)
    run_alerax(test_data, test_output, families_file, model, heterogeneity, cores)
    print("Test " + test_name + ": ok") 
  except Exception as inst:
    print(inst)
    print("Test " + test_name + ": FAILED") 
    return False
  return True


dataset_set = ["simulated_2"]
heterogeneity_set = ["global", "families", "species"]
model_set = ["UndatedDL", "UndatedDTL"]
cores_set = [1]
if (is_mpi_installed()):
  cores_set.append(3)

all_ok = True
for dataset in dataset_set:
  for model in model_set:
    for heterogeneity in heterogeneity_set:
      for cores in cores_set:
        all_ok = all_ok and run_test(dataset, model, heterogeneity, cores)

if (not all_ok):
  print("[Error] Some tests failed, please fix them!")
  exit(1)
else:
  print("All tests succeeded!")






