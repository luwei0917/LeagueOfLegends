import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import imp
import numpy as np
import fileinput
from itertools import product
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
from os import listdir
from scipy.interpolate import griddata
import matplotlib as mpl
from Bio.PDB.Polypeptide import one_to_three
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.PDBParser import PDBParser
from collections import defaultdict
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from helper_functions import *

parser = argparse.ArgumentParser(description="This is my playground for current project")
# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("-s", "--solution", type=str, default="test_nov16/GMAP_combined_nov06_post_modification_4.tsv")
parser.add_argument("-v", "--validation", type=str, default="dev_validation_set.tsv")
args = parser.parse_args()

with open('cmd_compare_results.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


# fileLocation = "dev_validation_set.tsv"
fileLocation = args.validation
dev_validation_set = read_result(fileLocation)
# fileLocation = "test/GMAP_combined_nov06_post_modification_3.tsv"
fileLocation = args.solution
mySolution = read_result(fileLocation)

# only consider sourceA is A for now.
# ground_truth = dev_validation_set.query("SourceA == 'B'")
ground_truth = dev_validation_set
print(len(ground_truth))
print(len(mySolution))
my_summary = {}

def compare_result(call, ground_truth, mySolution, my_summary):
    print(pad_with_dash(f"should be {call}"))
    count = 0
    correct_count = 0
    for i, line in ground_truth.query(f"Call == '{call}'").reset_index(drop=True).iterrows():
        # print(i)
        SourceA_Transcript_ID = line["SourceA_Transcript_ID"]
        call = line["Call"]
        solution = mySolution.query(f"SourceA_Transcript_ID == '{SourceA_Transcript_ID}'")
        if len(solution) == 0:
            ref = "\t".join(line.values.astype(str))
            sol = "None"
        elif len(solution) > 1:
            ref = "\t".join(line.values.astype(str))
            # print("ref:\n"+ref)
            sol = ""
            for j, one_sol in solution.iterrows():
                sol += "\t".join(one_sol.values.astype(str)) + '\n'
            # print("mySolution:\n"+sol)
        else:
            solution = solution.iloc[0]
            if np.alltrue(line.values.astype(str)[:8] == solution.values.astype(str)[:8]):
                # the all same.
                correct_count += 1
                continue
            else:
                ref = "\t".join(line.values.astype(str))
                sol = "\t".join(solution.values.astype(str))
        count += 1
        print("----------")
        print("ref:\n"+ref)
        print("mySolution:\n"+sol)
        print("----------")
    print("Correct: ", correct_count, "Wrong: ", count)
    my_summary[f"wrong_{call}"] = count
    my_summary[f"correct_{call}"] = correct_count

def compare_multiple_result(ground_truth, mySolution, my_summary):
    call = "multiple_transcript"
    print(pad_with_dash(f"should be {call}"))
    count = 0
    correct_count = 0
    sourceA_list = ground_truth.query(f"Call == '{call}'")["SourceA_Transcript_ID"].unique()
    for SourceA_Transcript_ID in sourceA_list:
        a = ground_truth.query(f"SourceA_Transcript_ID == '{SourceA_Transcript_ID}'")
        solution = mySolution.query(f"SourceA_Transcript_ID == '{SourceA_Transcript_ID}'")
        toTranscript_ground_truth = sorted(a["SourceB_Transcript_ID"])
        toTranscript_solution = sorted(solution["SourceB_Transcript_ID"])
        for toTranscript_solution_i in toTranscript_solution:
            if toTranscript_solution_i not in toTranscript_ground_truth:
                line = solution.query(f"SourceB_Transcript_ID == '{toTranscript_solution_i}'")
                assert len(line) == 1
                line = line.iloc[0]
                # print(toTranscript_solution_i, "Not exist in ground true")
                # print(line)
                sol = "\t".join(line.values.astype(str))
                ref = ""
                print("ref:\n"+ref)
                print("mySolution:\n"+sol)
                count += 1
        for toTranscript_ground_truth_i in toTranscript_ground_truth:
            if toTranscript_ground_truth_i not in toTranscript_solution:
                # print(toTranscript_ground_truth_i, "Not found in my solution")
                line = a.query(f"SourceB_Transcript_ID == '{toTranscript_ground_truth_i}'")
                assert len(line) == 1
                line = line.iloc[0]
                ref = "\t".join(line.values.astype(str))
                sol_line = mySolution.query(f"SourceA_Transcript_ID == '{toTranscript_ground_truth_i}'")
                sol = "\t".join(sol_line.values.astype(str))
                print("ref:\n"+ref)
                print("mySolution:\n"+sol)
                count += 1

    correct_count = len(ground_truth.query(f"Call == '{call}'")) - count
    print("Correct: ", correct_count, "Wrong: ", count)
    my_summary[f"wrong_{call}"] = count
    my_summary[f"correct_{call}"] = correct_count

compare_result("unique_transcript", ground_truth, mySolution, my_summary)
# compare_result("multiple_transcript", ground_truth, mySolution, my_summary)
compare_multiple_result(ground_truth, mySolution, my_summary)
compare_result("gene_fusion", ground_truth, mySolution, my_summary)
compare_result("absent_transcript", ground_truth, mySolution, my_summary)
compare_result("absent_gene", ground_truth, mySolution, my_summary)
compare_result("absent_genome", ground_truth, mySolution, my_summary)

print("---------------exist in reference but not in my solution-----------------")
count = 0
for i, line in ground_truth.iterrows():
    SourceA_Transcript_ID = line["SourceA_Transcript_ID"]
    call = line["Call"]
    solution = mySolution.query(f"SourceA_Transcript_ID == '{SourceA_Transcript_ID}'")
    if len(solution) == 0:
        # print(SourceA_Transcript_ID, call)
        ref = "\t".join(line.values.astype(str))
        print("ref:\n"+ref)
        print("mySolution:")
        # sol = "\t".join(solution.values.astype(str))
        print("None")
        count += 1
print("Occurrence: ", count)
my_summary["exist_in_ground_truth_but_not_in_my_solution"] = count

# print("---------------matched more than one-----------------")
# count = 0
# for i, line in ground_truth.iterrows():
#     SourceA_Transcript_ID = line["SourceA_Transcript_ID"]
#     call = line["Call"]
#     solution = mySolution.query(f"SourceA_Transcript_ID == '{SourceA_Transcript_ID}'")
#     if len(solution) > 1:
#         ref = "\t".join(line.values.astype(str))
#         print("ref:\n"+ref)
#         print("mySolution:")
#         for one_sol in solution.values:
#             sol = "\t".join(one_sol.astype(str))
#             print(sol)
#         count += 1
# print("Occurrence: ", count)
# my_summary["matched_more_than_one"] = count

# print("---------------in my solution but not in reference-----------------")
# count = 0
# for i, line in mySolution.iterrows():
#     SourceA_Transcript_ID = line["SourceA_Transcript_ID"]
#     call = line["Call"]
#     ref = ground_truth.query(f"SourceA_Transcript_ID == '{SourceA_Transcript_ID}'")
#     if len(ref) == 0:
#         ref = "None"
#         print("ref:\n"+ref)
#         sol = "\t".join(line.values.astype(str))
#         print("mySolution:\n"+sol)
#         count += 1
# print("Occurrence: ", count)
# my_summary["in_my_solution_but_not_in_ref"] = count
print(pad_with_dash("summary"))
total = 0
for item in my_summary:
    total += my_summary[item]
    print(item, my_summary[item])
print(total)