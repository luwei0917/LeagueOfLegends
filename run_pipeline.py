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


import gffutils
# dbA = gffutils.create_db("/Users/weilu/Dropbox/genome_algorithms_comp519/project/Challenge_9934185_A.chromosomes/A.gff3", "chromosomeA")
# dbB = gffutils.create_db("/Users/weilu/Dropbox/genome_algorithms_comp519/project/Challenge_9934185_B.chromosomes/B.gff3", "chromosomeB")
# db_A_to_B = gffutils.create_db("/Users/weilu/Dropbox/genome_algorithms_comp519/project/results/A_to_B/new_result/uniuqe_gff3_A_to_B.gff", "chromosome_A_to_B")
# db_A_to_B = gffutils.create_db("/Users/weilu/Dropbox/genome_algorithms_comp519/project/results/GMAP_A_to_B/AtoBmap_match_gene.gff3", "GMAP_chromosome_A_to_B")
# db_B_to_A = gffutils.create_db("/Users/weilu/Dropbox/genome_algorithms_comp519/project/results/GMAP_B_to_A/BtoAmap_match_gene.gff3", "GMAP_chromosome_B_to_A")

parser = argparse.ArgumentParser(description="This is my playground for current project")
# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("-f", "--fromChromosome", type=str, default="A")
parser.add_argument("-t", "--toChromosome", type=str, default="B")
parser.add_argument("--AtoBTracking", type=str, default="database/AtoBGMAP.tracking")
parser.add_argument("--BtoATracking", type=str, default="database/BtoAGMAP.tracking")
parser.add_argument("--dbA", type=str, default="database/chromosomeA")
parser.add_argument("--dbB", type=str, default="database/chromosomeB")
parser.add_argument("--dbAtoB", type=str, default="database/GMAP_chromosome_A_to_B")
parser.add_argument("--dbBtoA", type=str, default="database/GMAP_chromosome_B_to_A")
parser.add_argument("-o", "--out", type=str, default=".", help="where you store your results")
parser.add_argument("--validation", type=str, default="./dev_validation_set.tsv", help="focus on the one has ground truth. speed up everything")
parser.add_argument("-n", "--name", type=str, default="GMAP_combined_nov06")
parser.add_argument("--blastnAtoB", type=str, default="database/blastn_A_to_B.pkl")
parser.add_argument("--blastnBtoA", type=str, default="database/blastn_B_to_A.pkl")
parser.add_argument("-d", "--debug", action="store_true", default=False, help="print out additional information")
parser.add_argument("--trmap_AtoB", type=str, default="database/trmap_AtoBGMAP")
parser.add_argument("--trmap_BtoA", type=str, default="database/trmap_BtoAGMAP")
args = parser.parse_args()

with open('gg_cmd.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


# A = "A"
# B = "B"
A = args.fromChromosome
B = args.toChromosome

# dbA = gffutils.FeatureDB(f"chromosome{A}")
# dbB = gffutils.FeatureDB(f"chromosome{B}")
# db_A_to_B = gffutils.FeatureDB(f"GMAP_chromosome_{A}_to_{B}")
# db_B_to_A = gffutils.FeatureDB(f"GMAP_chromosome_{B}_to_{A}")
# fileLocation = "/Users/weilu/Dropbox/genome_algorithms_comp519/project/results/A_to_B/result_A_to_B.tracking"

dbA = gffutils.FeatureDB(args.dbA)
dbB = gffutils.FeatureDB(args.dbB)
db_A_to_B = gffutils.FeatureDB(args.dbAtoB)
db_B_to_A = gffutils.FeatureDB(args.dbBtoA)



# fileLocation = f"/Users/weilu/Dropbox/genome_algorithms_comp519/project/results/GMAP_{A}_to_{B}/{A}to{B}GMAP.tracking"
fileLocation = args.AtoBTracking
raw_data_A_to_B = read_tracking(fileLocation)


# fileLocation = f"/Users/weilu/Dropbox/genome_algorithms_comp519/project/results/GMAP_{B}_to_{A}/{B}to{A}GMAP.tracking"
fileLocation = args.BtoATracking
raw_data_B_to_A = read_tracking(fileLocation)

raw_data_A_to_B["fromGenome"] = A
raw_data_A_to_B["toGenome"] = B
raw_data_B_to_A["fromGenome"] = B
raw_data_B_to_A["toGenome"] = A

raw_data = pd.concat([raw_data_A_to_B, raw_data_B_to_A]).reset_index(drop=True)
database_dic = {A:dbA, B:dbB, A+"to"+B:db_A_to_B, B+"to"+A:db_B_to_A}

# pre = "/Users/weilu/Dropbox/genome_algorithms_comp519/project/Challenge_9934185_scoring/"
pre = args.out
os.system(f'mkdir -p {pre}')
write_to = f"{pre}/{args.name}.tsv"

if args.validation != "None":
    # fileLocation = "/Users/weilu/Dropbox/genome_algorithms_comp519/project/Challenge_9934185_scoring/dev_validation_set.tsv"
    fileLocation = args.validation
    dev_validation_set = read_result(fileLocation)
    # print(len(dev_validation_set))
    # print(len(dev_validation_set.query("SourceA == 'A'")))

    # keep transcriptA if it is inside SourceA_Transcript_ID or SourceB_Transcript_ID
    SourceA_Transcript_ID_list = dev_validation_set["SourceA_Transcript_ID"].tolist()
    SourceB_Transcript_ID_list = dev_validation_set["SourceB_Transcript_ID"].tolist()
    filtered_data = raw_data.query("fromTranscript in @SourceA_Transcript_ID_list or fromTranscript in @SourceB_Transcript_ID_list")
else:
    filtered_data = raw_data
# print("filtered_data", len(filtered_data))
# reverse_gffcompare_fileLocation = '/Users/weilu/Dropbox/genome_algorithms_comp519/project/results/A_to_B/fake_A_to_B.tracking'
# fake_A_to_B = read_tracking(reverse_gffcompare_fileLocation)
# filtered_fake_data = fake_A_to_B.query("fromTranscript in @SourceA_Transcript_ID_list or fromTranscript in @SourceB_Transcript_ID_list")
# # print("filtered_fake_data", len(filtered_fake_data))
# fake_cat_k_toTranscript_list = fake_A_to_B.query("category == 'k'")["toTranscript"].tolist()


data = filtered_data
out = ""
# fromGenome="A"
# toGenome="B"
# fromDB=dbA
# map_db = db_A_to_B
# toDB = dbB
# fromGenome="B"
# toGenome="A"
# fromDB=dbB
# map_db = db_B_to_A
# toDB = dbA
for i, line in data.iterrows():
    c = line["category"]
    fromTranscript = line["fromTranscript"]
    originalFromTranscript = line["originalFromTranscript"]
    toTranscript = line["toTranscript"]
    fromGenome = line["fromGenome"]
    fromDB = database_dic[fromGenome]
    toGenome = line["toGenome"]
    toDB = database_dic[toGenome]
    map_db = database_dic[fromGenome+"to"+toGenome]
    # get the fromGene information.
    fromGene = getFromGene(fromDB, fromTranscript)
    toGene = line["toGene"]
    extra = ""
    assign = ""
    assign = assign_unique(c, fromTranscript)
    # additional checks
    assign, extra = assign_absent_gene(c, assign, map_db, toDB, fromTranscript, toGene)
    assign, gene_fusion_list, toGene, extra = assign_gene_fusion(c, assign, fromGenome, toGenome, originalFromTranscript, toTranscript, fromGene, toGene, toDB, map_db)
    assign, extra = assign_absent_transcript(c, assign)

    out += assign_to_out(assign, fromGenome, toGenome, fromTranscript, toTranscript, fromGene, toGene)
    out += extra
    out += "\n"
with open(write_to, "w") as o:
    o.write(out)


print(pad_with_dash("Done initial processing"))
fileLocation = write_to
mySolution = read_result(fileLocation)

is_duplicated = mySolution.duplicated(subset="SourceA_Transcript_ID", keep=False)
no_duplication = mySolution[~is_duplicated]
with_duplication = mySolution[is_duplicated]

# print(len(with_duplication["SourceA_Transcript_ID"].unique()))


print(pad_with_dash("consider gene fusion"))
# if has gene fusion then use it.
# if not, the nothing change.
def only_keep_gene_fusion_if_it_exist(d):
    if "gene_fusion" in d["Call"].tolist():
        return d.query("Call == 'gene_fusion'")
    else:
        return d
# gene fusion
with_duplication = with_duplication.groupby("SourceA_Transcript_ID").apply(only_keep_gene_fusion_if_it_exist).reset_index(drop=True)


print(pad_with_dash("consider unique and absent transcript both exists"))
# if one transcript has both unique_transcript and absent_transcirpt match,
# only keep the unique one.
# for example
# weilu > grep "transcriptA128220" AtoBGMAP_sensitivity.tracking
# TCONS_00003300	XLOC_000853	geneB20367|transcriptB20365	k	q1:transcriptA128220.path2|transcriptA128220.mrna2|1|0.000000|0.000000|0.000000|3019
# TCONS_00145525	XLOC_040599	geneB11392|transcriptB11386	j	q1:transcriptA128220.path1|transcriptA128220.mrna1|11|0.000000|0.000000|0.000000|3426
def only_keep_unique(d):
    if "unique_transcript" in d["Call"].tolist():
        return d.query("Call == 'unique_transcript'")
    else:
        return d

with_duplication = with_duplication.groupby("SourceA_Transcript_ID").apply(only_keep_unique).reset_index(drop=True)

print(pad_with_dash("remove_duplicated_absent_transcript"))
# here, just keep the first one if we still have duplicated absent_transcript.
def remove_duplicated_absent_transcript(d):
    if "absent_transcript" in d["Call"].tolist():
        return d.head(1)
    else:
        return d
with_duplication = with_duplication.groupby("SourceA_Transcript_ID").apply(remove_duplicated_absent_transcript).reset_index(drop=True)

# with_duplication = with_duplication.drop_duplicates(subset="SourceA_Transcript_ID", keep="first")
# with_duplication
# fileLocation = "/Users/weilu/Dropbox/genome_algorithms_comp519/project/Challenge_9934185_scoring/cleaned_GMAP_oct30.tsv"
write_to = f"{pre}/{args.name}_post_modification_1.tsv"
df = pd.concat([no_duplication, with_duplication])
pandas_to_tsv(write_to, df)

print(pad_with_dash("consider multiple transcript"))
fileLocation = write_to
mySolution = read_result(fileLocation)
# AtoB_compare_filename = '/Users/weilu/Research/LeagueOfLegends/tracking/trmap_AtoBGMAP'
# BtoA_compare_filename = '/Users/weilu/Research/LeagueOfLegends/tracking/trmap_BtoAGMAP'
AtoB_compare_filename = args.trmap_AtoB
BtoA_compare_filename = args.trmap_BtoA
AtoB_bestmatch_list = AtoB_bestmatch(AtoB_compare_filename)
BtoA_bestmatch_list = BtoA_bestmatch(BtoA_compare_filename)
complete = AtoB_bestmatch_list
complete.update(BtoA_bestmatch_list)

newSolution = mySolution.copy()
for i, line in mySolution.iterrows():
    fromTranscript = line["SourceA_Transcript_ID"]
    toTranscript = line["SourceB_Transcript_ID"]
    toGenome = line["SourceB"]
    a = complete[fromTranscript]
    toDB = database_dic[toGenome]
    if len(a) == 0:
        # print("Not found", fromTranscript)
        pass
    elif len(a) == 1:
        # print("Unique transscript", fromTranscript)
        if toTranscript != a[0]:
            if line["Call"] == "gene_fusion":
                # ok
                pass
            else:
                print("why they differ", fromTranscript, toTranscript, a[0])
    else:
        # print("Multiple transcript", fromTranscript)
        toGenes = [getFromGene(toDB, x) for x in a]
        if len(set(toGenes)) != 1:
            print("They don't have same gene?", fromTranscript, a, toGenes)
        for x in a:
            new_line = line.copy()
            new_line["Call"] = "multiple_transcript"
            new_line["SourceB_Transcript_ID"] = x
            newSolution = newSolution.append(new_line)
    newSolution = newSolution.reset_index(drop=True)
# if one transcript has both multiple_transcript and others match,
# only keep the unique one.

def only_keep_multiple_transcript(d):
    if "multiple_transcript" in d["Call"].tolist():
        return d.query("Call == 'multiple_transcript'")
    else:
        return d

mySolution = newSolution.groupby("SourceA_Transcript_ID").apply(only_keep_multiple_transcript).reset_index(drop=True)

write_to = f"{pre}/{args.name}_post_modification_2.tsv"
pandas_to_tsv(write_to, mySolution)


print(pad_with_dash("blastn result further process"))
fileLocation = write_to
mySolution = read_result(fileLocation)
# fileLocation = "/Users/weilu/Dropbox/genome_algorithms_comp519/project/results/blastn/blastn_A_to_B.dms"
fileLocation = args.blastnAtoB
if args.blastnAtoB[-4:] == ".pkl":
    data1 = pd.read_pickle(fileLocation)
else:
    data1 = pd.read_csv(fileLocation, sep="\t", names=["qseqid", "sseqid",
                                            "pident", "length", "mismatch",
                                            "gapopen", "qstart", "qend", "sstart",
                                            "send", "evalue", "bitscore", "slen", "qlen"])
    # each query-subject pair, I only keep the first(the longest alignment)
    data1["query_gene"] = data1["qseqid"].apply(lambda x: getFromGene(dbA, x))
    data1["subject_gene"] = data1["sseqid"].apply(lambda x: getFromGene(dbB, x))
    data1.to_pickle("database/blastn_A_to_B.pkl")

fileLocation = args.blastnBtoA
if args.blastnBtoA[-4:] == ".pkl":
    data2 = pd.read_pickle(fileLocation)
else:
    data2 = pd.read_csv(fileLocation, sep="\t", names=["qseqid", "sseqid",
                                            "pident", "length", "mismatch",
                                            "gapopen", "qstart", "qend", "sstart",
                                            "send", "evalue", "bitscore", "slen", "qlen"])
    # each query-subject pair, I only keep the first(the longest alignment)
    data2["query_gene"] = data2["qseqid"].apply(lambda x: getFromGene(dbA, x))
    data2["subject_gene"] = data2["sseqid"].apply(lambda x: getFromGene(dbB, x))
    data2.to_pickle("database/blastn_B_to_A.pkl")

df1 = data1.groupby(["qseqid", "sseqid"]).head(1).reset_index(drop=True)
df2 = data2.groupby(["qseqid", "sseqid"]).head(1).reset_index(drop=True)


df = pd.concat([df1, df2]).reset_index(drop=True)


df["ratio"] = df["length"]**2/df["qlen"]
# filter out pident < 95
df = df.query("pident > 95").reset_index(drop=True)
# for those queries that target same subject. the one with largest ratio wins.
df = df.sort_values(["sseqid", "ratio"], ascending=[True, False]).reset_index(drop=True)
df["rank"] = df.groupby("sseqid")["ratio"].rank(ascending=False)
# for every subject transcript. only keep the one with highest ratio.
df2 = df.groupby("sseqid").head(1).reset_index(drop=True)
# fromDB = dbA
# toDB = dbB
# database_dic[x[10]]
# df2["query_gene"] = df2["qseqid"].apply(lambda x: getFromGene(database_dic[x[10]], x))
# df2["subject_gene"] = df2["sseqid"].apply(lambda x: getFromGene(database_dic[x[10]], x))

write_to = f"{pre}/{args.name}_post_modification_3.tsv"

with open(write_to, "w") as out:
    for i, line in mySolution.iterrows():
        call = line["Call"]
        if call == "absent_transcript":
            toGene = line["SourceB_Gene"]
            a = df2.query(f"subject_gene == '{toGene}'")
            if len(a) > 1:
                if args.debug:
                    print("---should not happen----")
                    print(a)
            elif len(a) == 1:
                # print(a)
                query_gene = a["query_gene"].iloc[0]
                fromGene = line["SourceA_Gene"]
                # print(i, fromGene)
                if fromGene == query_gene:
                    pass
                else:
                    # change to absent_gene
                    # print(i)
                    line["Call"] = "absent_gene"
                    line["SourceB_Gene"] = ""
        line = line.fillna('')
        line["Score"] = str(line["Score"])
        formated_line = "\t".join(line.values) + "\n"
        out.write(formated_line)

# notes:
# if you want to know which chromosome a transcript belongs, you can do the following.
# transcript = dbA["transcriptA1"]
# transcript.chrom

# if you have the following error
# ModuleNotFoundError: No module named 'pandas.core.internals.managers'; 'pandas.core.internals' is not a package
# weis-mbp-2:LeagueOfLegends weilu$ python3 -m pip install --upgrade pandas