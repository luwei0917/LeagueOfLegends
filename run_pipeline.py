#!/usr/bin/env python
"""
This pipeline generate the required output in one file for each pairwise mapping of genomes

to run the file:
    python run_pipeline.py --AtoBTracking <AtoBtrackingFilePath> --BtoATracking <BtoAtrackingFilePath> \
           -o <outFileDirectory> --validation <validationFileLocation>

example command:
    python run_pipeline.py --AtoBTracking ./AtoBGMAP.tracking --BtoATracking ./BtoAGMAP.tracking -o test --validation dev_validation_set.tsv

"""
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
import pickle
from scipy import stats
from sklearn.cluster import KMeans

from helper_functions import *


import gffutils
# dbA = gffutils.create_db("/Users/weilu/Dropbox/genome_algorithms_comp519/project/Challenge_9934185_A.chromosomes/A.gff3", "chromosomeA")
# dbB = gffutils.create_db("/Users/weilu/Dropbox/genome_algorithms_comp519/project/Challenge_9934185_B.chromosomes/B.gff3", "chromosomeB")
# db_A_to_B = gffutils.create_db("/Users/weilu/Dropbox/genome_algorithms_comp519/project/results/A_to_B/new_result/uniuqe_gff3_A_to_B.gff", "chromosome_A_to_B")
# db_A_to_B = gffutils.create_db("/Users/weilu/Dropbox/genome_algorithms_comp519/project/results/GMAP_A_to_B/AtoBmap_match_gene.gff3", "GMAP_chromosome_A_to_B")
# db_B_to_A = gffutils.create_db("/Users/weilu/Dropbox/genome_algorithms_comp519/project/results/GMAP_B_to_A/BtoAmap_match_gene.gff3", "GMAP_chromosome_B_to_A")

"""
parse the files command line argument

"""

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


'''

store and process the input files
- read input files into variables
- concatenate the tracking files
- when validating with dev_validation_set, only keep the rows of transcripts that are in the test file

'''

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


'''
initial processing:

- parse each line of the filtered database
- classify transcript in this sequence(may be have multiple calls for now):
    unique transcript, absent gene, gene fusion and absent transcript
- store the result of the assignment
'''

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


'''
additional processing:

- find the transcripts with multiple calls
- analyze all of its calls and keep one call only
- store the result of the assignment
'''

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
    '''
    If a transcript has a gene_fusion classification, even it may have other classifications,
    only classify the transcript as gene fusion

    input:
        d: Pandas DataFrame, result with duplicated classifications
    output:
        d: Pandas DataFrame, modified the calls of the transcripts that have gene fusions

    '''
    if "gene_fusion" in d["Call"].tolist():
        return d.query("Call == 'gene_fusion'")
    else:
        return d
# gene fusion
with_duplication = with_duplication.groupby("SourceA_Transcript_ID").apply(only_keep_gene_fusion_if_it_exist).reset_index(drop=True)



with_duplication = only_keep_unique_if_unique_and_absent_transcript_both_exist(with_duplication)
with_duplication = only_keep_first_one_if_has_more_than_one_absent_transcript(with_duplication)

# with_duplication = with_duplication.drop_duplicates(subset="SourceA_Transcript_ID", keep="first")
# with_duplication
# fileLocation = "/Users/weilu/Dropbox/genome_algorithms_comp519/project/Challenge_9934185_scoring/cleaned_GMAP_oct30.tsv"
write_to = f"{pre}/{args.name}_post_modification_1.tsv"
df = pd.concat([no_duplication, with_duplication])
pandas_to_tsv(write_to, df)


A_records = SeqIO.to_dict(SeqIO.parse("database/A.transcripts.fasta", "fasta"))
B_records = SeqIO.to_dict(SeqIO.parse("database/B.transcripts.fasta", "fasta"))
seq_records = {A:A_records, B:B_records}
print(pad_with_dash("consider multiple transcript"))
fileLocation = write_to
mySolution = read_result(fileLocation)
mySolution.Category = mySolution.Category.fillna('')
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
    # a is the list of transcripts that best match to fromTranscript.
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
        fromGenome = line["SourceA"]
        from_seq_len = len(seq_records[fromGenome][fromTranscript].seq)
        mrna = database_dic[fromGenome][fromTranscript]
        from_coverage_length = mrna.stop - mrna.start
        # consider if they all have similar length.
        # for x in a:
        #     transcript = seq_records[toGenome][x]
        #     seq_len = len(transcript.seq)
        # print("Multiple transcript", fromTranscript)
        toGenes = [getFromGene(toDB, x) for x in a]
        if len(set(toGenes)) != 1:
            print("They don't have same gene?", fromTranscript, a, toGenes)
        for i, x in enumerate(a):
            new_line = line.copy()
            new_line["Call"] = "multiple_transcript"
            new_line["SourceB_Transcript_ID"] = x
            new_line["SourceB_Gene"] = toGenes[i]
            reverse_call = mySolution.query(f"SourceA_Transcript_ID == '{x}'")["Call"].iloc[0]
            new_line["Category"] += f";reverse_call:{reverse_call};"

            transcript = seq_records[toGenome][x]
            seq_len = len(transcript.seq)
            mrna = database_dic[toGenome][x]
            coverage_length = mrna.stop - mrna.start
            new_line["Category"] += f"from_len={from_seq_len};to_len={seq_len};from_coverage_length={from_coverage_length};to_coverage_length={coverage_length}"
            newSolution = newSolution.append(new_line)
    newSolution = newSolution.reset_index(drop=True)
# if one transcript has both multiple_transcript and others match,
# only keep the unique one.

def only_keep_multiple_transcript(d):
    '''
           if one transcript has both unique_transcript and multiple_transcirpt match,
        only keep the multiple one.
    input:
        d: Pandas DataFrame, result with duplicated classifications
    output:
        d: Pandas DataFrame, modified the calls of the transcripts that have both unique_transcript and multiple_transcript calls
    '''
    if "multiple_transcript" in d["Call"].tolist():
        return d.query("Call == 'multiple_transcript'")
    else:
        return d

mySolution = newSolution.groupby("SourceA_Transcript_ID").apply(only_keep_multiple_transcript).reset_index(drop=True)

def multiple_gene_fusion(d):
    '''
        if transcript A1 is multiple transcript to B1, B2, B3
        if all transcriptB to A1 is gene fusion or absent gene:
            return d
        else:
            remove transcriptB in multiple transcript
        if transcript A1 match to only one transcriptB
        change A1's call to unique
    '''
    multiple_transcript = d.query("Call == 'multiple_transcript'")
    transcriptA = multiple_transcript["SourceA_Transcript_ID"].tolist()
    transcriptB = multiple_transcript["SourceB_Transcript_ID"].tolist()

    #form a dictionary for multiple transcript matchings
    A_to_B_dict = dict()
    for i in range(len(transcriptA)):
        if transcriptA[i] not in A_to_B_dict:
            A_to_B_dict[transcriptA[i]] = [transcriptB[i]]
        else:
            A_to_B_dict[transcriptA[i]].append(transcriptB[i])

    gene_fusion = d.query("Call == 'gene_fusion'")["SourceA_Transcript_ID"].tolist()
    absent_gene = d.query("Call == 'absent_gene'")["SourceA_Transcript_ID"].tolist()

    for transcriptA, transcriptBs in A_to_B_dict.items():
        new_list = []
        for transcriptB in transcriptBs:
            if transcriptB not in gene_fusion and transcriptB not in absent_gene:
                new_list.append(transcriptB)
        if new_list:
            A_to_B_dict[transcriptA] = new_list

    #add all non-multiple transcripts
    d_new = d.query("Call != 'multiple_transcript'")

    #add the gene fusion/absnet gene filtered multiple transcripts
    for transcriptA, transcriptBs in A_to_B_dict.items():
        if len(transcriptBs) == 1:
            transcriptB = transcriptBs[0]
            trAB = d.query("SourceA_Transcript_ID == '%s'"%transcriptA).query("SourceB_Transcript_ID == '%s'"%transcriptB)
            trAB["Call"] = 'unique_transcript'
            d_new = d_new.append(trAB)
        else:
            for transcriptB in transcriptBs:
                trAB = d.query("SourceA_Transcript_ID == '%s'"%transcriptA).query("SourceB_Transcript_ID == '%s'"%transcriptB)
                d_new = d_new.append(trAB)
    return d_new

mySolution = multiple_gene_fusion(mySolution)

write_to = f"{pre}/{args.name}_post_modification_2.tsv"
pandas_to_tsv(write_to, mySolution)


'''
further processing:

- read blast result(only the longest alignment is kept)
- filter out the result with percent identification < 95
- for those queries that target same subject. the one with largest ratio wins.
- for every subject transcript, only keep the one with highest ratio.

'''
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


# df["ratio"] = df["length"]**2/df["qlen"]
# df["ratio"] = df["length"]*(1+ df["length"]/df["qlen"]) * df["pident"]
# df["ratio"] = df["bitscore"] + df["length"]*(1+ df["length"]/df["qlen"])
# df["ratio"] = df["bitscore"]
# df["lengthScore"] = df["length"]*(1+ df["length"]/df["qlen"])
df["lengthScore2"] = df["length"]**2/df["qlen"]/df["slen"]  # df["slen"] part just used for normalization
# I could add an aditional penalty for mismatch.
df["ratio"] = df["bitscore"] * df["lengthScore2"]  # lengthScore2 disfavor larger qlen more than lengthScore.
# filter out pident < 95
df = df.query("pident > 95").reset_index(drop=True)
# for those queries that target same subject. the one with largest ratio wins.
df = df.sort_values(["sseqid", "ratio"], ascending=[True, False]).reset_index(drop=True)
df["rank"] = df.groupby("sseqid")["ratio"].rank(ascending=False)
# for every subject transcript. only keep the one with highest ratio.
df2 = df.groupby("sseqid").head(1).reset_index(drop=True)
# print(df.query("subject_gene == 'geneB24592'"))
# fromDB = dbA
# toDB = dbB
# database_dic[x[10]]
# df2["query_gene"] = df2["qseqid"].apply(lambda x: getFromGene(database_dic[x[10]], x))
# df2["subject_gene"] = df2["sseqid"].apply(lambda x: getFromGene(database_dic[x[10]], x))

write_to = f"{pre}/{args.name}_post_modification_3.tsv"


'''
additional info

'''

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


print(pad_with_dash("cluster multi transcirpt, further classification."))
'''
ref:
B	A	transcriptB32395	transcriptA52049	unique_transcript	0	geneB32395	geneA14875	target_contained
mySolution:
B	A	transcriptB32395	transcriptA125712	multiple_transcript	0	geneB32395	geneA14875	from_len=1809;to_len=3163;from_coverage_length=4574;to_coverage_length=7582
B	A	transcriptB32395	transcriptA86515	multiple_transcript	0	geneB32395	geneA14875	from_len=1809;to_len=4799;from_coverage_length=4574;to_coverage_length=8921
B	A	transcriptB32395	transcriptA99865	multiple_transcript	0	geneB32395	geneA14875	from_len=1809;to_len=4708;from_coverage_length=4574;to_coverage_length=8921
B	A	transcriptB32395	transcriptA7588	multiple_transcript	0	geneB32395	geneA14875	from_len=1809;to_len=5745;from_coverage_length=4574;to_coverage_length=8645
B	A	transcriptB32395	transcriptA139610	multiple_transcript	0	geneB32395	geneA14875	from_len=1809;to_len=5047;from_coverage_length=4574;to_coverage_length=8921
B	A	transcriptB32395	transcriptA126200	multiple_transcript	0	geneB32395	geneA14875	from_len=1809;to_len=5075;from_coverage_length=4574;to_coverage_length=8921
B	A	transcriptB32395	transcriptA74121	multiple_transcript	0	geneB32395	geneA14875	from_len=1809;to_len=4096;from_coverage_length=4574;to_coverage_length=7288
B	A	transcriptB32395	transcriptA52049	multiple_transcript	0	geneB32395	geneA14875	from_len=1809;to_len=1452;from_coverage_length=4574;to_coverage_length=4176

'''
fileLocation = write_to
mySolution = read_result(fileLocation)
call = "multiple_transcript"
subset = mySolution.query(f"Call == '{call}'").reset_index(drop=True)
subset_no_change = mySolution.query(f"Call != '{call}'")
consider_list = subset["SourceA_Transcript_ID"].unique()

new_subset_list = []

def get_coverage_length(Genome, Transcript, database_dic):
    mrna = database_dic[Genome][Transcript]
    return mrna.stop - mrna.start

subset["toCoverage"] = subset.apply(lambda x: get_coverage_length(x.SourceB, x.SourceB_Transcript_ID, database_dic), axis=1)
subset["fromCoverage"] = subset.apply(lambda x: get_coverage_length(x.SourceA, x.SourceA_Transcript_ID, database_dic), axis=1)
# subset["diff"] = abs(subset["fromCoverage"] - subset["toCoverage"])

for i in range(len(consider_list)):
    fromTranscript = consider_list[i]
    a = subset.query(f"SourceA_Transcript_ID=='{fromTranscript}'").reset_index(drop=True)
    z = np.abs(stats.zscore(a["toCoverage"]))
    X = np.append(a["fromCoverage"].iloc[0], a["toCoverage"].values).reshape(-1, 1)
    kmeans = KMeans(n_clusters=2)
    kmeans.fit(X)
    y_kmeans = kmeans.predict(X)
    if sum(y_kmeans == y_kmeans[0]) == 1:
        # only the source it self in one group. then we do nothing.
        new_subset_list.append(a)
        pass
    elif sum(y_kmeans == y_kmeans[0]) == 2:
        # we keep those in the same cluster of the source
        # if only one, then it is unique.
        b = a[y_kmeans[1:] == y_kmeans[0]].reset_index(drop=True)
        b["Call"] = "unique_transcript"
        new_subset_list.append(b)
    else:
        # if more than one, then keep those in the same cluster of the source
        b = a[y_kmeans[1:] == y_kmeans[0]]
        new_subset_list.append(b)
new_subset = pd.concat(new_subset_list).reset_index(drop=True)
mySolution = pd.concat([new_subset.iloc[:,:9], subset_no_change]).reset_index(drop=True)
write_to = f"{pre}/{args.name}_post_modification_4.tsv"
pandas_to_tsv(write_to, mySolution)

# fileLocation = write_to
# mySolution = read_result(fileLocation)
# A_transcript_length_dict = pickle.load(open(f"database/A_transcript_length_dict.pkl","rb"))
# B_transcript_length_dict = pickle.load(open(f"database/B_transcript_length_dict.pkl","rb"))
# transcript_length_dict = {A:A_transcript_length_dict, B:B_transcript_length_dict}

# call = "multiple_transcript"
# subset = mySolution.query(f"Call == '{call}'").reset_index(drop=True)
# subset_no_change = mySolution.query(f"Call != '{call}'")
# consider_list = subset["SourceA_Transcript_ID"].unique()
# subset["toTranscriptLength"] = subset.apply(lambda x: transcript_length_dict[x.SourceB][x.SourceB_Transcript_ID], axis=1)
# subset["fromTranscriptLength"] = subset.apply(lambda x: transcript_length_dict[x.SourceA][x.SourceA_Transcript_ID], axis=1)
# subset["diff"] = abs(subset["fromTranscriptLength"] - subset["toTranscriptLength"])
# new_subset_list = []
# for i in range(len(consider_list)):
#     fromTranscript = consider_list[i]
#     a = subset.query(f"SourceA_Transcript_ID=='{fromTranscript}'").reset_index(drop=True)
#     z = np.abs(stats.zscore(a["toTranscriptLength"]))
#     idx = a["diff"].idxmin()
#     if z[idx] > 2:
#         b = a.iloc[idx:idx+1]
#     else:
#         b = a[z<2]
#     new_subset_list.append(b)
# new_subset = pd.concat(new_subset_list).reset_index(drop=True)
# mySolution = pd.concat([new_subset.iloc[:,:9], subset_no_change]).reset_index(drop=True)
# write_to = f"{pre}/{args.name}_post_modification_4.tsv"
# pandas_to_tsv(write_to, mySolution)
# notes:
# if you want to know which chromosome a transcript belongs, you can do the following.
# transcript = dbA["transcriptA1"]
# transcript.chrom




# if you have the following error
# ModuleNotFoundError: No module named 'pandas.core.internals.managers'; 'pandas.core.internals' is not a package
# weis-mbp-2:LeagueOfLegends weilu$ python3 -m pip install --upgrade pandas