import argparse
import subprocess
import os
import gffutils
import pandas as pd
from helper_functions import *

def read_blastn(fileLocation, dbA, dbB):
    data = pd.read_csv(fileLocation, sep="\t", names=["qseqid", "sseqid",
                                            "pident", "length", "mismatch",
                                            "gapopen", "qstart", "qend", "sstart",
                                            "send", "evalue", "bitscore", "slen", "qlen"])
    # each query-subject pair, I only keep the first(the longest alignment)
    data["query_gene"] = data["qseqid"].apply(lambda x: getFromGene(dbA, x))
    data["subject_gene"] = data["sseqid"].apply(lambda x: getFromGene(dbB, x))
    data.to_pickle(f"{fileLocation}.pkl")


parser = argparse.ArgumentParser(description="Project pipeline")
parser.add_argument("-s", "--source", type=str)
parser.add_argument("-t", "--target", type=str)
parser.add_argument("-o", "--output", type=str)
args = parser.parse_args()

source = args.source
target = args.target
output = args.output


print("Building GMAP database...")
subprocess.run(["gmap_build", "-d", f"database/Genome{target}", f"database/{target}.chromosomes.fasta"])

print("Mapping with GMAP...")
command = f"gmap -d Genome{target} database/{target}.transcripts.fasta -D database/ -f 2 > database/{source}to{target}_match_gene.gff3"
os.system(command)

print("Generating tracking file using gffcompare...")
subprocess.run(["gffcompare", f"database/{source}to{target}_match_gene.gff3", "-r", f"database/{target}.gff3", "-o", f"result_{source}_to_{target}", "-T"])

print("Running blastn...")
command = f"blastn -subject {target}.transcripts.fasta> -query {source}.transcripts.fasta -outfmt 6 > blastn_{source}_to_{target}.dms"
os.system(command)

#print("Read and processing blastn results...")
#dbA = gffutils.FeatureDB("database/chromosomeA")
#dbC = gffutils.FeatureDB("database/chromosomeC")
#read_blastn("database/blastn_A_to_C", dbA, dbC)
#read_blastn("database/blastn_C_to_A", dbC, dbA)

print("Running trmap...")
subprocess.run(["trmap", "-o", f"database/{source}to{target}GMAP.gff3", f"database/{target}.gff3", f"{source}to{target}map_match_gene.gff3"])