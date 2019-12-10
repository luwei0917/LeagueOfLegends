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
subprocess.run(["gmap_build", "-d", f"database/Genome{source}", f"database/{source}.chromosomes.fasta"])

print("Mapping with GMAP...")
command = f"gmap -d Genome{target} database/{target}.transcripts.fasta -D database/ -f 2 > database/{source}to{target}_match_gene.gff3"
os.system(command)
command = f"gmap -d Genome{source} database/{source}.transcripts.fasta -D database/ -f 2 > database/{target}to{source}_match_gene.gff3"
os.system(command)

print("Generating tracking file using gffcompare...")
subprocess.run(["gffcompare", f"database/{source}to{target}_match_gene.gff3", "-r", f"database/{target}.gff3", "-o", f"database/{source}to{target}GMAP", "-T"])
subprocess.run(["gffcompare", f"database/{target}to{source}_match_gene.gff3", "-r", f"database/{source}.gff3", "-o", f"database/{target}to{source}GMAP", "-T"])

print("Running blastn...")
command = f"blastn -subject database/{target}.transcripts.fasta> -query database/{source}.transcripts.fasta -outfmt 6 > database/blastn_{source}_to_{target}"
os.system(command)
command = f"blastn -subject database/{source}.transcripts.fasta> -query database/{target}.transcripts.fasta -outfmt 6 > database/blastn_{target}_to_{source}"
os.system(command)

# creating database
dbA = gffutils.create_db(f"database/{source}.gff3", f"database/chromosome{source}")
dbB = gffutils.create_db(f"database/{target}.gff3", f"database/chromosome{target}")

print("Read and processing blastn results...")
read_blastn(f"database/blastn_{source}_to_{target}", dbA, dbB)
read_blastn(f"database/blastn_{target}_to_{source}", dbB, dbA)

print("Running trmap...")
subprocess.run(["trmap", "-o", f"database/{source}to{target}GMAP.gff3", f"database/{target}.gff3", f"database/{source}to{target}_match_gene.gff3"])
subprocess.run(["trmap", "-o", f"database/{target}to{source}GMAP.gff3", f"database/{source}.gff3", f"database/{target}to{source}_match_gene.gff3"])