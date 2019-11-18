import gffutils
import pandas as pd
from helper_functions import *
# creating database
# dbA = gffutils.create_db("database/A.gff3", "database/chromosomeA")
# dbB = gffutils.create_db("database/B.gff3", "database/chromosomeB")


def read_blastn(fileLocation, dbA, dbB):
    data = pd.read_csv(fileLocation, sep="\t", names=["qseqid", "sseqid",
                                            "pident", "length", "mismatch",
                                            "gapopen", "qstart", "qend", "sstart",
                                            "send", "evalue", "bitscore", "slen", "qlen"])
    # each query-subject pair, I only keep the first(the longest alignment)
    data["query_gene"] = data["qseqid"].apply(lambda x: getFromGene(dbA, x))
    data["subject_gene"] = data["sseqid"].apply(lambda x: getFromGene(dbB, x))
    data.to_pickle(f"{fileLocation}.pkl")

# print("creating database chromosomeC")
# dbC = gffutils.create_db("database/C_fixerror.gff3", "database/chromosomeC")
# print("creating database GMAP_chromosome_A_to_C")
# db_A_to_C = gffutils.create_db("database/AtoCmap_match_gene.gff3", "database/GMAP_chromosome_A_to_C")
# print("creating database GMAP_chromosome_C_to_A")
# db_C_to_A = gffutils.create_db("database/CtoAmap_match_gene.gff3", "database/GMAP_chromosome_C_to_A")

print("read and processing blastn results")
dbA = gffutils.FeatureDB("database/chromosomeA")
dbC = gffutils.FeatureDB("database/chromosomeC")
read_blastn("database/blastn_A_to_C", dbA, dbC)
read_blastn("database/blastn_C_to_A", dbC, dbA)


# dbA = gffutils.FeatureDB("database/chromosomeA")
# dbB = gffutils.FeatureDB("database/chromosomeB")
# read_blastn("database/blastn_A_to_B", dbA, dbB)
# read_blastn("database/blastn_B_to_A", dbB, dbA)
