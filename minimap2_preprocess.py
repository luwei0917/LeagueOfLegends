#Maize Pipeline

import subprocess
import os
import gffutils

def run(chromosomes_file, transcripts_file, reference_file):
    outprefix = transcripts_file[0] + "_to_" + chromosomes_file[0]
    sam_file = outprefix + ".sam"
    bed_file = outprefix + ".bed"
    gff_file = outprefix + ".gff"
    trmap_file = "Minimap_trmap_" + outprefix
    db_file = "Minimap_chromosome_" + outprefix

    # minimap2
    print("Running Minimap2")
    subprocess.run(["minimap2", "-ax", "splice", chromosomes_file, transcripts_file, "-o", sam_file]) #"-A3", "-B6", "-O2,16"

    # convert sam file to bed file
    print("Convert .sam file to .bed file")
    command = "paftools.js splice2bed " + sam_file + " > " + bed_file
    os.system(command)

    # convert bed file to gff file
    print("Convert .bed file to .gff file")
    subprocess.run(["python", "bed2gff.py", bed_file, gff_file])

    # gffcompare, generate tracking
    print("Generating tracking file")
    subprocess.run(["gffcompare", gff_file, "-r", reference_file, "-o", outprefix, "-T"])

    # trmap, generate trmap
    print("Generating trmap file")
    subprocess.run(["trmap", reference_file, gff_file, "-o", trmap_file])

    print("Creating gffutils database")
    gffutils.create_db(gff_file, db_file)

chromosomes_file = "A.chromosomes.fasta"
transcripts_file = "B.transcripts.fasta"
reference_file = "A.gff3"

run(chromosomes_file, transcripts_file, reference_file)

chromosomes_file = "B.chromosomes.fasta"
transcripts_file = "A.transcripts.fasta"
reference_file = "B.gff3"

run(chromosomes_file, transcripts_file, reference_file)

chromosomes_file = "A.chromosomes.fasta"
transcripts_file = "C.transcripts.fasta"
reference_file = "A.gff3"

run(chromosomes_file, transcripts_file, reference_file)

chromosomes_file = "C.chromosomes.fasta"
transcripts_file = "A.transcripts.fasta"
reference_file = "C.gff3"

run(chromosomes_file, transcripts_file, reference_file)

