import pandas as pd
from collections import defaultdict

def getFromTranscript(data):
    d = data.split(":")[1]
    d = d.split("|")[1]
    return d
def getToGene(data):
    try:
        gene = data.split("|")[0]
    except:
        gene = "None"
    return gene
def getToTranscript(data):
    try:
        transcript = data.split("|")[1]
    except:
        transcript = "None"
    return transcript

def read_tracking(fileLocation):
    # fileLocation = "/Users/weilu/Research/server/oct_2019/project/oct04/A_to_B/result_A_to_B.tracking"
    data = pd.read_csv(fileLocation, sep="\t", names=["tcons", "xloc", "reference", "category", "q1"])
    data["fromTranscript"] = data["q1"].apply(getFromTranscript)
    data["toTranscript"] = data["reference"].apply(getToTranscript)
    data["toGene"] = data["reference"].apply(getToGene)
    data = data.drop(["reference", "q1"], axis=1)
    # GMAP give a unique name to each transcirpt "transcriptA119577.mrna1"
    data["originalFromTranscript"] = data["fromTranscript"]
    data["fromTranscript"] = data["fromTranscript"].apply(lambda x:x.split(".")[0])
    return data

def read_result(fileLocation):
    return pd.read_csv(fileLocation, sep="\t", names=["SourceA", "SourceB", "SourceA_Transcript_ID", "SourceB_Transcript_ID", "Call",
               "Score", "SourceA_Gene", "SourceB_Gene", "Category"], comment="#", index_col=False)

# write pandas Dataframe to tsv format.
def pandas_to_tsv(fileLocation, df):
    df = df.fillna('')
    df["Score"] = df["Score"].astype(str)
    with open(fileLocation, "w") as out:
        for i in range(len(df)):
            formated_line = "\t".join(df.iloc[i].values) + "\n"
            out.write(formated_line)

def getFromGene(fromDB, fromTranscript):
    try:
        t = fromDB[fromTranscript]
        g = list(fromDB.parents(t))
        if len(g) != 1:
            print("Unexpected!!!!", g, fromTranscript)
        fromGene = g[0].id
    except:
        # print(fromTranscript, " no corresponding gff3 info")
        # return ""
        fromGene = fromTranscript    # use transcript as it's gene name since no corresponding gff3 info
    return fromGene

# def assign_unique(c, fromTranscript, fake_cat_k_toTranscript_list):
#     assign = ""
#     if c == "=" or c == "k":
#         assign = "unique"
#     elif c == "c":
#         # possible 'unique', 'multiple', 'gene_fusion'
#         if fromTranscript in fake_cat_k_toTranscript_list:
#             assign = "unique"
#         else:
#             assign = ""
#     return assign

def assign_unique(c, fromTranscript):
    assign = ""
    if c == "=" or c == "k" or c == "c":
        assign = "unique"
    return assign



def find_map_location(transcript, map_db):
    gene = map_db[transcript]
    return gene[0],gene[3],gene[4],gene[6]

def find_gene_annotation(genome_db,chromosome,start,end,strand):
    return list(genome_db.region(region=(chromosome, start, end), strand=strand, featuretype="gene"))


def assign_gene_fusion(c, assign, fromGenome, toGenome, fromTranscript, toTranscript, fromGene, toGene, genome_db, map_db):
    gene_fusion_list = []
    extra = ""
    if c not in ["=", "c", "k", "m", "n", "j", "e", "o"]:
        return assign, gene_fusion_list, toGene, extra

    chromosome,start,end,strand = find_map_location(fromTranscript, map_db)
    genes = find_gene_annotation(genome_db,chromosome,start,end,strand)
    if len(genes)>1:
        for gene in genes:
            gene_id = gene['ID'][0]
            gene_fusion_list.append(gene_id)
    if len(gene_fusion_list) > 0:
        assign = "gene_fusion"
        toGene = ";".join(gene_fusion_list)

    return assign, gene_fusion_list, toGene, extra

def assign_absent_gene(c, assign, map_db, toDB, fromTranscript, toGene):
    extra = ""
    # if c in ['i', 'j', 'n', 'p', 'm', 'o', 'e', 'y', 's']:
    #     t = map_db[fromTranscript]
    #     # print(fromTranscript, toGene, t.start, dbB[toGene].start)
    #     if (t.start - toDB[toGene].start) < -10000:
    #         # print(fromTranscript)
    #         assign = "absent_gene"
    #     extra = f"\t{t.start}, {toDB[toGene].start}, {t.start - toDB[toGene].start}\t"

    if c == 'u' or c == "x":
        # probably no match at all.
        assign = "absent_gene"
        extra = ""

    return assign, extra

# def assign_absent_gene(c, assign, db_A_to_B, toDB, fromTranscript, toGene):
#     extra = ""
#     if c in ['i', 'j', 'n', 'p', 'm', 'o', 'e', 'y', 's']:
#         t = db_A_to_B[fromTranscript]
#         # print(fromTranscript, toGene, t.start, dbB[toGene].start)
#         if (t.start - toDB[toGene].start) < -10000:
#             # print(fromTranscript)
#             assign = "absent_gene"
#         extra = f"\t{t.start}, {toDB[toGene].start}, {t.start - toDB[toGene].start}\t"

#     if c == 'u' or c == "x":
#         # probably no match at all.
#         assign = "absent_gene"
#         extra = ""

#     return assign, extra

def assign_absent_transcript(c, assign):
    extra = ""
    if assign == "":
        assign = 'absent_transcript'
    return assign, extra

def assign_to_out(assign, fromGenome, toGenome, fromTranscript, toTranscript, fromGene, toGene):
    out = ""
    if assign == 'unique':
        out += f"{fromGenome}\t{toGenome}\t{fromTranscript}\t{toTranscript}\tunique_transcript\t0\t{fromGene}\t{toGene}"
    elif assign == 'absent_transcript':
        out += f"{fromGenome}\t{toGenome}\t{fromTranscript}\t\tabsent_transcript\t0\t{fromGene}\t{toGene}"
    elif assign == 'absent_gene':
        out += f"{fromGenome}\t{toGenome}\t{fromTranscript}\t\tabsent_gene\t0\t{fromGene}\t"
    elif assign == 'gene_fusion':
        out += f"{fromGenome}\t{toGenome}\t{fromTranscript}\t\tgene_fusion\t0\t{fromGene}\t{toGene}"
    else:
        print("unknown assign", assign, fromGenome, toGenome, fromTranscript, toTranscript, fromGene, toGene)
    return out


def pad_with_dash(s):
    return '{s:{c}^{n}}'.format(s=s,n=80,c='-')


def clean_query(string):
    words = string.split(".")
    return words[0][1:]

def clean_ref(string):
    words = string.split()
    return words[5]




def extract_map(lines):
    if not lines:#first lines
        return []
    qtranscript = clean_query(lines[0])
    transcripts = [qtranscript]
    for line in lines[1:]:
        if line[0] in ["c","k","="]:
            rtranscript = clean_ref(line)
            transcripts.append(rtranscript)
    return transcripts

def read_tracking_file(filename):

    tracking_transcript_list = defaultdict(list)
    fi = open(filename,'r')
    line = fi.readline()
    lines = []
    while(line):
        if line[0] == '>':
            transcripts = extract_map(lines)
            lines = [line]
            if len(transcripts)>1:
                tracking_transcript_list[transcripts[0]] = transcripts[1:]

        else:
            lines.append(line)

        line = fi.readline()

    # last transcript
    transcripts = extract_map(lines)
    if len(transcripts)>1:
        tracking_transcript_list[transcripts[0]] = transcripts[1:]


    fi.close()
    return tracking_transcript_list

def read_validation_file(filename = "dev_validation_set.tsv", source = 'A', \
                         calls = ['unique_transcript','multiple_transcript']):
    vali_list = pd.read_csv(filename, skiprows=1, sep="\t",\
        names=["SourceA", "SourceB", "SourceA_Transcript_ID", "SourceB_Transcript_ID", "Call",\
               "Score", "SourceA_Gene", "SourceB_Gene", "Category"])
    matched_gene_list = vali_list.query("("+f"or".join(" Call == '%s' "%i for i in calls) +") "+\
    f"and SourceA == '{source}'")

    return matched_gene_list

def AtoB_bestmatch(AtoB_compare_filename):

    #faketracking_transcript_list = read_fake_file(AtoB_compare_dir+fakeAtoB_compare_filename)
    tracking_transcript_list = read_tracking_file(AtoB_compare_filename)
    #tracking_transcript_list.update(faketracking_transcript_list)

    return tracking_transcript_list

def BtoA_bestmatch(BtoA_compare_filename):
    #faketracking_transcript_list = read_fake_file(BtoA_compare_dir+fakeBtoA_compare_filename)

    tracking_transcript_list = read_tracking_file(BtoA_compare_filename)
    #tracking_transcript_list.update(faketracking_transcript_list)
    return tracking_transcript_list

def match_sen_pre(tracking_transcript_list, matched_gene_list):

    count = 0
    NoMatchList = []
    WrongMatchList = []
    WrondDict = {}
    NoDict = {}
    for index, row in matched_gene_list.iterrows():
        if row['SourceA_Transcript_ID'] in tracking_transcript_list:
            if row['SourceB_Transcript_ID'] in tracking_transcript_list[row['SourceA_Transcript_ID']]:
                count += 1
            else:
                WrongMatchList.append(row['SourceA_Transcript_ID'])
                WrondDict[row['SourceA_Transcript_ID']] = \
                [row['SourceB_Transcript_ID'], tracking_transcript_list[row['SourceA_Transcript_ID']]]
        else:
            NoMatchList.append(row['SourceA_Transcript_ID'])
            NoDict[row['SourceA_Transcript_ID']] = row['SourceB_Transcript_ID']
    print("%d queries"%len(matched_gene_list))
    print("Sensitivity %.2f " % (count/(len(matched_gene_list)-len(WrongMatchList))))
    print("Presition %.2f " % (count/(len(matched_gene_list)-len(NoMatchList))))
    print("Transcripts with no matched gene", NoMatchList)
    print("Transcripts with wrong matched gene", WrongMatchList)
    print("Wrong matches:")
    print("Transcript \t\tGround_Truth \tGMAP")
    for i in WrongMatchList:
        print("%-20s%-15s\t%s"%(i,WrondDict[i][0],WrondDict[i][1]))
    print("No matches:")
    print("Transcript \t\tGround_Truth \tGMAP")
    for i in NoMatchList:
        print("%-20s%-15s\t%s"%(i,NoDict[i],"Not found"))
    return NoMatchList, WrongMatchList

def BtoA_match_sen_pre(tracking_transcript_list, matched_gene_list):

    count = 0
    NoMatchList = []
    WrongMatchList = []
    WrondDict = {}
    NoDict = {}
    print("Wrong matches:")
    print("Transcript \t\tGround_Truth \tGMAP")
    for index, row in matched_gene_list.iterrows():
        if row['SourceA_Transcript_ID'] in tracking_transcript_list:
            if row['SourceB_Transcript_ID'] in tracking_transcript_list[row['SourceA_Transcript_ID']]:
                count += 1
            else:
                WrongMatchList.append(row['SourceA_Transcript_ID'])
                WrondDict[row['SourceA_Transcript_ID']] = \
                [row['SourceB_Transcript_ID'], tracking_transcript_list[row['SourceA_Transcript_ID']]]
                print(row['SourceA_Transcript_ID'],\
                      row['SourceB_Transcript_ID'], tracking_transcript_list[row['SourceA_Transcript_ID']])
        else:
            NoMatchList.append(row['SourceA_Transcript_ID'])
            NoDict[row['SourceA_Transcript_ID']] = row['SourceB_Transcript_ID']
    print("%d queries"%len(matched_gene_list))
    print("Sensitivity %.2f " % (count/(len(matched_gene_list)-len(WrongMatchList))))
    print("Presition %.2f " % (count/(len(matched_gene_list)-len(NoMatchList))))
    print("Transcripts with no matched gene", NoMatchList)
    print("Transcripts with wrong matched gene", WrongMatchList)


    print("No matches:")
    print("Transcript \t\tGround_Truth \tGMAP")
    for i in NoMatchList:
        print("%-20s%-15s\t%s"%(i,NoDict[i],"Not found"))
    return NoMatchList, WrongMatchList


