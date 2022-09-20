#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 13:00:52 2022

@author: hbouzek
"""
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import urllib.error
import pandas as pd
import os, re, sys, getopt
import shutil


# user functions

## the csv must have columns named genome_id, genome_name, contig_name, contig_start, contig_end
### the csv is output from OCTAPUS
### columns and column names: 
### genome_id - the NCBI identifier for the organism's sequence
### genome_name - Genus species name of the organism
### contig_name - Genbank identifier used to locate the genbank file
### contig_start - start of the contig provided to determine neighborhood of gene
### contig_end - end of the contig provided to determine neighborhood of the gene

# Function that allows the python code to accept arguments
# For testing:
#print(input_csv,window,gene_folder,output_gbk_path,Entrez_email)
#input_csv = "odd_gbk.csv"
#window = 1000
#output_gbk_path = "genbank_files"
#fasta_folder = "test01"
#Entrez_email = "hbouzek@fredhutch.org"

def dataset_input(argv):
    #get inputs
    input_csv = ''
    fasta_folder = ''
    Entrez_email = ''
    output_gbk_path = ''
    window = ''
    try:
        opts, args = getopt.getopt(argv, "hi:w:f:o:m:", ["help", "input_csv=","window=", "fasta_folder=","output_gbk_path=", "entrez_email="])
    except getopt.GetoptError:
        print('Error! Usage: python pullrm.py -i <input OCTAPUS csv> -w <contig_search_window> -f <fasta_folder> -o <output gbk path> -m <email for ncbi search>' )
        print('   or: python pullrm.py --input_csv <OTU csv> --window <contig_search_window> --fasta_folder <fasta output folder> --output <output gbk path> --email <email for ncbi search>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('Usage: python pullrm.py -i <input OCTAPUS csv> -w <contig_search_window> -f <fasta output path > -o <output gbk path> -m <email for ncbi search>' )
            print('   or: python pullrm.py --input_csv <OTU csv> --window <contig_search_window> --fasta_folder <fasta output folder> --output <output gbk path> --email <email for ncbi search>')
            sys.exit()
        elif opt in ("-i", "--input_csv"):
            input_csv = arg
        elif opt in ("-w", "--window"):
            window = arg
        elif opt in ("-f", "--fasta_folder"):
            fasta_folder = arg
        elif opt in ("-o", "--output"):
            output_gbk_path = arg
        elif opt in ("-m", "--email"):
            Entrez_email = arg
    return input_csv,int(window),fasta_folder,output_gbk_path,Entrez_email

# take inputs
input_csv,window,fasta_folder,output_gbk_path,Entrez_email = dataset_input(sys.argv[1:])


########

# 1. Take in a file of OCTAPUS results
## csv upload

# 1a. Create file folders if they don't already exist
rp = os.getcwd()
# Create a new folder for Genbank file
grp = rp + "/" + output_gbk_path
if not os.path.exists(grp):
    os.mkdir(grp)
# Create a new folder for Fasta file
orp = rp + "/" + fasta_folder
if not os.path.exists(orp):
    os.mkdir(orp)
prp = rp + "/" + "for_annotation"
if not os.path.exists(prp):
    os.mkdir(prp)
    
# 1b. Grab file names, locations, and strand, from the OCTAPUS csv file
    xldoc = pd.read_csv(input_csv)
    # Contig name has the genbank file names
    df = xldoc[['genome_id', 'genome_name', 'strand', 'contig_name', 'contig_start', 'contig_end']]
    fn = df.loc[:,('contig_name')].map(str)+ '.gbk'
    df = df.assign(contig_filename = fn)
    dfs = df.loc[:,('contig_start', 'contig_end')].min(axis = 1)
    df = df.assign(search_start = dfs)
    dfe = df.loc[:,('contig_start', 'contig_end')].max(axis = 1)
    df = df.assign(search_end = dfe)

##########

# 2. Check for genebank file in genebank folder
## if none, download from NCBI

# 2a. Check to see if files are in the folder
fileList = df['contig_filename']
contents = os.listdir(grp)
# make sure fileList from Octapus and contents from genebank folder are the same
indices = [i for i, element in enumerate(fileList) if element in contents]
indices2 = [i for i, element in enumerate(fileList) if element not in contents]
df2 = df.loc[indices]
df2 = df2.reset_index(drop = True)
print("Found %d files in folder" % (len(df2)))

# 2b. for files not found within the folder - grab them from NCBI

# Function to download gbk files listed in Octapus
def gbk_download(query,output_path):
    filename = output_path + '/'+ query + '.gbk'
    # Downloading...
    print("Downloading :",query)
    net_handle = Entrez.efetch(db="nuccore",id=query,rettype="gb", retmode="text")
    out_handle = open(filename, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    print ("Saved")
    return

# Query Entrez
df3 = df.loc[indices2]

try:
    if len(df3) > 0:
        query_list = df3.loc[:, 'contig_name'].to_list()
        Entrez.email = Entrez_email
        n=1
        for each_query in query_list:
            gbk_download(each_query, grp)
            print ("count ",n)
            print (' ')
            n+=1
            contents = os.listdir(grp)
            indices = [i for i, element in enumerate(fileList) if element in contents]
            df2 = df.loc[indices]
            df2 = df2.reset_index(drop = True)
except urllib.error.HTTPError as err:
    print(err.code)
    print("Server error, try again later.")
    
############

# 3. From each genbank file, pull CDS
## if no CDS - annotate genome

# Extract CDS
def extract_cds (dfrow, gb, gbk_path, shifted_contig_path):
    # dfrow = row in the data frame containing contig filenames
    # gb = 
    # gbk_path = genebank file path
    # shifted_contig_path = shifted contig path
    op = gbk_path + "/" + dfrow['contig_filename']
    np = shifted_contig_path + "/" + dfrow['contig_filename']
    new_cds = pd.DataFrame()
    # setting col width any lower will cause sequence to be truncated
    pd.set_option("max_colwidth", 9999)
    cds1 = [feature for feature in gb.features if feature.type == "CDS"]
    gene1 = [feature for feature in gb.features if feature.type == "gene"]
    if len(cds1) == 0:
        shutil.copyfile(op, np)
    else:
        for feature in cds1:
            position = feature.location
            start = position.start.real
            stop = position.end.real
            strand = position.strand.real
            try:
                cds = feature.qualifiers["translation"][0]
            except KeyError:
                cds = ""
            try:
                product = feature.qualifiers['product']
            except KeyError:
                product = ""
            new_val = pd.DataFrame({'description':[dfrow['genome_name']],'gb':[dfrow['genome_id']],'emb':[gb.id], 'start':[start], \
            'stop':[stop], 'strand':[strand], 'product':[product], 'cds':[cds], 'contig_filename':[dfrow['contig_filename']]})
            new_cds = pd.concat([new_cds, new_val], axis = 1, ignore_index = True)
    return new_cds

# find features for genes and cds
def feature_extract (dfrow, gbk_path, working_path, annotation_folder_path):
    # change to the directory with genbank files    
    # for testing rows in dfrow `dfrow = df2.iloc[[0]]`
    os.chdir(gbk_path)
    np = annotation_folder_path + "/" + dfrow['contig_filename']
    fn = dfrow['contig_filename']
    fn = fn.item().strip()
    descr = dfrow['genome_name']
    descr = descr.item().strip()
    gbid = dfrow['genome_id']
    gbid = gbid.item().strip()
    op = "./" + fn
    new_cds = pd.DataFrame()
    new_gene = pd.DataFrame()
    # setting col width any lower will cause sequence to be truncated
    pd.set_option("max_colwidth", 9999)
    try:
        gb = SeqIO.read(fn, 'genbank') #open genbank file
        print("Genbank file read for " + fn)
    except ValueError:
        gb = None
        print("No records found for " + fn)
        exit
    # Parse features
    if len(gb.features) < 1:
        print("No CDS found for %s." % gb.id)
        # Send files with no CDS to "for_annotation" folder
        shutil.copyfile(op, np)
        exit
    else:
        cds1 = [feature for feature in gb.features if feature.type == "CDS"]
        print(cds1)
        gene1 = [feature for feature in gb.features if feature.type == "gene"]
        if len(cds1) == 0:
            print("No CDS found for %s." % gb.id)
            shutil.copyfile(op, np)
            return
        else:
            for feature in gene1:
                start = feature.location.start.real
                stop = feature.location.end.real
                seq = gb.seq[start:stop]
                seq = ''.join(seq)
                new_val = pd.DataFrame({'start':[start],'seq':[seq]})
                new_gene = pd.concat([new_gene, new_val], axis = 0, ignore_index = True)
                print(new_gene)
            for feature in cds1:
                start = feature.location.start.real
                print(start)
                stop = feature.location.end.real
                strand = feature.location.strand.real
                try:
                    cds = feature.qualifiers["translation"][0]
                except KeyError:
                    record = ftx[ftx['start'] == start]['seq']
                    record = Seq(record.item())
                    if len(record) %3 ==0:
                        cds = record
                    elif (len(record)+1) %3 ==0:
                        cds = record + Seq('N')
                    else:
                        cds = record + Seq('NN')
                    cds = cds.reverse_complement().translate(table='Bacterial', to_stop=True) 
                    cds = ''.join(cds)
                product = feature.qualifiers['product'][0]
                new_val = pd.DataFrame({\
                                        'description':[descr], \
                                        'gb':[gbid],
                                        'emb':[gb.id], \
                                        'start':[start], 'stop':[stop], \
                                        'strand':[strand], 'product':[product], \
                                        'cds':[cds],\
                                        'contig_filename':[fn]
                                        })
                new_cds = pd.concat([new_cds, new_val], axis = 0, ignore_index = True)
                print(new_cds)
            
        gene_cds = new_gene.merge(new_cds, on='start', how='left')
        count = len(gene_cds)
        print('%d CDS features collected for %s.' % (count, gb.id))
    os.chdir(working_path)
    return gene_cds

def clean_seq(found_gene_row):
    seq1 = found_gene_row['cds']
    regex = re.compile('[^A-Z]')
    seq1 = regex.sub('', seq1)
    gb1 = re.sub(r"\d+ ","", found_gene_row['gb'])
    emb1 = re.sub(r"\d+ ","", found_gene_row['emb'])
    desc1 = re.sub(r"\d+ ","", found_gene_row['description'])
    try:
        prod1 = re.sub(r"\d+ ","", found_gene_row['product'][0])
    except KeyError:
        prod1 = ''
    start1 = found_gene_row['start']
    stop1 = found_gene_row['stop']
    id1 = desc1.strip().replace(" ", "_")
    name1 = 'gb|' + gb1.strip() + '|emb|' + emb1.strip()
    desc2 = prod1 + '_' + str(start1) + '_' + str(stop1)
    sr = SeqRecord(Seq(seq1), id = id1, description = desc2, name = name1)
    return sr

# Create a blank fasta file and a blank data frame
gene_rm_seq = orp + "/" + input_csv[:-4] + "gene_rm.fasta"
cds_rm_seq = orp + "/" + input_csv[:-4] + "cds_rm.fasta"
gene_sens_seq = orp + "/" + input_csv[:-4] + "gene_sens.fasta"
cds_sens_seq = orp + "/" + input_csv[:-4] + "cds_sens.fasta"
phylo_only = pd.DataFrame()
noseq = pd.DataFrame()
with open(cds_rm_seq, "w") as output_handle: # send output to a fasta file
    for index, row in df2.iterrows():
        ldf = len(df2)
        print("RMtase CDS File " + str(index) + " of " + str(ldf))
        ftx = feature_extract(row, grp, rp, prp)
        # find RM - reduce list to putative start and stop -/+ window respectively
        if ftx is None or len(ftx) == 0:
            print("Cannot find " + row['contig_filename'])
            # want a more specific error recorded
            noseq = pd.concat([noseq, row], axis = 1, ignore_index = True)
        else: 
           # find RM
           ftx = ftx.assign(close_start = abs(ftx['start']-int(row['search_start'])))
           minstart = min(ftx['close_start'])
           ftx = ftx[ftx['close_start'] == minstart]
           sr = clean_seq(ftx['cds'])
           # track rows added to fasta file
           phylo_only = pd.concat([phylo_only, row], axis = 1, ignore_index = True)
           SeqIO.write(sr, output_handle, 'fasta')
           
with open(gene_rm_seq, "w") as output_handle: # send output to a fasta file
    for index, row in df2.iterrows():
        ldf = len(df2)
        print("RMtase Gene File " + str(index) + " of " + str(ldf))
        ftx = feature_extract(row, grp, rp, prp)
        # find RM - reduce list to putative start and stop -/+ window respectively
        if ftx is None or len(ftx) == 0:
            print("Cannot find " + row['contig_filename'])
        else: 
           # find RM
           ftx = ftx.assign(close_start = abs(ftx['start']-int(row['search_start'])))
           minstart = min(ftx['close_start'])
           ftx = ftx[ftx['close_start'] == minstart]
           sr = clean_seq(ftx['gene'])
           # track rows added to fasta file
           SeqIO.write(sr, output_handle, 'fasta')
       
with open(gene_sens_seq, "w") as output_handle: # send output to a fasta file
    for index, row in df2.iterrows():
        ldf = len(df2)
        print("Sensitivity Subunit Gene File " + str(index) + " of " + str(ldf))
        ftx = feature_extract(row, grp, rp, prp)
        # find RM - reduce list to putative start and stop -/+ window respectively
        if ftx is None or len(ftx) == 0:
            print("Cannot find " + row['contig_filename'])
        else: 
           # find Sens
           ftx = ftx.assign(close_start = abs(ftx['start']-int(row['search_start'])))
           minstart = min(ftx['close_start'])
           ftx['closest'] = ftx['close_start'] == minstart
           
           sr = clean_seq(ftx['gene'])
           # track rows added to fasta file
           SeqIO.write(sr, output_handle, 'fasta')         
            


# 4. Pull Restriction-Methyltransferase using OCTAPUS coordinates
## OCTAPUS coords are not always correct, may be shifted
## Directionality can be inferred from start/stop

# 5. Pull Specificity subunit using directionality and rmtase coords
## Usually within 4-7 base pairs
## Remove transposon hits

# 6. Create a FASTA file of DNA and AA for results
## One for RMtase and one for Specificity

#7. retain name, taxonomy, length of sequence, and accession in a csv file
## Can check for outliers