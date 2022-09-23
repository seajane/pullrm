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
#input_csv = "Odd_ones2.csv"
#window = 500
#output_gbk_path = "genbank_files"
#fasta_folder = "test01"
#Entrez_email = "hbouzek@fredhutch.org"

# for running on cmdline:
    ## python3 pullrms.py -i 2022_16_02_NRM_Octapus_full.csv -w 500 -f 'test02' -o 'gbk_files' -m 'hbouzek@fredhutch.org'

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
if len(xldoc) > 0:
    print("Imported OCTAPUS file.")
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


# find features for genes and cds
def feature_extract (dfrow, gbk_path, working_path, annotation_folder_path):
    # change to the directory with genbank files    
    # for testing rows in dfrow `dfrow = df2.iloc[[0]]`
    os.chdir(gbk_path)
    np = annotation_folder_path + "/" + dfrow['contig_filename']
    fn = dfrow['contig_filename']
    descr = dfrow['genome_name']
    gbid = dfrow['genome_id']
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
        gene1 = [feature for feature in gb.features if feature.type == "gene"]
        if len(cds1) == 0:
            print("No CDS found for %s." % gb.id)
            shutil.copyfile(op, np)
            return
        else:
            for feature in gene1:
                start = feature.location.start.real
                stop = feature.location.end.real
                try:
                    seq = gb.seq[start:stop]
                    gene_seq = ''.join(seq)
                except: 
                    gene_seq = ''
                new_val = pd.DataFrame({'start':[start],'seq':[gene_seq]})
                new_gene = pd.concat([new_gene, new_val], axis = 0, ignore_index = True)
            for feature in cds1:
                start = feature.location.start.real
                stop = feature.location.end.real
                strand = feature.location.strand.real
                if 'product' in feature.qualifiers:
                    product = feature.qualifiers['product'][0]
                elif 'note' in feature.qualifiers:
                    product = feature.qualifiers['note'][0]
                else:
                    product = ''
                if 'translation' in feature.qualifiers:
                    cdsb = feature.qualifiers["translation"][0]
                else:
                    cdsb = None
                    if gene1 is not None:
                        record = Seq(gene_seq)
                        if len(record) %3 ==0:
                            cds2 = record
                        elif (len(record)+1) % 3 == 0:
                            cds2 = record + Seq('N')
                        else:
                            cds2 = record + Seq('NN')
                        if isinstance(cds2, str):
                            cds2 = Seq(cds2)
                        cds2 = cds2.reverse_complement().translate(table = 'Bacterial', to_stop=True) 
                        cdsb = ''.join(cds2)
                    else: 
                        cdsb = Seq('NNNN')
                new_val = pd.DataFrame({\
                                        'description':[descr], \
                                        'gb':[gbid],
                                        'emb':[gb.id], \
                                        'start':[start], 'stop':[stop], \
                                        'strand':[strand], 'product':[product], \
                                        'AA':[cdsb],'contig_filename':[fn]
                                        })
                new_cds = pd.concat([new_cds, new_val], axis = 0, ignore_index = True)
        gene_cds = new_gene.merge(new_cds, on='start', how='left')
        count = len(gene_cds)
        print('%d CDS features collected for %s.' % (count, gb.id))
    os.chdir(working_path)
    return gene_cds

def clean_seq(gene_cds_row, seq_type):
    if isinstance(gene_cds_row, pd.DataFrame) or  isinstance(gene_cds_row, pd.Series):
        if seq_type == "DNA":
            seq1 = gene_cds_row["seq"]
        else:
            seq1 = gene_cds_row["AA"]
        try:
            prod1 = gene_cds_row["product"]
        except KeyError:
            prod1 = ''
        desc1 = gene_cds_row['description']
        gb1 = gene_cds_row['gb']
        emb1 = gene_cds_row['emb']
        start1 = gene_cds_row['start']
        stop1 = gene_cds_row['stop']
    elif isinstance(gene_cds_row, list):
        # use numbers
        if seq_type == "DNA":
            seq1 = gene_cds_row[1]
        else:
            seq1 = gene_cds_row[8]
        try:
            prod1 = gene_cds_row[7]
        except KeyError:
            prod1 = ''
        desc1 = gene_cds_row[2]
        gb1 = gene_cds_row[3]
        emb1 = gene_cds_row[4]
        start1 = gene_cds_row[0]
        stop1 = gene_cds_row[5]
    elif isinstance(gene_cds_row, tuple):
        gene_cds_row = pd.Series(gene_cds_row)[1]
        # use numbers
        if seq_type == "DNA":
            seq1 = gene_cds_row.iloc[1]
        else:
           seq1 = gene_cds_row.iloc[8]
        #print(seq1)
        try:
           prod1 = gene_cds_row[7]
        except KeyError:
           prod1 = ''
        desc1 = gene_cds_row[2]
        gb1 = gene_cds_row[3]
        emb1 = gene_cds_row[4]
        start1 = gene_cds_row[0]
        stop1 = gene_cds_row[5]
    desc1 = re.sub(r"\d+ ","", desc1)
    id1 = desc1.strip().replace(" ", "_")
    name1 = 'gb|' + gb1 + '|emb|' + emb1
    desc2 = prod1 + '_' + str(start1) + '_' + str(stop1).replace(" ", "_")
    sr = SeqRecord(Seq(seq1), id = id1, description = desc2, name = name1)
    return sr

# 4. Pull Restriction-Methyltransferase using OCTAPUS coordinates
## OCTAPUS coords are not always correct, may be shifted
## Directionality can be inferred from start/stop

def store_seq(dataframe, seq_type, output_csv, output_fasta, output_fail_csv, \
              genbank_folder_path = grp, home_path = rp, \
              annotation_folder_path = prp):
    #colnms = list(dataframe.columns.values)
    with open(output_fasta, "w") as output_handle: # send output to a fasta file
        for index, row in dataframe.iterrows():
            ldf = len(dataframe)
            print("Feature File " + str(index + 1) + " of " + str(ldf))
            ftx = feature_extract(row, genbank_folder_path, home_path, annotation_folder_path)
            # find RM - reduce list to putative start and stop -/+ window respectively
            if ftx is None or len(ftx) == 0:
                print("Storing information for " + row['contig_filename'])
                # want a more specific error recorded
                if seq_type == "AA":
                    output_fail_csv = pd.concat([output_fail_csv, pd.DataFrame(row).T], axis = 0, ignore_index = True)
            else: 
                # find RM
                ftx = ftx.assign(close_start = abs(ftx['start']-int(row['search_start'])))
                mask = ftx.strand.isnull()
                column_name = 'close_start'
                ftx.loc[mask, column_name] = None
                ftx.loc[ftx.strand.isnull(), 'close_start'] = None
                minstart = min(ftx['close_start'])
                ftx = ftx[ftx['close_start'] == minstart]
                for row1 in ftx.iterrows():
                    # clean up an annotate
                    if seq_type == "DNA":
                        sr = clean_seq(row1, seq_type = "DNA")
                    else:
                        sr = clean_seq(row1, seq_type = "AA")
                    # track rows added to fasta file
                    if seq_type == "AA":
                        output_csv = pd.concat([output_csv, ftx], axis = 0, ignore_index = True)
                    SeqIO.write(sr, output_handle, 'fasta')
        #write out csv
        if seq_type == "AA":
            nm_pass = annotation_folder_path + "/" + "RM_pass.csv"
            nm_fail = annotation_folder_path + "/" + "RM_fail.csv"
            output_csv.to_csv(nm_pass)
            output_fail_csv.to_csv(nm_fail)
    return


# Create a blank fasta file and a blank data frame
gene_rm_seq = orp + "/" + input_csv[:-4] + "_gene_rm.fasta"
cds_rm_seq = orp + "/" + input_csv[:-4] + "_cds_rm.fasta"
gene_sens_seq = orp + "/" + input_csv[:-4] + "_gene_sens.fasta"
cds_sens_seq = orp + "/" + input_csv[:-4] + "_cds_sens.fasta"
phylo_only_aa = pd.DataFrame()
phylo_only_dna = pd.DataFrame()
sens_aa = pd.DataFrame()
sens_dna = pd.DataFrame()
noseq_aa_rm = pd.DataFrame()
noseq_dna_rm = pd.DataFrame()
noseq_aa_sens= pd.DataFrame()
noseq_dna_sens = pd.DataFrame()
specificity_lookup = pd.DataFrame()

store_seq(dataframe = df2, seq_type = "DNA", output_csv = phylo_only_dna, output_fail_csv = noseq_dna_rm, output_fasta = gene_rm_seq)
store_seq(dataframe = df2, seq_type = "AA", output_csv = phylo_only_aa, output_fail_csv = noseq_aa_rm, output_fasta = cds_rm_seq)
           
# 5. Pull Specificity subunit using directionality and rmtase coords
## Usually within 4-7 base pairs
## Remove transposon hits
        
def store_sens(dataframe, seq_type, output_csv, output_fasta, output_fail_csv, \
              genbank_folder_path = grp, home_path = rp, \
              annotation_folder_path = prp):
    with open(output_fasta, "w") as output_handle: # send output to a fasta file
        for index, row in dataframe.iterrows():
            ldf = len(dataframe)
            print("Feature File " + str(index + 1) + " of " + str(ldf))
            ftx = feature_extract(row, genbank_folder_path, home_path, annotation_folder_path)
            # find RM - reduce list to putative start and stop -/+ window respectively
            if ftx is None or len(ftx) == 0:
                print("Storing information for " + row['contig_filename'])
                # want a more specific error recorded
                if seq_type == "AA":
                    output_fail_csv = pd.concat([output_fail_csv, pd.DataFrame(row).T], axis = 0, ignore_index = True)
            else: 
                # find Sens
                lft = len(ftx)
                ftx = ftx.assign(close_start = abs(ftx['start']-int(row['search_start'])))
                ftx = ftx.assign(indx = list(range(0,lft)))
                minstart = min(ftx['close_start'])
                ftx['closest'] = ftx['close_start'] == minstart
                f2 = ftx[ftx['close_start'] == minstart]
                f2sens = int(f2['indx']) + int(f2['strand'])
                try:
                    ftx  = ftx[ftx['indx']== f2sens]
                except:
                    print("Contig doesn't contain Specificity subunit.")
                    output_csv = pd.concat([specificity_lookup, ftx], axis = 0, ignore_index = True)
                    return output_csv
                # clean up an annotate
                for row1 in ftx.iterrows():
                    # clean up an annotate
                    if seq_type == "DNA":
                        sr = clean_seq(row1, seq_type = "DNA")
                    else:
                        sr = clean_seq(row1, seq_type = "AA")
                    # track rows added to fasta file
                    if seq_type == "AA":
                        output_csv = pd.concat([output_csv, ftx], axis = 0, ignore_index = True)
                    SeqIO.write(sr, output_handle, 'fasta')
    #write out csv
    nm_pass = annotation_folder_path + "/" + "Sens_pass.csv"
    nm_fail = annotation_folder_path + "/" + "Sens_fail.csv"
    output_csv.to_csv(nm_pass)
    output_fail_csv.to_csv(nm_fail)

    
            
store_sens(dataframe = df2, seq_type = "DNA", output_csv = sens_dna, output_fail_csv = noseq_dna_sens, output_fasta = gene_sens_seq)
store_sens(dataframe = df2, seq_type = "AA", output_csv = sens_aa, output_fail_csv = noseq_aa_sens, output_fasta = cds_sens_seq)


# 6. Create a FASTA file of DNA and AA for results
## One for RMtase and one for Specificity

#7. retain name, taxonomy, length of sequence, and accession in a csv file
## Can check for outliers

