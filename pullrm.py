#!/bin/python
import os
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqFeature import BeforePosition, AfterPosition
import pandas as pd
import re

# first argument is the name of the csv file containing the file names and locations

## the csv must have columns named genome_id, genome_name, contig_name, contig_start, contig_end
### genome_id - the NCBI identifier for the organism's sequence
### genome_name - Genus species name of the organism
### contig_name - Genbank identifier used to locate the genbank file
### contig_start - start of the contig provided to determine neighborhood of gene
### contig_end - end of the contig provided to determine neighborhood of the gene
# another argument will be the folder containing genbank files
# another argument will be the window for searching from the contig_start and contig_end

import sys,getopt
def dataset_input(argv):
    #get inputs
    input_csv = ''
    gene_folder = ''
    Entrez_email = ''
    output_gbk_path = ''
    window = ''
    try:
        opts, args = getopt.getopt(argv, "hi:w:g:o:m:", ["help", "input_csv=","window=", "gene_folder=","output_gbk_path=", "entrez_email="])
    except getopt.GetoptError:
        print('Error! Usage: python pullrm.py -i <input OCTAPUS csv> -w <contig_search_window> -g <gene_folder> -o <output gbk path> -m <email for ncbi search>' )
        print('   or: python pullrm.py --input_csv <OTU csv> --window <contig_search_window> --gene_folder <gene_folder> --output <output gbk path> --email <email for ncbi search>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('Usage: python pullrm.py -i <input OCTAPUS csv> -w <contig_search_window> -g <gene_folder> -o <output gbk path> -m <email for ncbi search>' )
            print('   or: python pullrm.py --input_csv <OTU csv> --window <contig_search_window> --gene_folder <gene_folder> --output <output gbk path> --email <email for ncbi search>')
            sys.exit()
        elif opt in ("-i", "--input_csv"):
            input_csv = arg
        elif opt in ("-w", "--window"):
            window = arg       
        elif opt in ("-g", "--gene_folder"):
            gene_folder = arg            
        elif opt in ("-o", "--output"):
            output_gbk_path = arg
        elif opt in ("-m", "--email"):
            Entrez_email = arg            
    return input_csv,window,gene_folder,output_gbk_path,Entrez_email

# take inputs 
input_csv,window,gene_folder,output_gbk_path,Entrez_email = dataset_input(sys.argv[1:])
#print(input_csv,window,gene_folder,output_gbk_path,Entrez_email)
# Grab file names and locations
xldoc = pd.read_csv(input_csv)

# contig name has the genbank file names
df = xldoc[['genome_id', 'genome_name', 'contig_name', 'contig_start', 'contig_end']]
fn = df.loc[:,('contig_name')].map(str)+ '.gbk'
df = df.assign(contig_filename = fn)
dfs = df.loc[:,('contig_start', 'contig_end')].min(axis =1)
df = df.assign(search_start = dfs)
dfe = df.loc[:,('contig_start', 'contig_end')].max(axis =1)
df = df.assign(search_end = dfe)

# create a new genbank file and folder
gene_folder = gene_folder
if not os.path.exists(gene_folder):
    os.mkdir(gene_folder)

# Check to see if files are in the folder
fileList = df['contig_filename']
contents = os.listdir('./genebank_files')
indices = [i for i, element in enumerate(fileList) if element in contents]
indices2 = [i for i, element in enumerate(fileList) if element not in contents]
df2 = df.loc[indices]
df2 = df2.reset_index(drop = True)
print("Found %d files in folder" % (len(df2)))
df3 = df.loc[indices2]

if len(df3) > 0:
    query_list = df3.loc[:, 'contig_name'].to_list()
    # Hanrui's code for downloading genbank files
    #############################################################################
    Entrez.email = Entrez_email

    def gbk_download(query,output_path):
        filename = output_path+'/'+query+'.gbk'
        print("query :",query)
        # Downloading...
        print("Downloading ",query)
        net_handle = Entrez.efetch(db="nucleotide",id=query,rettype="gb", retmode="text")
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print ("Saved")
        return

    output_path = output_gbk_path
    n=1
    for each_query in query_list:
        gbk_download(each_query,output_path)
        print ("count ",n)
        print (' ')
        n+=1
    ###############################################################################
    # now remake df2 with new files
    contents = os.listdir('./genebank_files')
    indices = [i for i, element in enumerate(fileList) if element in contents]
    df2 = df.loc[indices]
    df2 = df2.reset_index(drop = True)
# extract CDS

def feature_extract (dfrow):
    gb = SeqIO.read(dfrow['contig_filename'], 'genbank')
    lgb = len(gb.features)
    if lgb < 2:
        print("No CDS found for %s." % gb.id)
        new_cds = None
        # put files into a new folder
        no_cds_folder = "./for_prokka"
        if not os.path.exists(no_cds_folder):
            os.mkdir(no_cds_folder)
        op = "./" + dfrow['contig_filename']
        np = "./" + no_cds_folder + "/" + dfrow['contig_filename']
        os.rename(op, np)
    else:
        count = 1
        new_cds = pd.DataFrame()
        pd.set_option("max_colwidth", 9999)
        for feature in gb.features:
            if feature.type=='CDS':
                position = feature.location
                start = position.start.real
                stop = position.end.real
                strand = position.strand.real
                cds = feature.qualifiers["translation"][0]
                product = feature.qualifiers['product']
                new_val = pd.DataFrame({'description':[dfrow['genome_name']],\
                'gb':[dfrow['genome_id']],'emb':[gb.id], 'start':[start], \
                'stop':[stop], 'strand':[strand], 'product':[product], \
                'cds':[cds], 'contig_filename':[dfrow['contig_filename']]})
                new_cds = new_cds.append(new_val, ignore_index = True)
                count+=1
        print('%d CDS features collected for %s' % (count, gb.id))
    return new_cds


# use feature_extract parse through df2
rp = os.getcwd()
os.chdir(rp+'/genebank_files')
with open("./01_Subset_genbank/fullcds.fasta", "w") as output_handle:
    for index, row in df2.iterrows():
        ftx = feature_extract(row)
        if ftx is not None:
            # find the region within 1000 bases of start and/or 1000 bases of stop
            # HW: change 1000 to outside argument 'window'
            find_gene1 = ftx[ftx['stop'] < (row['search_end'] + window)]
            find_gene = find_gene1[find_gene1['start'] > (row['search_start'] - window)]
            seq1 = str(find_gene['cds'])
            regex = re.compile('[^A-Z]')
            seq1 = regex.sub('', seq1)
            gb1 = re.sub(r"\d+ ","", find_gene['gb'].to_string())
            emb1 = re.sub(r"\d+ ","", find_gene['emb'].to_string())
            desc1 = re.sub(r"\d+ ","", find_gene['description'].to_string())
            prod1 = re.sub(r"\d+ ","", find_gene['product'].to_string())
            start1 = re.sub(r"\d+ ","", find_gene['start'].to_string())
            stop1 = re.sub(r"\d+ ","", find_gene['stop'].to_string())
            id1 = desc1.strip() + "_" + ftx['contig_filename']
            name1 = 'gb|'+ gb1.strip() + '|emb|' + emb1.strip()
            desc2 = prod1.strip() + '_' + start1.strip() + '_' + stop1.strip()
            sr = SeqRecord(Seq(seq1), id = id1, name = name1, description = desc2)
            print(sr)
            SeqIO.write(sr, output_handle, 'fasta')

os.chdir(rp)

