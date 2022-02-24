#!/bin/python

from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import pandas as pd
import os, re, sys, getopt
from os.path import exists

## the csv must have columns named genome_id, genome_name, contig_name, contig_start, contig_end
### genome_id - the NCBI identifier for the organism's sequence
### genome_name - Genus species name of the organism
### contig_name - Genbank identifier used to locate the genbank file
### contig_start - start of the contig provided to determine neighborhood of gene
### contig_end - end of the contig provided to determine neighborhood of the gene

# Function that allows the python code to accept arguments
def dataset_input(argv):
    #get inputs
    input_csv = ''
    fasta_folder = ''
    Entrez_email = ''
    output_gbk_path = ''
    window = ''
    try:
        opts, args = getopt.getopt(argv, "hi:w:g:o:m:", ["help", "input_csv=","window=", "fasta_folder=","output_gbk_path=", "entrez_email="])
    except getopt.GetoptError:
        print('Error! Usage: python pullrm.py -i <input OCTAPUS csv> -w <contig_search_window> -g <gene_folder> -o <output gbk path> -m <email for ncbi search>' )
        print('   or: python pullrm.py --input_csv <OTU csv> --window <contig_search_window> --gene_folder <gene_folder> --output <output gbk path> --email <email for ncbi search>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('Usage: python pullrm.py -i <input OCTAPUS csv> -w <contig_search_window> -g <fasta output path > -o <output gbk path> -m <email for ncbi search>' )
            print('   or: python pullrm.py --input_csv <OTU csv> --window <contig_search_window> --gene_folder <gene_folder> --output <output gbk path> --email <email for ncbi search>')
            sys.exit()
        elif opt in ("-i", "--input_csv"):
            input_csv = arg
        elif opt in ("-w", "--window"):
            window = arg
        elif opt in ("-f", "--fasta_folder"):
            gene_folder = arg
        elif opt in ("-o", "--output"):
            output_gbk_path = arg
        elif opt in ("-m", "--email"):
            Entrez_email = arg
    return input_csv,int(window),fasta_folder,output_gbk_path,Entrez_email

# take inputs
input_csv,window,fasta_folder,output_gbk_path,Entrez_email = dataset_input(sys.argv[1:])
# For testing:
#print(input_csv,window,gene_folder,output_gbk_path,Entrez_email)
#input_csv = "odd_gbk.csv"
#window = 1000
#output_gbk_path = "genebank_files"
#fasta_folder = "test01"
#Entrez_email = "hbouzek@fredhutch.org"



# if on mac and get this error: "urllib.error.URLError: <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local issuer certificate (_ssl.c:1108)":
# 1. Open Finder and navigated to Applications, Python 3.XX
# 2. Double-click on Install Certificates.command file

# Create file folders
rp = os.getcwd()
# Create a new folder for Genbank file
grp = rp + "/" + output_gbk_path
if not os.path.exists(grp):
    os.mkdir(grp)

# Create a new folder for Fasta file
orp = rp + "/" + fasta_folder
if not os.path.exists(orp):
    os.mkdir(orp)

prp = rp + "/" + "for_prokka"
if not os.path.exists(prp):
    os.mkdir(prp)

crp = rp + "/" + "contig_shift"
if not os.path.exists(crp):
    os.mkdir(crp)

# Function to download gbk files
def gbk_download(query,output_path):
    filename = output_path + '/'+ query + '.gbk'
    print("Query :",query)
    # Downloading...
    print("Downloading :",query)
    net_handle = Entrez.efetch(db="nucleotide",id=query,rettype="gb", retmode="text")
    out_handle = open(filename, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    print ("Saved")
    return

# Grab file names and locations from the csv file
xldoc = pd.read_csv(input_csv)

# Contig name has the genbank file names
df = xldoc[['genome_id', 'genome_name', 'contig_name', 'contig_start', 'contig_end']]
fn = df.loc[:,('contig_name')].map(str)+ '.gbk'
df = df.assign(contig_filename = fn)
dfs = df.loc[:,('contig_start', 'contig_end')].min(axis =1)
df = df.assign(search_start = dfs)
dfe = df.loc[:,('contig_start', 'contig_end')].max(axis =1)
df = df.assign(search_end = dfe)

# Check to see if files are in the folder
fileList = df['contig_filename']
contents = os.listdir(grp)
indices = [i for i, element in enumerate(fileList) if element in contents]
indices2 = [i for i, element in enumerate(fileList) if element not in contents]
df2 = df.loc[indices]
df2 = df2.reset_index(drop = True)
print("Found %d files in folder" % (len(df2)))
df3 = df.loc[indices2]

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

def extract_cds (dfrow, gb, gbk_path, shifted_contig_path):
    op = gbk_path + "/" + dfrow['contig_filename']
    np = shifted_contig_path + "/" + dfrow['contig_filename']
    new_cds = pd.DataFrame()
    # setting col width any lower will cause sequence to be truncated
    pd.set_option("max_colwidth", 9999)
    cds1 = [feature for feature in gb.features if feature.type == "CDS"]
    if len(cds1) == 0:
        os.rename(op, np)
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
            product = feature.qualifiers['product']
            new_val = pd.DataFrame({'description':[dfrow['genome_name']],'gb':[dfrow['genome_id']],'emb':[gb.id], 'start':[start], \
            'stop':[stop], 'strand':[strand], 'product':[product], 'cds':[cds], 'contig_filename':[dfrow['contig_filename']]})
            new_cds = new_cds.append(new_val, ignore_index = True)
    return new_cds

# extract CDS
def feature_extract (dfrow, gbk_path, working_path, prokka_folder_path, shifted_contig_folder_path):
    os.chdir(gbk_path)
    np = prokka_folder_path + "/" + dfrow['contig_filename']
    fn = dfrow['contig_filename'].strip()
    op = "./" + fn
    gb = SeqIO.read(fn, 'genbank')
    # Parse features
    lgb = len(gb.features)
    if lgb < 1:
        print("No CDS found for %s." % gb.id)
        new_cds = None
        # Send files with no CDS to "for_prokka" folder
        os.rename(op, np)
    else:
        count = 1
        cds = extract_cds(dfrow, gb, grp, crp)
        count += 1
        print('%d CDS features collected for %s' % (count, gb.id))
    os.chdir(working_path)
    return cds

def yousendme (features_extracted, rowinfo, window, genbank_path, cwd):
    fn = rowinfo['contig_filename']
    ingbk = genbank_path + "/" + fn
    np = cwd + "/" + "for_prokka" + "/" + fn
    if len(features_extracted) == 0 and exists(ingbk):
        os.rename(ingbk, np)
    elif len(features_extracted) == 0 and not exists(ingbk):
        return
    # find the region within 1000 bases of start and/or 1000 bases of stop
    else:
        s_end = rowinfo['search_end'] + window
        s_start = rowinfo['search_start'] - window
        find_gene1 = features_extracted[features_extracted['stop'] < s_end]
        find_gene = find_gene1[find_gene1['start'] > s_start]
        # put files into a misfit folder
        if len(find_gene) == 0 or len(find_gene['cds']) == 0:
            os.rename(ingbk, np)
        else:
            return find_gene

def clean_seq(found_gene_row):
       print(found_gene_row)
       seq1 = found_gene_row['cds']
       if len(seq1) < 1:
           return None
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
       print(sr)
       return sr

# use feature_extract parse through df
ofn = orp + "/" + input_csv[:-4] + ".fasta"
with open(ofn, "w") as output_handle:
    for index, row in df2.iterrows():
        print(row)
        ftx = feature_extract(row, grp, rp, prp, crp)
        op = grp + "/" + row['contig_filename']
        find_gene = yousendme(ftx, row, window, grp, rp)
        if find_gene is None:
            continue
        else:
            find_gene.reset_index()
            for index, newrow in find_gene.iterrows():
                sr = clean_seq(newrow)
                SeqIO.write(sr, output_handle, 'fasta')

# Find file contents of different bins and save as a data frame
os.chdir(rp)
cnts_gbkfolder = os.listdir(grp)
ndf1 = pd.DataFrame({'files':cnts_gbkfolder, 'code':"CDS"})
cnts_prokka = os.listdir(prp)
ndf2 = pd.DataFrame({'files':cnts_prokka, 'code':"Prokka"})
ndf = ndf1.append(ndf2, ignore_index = True)
cnts_geneshift = os.listdir(crp)
ndf3 = pd.DataFrame({'files':cnts_geneshift, 'code':"Gene Shift"})
ndf1 = ndf.append(ndf3, ignore_index = True)
# Print out to csv
ndf1.to_csv(rp +"/" + "results.csv")
