#!/bin/python

from Bio import SeqIO, Entrez
import requests
import os
import pandas as pd

# Functions
## import arguements
def dataset_input(argv):
    #get inputs
    input_csv = ''
    prokka_folder = ''
    fasta_folder = ''
    Entrez_email = ''

    try:
        opts, args = getopt.getopt(argv, "hi:p:f:m:", ["help", "input_csv=",  "prokka_folder=", "fasta_folder=", "entrez_email="])
    except getopt.GetoptError:
        print('Error! Usage: python fastadown.py -i <pullrm issues csv> -p <for_prokka folder> -f <fasta_folder> -m <email for ncbi search>' )
        print('   or: python fastadown.py --input_csv <issues csv from pullrm> --prokka_folder <prokka folder of gbk files> --fasta_folder <fasta output folder> --email <email for ncbi search>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('Error! Usage: python fastadown.py -i <pullrm issues csv> -p <for_prokka folder> -f <fasta_folder> -m <email for ncbi search>' )
            print('   or: python fastadown.py --issues_csv <issues csv from pullrm> --prokka_folder <prokka folder of gbk files> --fasta_folder <fasta output folder> --email <email for ncbi search>')
            sys.exit()
        elif opt in ("-i", "--issues_csv"):
            input_csv = arg
        elif opt in ("-p", "--prokka_folder"):
            prokka_folder = arg
        elif opt in ("-f", "--fasta_folder"):
            fasta_folder = arg
        elif opt in ("-m", "--email"):
            Entrez_email = arg
    return input_csv,int(window),fasta_folder,output_gbk_path,Entrez_email

# take inputs
issues_csv,prokka_folder,fasta_folder,Entrez_email = dataset_input(sys.argv[1:])

#Entrez_email = "hbouzek@fredhutch.org"
rp = os.getcwd()
# Create a new folder for FASTA files
frp = rp + "/" + fasta_folder
if not os.path.exists(frp):
    os.mkdir(frp)

uid_list = os.listdir("./" + prokka_folder)
Entrez.email = Entrez_email

def remove_end(s):
    return s[0:-4]

# Find all sequences that were troublesome
uid_list = [remove_end(s) for s in uid_list]
uid_list2 = pd.read_csv(issues_csv)['contig_name']
uid_list = uid_list + uid_list2.tolist()
uid_list = list(set(uid_list))

def fasta_download(query,output_path):
    filename = output_path + '/'+ query + '.fasta'
    # Downloading...
    print("Downloading :", query)
    net_handle = Entrez.efetch(db="nuccore", id=query, rettype="fasta", retmode="text")
    out_handle = open(filename, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    print ("Saved")
    return

def process_file(filename, file_type):
    for seq_record in SeqIO.parse(filename, file_type):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))

if len(uid_list) > 0:
    Entrez.email = Entrez_email
    n=1
    for each_query in uid_list:
        fasta_download(each_query, frp)
        print ("count ",n)
        print (' ')
        n+=1
