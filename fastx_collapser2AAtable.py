#!/usr/bin/env python

#######################

import os.path
import sys
import argparse
from Bio import SeqIO

#######################

#Globals

#######################

    
#######################

#parse command-line arguments
parser = argparse.ArgumentParser(description='Convert USEARCH unique FASTA file to table of counts with AA translation')
parser.add_argument('-i', '--inputFile', help='Path to USEARCH unique FASTA file (with counts in sequence name)', required=True)
parser.add_argument('-o', '--outputFile', help='Path to output file', required=True)

args = parser.parse_args()

input_FASTA = args.inputFile
output_tab = args.outputFile

#######################

#Check whether output files already exist
if os.path.isfile(output_tab):
    print("file output path exists")
    sys.exit()
else:
	output_file = open(output_tab, 'w')
	output_file.write('\t'.join(['nt_seq', 'count', 'aa_seq'])+'\n')
	for seq_record in SeqIO.parse(input_FASTA, "fasta"):
		#Initialise output dictionary
		temp = {}
		#Nucleotide sequence
		temp['nt_seq'] = str(seq_record.seq)
		temp['count'] = seq_record.name.split('-')[1]
		#Convert nucleotide sequence to lower case
		temp['nt_seq'] = temp['nt_seq'].lower()
		#Translate nucleotide sequence
		temp['aa_seq'] = str(seq_record.translate().seq)
		#Write to output file
		output_file.write('\t'.join([temp['nt_seq'], temp['count'], temp['aa_seq']])+'\n')

#Exit
sys.exit()
