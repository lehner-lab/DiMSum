#!/usr/bin/env python

#######################

import os.path
import sys
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import math
import statistics
import numpy

#######################

#Globals
#Illumina 1.8+ Phred+33 read score dictionary
phred_string = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
phred_dict = {i:phred_string.index(i) for i in list(phred_string)}
qual_dict = {j:i for (i,j) in phred_dict.items()}
eprob_dict = {i:math.pow(10, phred_string.index(i)/(-10.0)) for i in list(phred_string)}

#######################

#Functions
def readnlines(file_handle, n):
    lines = []
    for i in range(n):
        try:
            lines.append(next(file_handle))
        except StopIteration:
            return None
    return(lines)

def next_FASTQ_record(file_handle):
    temp_record = readnlines(file_handle, 4)
    if temp_record == None:
        return(None)
    else:
        record_list = [i.rstrip() for i in temp_record]
        fq_record = {}
        fq_record['read_name'] = record_list[0]
        fq_record['read_seq'] = record_list[1]
        fq_record['read_name_qual'] = record_list[2]
        fq_record['read_qual'] = record_list[3]
        return(fq_record)

def write_FASTQ_record(file_handle, fq_record):
    #Add UMI and corrected length to read name
    file_handle.write(fq_record['read_name']+'\n'+
                      fq_record['read_seq']+'\n'+
                      fq_record['read_name_qual']+'\n'+
                      fq_record['read_qual']+'\n')
    
#######################

#parse command-line arguments
parser = argparse.ArgumentParser(description='Attempt to merge paired end reads according to a specified alignment length')
parser.add_argument('-i1', '--inputFile1', help='Path to first read FASTQ file ', required=True)
parser.add_argument('-i2', '--inputFile2', help='Path to second read FASTQ file', required=True)
parser.add_argument('-o', '--outputFile', help='Path to output FASTQ file', required=True)
parser.add_argument('-r', '--reportFile', help='Path to output report file', required=True)
parser.add_argument('-n', '--numNuc', type=int, help='Enforced alignment length in base pairs', required=True)
parser.add_argument('-m', '--minqual', type=int,  help='Minimum observed base quality to retain read pair', required=True)
parser.add_argument('-e', '--maxee', type=float, help='Maximum number of expected errors to retain read pair', required=True)

args = parser.parse_args()

input_FASTQ1 = args.inputFile1
input_FASTQ2 = args.inputFile2
output_FASTQ = args.outputFile
output_REPORT = args.reportFile
num_nuc = args.numNuc
min_qual = args.minqual
max_ee = args.maxee

#######################

#Alignment statistics
a_stats = {}
a_stats['Pairs'] = 0
a_stats['Merged'] = 0
a_stats['Alignments with zero diffs'] = 0
a_stats['Too many diffs (> ' + str(num_nuc) + ')'] = 0
a_stats['Exp.errs. too high (max=' + str(max_ee) + ')'] = 0
a_stats['Min Q too low (' + str(min_qual) + ')'] = 0
a_stats['Fwd too short (<64)'] = 0
a_stats['Rev too short (<64)'] = 0
a_stats['merged_lengths'] = []
#Open FASTQ file
input_file1 = open(input_FASTQ1, 'r')
input_file2 = open(input_FASTQ2, 'r')
#Check whether output files already exist
if os.path.isfile(output_FASTQ) or os.path.isfile(output_REPORT):
    print("file output path exists")
    sys.exit()
else:
    output_file = open(output_FASTQ, 'w')
    output_report = open(output_REPORT, 'w')
    #Loop over FASTQ records in files
    while True:
        fqr1 = next_FASTQ_record(input_file1)
        fqr2 = next_FASTQ_record(input_file2)
        fqo = {}
        if fqr1 == None:
            break

        #Update statistics
        a_stats['Pairs'] += 1
        #Both sequences at least 64 bp long
        if len(fqr1['read_seq']) >= 64 and len(fqr2['read_seq']) >= 64:
            #Reverse complement read2
            temp_seq = Seq(fqr2['read_seq'], generic_dna)
            fqr2['read_seq'] = str(temp_seq.reverse_complement())
            #Do reads align exactly?
            if fqr1['read_seq'][-num_nuc:] == fqr2['read_seq'][:num_nuc]:
                #Update statistics
                a_stats['Alignments with zero diffs'] += 1
                #Reverse qualities
                fqr2['read_qual'] = fqr2['read_qual'][::-1]
                #Read error probabilities
                fqr1['read_eprob'] = [eprob_dict[i] for i in list(fqr1['read_qual'])]
                fqr2['read_eprob'] = [eprob_dict[i] for i in list(fqr2['read_qual'])]
                #Merged region error probabilities
                px = fqr1['read_eprob'][-num_nuc:]
                py = fqr2['read_eprob'][:num_nuc]
                merged_region_eprob = [((i[0]*i[1] / 3.0) / (1.0 - i[0] - i[1] + 4.0*i[0]*i[1] / 3.0)) for i in list(zip(px, py))]
                #Merged region qualities
                merged_region_qualscore = [(-10) * math.log10(i) for i in merged_region_eprob]
                #Are merged region qualities less than specified?
                merge_num_bases_too_low_qual = sum([i<min_qual for i in merged_region_qualscore])
                if merge_num_bases_too_low_qual == 0:
                    #Are non-merged region qualities less than specified?
                    non_merge_num_bases_too_low_qual = sum([phred_dict[i]<min_qual for i in fqr1['read_qual'][:-num_nuc] + fqr2['read_qual'][num_nuc:]])
                    if non_merge_num_bases_too_low_qual == 0:
                        fqo['read_eprob'] = fqr1['read_eprob'][:-num_nuc] + merged_region_eprob + fqr2['read_eprob'][num_nuc:]
                        #Are expected number of read errors less than specified?
                        exp_num_read_errors = sum(fqo['read_eprob'])
                        if exp_num_read_errors < max_ee:
                            fqo['read_seq'] = fqr1['read_seq'][:-num_nuc] + fqr2['read_seq']
                            a_stats['merged_lengths'].append(len(fqo['read_seq']))
                            merged_region_qual = [qual_dict[int(i)] if int(i) <= 41 else 'J' for i in merged_region_qualscore]
                            fqo['read_qual'] = fqr1['read_qual'][:-num_nuc] + ''.join(merged_region_qual) + fqr2['read_qual'][num_nuc:]
                            fqo['read_name'] = fqr1['read_name']
                            fqo['read_name_qual'] = '+'
                            write_FASTQ_record(output_file, fqo)
                            a_stats['Merged'] += 1
                        else:
                            #Update statistics
                            a_stats['Exp.errs. too high (max=' + str(max_ee) + ')'] += 1
                    else:
                        #Update statistics
                        a_stats['Min Q too low (' + str(min_qual) + ')'] += 1
                else:
                    #Update statistics
                    a_stats['Min Q too low (' + str(min_qual) + ')'] += 1
            else:
                #Update statistics
                a_stats['Too many diffs (> ' + str(num_nuc) + ')'] += 1
        else:
            #Update statistics
            if len(fqr1['read_seq']) < 64:
                a_stats['Fwd too short (<64)'] += 1
            if len(fqr2['read_seq']) < 64:
                a_stats['Rev too short (<64)'] += 1
    #Report
    output_report.write('\nMerge\n\tFwd ' + input_FASTQ1 + '\n\tRev ' + input_FASTQ2 + '\n\tKeep read labels\n\t' + str(a_stats['Merged']) + ' / ' + str(a_stats['Pairs']) + ' pairs merged\n\n')
    output_report.write('Merged length distribution:\n')
    output_report.write('\t ' + str(min(a_stats['merged_lengths'])) + '  Min\n')
    output_report.write('\t ' + str(int(numpy.percentile(a_stats['merged_lengths'], 25))) + '  Low quartile\n')
    output_report.write('\t ' + str(int(numpy.percentile(a_stats['merged_lengths'], 50))) + '  Median\n')
    output_report.write('\t ' + str(int(numpy.percentile(a_stats['merged_lengths'], 75))) + '  High quartile\n')
    output_report.write('\t ' + str(max(a_stats['merged_lengths'])) + '  Max\n\nTotals:\n')
    output_report.write('\t ' + str(a_stats['Pairs']) + '  Pairs ()\n')
    output_report.write('\t ' + str(a_stats['Merged']) + '  Merged (, )\n')
    output_report.write('\t ' + str(a_stats['Alignments with zero diffs']) + '  Alignments with zero diffs ()\n')
    output_report.write('\t ' + str(a_stats['Too many diffs (> ' + str(num_nuc) + ')']) + '  Too many diffs (> 1) ()\n')
    output_report.write('\t ' + str(0) + '  Fwd tails Q <= 2 trimmed ()\n')
    output_report.write('\t ' + str(0) + '  Rev tails Q <= 2 trimmed ()\n')
    output_report.write('\t ' + str(a_stats['Fwd too short (<64)']) + '  Fwd too short (< 64) after tail trimming ()\n')
    output_report.write('\t ' + str(a_stats['Rev too short (<64)']) + '  Rev too short (< 64) after tail trimming ()\n')
    output_report.write('\t ' + str(0) + '  No alignment found ()\n')
    output_report.write('\t ' + str(0) + '  Alignment too short ( ) ()\n')
    output_report.write('\t ' + str(a_stats['Exp.errs. too high (max=' + str(max_ee) + ')']) + '  Exp.errs. too high (max=' + str(max_ee) + ') ()\n')
    output_report.write('\t ' + str(a_stats['Min Q too low (' + str(min_qual) + ')']) + '  Min Q too low (<' + str(min_qual) + ') ()\n')
    #Close files
    input_file1.close()
    input_file2.close()
    output_file.close()
    output_report.close()





