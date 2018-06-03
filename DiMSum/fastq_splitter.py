#!/usr/bin/env python

#######################

import os.path
import sys
import argparse

#######################

#Globals

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
    file_handle.write(fq_record['read_name']+'\n'+
                      fq_record['read_seq']+'\n'+
                      fq_record['read_name_qual']+'\n'+
                      fq_record['read_qual']+'\n')
    
#######################

#parse command-line arguments
parser = argparse.ArgumentParser(description='Split FASTQ file into roughly equally sized chunks (bytes)')
parser.add_argument('-i', '--inputFile', help='Path to FASTQ file', required=True)
parser.add_argument('-o', '--outputFilePrefix', help='Prefix to output FASTQ file', required=True)
parser.add_argument('-c', '--chunkSize', type=int, help='Chunk size in bytes', required=False)
parser.add_argument('-n', '--numRecords', type=int, help='Number of FASTQ records per file', required=False)

args = parser.parse_args()

input_FASTQ = args.inputFile
output_prefix = args.outputFilePrefix
chunk_size = args.chunkSize

#######################

#Determine number of fastq records to write per output file
num_fastqrecords = ""
if args.numRecords != None:
	#Number of records provided
	num_fastqrecords = args.numRecords
elif args.chunkSize != None:
	#Number of bytes per file provided
	#Open FASTQ file
	input_file = open(input_FASTQ, 'r')
	#Number of bytes in first FASTQ record
	bytesperrecord = len(''.join(readnlines(input_file, 4)))
	#Estimated lines per file
	num_fastqrecords = int(chunk_size/bytesperrecord)
	#Close file
	input_file.close()
else:
	print("Either --chunkSize or --numRecords arguments need to be specified")
	sys.exit()

#######################

#Split all FASTQ records into numbered output files with max num_fastqrecords per file
count = 0
records = 0
output_file = ""
#Open FASTQ file
input_file = open(input_FASTQ, 'r')
while True:
	this_record = next_FASTQ_record(input_file)
	if this_record == None:
		break
	if count == 0:
		count += 1
		output_FASTQ = output_prefix + str(count) + '.fastq'
		#Check whether output files already exist
		if os.path.isfile(output_FASTQ):
		    print("file output path exists")
		    sys.exit()
		else:
			output_file = open(output_FASTQ, 'w')
	#This file is full
	elif records == num_fastqrecords:
		#Close file
		output_file.close()
		#Increment number of files
		count += 1
		#Reset record cout
		records = 0
		#Open new file
		output_FASTQ = output_prefix + str(count) + '.fastq'
		#Check whether output files already exist
		if os.path.isfile(output_FASTQ):
		    print("file output path exists")
		    sys.exit()
		else:
			output_file = open(output_FASTQ, 'w')
	#Write FASTQ record
	write_FASTQ_record(output_file, this_record)
	records += 1

#Close last file (if not closed already)
output_file.close()

#Return number of fastq records per file
print(num_fastqrecords)

#Exit
sys.exit()


