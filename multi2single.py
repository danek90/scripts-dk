#WHOLE GENOME SEQUENCE MULTIPLE LOCUS VARIABLE NUMBER TANDEM REPEAT SEARCH 
#last edit: 03/13/2018

#Imported Modules.
import sys
import re
import os
import argparse
from Bio import SeqIO

#Program requires 2 arguments to be entered and will autormatically
#load in the VNTR chart data.

parser = argparse.ArgumentParser(prog = " MLVA in silico Program", \
								usage = 'Find MLVA code from closed \
										fasta sequence.')
parser.add_argument('-in', '--inputFile', required = True , \
					help = 'input fasta file. Required.')
parser.add_argument('-n', '--name', required = False, \
					help = 'Assign file name. If not assigned it will use \
							the file name of the fasta file.')
args = parser.parse_args()
inFile = args.inputFile
sampleID = args.name
outFile = sampleID +"out.fasta"


#Reads sequence file list and stores it as a string object. Close file.
#Splits string at the start of a line.
#First fasta in the file is split into an empty element and the first fasta.
#"del" removes this empty element.
with open(inFile,"r") as newFile:
	sequences = newFile.read()
	sequences = re.split("^>", sequences, flags = re.MULTILINE)
	del sequences[0]
	newFile.close()
#Conversts multiline fasta to single line. Writes new fasta to file: "X_out".
#1. Split each fasta into header and sequence.
#2. Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
#3. Replace newlines in sequence, remembering to add one to the end.
with open(outFile,"w") as newFasta:
	for fasta in sequences:
		try:
			header, sequence = fasta.split("\n", 1)
		except ValueError:
			print fasta
		header = ">" + header + "\n"
		sequence = sequence.replace("\n","") + "\n"
		newFasta.write(header + sequence)
	newFasta.close()