import os
import sys
import argparse
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(prog = "dna2AA", \
								usage = 'translates dna fragment to amino acid.')
parser.add_argument('-i', '--inputFile', required = True , \
					help = 'input nucleotide sequence file.')


#records = list(SeqIO.parse("example.fasta", "fasta"))
args = parser.parse_args()
inFile = args.inputFile
outFile = inFile + ".out"
#outFile = args.outputFile

with open(inFile, "r") as f:
	with open(outFile, "w") as output:
		for record in SeqIO.parse(f, "fasta"):
			data = record.seq
			name = record.id
			AAdata = name, "\n", data.translate(table=11)
			print >>output, AAdata
			print "translation complete."



