#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my$input = $ARGV[0];
($ARGV[0]) or &HELP_MESSAGE;

my$seqI = Bio::SeqIO->new(-file=> "$input", -format=> 'fasta');

while(my$seq = $seqI->next_seq){
	print ">". $seq->display_id()."\n";
	my$prot = $seq->translate;
	my$out = $prot->seq;
	for( my$o=0; $o <= length($out); $o+=60){
		print substr($out,$o,60) . "\n";
	}
}

sub HELP_MESSAGE { die "
.Description:
   Translates FASTA of nucleotide sequences into amino acid sequences using BioPerl.

   WARNING: nucleotide sequences are assummed to be genes. Currently uses codon table 11.

.Usage: $0 input.fna > output.faa

   [mandatory]
   <input.fna>	Input fasta file of gene sequences.


   [dependencies]
   1. BioPerl		Bio::SeqIO

" }
