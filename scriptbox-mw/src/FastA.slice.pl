#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
#usage: FastA.slice.pl input.fasta start stop > output.fasta

my$input = $ARGV[0];

($ARGV[0] and $ARGV[1] and $ARGV[2]) or &HELP_MESSAGE;
($ARGV[1] < $ARGV[2]) or &HELP_MESSAGE;


my$seqI = Bio::SeqIO->new(-file=> "$input", -format=> 'fasta');

while(my$seq = $seqI->next_seq){
	print ">". $seq->display_id()."_".$ARGV[1]."-".$ARGV[2]."\n";
	my$out = $seq->subseq($ARGV[1],$ARGV[2]);
	for( my$o=0; $o <= length($out); $o+=60){
		print substr($out,$o,60) . "\n";	
	}
}

sub HELP_MESSAGE { die "
.Description:	
   Takes a slice of sequence from a FASTA file between two given coordinates (inclusive) using BioPerl.
   
   WARNING: only works on files containing a single sequence/header.

.Usage: $0 input.fasta start stop > output.fasta
   
   [mandatory]
   <input.fasta>	Input fasta file.
   <start>		Start coordinate, must be >1 and <stop.
   <stop>		End coordinate, must be >start and <length.
   
   [dependencies]
   1. BioPerl		Bio::SeqIO
 
" }
