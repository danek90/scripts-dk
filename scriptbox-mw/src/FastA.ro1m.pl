#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
#usage: FastA.ro1m.pl input.fasta > output.fasta

($ARGV[0]) or &HELP_MESSAGE;
($ARGV[0] =~ /fasta|fa|fna/) or &HELP_MESSAGE;

my$seqI = Bio::SeqIO->new(-file=> "$ARGV[0]", -format=> 'fasta');

while(my$seq = $seqI->next_seq){
	print ">". $seq->display_id().".ro1m\n";
	
	my$l=$seq->length();
	my$out = $seq->subseq($l-999998, $l) . $seq->subseq(1, $l-999999);
	
	for( my$o=0; $o <= length($out); $o+=60){
		print substr($out,$o,60) . "\n";	
	}
}

sub HELP_MESSAGE { die "
.Description:	
   Takes a sequence in FASTA format, presumably a circular genome, and re-orients the ends to position 1,000,000 using BioPerl.
   
   WARNING: only works on files containing a single sequence/header.

.Usage: $0 input.fasta > output.fasta
   
   [mandatory]
   <input.fasta>	Input fasta file, expects .fasta, .fa, or .fna extension.
   
   [dependencies]
   1. BioPerl		Bio::SeqIO
 
" }
