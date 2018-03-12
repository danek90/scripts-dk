#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Getopt::Long;

&GetOptions( 'h' => \my$h,
		'in=s' => \my$in,
		'help' => \my$help);

if($h || $help){&HELP_MESSAGE};
unless($in){&HELP_MESSAGE};
unless($in =~ /gbk$/ || $in =~ /gb$/){&HELP_MESSAGE};

my$gbk = Bio::SeqIO->new(-file => "$in", -format => 'genbank');
my$out = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);

while(my$seq = $gbk->next_seq){
	for my$feat ($seq->get_SeqFeatures){
		$out->write_feature($feat);
	}
}

sub HELP_MESSAGE { die "
.Description:
   Simply converts Genbank file to GFF3 format using BioPerl.

.Usage: $0 -in input.gbk > output.gff

   [mandatory]
   <input.gbk>	Input file in Genbank format (must be *.gbk or *.gb).
   <output.gff>		Output file in GFF3 format.

   [dependencies]
   1. BioPerl		Bio::SeqIO, Bio::Tools::GFF

" }
