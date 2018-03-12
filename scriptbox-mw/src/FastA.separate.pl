#!/usr/bin/perl -w
use strict;
use Getopt::Long;

&GetOptions(	'in=s' => \my$fasta,		# input multi-fasta
		'q' => \my$q,
		'dir=s' => \my$dir);		# output directory for individual files
($fasta and $dir) or &HELP_MESSAGE;

unless( -d "$dir") {
	system("mkdir $dir");
}

open IN, "$fasta";
while(my$i = <IN>){
	if($i =~ /^>/){
		chomp$i;
		(my$name = $i) =~ s/>//;
		$name =~ s/[\W]+/_/g;
		$name =~ s/^lcl\|//;
		$name =~ s/"?"//;
		unless($q){ print $name."\n" };
		my$out = $name.".fasta";
		open OUT, ">$dir/$out";
		print OUT ">" . $name."\n";
	}else{ print OUT $i };

}

sub HELP_MESSAGE { die "
.Description:
   Separates a multi-sequence fasta in set of single-sequence files. New files are named with existed fasta deflines stripped of non-word characters.

.Usage: $0 -in in.fasta -dir out

   [mandatory options]
   -in	<in.fasta>	Input multi-fasta.
   -dir <path>		Output directory for individual files.

   [optional]
   -q		quiet

" }
