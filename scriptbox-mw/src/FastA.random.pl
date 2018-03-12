#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my$sample=1;
&GetOptions( 'n=s'=>\$sample);
($ARGV[0]) or &HELP_MESSAGE;

my$rand = qx( grep ">" $ARGV[0] | sed 's/>//g' | shuf );
chomp$rand;
my%li = ();
my@random = split("\n",$rand);
for(my$n=0; $n<$sample; $n++){
  $li{ $random[$n] } =1;
}

open IN, "$ARGV[0]";
my$good=0;
while(my$ln = <IN>){
   next if $ln =~ /^;/;
   chomp $ln;
   #adapted from FastA.filter.pl (enveomics)
   if($ln =~ m/^>((\S+).*)/){ $good = (exists $li{$1} or exists $li{">$1"} or exists $li{$2} or exists $li{$ln}) }
   print "$ln\n" if ($good);
}
close IN;


sub HELP_MESSAGE { die "
.Description:
   Randomly outputs sequence(s) subset from a Fasta file.

.Usage: $0 [options] input.fasta > output.fasta

   [mandatory]
   input.fasta		Input fasta sequence.
   output.fasta		Filename of new output fasta containing random subset.

   [options]
   -n   <int>     Number of sequences to sample (default = 1).

   [dependencies]


" }
