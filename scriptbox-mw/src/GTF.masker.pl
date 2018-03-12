#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;

&GetOptions(
        'in=s'=>\my$input,      #
        'fasta=s'=>\my$fasta,       #
        'query=s'=>\my$query,   #fasta file of sequence matches to be masked
        'q'=>\my$q,
        'r'=>\my$rev,
        'f=s'=>\my$flank2,
        'out=s'=>\my$output);      #new output FASTA

($input and $fasta and $query and $output) or &HELP_MESSAGE;

my%blast=();
my@hits = &doblast($query,$fasta);
foreach my$h (@hits){
  if($flank2){
    $blast{ @{$h}[0] - $flank2 } = @{$h}[1] + $flank2; #add specified flanking sequence
  }else{
    $blast{ @{$h}[0] } = @{$h}[1];
  }
}
unless($q){ print "\n\tFound " .scalar(keys%blast)." BLASTn hits in '$fasta'.\n\n" };

open IN, "$input";
open OUT, ">$output";
while(my$i = <IN>){
  if($i =~ /^\#/){
    print OUT $i;
  }else{
    chomp$i;
    my@si = split("\t",$i);
    my$mask = &hitcheck( $si[3], $si[4], \%blast);
    if($mask == 1){
      if($rev){ print OUT $i."\n" };
      #print join("-",($si[3], $si[4]))."\n";
    }else{
      unless($rev){ print OUT $i."\n" };
    }
  }
}


########################
sub doblast {
  my$blast = qx( blastn -query $_[0] -subject $_[1] -outfmt \'10 sstart send\' -qcov_hsp_perc 75 -perc_identity 80 );
  chomp$blast;
  my@hits = split("\n",$blast);
  return( &hitsorter( \@hits ) );
}

sub hitsorter {
  my@hits = @{ $_[0] };
  for(my$h=0; $h < scalar@hits;$h++){
    #print $h."\t".$hits[$h]."\n";
    my@s = split(",",$hits[$h]);
    my@ss = (sort {$a <=> $b}(@s));
    my$sss = join(",",@ss);
    #print "\t".join(",",@ss)."\n";
    splice(@hits, $h,1,\@ss); #replaces each hit with array ref of sorted coordinates
  }
  return( @hits );
}

sub midpoint {
  my$mid = ($_[1] - $_[0])/2 + $_[0];
  return( $mid );
}

sub hitcheck {
  my$start = $_[0];
  my$end = $_[1];
  my%blast = %{$_[2]};
  my$fmid = &midpoint($start, $end);
  my$match=0;

  foreach my$c1 (keys%blast){
    my$bmid = &midpoint( $c1, $blast{$c1});
    if( ($c1 <= $start) && ($fmid <= $blast{$c1}) ){
      $match=1;
      last;
    }elsif( ($end <= $blast{$c1}) && ($fmid >= $c1) ){
      $match=1;
      last;
    }elsif( ($start <= $c1) && ($bmid <= $end) ){
      $match=1;
      last;
    }elsif( ($blast{$c1} <= $end) && ($bmid >= $start) ){
      $match=1;
      last;
    }
  }
  return( $match );
}

sub ERROR {
  print "\nERROR: Looks like your input file has >1 sequence. Sorry, this script can only handle single-contig genomes.\n\n";
  &HELP_MESSAGE;
}

sub HELP_MESSAGE { die "
.Description:
   Filters annotated features in a *.gtf file whose coordinates overlap with BLASTn alignment of a query fasta.

.Usage: $0 -in <input.gtf> -fasta <input.fasta> -query <q.fasta> -out <output.gtf>

   [mandatory]
   -in      Input file in gtf format.
   -fasta   Cooresponding genome sequence fasta.
   -query   Fasta for BLASTn search.
   -out     Filename of new gtf file.

   [options]
   -f       Mask additional bp flanking coordinates of BLASTn matches. (eg. -f 6)
   -r       Reverse, outputs features that overlap BLASTn matches.
   -q       Runs quietly.

   [dependencies]
   Bioperl
   BLASTn (-qcov_hsp_perc 75 -perc_identity 80)

" }
