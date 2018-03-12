#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;

&GetOptions(
        'in=s'=>\my$input,      #input FASTA
        'bed=s'=>\my$bed,       #tabular bed file of coordinates to be masked
        'query=s'=>\my$query,   #fasta file of sequence matches to be masked
        'header=s'=>\my$header, #header for output FASTA
        'q'=>\my$q,
        'b=s'=>\my$flank,
        'f=s'=>\my$flank2,
        'out=s'=>\my$out);      #new output FASTA

($input and $header and $out) or &HELP_MESSAGE;
($bed or $query) or &HELP_MESSAGE;

my$seqI = Bio::SeqIO->new(-file=> "$input", -format=> 'fasta');
my$original = '';
my$count=0;
while(my$seq = $seqI->next_seq){
  if($count == 1){ &ERROR };
  $original .= $seq->seq();
  $count++;
}
my$seqIN = Bio::Seq->new(-seq => $original, -display_id=>"original", -format=>'fasta');
my%hash=();
my%blast=();

if($bed){
  open BED, "$bed";
  while(my$i = <BED>){
  	chomp$i;
  	my@si=split("\t",$i);
  	my@coords = sort { $a <=> $b}($si[-1],$si[-2]);
  	if($flank){
      $hash{ $coords[0] - $flank } = $coords[1] + $flank; #add specified flanking sequence
    }else{
      $hash{ $coords[0] } = $coords[1];
    }
  }
  close BED;
}

if($query){
  my@hits = &doblast($query,$input);
  foreach my$h (@hits){
    if($flank2){
      $blast{ @{$h}[0] - $flank2 } = @{$h}[1] + $flank2; #add specified flanking sequence
    }else{
      $blast{ @{$h}[0] } = @{$h}[1];
    }
  }
}

my$bedmasked = &nmasker( $seqIN, \%hash);
my$blastmasked = &nmasker( $bedmasked, \%blast);

=head
while(my$seq = $seqI->next_seq){
	print OUT ">". $header ."\n"; # $seq->display_id . "-masked\n";
	$end = $seq->length;

	for (my$j=0; $j <= @sorted-1; $j++){
		if($j==0){
			$outseq .= $seq->subseq(1,$sorted[$j]-1);
			$outseq .= "N"x( $hash{$sorted[$j]}-$sorted[$j]+1);
			unless($q){print "\t"."1--". ($sorted[$j]-1) ."\n"}; #[".$sorted[$j]."-".$hash{$sorted[$j]}."]\n";
			unless($q){print "\t". (($hash{$sorted[$j]}) - $sorted[$j] + 1) . " Ns\t\t[".$sorted[$j]."-".$hash{$sorted[$j]}."]\n"};
		}else{
			if( $hash{$sorted[$j-1]} >= $sorted[$j]){
				$outseq .= "N"x($hash{$sorted[$j]}-$hash{$sorted[$j-1]});
				unless($q){print "\t". ($hash{$sorted[$j]}-$hash{$sorted[$j-1]})." Ns\t\t[".($hash{$sorted[$j-1]}+1)."-".$hash{$sorted[$j]}."]\n"};
			}else{
				$outseq .= $seq->subseq( $hash{ $sorted[$j-1]}+1, $sorted[$j]-1);
				unless($q){print "\t". ($hash{ $sorted[$j-1]}+1) ."--". ($sorted[$j]-1) ."\n"};
				$outseq .= "N"x( $hash{$sorted[$j]}-$sorted[$j]+1);
		        	unless($q){print "\t". ($hash{$sorted[$j]}-$sorted[$j]+1)." Ns\t\t[".$sorted[$j]."-".$hash{$sorted[$j]}."]\n"};
			}
		}
	}
	$outseq .= $seq->subseq( $hash{$sorted[-1]}+1,$end);
	unless($q){print "\t".($hash{$sorted[-1]}+1)."--".$end."\n"};
}
=cut
unless($q){print "Original length = ".$seqIN->length."\n"};
unless($q){print "Masked length = ".$blastmasked->length ."\n"};

open OUT, ">$out";
print OUT ">".$header."\n";
my$outseq = $blastmasked->seq;
for( my$o=0; $o <= length($outseq); $o+=60){
        print OUT substr($outseq,$o,60) . "\n";
}

########################
sub doblast {
  my$blast = qx( blastn -query $_[0] -subject $_[1] -outfmt \'10 sstart send\' -qcov_hsp_perc 85 -perc_identity 90 );
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

sub nmasker {
  my$bioseq = $_[0];
  my%hash = %{$_[1]};
  if(scalar(keys%hash) == 0){
    return ($bioseq);
  }else{
    my@sorted = (sort {$a <=> $b}(keys%hash));
    my$outseq = '';
    my$end = $bioseq->length;
    for (my$j=0; $j <= @sorted-1; $j++){
  		if($j==0){
  			$outseq .= $bioseq->subseq(1,$sorted[$j]-1);
  			$outseq .= "N"x( $hash{$sorted[$j]}-$sorted[$j]+1);
  			#unless($q){print "\t"."1--". ($sorted[$j]-1) ."\n"}; #[".$sorted[$j]."-".$hash{$sorted[$j]}."]\n";
  			#unless($q){print "\t". (($hash{$sorted[$j]}) - $sorted[$j] + 1) . " Ns\t\t[".$sorted[$j]."-".$hash{$sorted[$j]}."]\n"};
      }else{
  			if( $hash{$sorted[$j-1]} >= $sorted[$j]){
  				$outseq .= "N"x($hash{$sorted[$j]}-$hash{$sorted[$j-1]});
  				#unless($q){print "\t". ($hash{$sorted[$j]}-$hash{$sorted[$j-1]})." Ns\t\t[".($hash{$sorted[$j-1]}+1)."-".$hash{$sorted[$j]}."]\n"};
  			}else{
  				$outseq .= $bioseq->subseq( $hash{ $sorted[$j-1]}+1, $sorted[$j]-1);
  				#unless($q){print "\t". ($hash{ $sorted[$j-1]}+1) ."--". ($sorted[$j]-1) ."\n"};
  				$outseq .= "N"x( $hash{$sorted[$j]}-$sorted[$j]+1);
  		    #unless($q){print "\t". ($hash{$sorted[$j]}-$sorted[$j]+1)." Ns\t\t[".$sorted[$j]."-".$hash{$sorted[$j]}."]\n"};
  			}
  		}
  	}
    $outseq .= $bioseq->subseq( $hash{$sorted[-1]}+1,$end);
    #unless($q){print "\t".($hash{$sorted[-1]}+1)."--".$end."\n"};
    my$seqOUT=Bio::Seq->new(-seq=>$outseq, -display_id=>'masked',-format=>'fasta');
    return ( $seqOUT );
  }
}

sub ERROR {
  print "\nERROR: Looks like your input file has >1 sequence. Sorry, this script can only handle single-contig genomes.\n\n";
  &HELP_MESSAGE;
}

sub HELP_MESSAGE { die "
.Description:
   Masks specified genomic coordinates with 'Ns' in new fasta. Coordinates can be supplied as a bed file and/or identified by BLASTn alignment with a query fasta.
   At least one .bed or .fasta file is required.

.Usage: $0 -in <input.fasta> -bed <coords.bed> -query <q.fasta> -header <header> -out <output.fasta>

   [mandatory]
   -in		Input fasta sequence.
   -header	Header for output fasta.
   -out		Filename of new output fasta.

   [options]
   -bed     List of coordinates to mask in bed format.
   -query   Fasta file of sequence(s) to mask according to BLASTn alignment.
            (-qcov_hsp_perc 85 -perc_identity 90)
   -b       Mask additional bp flanking coordinates in bed file. (eg. -b 6)
   -f       Mask additional bp flanking coordinates of BLASTn matches. (eg. -f 6)
   -q       Runs quietly.

   [dependencies]
   Bioperl
   BLASTn

" }
