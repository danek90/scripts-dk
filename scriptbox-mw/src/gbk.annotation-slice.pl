#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
my$coord2=0;
&GetOptions(	'genome=s' => \my$genome,	# input genome
		'c1=s' => \my$coord1,		# first coordinate
		'c2=s' => \$coord2,		# second coordinate
		'w=s' => \my$window,		# number of basepairs around central coordinate
		'q' => \my$q,
		'f' => \my$f,
		'out=s' => \my$out);		# output table

($genome and $coord1 and $out) or &HELP_MESSAGE;
($coord2 or $window) or &HELP_MESSAGE;
if($coord2 and $window){ print "\n# ERROR: You cannot define -c2 AND -w. Please only pick one.\n" and &HELP_MESSAGE };

my@coords=();
if($coord2){ @coords = sort{$a <=> $b}($coord1,$coord2) };
if($window){ @coords = (($coord1 - $window), ($coord1 + $window)) };

open OUT, ">$out";
unless($f){ print OUT "Locus_tag\tStrand\tCoords\tPGAP annotation\tTranslation\n"}

my$count=0;
my$origins=0;
my$gbk = Bio::SeqIO->new(-file => $genome, -format => 'genbank');
while(my$seq = $gbk->next_seq){
	if($origins == 1){ &ERROR };
	$origins++;

	my($contig)=$seq->id();
	my$l=$seq->length();
	if(($coord1 > $l) || ($coord2 > $l)){ print "\n#ERROR: please enter coordinates between 1 and $l\n" and &HELP_MESSAGE };

	if($coords[0] < 1){
		unshift(@coords, ( ($l + $coords[0]), $l));
		$coords[2] = 1;
		unless($q){ print "\n#WARNING: Defined window extends beyond sequence length. Hopefully your sequence is circular.\n"}
	}elsif($coords[1] > $l){
		push(@coords, (1, ($coords[1] - $l)) );
		$coords[1] = $l;
		unless($q){ print "\n#WARNING: Defined window extends beyond sequence length. Hopefully your sequence is circular.\n"}
	}

	unless($q){
		print "\n\tSearching within range(s): ".join("-",($coords[0],$coords[1]));
		if(scalar(@coords) == 4){ print " and ".join("-",($coords[2],$coords[3]))."\n";
		}else{ print "\n"};
	}

	for my$feat ($seq->get_SeqFeatures){
		if(($feat->primary_tag eq "CDS") && ($feat->has_tag("protein_id"))){
			my($s1) = $feat->start;
			my($s2) = $feat->end;
			my($sd) = $feat->strand;
			my($id) = $feat->get_tag_values("locus_tag");
			my($prot) = $feat->get_tag_values("translation");
			my($prod) =  $feat->get_tag_values("product");

			for(my$o=0; $o < scalar(@coords); $o+=2){
				if(($s2 >= $coords[$o]) && ($s1 <= $coords[$o+1])){
					if($f){ print OUT ">".$id."\n".$prot."\n";
					}else{ print OUT $id."\t".$sd."\t".$s1."-".$s2."\t".$prod."\t".$prot."\n"}
					$count++;
				}
			}
		}
	}
}

unless($q){
	print "\tFound ".$count." protein-coding genes in: ".$genome."\n\n";
}

sub ERROR {
  print "\nERROR: Looks like your input file has >1 sequence. Sorry, this script can only handle single-contig genomes.\n\n";
  &HELP_MESSAGE;
}

sub HELP_MESSAGE { die "
.Description:
   Returns a table of annotated genes within a given range of genomic coordinates, including those which extend beyond range.
	 Range must be defined as either [1] start - stop, OR [2] start +/- window.

.Usage: $0 -genome in.gbk -c1 start -c2 stop -w window -out out.txt

   [mandatory]
   -genome	<in.gbk>	Input full genbank file. Assumes single, circular contig.
   -c1		<number>	First coordinate.
   -c2		<number>	Second coordinate.
   -w		<number>	Window size in bp +/- around coordinate 1 (eg. 10000).
   -out		<out.txt>	New output table.

   [optional]
   -f				Output as FASTA instead of table.
   -q				Run quietly.

   [future]
   -cir				Input sequence is circular, extend window beyond end.
   -o		<fasta>		Assign ortholog id's based on fasta sequence using rmb.rb (nt or aa). NOT AVAILABLE
   -nt				Output nt gene sequences instead of protein (default).

   [dependencies]
   BioPerl

" }
