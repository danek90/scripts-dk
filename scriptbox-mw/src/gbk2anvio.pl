#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;

&GetOptions( 'h' => \my$h,
		'in=s' => \my$in,
		'fasta=s' => \my$fasta,
		'genes=s' => \my$genes,
		'func=s' => \my$func,
		'help' => \my$help);

if($h || $help){&HELP_MESSAGE};
($in and $fasta and $genes and $func) or &HELP_MESSAGE;
unless($in =~ /gbk$/ || $in =~ /gb$/){&HELP_MESSAGE};


my$pgap = qx(grep "Annotation Software revision" $in);
my($pgapv) = $pgap =~ m/([\d.]+)/;
#print $pgapv."\n";
my$contig = qx( grep ">" "$fasta" | sed 's/>//');
chomp$contig;


my$gbk = Bio::SeqIO->new(-file => "$in", -format => 'genbank');
open GENES, ">$genes";
open FUNC, ">$func";

print GENES "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tsource\tversion\n";
print FUNC "gene_callers_id\tsource\taccession\tfunction\te_value\n";

while(my$seq = $gbk->next_seq){
	#my($contig)=$seq->id();
	my$chrom = $seq->length;
	for my$feat ($seq->get_SeqFeatures){
		my@cds = ();

		if(($feat->primary_tag eq "CDS") && ($feat->has_tag("protein_id"))){
			@cds = &cds_values( $feat );

		# }elsif(($feat->primary_tag eq "gene") && ($feat->has_tag("pseudo")) && $ps ){
		# 	@cds = &pseudo_values( $feat );
		# 	$cds[0] .= "_ps";
		}else{
			next;
		}
		my$locus=$cds[0];
		$cds[0] =~ s/\w+_//;

		if($cds[3] > 0){ $cds[3] = "f";
		}else{ $cds[3] = "r"};

		#checks for 'partial' genes not evenly divisible by 3
		my$partial = 0;
		if( ($cds[2] - $cds[1] + 1) % 3){ $partial = 1 };

		#because anvio is 0-indexed
		$cds[1] = $cds[1] - 1;
		#$cds[2] = $cds[2] - 1;

		print GENES join("\t",($cds[0],$contig,$cds[1],$cds[2],$cds[3],$partial,"PGAP","v".$pgapv))."\n";
		print FUNC join("\t",($cds[0],"PGAP",$locus,$cds[-1],"0"))."\n";
	}
}




################
sub cds_values {
	my$feat = $_[0];
	my($id) = $feat->get_tag_values("locus_tag");
	my($s1) = $feat->start;
	my($s2) = $feat->end;
	my($sd) = $feat->strand;
	my($prod) = 'unknown';
	if($feat->has_tag("product")){
		($prod) = $feat->get_tag_values("product");
	}
	my$gene = 'null';
	# if($s2 > $s1){
	# 	($gene) = $feat->seq->seq;
 # 	}
	# my($prot) = $feat->get_tag_values("translation");
	return( ($id,$s1,$s2,$sd,$prod) ); #,$gene,$prot) );
}


# sub pseudo_values {
# 	my$feat = $_[0];
# 	my($id) = $feat->get_tag_values("locus_tag");
# 	my($s1) = $feat->start;
# 	my($s2) = $feat->end;
# 	my($sd) = $feat->strand;
# 	#my($gene) = $feat->seq->seq;
# 	return( ($id,$s1,$s2,$sd)); #,$gene));
# }

sub HELP_MESSAGE { die "
.Description:
   Extracts gene call info from Genbank file into table formated for Anvi'o using BioPerl.

	 **Please note - Anvi'o using 0-indexing for gene/genome sequences**

.Usage: $0 [options] -in input.gbk -fasta genome.fasta -genes genes.tsv -func func.tsv

   [mandatory]
   <input.gbk>		Input file in Genbank format (must be *.gbk or *.gb).
   <genome.fasta>	Corresponding genome sequence in FASTA format.
   <genes.tsv>		Output tab-separated table of gene calls.
   <func.tsv>			Output tab-separated table of functional annotations.

   [optional]

   [dependencies]
   1. BioPerl		Bio::SeqIO

" }
