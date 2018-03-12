#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Perl;
use Getopt::Long;

&GetOptions(
		'in=s' => \my$in,
		's' => \my$s,
		'p' => \my$p,
		'c' => \my$c,
		'pseudo' => \my$ps,
		'h' => \my$h,
		'out=s' => \my$outfile,
		);

if($h){ &HELP_MESSAGE };
if($p and $ps){
	print "WARNING: option '-p' masks option '-pseudo'.\n";
}
($in) or &HELP_MESSAGE;

my$gbk = Bio::SeqIO->new(-file => $in, -format => 'genbank');

while(my$seq = $gbk->next_seq){
	my($contig)=$seq->id();
	my$chrom = $seq->length;
	for my$feat ($seq->get_SeqFeatures){
		if(($feat->primary_tag eq "CDS") && ($feat->has_tag("protein_id"))){
			my@cds = &cds_values( $feat );

			#print fasta header
			print ">";
			if($c){print $contig."--"};
			print $cds[0];
			unless($s){
				print " [protein=".$cds[4]."] [location=";
				if($cds[3] == 1){ print $cds[1]."..".$cds[2]."]";
			}else{ print "complement(".$cds[1]."..".$cds[2].")]"};
			}
			print "\n";

			#print fastq sequence
			if($p){ &print_seq( $cds[-1] );
			}else{
				if($cds[-2] eq 'null'){
					$cds[-2] = &boundary_gene( $seq, $cds[1],$cds[2],$chrom,$cds[3]);
				}
				&print_seq( $cds[-2] );
			}

		}elsif(($feat->primary_tag eq "gene") && ($feat->has_tag("pseudo")) && $ps){
			my@ps_feat = &pseudo_values( $feat );
			unless($p){
				print ">";
				if($c){print $contig."--"};
				print $ps_feat[0];
				unless($s){
					print " [protein=pseudogene] [location=";
					if($ps_feat[3] == 1){ print $ps_feat[1]."..".$ps_feat[2]."]";
				}else{ print "complement(".$ps_feat[1]."..".$ps_feat[2].")]"};
				}
				print "\n";
				&print_seq( $ps_feat[-1]);
			}
		}
	}
}

###############################################
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
	if($s2 > $s1){
		($gene) = $feat->seq->seq;
 	}
	my($prot) = $feat->get_tag_values("translation");
	return( ($id,$s1,$s2,$sd,$prod,$gene,$prot) );
}

sub boundary_gene {
	my($seq,$start,$end,$total,$strand) = @_;
	my$gene = '';
	if($strand > 0){
		$gene = $seq->subseq($start,$total);
		$gene .= $seq->subseq(1,$end);
	}else{
		my$tmpgene = $seq->subseq($start,$total);
		$tmpgene .= $seq->subseq(1,$end);
		$gene = revcom($tmpgene)->seq;
	}
	return( $gene );
}

sub pseudo_values {
	my$feat = $_[0];
	my($id) = $feat->get_tag_values("locus_tag");
	my($s1) = $feat->start;
	my($s2) = $feat->end;
	my($sd) = $feat->strand;
	my($gene) = $feat->seq->seq;
	return( ($id,$s1,$s2,$sd,$gene));
}


sub print_seq {
	for( my$o=0; $o < length($_[0]); $o+=60){
  	print substr($_[0],$o,60) . "\n";
	}
}

sub HELP_MESSAGE { die "
.Description:
   Extracts gene nucleotide (default) or amino acid sequences from a genbank file.

.Usage: $0 -in in.gbk > out.fasta

   [mandatory]
   -in	<in.gbk>	Input full genbank file.

   [optional]
   -s			Output simple headers, without annotation or coordinate info.
   -c			Include contig accession in definition line.
   -p			Output translated amino acid sequences.
   -pseudo		Include pseudogenes (/pseudo).
   -h			This help message.

   [dependencies]
   BioPerl

" }
