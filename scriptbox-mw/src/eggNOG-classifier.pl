#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

&GetOptions(	'in=s' => \my$fasta,		# input multi-fasta of proteins
		'dir=s' => \my$dir,		# output directory for intermediate files
		'nocat' => \my$nocat,		# option to NOT produce category count table
		'anvio' => \my$anvio,
		'out=s' => \my$out);		# basename for output annotation table and category counts
($fasta and $dir and $out) or &HELP_MESSAGE;

unless( -d "$dir") {
	system("mkdir $dir");
}
#Step 1: Load eggNOG4 data.

#Load bproNOG annoations
my%bproNOG = ();
open BPRO, "/home/yrh8/Tools/eggNOG/bproNOG.annotations.tsv";
while(my$bp = <BPRO>){
	chomp$bp;
	my@bpro=split(/\s+/,$bp);
	#print $bpro[1]."\t".substr($bp,0,100)."\n";
	$bproNOG{ $bpro[1] } = substr($bp,0,100);
}

#Load eggNOG4 functional categories
my%eggnog4 = ();
open EGG4, "/home/yrh8/Tools/eggNOG/eggnog4.functional_categories.txt";
while(my$e = <EGG4>){
	chomp$e;
	if($e =~ /\[[JAKLB]/){
		(my$cat) = $e =~ /(?<=\[)(\w)/;
		$eggnog4{$cat}{desc}= substr($e,1) . " (INFORMATION STORAGE AND PROCESSING)";
		$eggnog4{$cat}{count}=0;
	}elsif($e =~ /\[[DYVTMNZWUO]/){
		(my$cat) = $e =~ /(?<=\[)(\w)/;
		$eggnog4{$cat}{desc} = substr($e,1) . " (CELLULAR PROCESSES AND SIGNALING)";
		$eggnog4{$cat}{count}=0;
	}elsif($e =~ /\[[CGEFHIPQ]/){
		(my$cat) = $e =~ /(?<=\[)(\w)/;
		$eggnog4{$cat}{desc} = substr($e,1) . " (METABOLISM)";
		$eggnog4{$cat}{count}=0;
	}elsif($e =~ /\[[RS]/){
		(my$cat) = $e =~ /(?<=\[)(\w)/;
		$eggnog4{$cat}{desc} = substr($e,1) . " (POORLY CHARACTERIZED)";
		$eggnog4{$cat}{count}=0;
	}
}



#Step 2: Split multi-fasta into individual fastas.

system("/home/yrh8/Scriptbox/src/FastA.separate.pl -q -in $fasta -dir $dir");

#Step 3: Search each fasta against bproNOG using hmmscan.

my$table = $out;
if($anvio){ $table .= "-func4anvio" };
$table .= ".table.txt";

open OUT, ">$table";

if($anvio){ print OUT join("\t",("gene_callers_id","source","accession","function","e_value"))."\n" };

my$nomatch=0;
my$count=1;
my$numfasta = qx( find $dir/*fasta | wc -l);
chomp$numfasta;
print "\nPerforming HMM searches...\n\n";
foreach my$f (<$dir/*fasta>){
	#print $f . "\n";
	(my$hmm = $f) =~ s/\.fasta/\.hmmscan/;
	(my$name = basename($f)) =~ s/\.fasta//;
	#print $name."\n";
	{
		local $\;
		print "\r\tScanning sequence: $count / $numfasta\t";
	}
	$count++;

	system("/home/yrh8/Tools/hmmer-3.1b2/src/hmmscan -o $hmm.out --tblout $hmm -T 10 /home/yrh8/Tools/eggNOG/bproNOG_hmm-db $f");

	open HMM, "$hmm";
	my$counter=0;
	while(my$h = <HMM>){
		unless($h =~ /^#/){
			my@sh = split(/\s+/,$h);
			(my$enog) = $sh[0] =~ /(?<=bproNOG\.)(\w+)/; #(?=\.meta_raw)/;
			my@senog = split("\t", $bproNOG{ $enog });

			if($anvio){
				$name =~ s/\w+_//;
				print OUT join("\t",($name, "eggNOG4",$enog,$senog[5],$sh[4]))."\n";
				print OUT join("\t",($name, "COG",substr($senog[-2],0,1), $eggnog4{ substr($senog[-2],0,1) }{desc}, $sh[4] ) )."\n";
			}else{
				print OUT join("\t",$name,$enog,@sh[4,5,6])."\t". $bproNOG{ $enog }."\t". $eggnog4{ substr($senog[-2],0,1) }{desc} . "\n";
			}

			$eggnog4{ substr($senog[-2],0,1) }{count}++;
			last;
		}else{
			$counter++;
		}
		if($counter >= 5){
			print OUT $name."\tNo matches found in bproNOG\t".join("\t",(("-")x10))."\n";;
			$nomatch++;
			last;
		}
	}
}
print "\n\n";

#Step 4: Print out eggNOG category counts.
unless($nocat || $anvio){
	my$cout = $out.".counts.txt";
	open OUT2, ">$cout";
	print OUT2 "INFORMATION STORAGE AND PROCESSING\n";
	foreach (split(//,"JAKLB")){
		print OUT2 $eggnog4{ $_ }{count}."\t". $eggnog4{ $_ }{desc}."\n";
	}
	print OUT2 "CELLULAR PROCESSES AND SIGNALING\n";
	foreach (split(//,"DYVTMNZWUO")){
		print OUT2 $eggnog4{ $_ }{count}."\t". $eggnog4{ $_ }{desc}."\n";
	}
	print OUT2 "METABOLISM\n";
	foreach (split(//,"CGEFHIPQ")){
		print OUT2 $eggnog4{ $_ }{count}."\t". $eggnog4{ $_ }{desc}."\n";
	}
	print OUT2 "POORLY CHARACTERIZED\n";
	foreach (split(//,"RS")){
		print OUT2 $eggnog4{ $_ }{count}."\t". $eggnog4{ $_ }{desc}."\n";
	}
	print OUT2 $nomatch."\t[ ] No match in pbroNOG database\n";
}

#####################
sub HELP_MESSAGE { die "
.Description:
   Annotates a set of protein sequences according to bproNOG functional database using hmmer.

.Usage: $0 -in -dir -out

   [mandatory]
   -in	<in.fasta>	Input multi-fasta of proteins.
   -dir	<dir>		Output directory for intermediate files.
   -out	<name>		Basename for output annotation table and category counts.

   [optional]
   -nocat		Option to supress output of category counts.
   -anvio		Output table in format for anvi'o.

" }
