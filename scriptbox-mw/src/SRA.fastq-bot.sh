#!/bin/bash

if [[ "$1" == "" || "$1" == "-h" ]] ; then
   echo "
   This is a simple bot for downloading Illumina fastq for a list of SRA run accessions.

   " >&2 ;
   exit 1 ;
fi ;

runlist="$1";

function downloadit {
	echo $1;
	~/Tools/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump \
	--split-spot --split-files --gzip \
	--defline-seq '@$ac.$si.$rl.$ri' --defline-qual '+$ac.$si.$rl.$ri' $1;
	echo

}
export -f downloadit

mapfile -t acc < $runlist
parallel -j 4 -k downloadit {} ::: ${acc[@]};
