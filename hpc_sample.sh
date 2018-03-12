#!/bin/bash -l
#
# Embedded Grid Engine parameters
#
#$ -N snippy-20170105-THISRUN
#$ -l mem_free=20G
#$ -q dbd.q
#
# Refer all file reference to work the current working directory which is
# the directory from which the script was qsubbed
#$ -cwd

#Send email
#$ -M yrh8@cdc.gov
#$ -m be

# Load the module
module load Python/2.7.3
module load htslib/1.3.1
module load perl/5.12.3
module load bwa/0.7.12
module load samtools/1.3.1
module load freebayes/1.0.2
module load GNUScience/1.9
module load vcflib/2015.10.19
module load vcftools/0.1.11
module load snpeff/4.0e
module load snippy/3.1
module load parallel/20131122
export TMPDIR=$HOME/tmp

echo "snippy vs NC_002929" 1>&2
echo "THISRUN" 1>&2
echo "READS1  READS2" 1>&2
date 1>&2

snippy --outdir /scicomp/home/yrh8/Hypermutation/Tohama-snippy/OUTDIR \
--mincov 10 --minfrac 0.75 \
--reference /scicomp/home/yrh8/Hypermutation/TohamaI_NC_002929.fasta \
--R1 READS1 \
--R2 READS2;


#--peil READS1




date 1>&2



