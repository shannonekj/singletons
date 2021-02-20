#!/bin/bash -l
#SBATCH --mail-user=sejoslin@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J recalGATK
#SBATCH -e GATK_recalibrate.%j.err
#SBATCH -o GATK_recalibrate.%j.out
#SBATCH -c 20
#SBATCH --time=03:03:00
#SBATCH --mem=60G
#SBATCH -p bigmemh

set -e # exits upon failing command
set -v # verbose -- all lines
hostname
start=`date +%s`
echo "My SLURM_JOB_ID: $SLURM_JOB_ID"
THREADS=${SLURM_NTASKS}
echo "Threads: $THREADS"
MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024)
echo "Mem: $MEM"

#### NOTES: ####
# This script is to try to figure out why we are getting so many singletons in SFS.
# GATK Base Quality Score Recalibration


#### SETUP ####
module load bio
module load R
module load java
module load maven
module load GATK/4.0
module load samtools
module load picardtools/2.7.1

# directories
proj="singletons"
data_dir="/home/sejoslin/projects/${proj}/data"
raw_dir="${data_dir}/bams_raw"
out_dir="${data_dir}/bams_recalibrated"
bamlist="Rand20.bamlist"
ref="DS_history_contigs_250.fasta"

mkdir -p ${out_dir}

## files
cd ${raw_dir}
# index fasta
samtools faidx ${raw_dir}/${ref}

# create dictionary
ref_base=$(basename ${ref} .fasta)
picard-tools CreateSequenceDictionary \
        R=${raw_dir}/${ref} \
        O=${raw_dir}/${ref_base}.dict

cd ${out_dir}
cat << EOT >> README
This directory contains recalibrated bam files from ${raw_dir}.

Files:
EOT

# make bamlist
cd ${raw_dir}
ls -alh *.bam | awk '{print $9}' >> ${bamlist}



#######################
###   known sites   ###
#######################
# need to generate a list of polymorphic loci that are 'good'
# 	I'm going to try with sites that have MAF >= 0.05
#	-doVcf
angsd --version
call_snps="angsd -b ${bamlist} \
	-dovcf 1 \
	-gl 1 \
	-dopost 1 \
	-domajorminor 1 \
	-domaf 1 \
	-snp_pval 1e-6 \
	-minMaf 0.05"
echo ${call_snps}
eval ${call_snps}

echo "    angsdput.vcf.gz == output from creating vcf file from raw bams" >> README



############################
###   BaseRecalibrator   ###
############################
# generate recalibration table based on various covariates

cd ${raw_dir}
x=$(wc -l ${raw_dir}/${bamlist} | awk '{print $1}')

for i in $(seq 1 $x)
do
	# pull bam file name
	bam=$(awk -v LINE="$i" '{if (NR==LINE) print $0}' ${raw_dir}/${bamlist})
	prefix=$(echo ${bam} | cut -d. -f1)
	call_BR="gatk BaseRecalibrator \
		-I ${raw_dir}/${bam} \
		-O ${out_dir}/recal_${prefix}.table
		-R ${raw_dir}/${ref}
		--known-sites angsdput.vcf.gz"
	echo ${call_BR}
	eval ${call_BR}
done

echo "    recal_*.table == recalibration table with default recalibration settings" >> README



#####################
###   ApplyBQSR   ###
#####################
# adjust quality score

cd ${out_dir}
x=$(wc -l ${raw_dir}/${bamlist} | awk '{print $1}')

for i in $(seq 1 ${x})
do
	bam=$(awk -v LINE="$i" '{if (NR==LINE) print $0}' ${raw_dir}/${bamlist})
	prefix=$(echo ${bam} | cut -d/ -f8 | cut -d. -f1)
	call_AB="gatk ApplyBQSR \
		-R ${raw_dir}/${ref} \
		-I ${raw_dir}/${bam} \
		--bqsr-recal-file recal_${prefix}.table \
		-O ${prefix}.recal.bam"
	echo ${call_AB}
	eval ${call_AB}
done

echo "    *.recal.bam == recalibrated base calls" >> README



end=`date +%s`
runtime=$((end-start))

echo Runtime: $runtime
