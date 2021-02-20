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
#module load bio
module load R
module load java
module load maven
module load GATK/4.0
module load samtools
module load picardtools/2.7.1

## directories
proj="singletons"
data_dir="/home/sejoslin/projects/${proj}/data"
raw_dir="${data_dir}/bams_raw"
out_dir="${data_dir}/bams_recalibrated"
bamlist="Rand20.bamlist"
ref="DS_history_contigs_250.fasta"

mkdir -p ${raw_dir}/test
mkdir -p ${raw_dir}/badReads

## files
cd ${raw_dir}

# index fasta
idx_ref="samtools faidx ${raw_dir}/${ref}"
echo ${idx_ref}
eval ${idx_ref}

# create dictionary
ref_base=$(basename ${ref} .fasta)
mk_dict="picard-tools CreateSequenceDictionary \
	R=${raw_dir}/${ref} \
	O=${raw_dir}/${ref_base}.dict"
echo ${mk_dict}
eval ${mk_dict}

# change header in angsd.vcf file
#gunzip angsdput.vcf.gz
#sed -i 's/(angsd version)//g' angsdput.vcf

# make bamlist
#	this list will just be file names
for i in *.bam
do
	echo ${i} >> ${bamlist}
done


########################
###  Recalibration!  ###
########################

x=$(wc -l ${bamlist} | awk '{print $1}')

for i in $(seq 1 $x)
do
	cd ${raw_dir}
	# pull bam file name
        bam=$(awk -v LINE="$i" '{if (NR==LINE) print $0}' ${raw_dir}/${bamlist})
        prefix=$(echo ${bam} | cut -d. -f1)
	echo ${prefix}
	echo "Recalibrating bam file $i out of $x"

###   Picard Work   ###
#######################
# test BAM files
	call_val1="picard-tools ValidateSamFile \
		I=${bam} \
		MODE=SUMMARY"
        echo ${call_val1}
	eval ${call_val1}
# add read groups
        call_RG="picard-tools AddOrReplaceReadGroups \
		I=${prefix}.sort.proper.rmdup.bam \
		O=${prefix}.sort.proper.rmdup.RG.bam \
		RGID=1 \
		RGLB=lib1 \
		RGPL=illumina \
		RGPU=run \
		RGSM=${prefix}"
	echo ${call_RG}
	eval ${call_RG}
# test BAM files
	call_val2="picard-tools ValidateSamFile \
		I=${prefix}.sort.proper.rmdup.RG.bam \
		O=test/${prefix}.test \
		MODE=VERBOSE  \
		IGNORE_WARNINGS=true"
	echo ${call_val2}
	eval ${call_val2}
# grab bad mates
	echo "Making list of bad reads in ${prefix}.sort.proper.rmdup.RG.bam"
	y=$(wc -l test/${prefix}.test | awk '{print $1}')
	for i in $(seq 1 $y)
	do 
		awk -v LINE="$i" '{if (NR==LINE) print $0}' test/${prefix}.test | cut -d "," -f1 | cut -d " " -f4 >> badReads/${prefix}.badReads.list
	done
# filter out bad reads
	call_filt="picard-tools FilterSamReads \
		I=${prefix}.sort.proper.rmdup.RG.bam \
		O=${prefix}.sort.proper.rmdup.FILT.RG.bam \
		READ_LIST_FILE=badReads/${prefix}.badReads.list \
		FILTER=excludeReadList"
	echo ${call_filt}
	eval ${call_filt}
# validate BAMs
	call_val3="picard-tools ValidateSamFile \
		I=${prefix}.sort.proper.rmdup.FILT.RG.bam \
		MODE=VERBOSE  \
		IGNORE_WARNINGS=true"
	echo ${call_val3}


###   GATK   ###
################
# generate recalibration table based on various covariates
	cd ${raw_dir}
	call_BR="gatk BaseRecalibrator \
		-I ${raw_dir}/${prefix}.sort.proper.rmdup.FILT.RG.bam \
		-O ${out_dir}/${prefix}.table
		-R ${raw_dir}/${ref}
		--known-sites ${raw_dir}/angsdput.vcf"
	echo ${call_BR}
	eval ${call_BR}

# adjust quality score
	cd ${out_dir}
	call_AB="gatk ApplyBQSR \
		-R ${raw_dir}/${ref} \
		-I ${raw_dir}/${prefix}.sort.proper.rmdup.FILT.RG.bam \
		--bqsr-recal-file ${out_dir}/${prefix}.table \
		-O ${out_dir}/${prefix}.recal.bam"
	echo ${call_AB}
	eval ${call_AB}

done

###   README   ###
##################
cat << EOT >> README
    *.table == recalibration table with default recalibration settings
    *.recal.bam == bam files with recalibrated base calls (previously sorted, removed dups, and properly paired)
EOT


end=`date +%s`
runtime=$((end-start))

echo Runtime: $runtime
