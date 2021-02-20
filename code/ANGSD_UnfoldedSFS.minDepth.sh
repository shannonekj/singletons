#!/bin/bash -l
#SBATCH --mail-user=sejoslin@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J SFSsing
#SBATCH -e ANGSD_UnfoldedSFS.%j.err
#SBATCH -o ANGSD_UnfoldedSFS.%j.out
#SBATCH -c 20
#SBATCH --time=10:03:00
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
# This script is to try to figure out why we are getting so many singletons in out SFS.
#	Altered flags:	-baq
#			-C
#			-minMapQ

#### SETUP ####
# directories
proj="singletons"
data_dir="/home/sejoslin/projects/${proj}/data"
raw_dir="${data_dir}/bams_raw"
bams_dir="${data_dir}/bams_recalibrated"
SFS_dir="${bams_dir}/results_SFS_unfold_wParalogs"

ref="${raw_dir}/DS_history_contigs_250.fasta"
raw_bams="${raw_dir}/Rand20.bamlist"
recal_bams="${bams_dir}/Rand20.bamlist"

mkdir -p ${SFS_dir}
cd ${SFS_dir}

cat << EOT >> README
This directory's purpose is to house SFSs to try to figure out why we are seeing so many singletons in our SFS generated from RAD sequencing data. I will be generating and comparing unfolded SFS

$(date)
Bamlist for recalibrated SFS can be found at ${recal_bams}
Bamlist for raw SFS can be found at ${raw_bams}
Reference partial genome assembly is located at ${ref}
All bams were processed with GATK Base Recalibration (BQSR) prior to generating SFS.

EOT

# make bamlist
cd ${bams_dir}
ls -alh *.bam | awk '{print $9}' >> ${recal_bams}

###################
###   get sfs   ###
###################

# ALREADY COMPLETED #echo "Indexing reference."
# ALREADY COMPLETED #samtools faidx ../${ref}
# ALREADY COMPLETED #sleep 1m
# ALREADY COMPLETED #touch ../${ref}.fai

module load bio 
angsd --version

SFS_1a="angsd -bam ${recal_bams} \
	-out ${SFS_dir}/SFS_minMapQ10 \
	-anc ${ref} \
	-GL 2 -doSaf 1 \
	-minMapQ 10 -minQ 20"

SFS_1b="angsd -bam ${raw_bams} \
        -out ${SFS_dir}/SFS_noRecal_minMapQ10 \
        -anc ${ref} \
        -GL 2 -doSaf 1 \
        -minMapQ 10 -minQ 20"

SFS_2a="angsd -bam ${recal_bams} \
	-out ${SFS_dir}/SFS_minMapQ30 \
	-anc ${ref} \
	-GL 2 -doSaf 1 \
	-minMapQ 30 -minQ 20"

SFS_2b="angsd -bam ${raw_bams} \
        -out ${SFS_dir}/SFS_noRecal_minMapQ30 \
        -anc ${ref} \
        -GL 2 -doSaf 1 \
        -minMapQ 30 -minQ 20"

SFS_3a="angsd -bam ${recal_bams} \
        -out ${SFS_dir}/SFS_minMapQ10_baq1_C50 \
        -anc ${ref} \
        -GL 2 -doSaf 1 \
        -minMapQ 10 -minQ 20 \
	-ref ${ref} \
        -baq 1 -C 50"

SFS_3b="angsd -bam ${raw_bams} \
        -out ${SFS_dir}/SFS_noRecal_minMapQ10_baq1_C50 \
        -anc ${ref} \
        -GL 2 -doSaf 1 \
        -minMapQ 10 -minQ 20 \
        -ref ${ref} \
        -baq 1 -C 50"

SFS_4a="angsd -bam ${recal_bams} \
        -out ${SFS_dir}/SFS_minMapQ30_baq1_C50 \
        -anc ${ref} \
        -GL 2 -doSaf 1 \
        -minMapQ 30 -minQ 20 \
	-ref ${ref} \
	-baq 1 -C 50"

SFS_4b="angsd -bam ${raw_bams} \
        -out ${SFS_dir}/SFS_noRecal_minMapQ30_baq1_C50 \
        -anc ${ref} \
        -GL 2 -doSaf 1 \
        -minMapQ 30 -minQ 20 \
        -ref ${ref} \
        -baq 1 -C 50"

SFS_5a="angsd -bam ${recal_bams} \
        -out ${SFS_dir}/SFS_minMapQ10_baq2_C50 \
        -anc ${ref} \
        -GL 2 -doSaf 1 \
        -minMapQ 10 -minQ 20 \
        -ref ${ref} \
        -baq 2 -C 50"

SFS_5b="angsd -bam ${raw_bams} \
        -out ${SFS_dir}/SFS_noRecal_minMapQ10_baq2_C50 \
        -anc ${ref} \
        -GL 2 -doSaf 1 \
        -minMapQ 10 -minQ 20 \
        -ref ${ref} \
        -baq 2 -C 50"

SFS_6a="angsd -bam ${recal_bams} \
        -out ${SFS_dir}/SFS_minMapQ30_baq2_C50 \
        -anc ${ref} \
        -GL 2 -doSaf 1 \
        -minMapQ 30 -minQ 20 \
        -ref ${ref} \
        -baq 2 -C 50"

SFS_6b="angsd -bam ${raw_bams} \
        -out ${SFS_dir}/SFS_noRecal_minMapQ30_baq2_C50 \
        -anc ${ref} \
        -GL 2 -doSaf 1 \
        -minMapQ 30 -minQ 20 \
        -ref ${ref} \
        -baq 2 -C 50"



cd ${SFS_dir}
cat << EOT2 >> README
$(date) : Creating site allele frequency likelihood based on genotype likelihoods assuming HWE.
All SFS were run with the following base commands:
    angsd -bam ${recal_bams} -anc ${ref} -GL 2 -doSaf 1 -minQ 20

And were run with additional commands stated in each file name
EOT2

echo "$(date) : Running ANGSD for:"
# recalibrated bams SFS
echo "Recalibrated bams"
cd ${bams_dir}

echo ${SFS_1a}
eval ${SFS_1a}

echo ${SFS_2a}
eval ${SFS_2a}

echo ${SFS_3a}
eval ${SFS_3a}

echo ${SFS_4a}
eval ${SFS_4a}
echo ""

# raw bams SFS
echo "Raw bams"
cd ${raw_dir}

echo ${SFS_1b}
eval ${SFS_1b}

echo ${SFS_2b}
eval ${SFS_2b}

echo ${SFS_3b}
eval ${SFS_3b}

echo ${SFS_4b}
eval ${SFS_4b}

cd ${SFS_dir}
cat << EOT3 >> README
    *.saf.gz == bzgf compressed flat floats (see angsd/doc/formats.pdf documentation) resulting from -doSaf
    *.saf.pos.gz == compressed positions resulting from -doSaf
    *.idx == uncompressed binary index file resulting from -doSaf
EOT3
echo ""


echo "$(date) : Generating unfolded site frequency spectrum for:" 
cd ${SFS_dir}

#recalibrated SFSs
echo "Recalibrated bams"
realSFS SFS_minMapQ10.saf.idx -maxIter 100 > SFS_minMapQ10.sfs
echo "    Completed SFS_minMapQ10.sfs"
realSFS SFS_minMapQ30.saf.idx -maxIter 100 > SFS_minMapQ30.sfs
echo "    Completed SFS_minMapQ30.sfs"
realSFS SFS_minMapQ10_baq1_C50.saf.idx -maxiter 100 > SFS_minMapQ10_baq1_C50.sfs
echo "    Completed SFS_minMapQ10_baq1_C50.sfs"
realSFS SFS_minMapQ30_baq1_C50.saf.idx -maxiter 100 > SFS_minMapQ30_baq1_C50.sfs
echo "    Completed SFS_minMapQ30_baq1_C50.sfs"
echo ""

# raw SFSs
echo "Raw bams"
realSFS SFS_noRecal_minMapQ10.saf.idx -maxIter 100 > SFS_noRecal_minMapQ10.sfs
echo "    Completed SFS_noRecal_minMapQ10.sfs"
realSFS SFS_noRecal_minMapQ30.saf.idx -maxIter 100 > SFS_noRecal_minMapQ30.sfs
echo "    Completed SFS_noRecal_minMapQ30.sfs"
realSFS SFS_noRecal_minMapQ10_baq1_C50.saf.idx -maxiter 100 > SFS_noRecal_minMapQ10_baq1_C50.sfs
echo "    Completed SFS_noRecal_minMapQ10_baq1_C50.sfs"
realSFS SFS_noRecal_minMapQ30_baq1_C50.saf.idx -maxiter 100 > SFS_noRecal_minMapQ30_baq1_C50.sfs
echo "    Completed SFS_noRecal_minMapQ30_baq1_C50.sfs"
echo ""
echo ""

echo "    *.sfs == a maximum likelihood estimate of the SFS" >> README


echo "$(date) : plotting SFS for:"
# recalibrated
echo "Recalibrated bams"
~/scripts/plotSFS.R SFS_minMapQ10.sfs
echo "    Completed SFS_minMapQ10.sfs"
~/scripts/plotSFS.R SFS_minMapQ30.sfs
echo "    Completed SFS_minMapQ30.sfs"
~/scripts/plotSFS.R SFS_minMapQ10_baq1_C50.sfs
echo "    Completed SFS_minMapQ10_baq1_C50.sfs"
~/scripts/plotSFS.R SFS_minMapQ30_baq1_C50.sfs
echo "    Completed SFS_minMapQ30_baq1_C50.sfs"
echo ""

# raw
echo "Raw bams"
~/scripts/plotSFS.R SFS_noRecal_minMapQ10.sfs
echo "    Completed SFS_noRecal_minMapQ10.sfs"
~/scripts/plotSFS.R SFS_noRecal_minMapQ30.sfs
echo "    Completed SFS_noRecal_minMapQ30.sfs"
~/scripts/plotSFS.R SFS_noRecal_minMapQ10_baq1_C50.sfs
echo "    Completed SFS_noRecal_minMapQ10_baq1_C50.sfs"
~/scripts/plotSFS.R SFS_noRecal_minMapQ30_baq1_C50.sfs
echo "    Completed SFS_noRecal_minMapQ30_baq1_C50.sfs"


echo "    *pdf == SFS plot resulting from ~/scripts/plotSFS.R <file>.sfs" >> README


end=`date +%s`
runtime=$((end-start))

echo Runtime: $runtime
