#!/bin/bash -l
#SBATCH --mail-user=sejoslin@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J SFSsing
#SBATCH -e ANGSD_UnfoldedSFS.baq2.%j.err
#SBATCH -o ANGSD_UnfoldedSFS.baq2.%j.out
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


###################
###   get sfs   ###
###################


module load bio 
angsd --version

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

echo "$(date) : Running ANGSD for:"
# recalibrated bams SFS
echo "Recalibrated bams"
cd ${bams_dir}

echo ${SFS_5a}
eval ${SFS_5a}

echo ${SFS_6a}
eval ${SFS_6a}

echo ""

# raw bams SFS
echo "Raw bams"
cd ${raw_dir}

echo ${SFS_5b}
eval ${SFS_5b}

echo ${SFS_6b}
eval ${SFS_6b}


cd ${SFS_dir}

echo "$(date) : Generating unfolded site frequency spectrum for:" 
# recalibrated
echo "Recalibrated bams"
realSFS SFS_minMapQ10_baq2_C50.saf.idx -maxIter 100 > SFS_minMapQ10_baq2_C50.sfs
echo "    Completed SFS_minMapQ10_baq2_C50.sfs"
realSFS SFS_minMapQ30_baq2_C50.saf.idx -maxIter 100 > SFS_minMapQ30_baq2_C50.sfs
echo "    Completed SFS_minMapQ30_baq2_C50.sfs"
# raw
echo "Raw bams"
realSFS SFS_noRecal_minMapQ10_baq2_C50.saf.idx -maxIter 100 > SFS_noRecal_minMapQ10_baq2_C50.sfs
echo "    Completed SFS_noRecal_minMapQ10_baq2_C50.sfs"
realSFS SFS_noRecal_minMapQ30_baq2_C50.saf.idx -maxIter 100 > SFS_noRecal_minMapQ30_baq2_C50.sfs
echo "    Completed SFS_noRecal_minMapQ30_baq2_C50.sfs"
echo ""
echo ""



echo "$(date) : plotting SFS for:"
# recalibrated
echo "Recalibrated bams"
~/scripts/plotSFS.R SFS_minMapQ10_baq2_C50.sfs
echo "    Completed SFS_minMapQ10_baq2_C50.sfs"
~/scripts/plotSFS.R SFS_minMapQ30_baq2_C50.sfs
echo "    Completed SFS_minMapQ30_baq2_C50.sfs"
echo ""

# raw
echo "Raw bams"
~/scripts/plotSFS.R SFS_noRecal_minMapQ10_baq2_C50.sfs
echo "    Completed SFS_noRecal_minMapQ10_baq2_C50.sfs"
~/scripts/plotSFS.R SFS_noRecal_minMapQ30_baq2_C50.sfs
echo "    Completed SFS_noRecal_minMapQ30_baq2_C50.sfs"
echo ""



end=`date +%s`
runtime=$((end-start))

echo Runtime: $runtime
