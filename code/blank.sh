#!/bin/bash -l
#SBATCH --mail-user=sejoslin@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J BLANK
#SBATCH -e BLANK.%j.err
#SBATCH -o BLANK.%j.out
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


#### SETUP ####
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


############################
###   BaseRecalibrator   ###
############################
# generate recalibration table based on various covariates




end=`date +%s`
runtime=$((end-start))

echo Runtime: $runtime
