#!/bin/bash -l
#SBATCH --mail-user=sejoslin@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J SFSsing
#SBATCH -e ANGSD_UnfoldedSFS.untrimmed.j%j.err
#SBATCH -o ANGSD_UnfoldedSFS.untrimmed.j%j.out
#SBATCH -c 20
#SBATCH --time=01-02:03:04
#SBATCH --mem=MaxMemPerNode
#SBATCH -p bigmemh

set -e # exits upon failing command
set -v # verbose -- all lines

hostname
start=`date +%s`
echo "My SLURM_JOB_ID: $SLURM_JOB_ID"

#### NOTES: ####
# This script is to try to figure out why we are getting so many singletons in out SFS.
#	Altered flags:	-baq [0,1,2]
#			-C [0,50]
#			-minMapQ [10, 30]
#	Bam types:	raw
#			recalibrated

#### SETUP ####
# directories
proj="singletons"
data_dir="/home/sejoslin/projects/${proj}/data"
raw_dir="${data_dir}/bams_raw"
rec_dir="${data_dir}/bams_recalibrated"
SFS_dir="${data_dir}/results_SFS_untrimmed_wParalogs"

ref="${raw_dir}/DS_history_contigs_250.fasta"
raw_bams="${data_dir}/Rand20.raw.bamlist"
recal_bams="${data_dir}/Rand20.recal.bamlist"
tag="untrim"

mkdir -p ${SFS_dir}
cd ${SFS_dir}

cat << EOT >> README
This directory's purpose is to house SFSs to try to figure out why we are seeing so many singletons in our SFS generated from RAD sequencing data. I will be generating and comparing unfolded SFS

$(date +%D' '%T)
Bamlist for recalibrated SFS can be found at ${recal_bams}
Bamlist for raw SFS can be found at ${raw_bams}
Reference partial genome assembly is located at ${ref}
All bams were processed with GATK Base Recalibration (BQSR) prior to generating SFS.

EOT

## make bamlists
#  recal
cd ${rec_dir}
ls -dalh "$PWD"/*.bam | awk '{print $9}' >> ${recal_bams}
# raw
cd ${raw_dir}
ls -dalh "$PWD"/*.sort.proper.rmdup.bam | awk '{print $9}' >> ${raw_bams}

###################
###   GET SFS   ###
###################

## Index Reference ##
echo "Indexing reference."
samtools faidx ${ref}
sleep 1m
touch ${ref}.fai


## SETUP SAF COMMANDS ##
module load bio 
angsd --version
for i in $(seq 1 2)
do
# minMap10
SFS_1a="angsd -bam ${recal_bams} \
	-out ${SFS_dir}/SFS_recal_${tag}_GL${i}_minMapQ10 \
	-anc ${ref} \
	-GL ${i} \
	-doSaf 1 \
	-minMapQ 10 \
	-minQ 20"
SFS_1b="angsd -bam ${raw_bams} \
        -out ${SFS_dir}/SFS_raw_${tag}_GL${i}_minMapQ10 \
        -anc ${ref} \
        -GL ${i} \
	-doSaf 1 \
        -minMapQ 10 \
	-minQ 20"

# minMap30
SFS_2a="angsd -bam ${recal_bams} \
	-out ${SFS_dir}/SFS_recal_${tag}_GL${i}_minMapQ30 \
	-anc ${ref} \
	-GL ${i} \
	-doSaf 1 \
	-minMapQ 30 \
	-minQ 20"
SFS_2b="angsd -bam ${raw_bams} \
        -out ${SFS_dir}/SFS_raw_${tag}_GL${i}_minMapQ30 \
        -anc ${ref} \
        -GL ${i} \
	-doSaf 1 \
        -minMapQ 30 \
	-minQ 20"

# minMap10 baq1 C50
SFS_3a="angsd -bam ${recal_bams} \
        -out ${SFS_dir}/SFS_recal_${tag}_GL${i}_minMapQ10_baq1_C50 \
        -anc ${ref} \
        -GL ${i} \
	-doSaf 1 \
        -minMapQ 10 \
	-minQ 20 \
	-ref ${ref} \
        -baq 1 \
	-C 50"
SFS_3b="angsd -bam ${raw_bams} \
        -out ${SFS_dir}/SFS_raw_${tag}_GL${i}_minMapQ10_baq1_C50 \
        -anc ${ref} \
        -GL ${i} \
	-doSaf 1 \
        -minMapQ 10 \
	-minQ 20 \
        -ref ${ref} \
        -baq 1 \
	-C 50"

# minMap30 baq1 C50
SFS_4a="angsd -bam ${recal_bams} \
        -out ${SFS_dir}/SFS_recal_${tag}_GL${i}_minMapQ30_baq1_C50 \
        -anc ${ref} \
        -GL ${i} \
	-doSaf 1 \
        -minMapQ 30 \
	-minQ 20 \
	-ref ${ref} \
	-baq 1 \
	-C 50"
SFS_4b="angsd -bam ${raw_bams} \
        -out ${SFS_dir}/SFS_raw_${tag}_GL${i}_minMapQ30_baq1_C50 \
        -anc ${ref} \
        -GL ${i} \
	-doSaf 1 \
        -minMapQ 30 \
	-minQ 20 \
        -ref ${ref} \
        -baq 1 \
	-C 50"

# minMap10 baq2 C50
SFS_5a="angsd -bam ${recal_bams} \
        -out ${SFS_dir}/SFS_recal_${tag}_GL${i}_minMapQ10_baq2_C50 \
        -anc ${ref} \
        -GL ${i} \
	-doSaf 1 \
        -minMapQ 10 \
	-minQ 20 \
        -ref ${ref} \
        -baq 2 \
	-C 50"
SFS_5b="angsd -bam ${raw_bams} \
        -out ${SFS_dir}/SFS_raw_${tag}_GL${i}_minMapQ10_baq2_C50 \
        -anc ${ref} \
        -GL ${i} \
	-doSaf 1 \
        -minMapQ 10 \
	-minQ 20 \
        -ref ${ref} \
        -baq 2 \
	-C 50"

# minMap30 baq2 C50
SFS_6a="angsd -bam ${recal_bams} \
        -out ${SFS_dir}/SFS_recal_${tag}_GL${i}_minMapQ30_baq2_C50 \
        -anc ${ref} \
        -GL ${i} \
	-doSaf 1 \
        -minMapQ 30 \
	-minQ 20 \
        -ref ${ref} \
        -baq 2 \
	-C 50"
SFS_6b="angsd -bam ${raw_bams} \
        -out ${SFS_dir}/SFS_raw_${tag}_GL${i}_minMapQ30_baq2_C50 \
        -anc ${ref} \
        -GL ${i} \
	-doSaf 1 \
        -minMapQ 30 \
	-minQ 20 \
        -ref ${ref} \
        -baq 2 \
	-C 50"

echo "$(date +%D' '%T) : Running ANGSD for:"
# recalibrated bams SFS
echo "Recalibrated bams and GL${i}"

echo ${SFS_1a}
eval ${SFS_1a}

echo ${SFS_2a}
eval ${SFS_2a}

echo ${SFS_3a}
eval ${SFS_3a}

echo ${SFS_4a}
eval ${SFS_4a}

echo ${SFS_5a}
eval ${SFS_5a}

echo ${SFS_6a}
eval ${SFS_6a}
echo ""

# raw bams SFS
echo "Raw bams and GL${i}"

echo ${SFS_1b}
eval ${SFS_1b}

echo ${SFS_2b}
eval ${SFS_2b}

echo ${SFS_3b}
eval ${SFS_3b}

echo ${SFS_4b}
eval ${SFS_4b}

echo ${SFS_5b}
eval ${SFS_5b}

echo ${SFS_6b}
eval ${SFS_6b}
echo ""

## CREATE SFS ##
echo "$(date +%D' '%T) : Generating unfolded site frequency spectrum for:" 
cd ${SFS_dir}

#recalibrated SFSs
echo "Recalibrated bams"
realSFS SFS_recal_${tag}_GL${i}_minMapQ10.saf.idx -maxIter 100 > SFS_recal_${tag}_GL${i}_minMapQ10.sfs
echo "    Completed SFS_recal_${tag}_GL${i}_minMapQ10.sfs"
realSFS SFS_recal_${tag}_GL${i}_minMapQ30.saf.idx -maxIter 100 > SFS_recal_${tag}_GL${i}_minMapQ30.sfs
echo "    Completed SFS_recal_${tag}_GL${i}_minMapQ30.sfs"
realSFS SFS_recal_${tag}_GL${i}_minMapQ10_baq1_C50.saf.idx -maxiter 100 > SFS_recal_${tag}_GL${i}_minMapQ10_baq1_C50.sfs
echo "    Completed SFS_recal_${tag}_GL${i}_minMapQ10_baq1_C50.sfs"
realSFS SFS_recal_${tag}_GL${i}_minMapQ30_baq1_C50.saf.idx -maxiter 100 > SFS_recal_${tag}_GL${i}_minMapQ30_baq1_C50.sfs
echo "    Completed SFS_recal_${tag}_GL${i}_minMapQ30_baq1_C50.sfs"
realSFS SFS_recal_${tag}_GL${i}_minMapQ10_baq2_C50.saf.idx -maxiter 100 > SFS_recal_${tag}_GL${i}_minMapQ10_baq2_C50.sfs
echo "    Completed SFS_recal_${tag}_GL${i}_minMapQ10_baq2_C50.sfs"
realSFS SFS_recal_${tag}_GL${i}_minMapQ30_baq2_C50.saf.idx -maxiter 100 > SFS_recal_${tag}_GL${i}_minMapQ30_baq2_C50.sfs
echo "    Completed SFS_recal_${tag}_GL${i}_minMapQ30_baq2_C50.sfs"
echo ""

# raw SFSs
echo "Raw bams"
realSFS SFS_raw_${tag}_GL${i}_minMapQ10.saf.idx -maxIter 100 > SFS_raw_${tag}_GL${i}_minMapQ10.sfs
echo "    Completed SFS_raw_${tag}_GL${i}_minMapQ10.sfs"
realSFS SFS_raw_${tag}_GL${i}_minMapQ30.saf.idx -maxIter 100 > SFS_raw_${tag}_GL${i}_minMapQ30.sfs
echo "    Completed SFS_raw_${tag}_GL${i}l_minMapQ30.sfs"
realSFS SFS_raw_${tag}_GL${i}_minMapQ10_baq1_C50.saf.idx -maxiter 100 > SFS_raw_${tag}_GL${i}_minMapQ10_baq1_C50.sfs
echo "    Completed SFS_raw_${tag}_GL${i}_minMapQ10_baq1_C50.sfs"
realSFS SFS_raw_${tag}_GL${i}_minMapQ30_baq1_C50.saf.idx -maxiter 100 > SFS_raw_${tag}_GL${i}_minMapQ30_baq1_C50.sfs
echo "    Completed SFS_raw_${tag}_GL${i}_minMapQ30_baq1_C50.sfs"
realSFS SFS_raw_${tag}_GL${i}_minMapQ10_baq2_C50.saf.idx -maxiter 100 > SFS_raw_${tag}_GL${i}_minMapQ10_baq2_C50.sfs
echo "    Completed SFS_raw_${tag}_GL${i}_minMapQ10_baq2_C50.sfs"
realSFS SFS_raw_${tag}_GL${i}_minMapQ30_baq2_C50.saf.idx -maxiter 100 > SFS_raw_${tag}_GL${i}_minMapQ30_baq2_C50.sfs
echo "    Completed SFS_raw_${tag}_GL${i}_minMapQ30_baq2_C50.sfs"
echo ""
echo ""


## PLOT SFS ##
echo "$(date +%D' '%T) : plotting SFS for:"
# recalibrated
echo "Recalibrated bams"
~/scripts/plotSFS.R SFS_recal_${tag}_GL${i}_minMapQ10.sfs
echo "    Completed SFS_recal_${tag}_GL${i}_minMapQ10.sfs"
~/scripts/plotSFS.R SFS_recal_${tag}_GL${i}_minMapQ30.sfs
echo "    Completed SFS_recal_${tag}_GL${i}_minMapQ30.sfs"
~/scripts/plotSFS.R SFS_recal_${tag}_GL${i}_minMapQ10_baq1_C50.sfs
echo "    Completed SFS_recal_${tag}_GL${i}_minMapQ10_baq1_C50.sfs"
~/scripts/plotSFS.R SFS_recal_${tag}_GL${i}_minMapQ30_baq1_C50.sfs
echo "    Completed SFS_recal_${tag}_GL${i}_minMapQ30_baq1_C50.sfs"
~/scripts/plotSFS.R SFS_recal_${tag}_GL${i}_minMapQ10_baq2_C50.sfs
echo "    Completed SFS_recal_${tag}_GL${i}_minMapQ10_baq2_C50.sfs"
~/scripts/plotSFS.R SFS_recal_${tag}_GL${i}_minMapQ30_baq2_C50.sfs
echo "    Completed SFS_recal_${tag}_GL${i}_minMapQ30_baq2_C50.sfs"
echo ""

# raw
echo "Raw bams"
~/scripts/plotSFS.R SFS_raw_${tag}_GL${i}_minMapQ10.sfs
echo "    Completed SFS_raw_${tag}_GL${i}_minMapQ10.sfs"
~/scripts/plotSFS.R SFS_raw_${tag}_GL${i}_minMapQ30.sfs
echo "    Completed SFS_raw_${tag}_GL${i}_minMapQ30.sfs"
~/scripts/plotSFS.R SFS_raw_${tag}_GL${i}_minMapQ10_baq1_C50.sfs
echo "    Completed SFS_raw_${tag}_GL${i}_minMapQ10_baq1_C50.sfs"
~/scripts/plotSFS.R SFS_raw_${tag}_GL${i}_minMapQ30_baq1_C50.sfs
echo "    Completed SFS_raw_${tag}_GL${i}_minMapQ30_baq1_C50.sfs"
~/scripts/plotSFS.R SFS_raw_${tag}_GL${i}_minMapQ10_baq2_C50.sfs
echo "    Completed SFS_raw_${tag}_GL${i}_minMapQ10_baq2_C50.sfs"
~/scripts/plotSFS.R SFS_raw_${tag}_GL${i}_minMapQ30_baq2_C50.sfs
echo "    Completed SFS_raw_${tag}_GL${i}_minMapQ30_baq2_C50.sfs"
echo ""
done


## UPDATE README ## 
echo "Updating README and then you will be done!"
cat << EOT2 >> ${SFS_dir}/README
$(date +%D' '%T) : Creating site allele frequency likelihoods for untrimmed bams based on genotype likelihoods assuming HWE.
    angsd -anc ${ref}  -doSaf 1 -minQ 20 (and run with additional commands stated in each file name).
	Parameters altered:
		-bam [raw,recalibrated]
		-GL [1,2]
		-baq [0,1,2]
		-C [0,50]
		-minMapQ [10,30]
    *.saf.gz == bzgf compressed flat floats (see angsd/doc/formats.pdf documentation) resulting from -doSaf
    *.saf.pos.gz == compressed positions resulting from -doSaf
    *.idx == uncompressed binary index file resulting from -doSaf
    *.sfs == a maximum likelihood estimate of the SFS
    *.pdf == SFS plot resulting from ~/scripts/plotSFS.R <file>.sfs
EOT2
echo ""

end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime
