#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"
export NXF_SYNTAX_PARSER=v1

WORK_DIR=/mnt/meister.data/jsemple/ChIPseq_H3K9me2_20260405
#CONFIG_FILE=/mnt/external.data/MeisterLab/nf-core/unibe_izb_noMACSmodel.config
CONFIG_FILE=/mnt/meister.data/nf-core/unibe_izb.config

nextflow run nf-core/chipseq -profile singularity --input ${WORK_DIR}/sampleSheet.csv --outdir $WORK_DIR -c $CONFIG_FILE --genome WBcel235 --read_length 150 --min_reps_consensus 3 --macs_fdr 0.05 -r 2.1.0

#nextflow run nf-core/chipseq -profile singularity --input sampleSheet.csv --outdir $WORK_DIR -c $CONFIG_FILE --genome WBcel235 --read_length 150 --min_reps_consensus 3 --skip_fastqc --skip_picard_metrics --skip_preseq --skip_plot_profile --skip_plot_fingerprint --skip_spp --skip_igv --macs_fdr 0.05 -r 2.1.0 
