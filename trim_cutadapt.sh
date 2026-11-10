#!/bin/bash
#SBATCH --job-name=trim_fastqc
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-24:00:00
#SBATCH --output=trim_fastqc_%A_%a.out
#SBATCH --error=trim_fastqc_%A_%a.err
#SBATCH --array=0-7
#SBATCH --mail-type=ALL

# Load modules
module load load_rhe17_apps
module load cutadapt/4.5
module load fastqc/0.11.8

# Base directory containing all sample folders
BASE_DIR="2026_3_30_CCL4_Nucleosome_Seq/01.RawData"

# Make a sample list (folder names)
SAMPLES=($(cat "2026_3_30_CCL4_Nucleosome_Seq/sample_list.txt"))
SAMPLE_NAME=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Adapter sequence for Illumina
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER1="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
# Go to sample folder
cd "$BASE_DIR/$SAMPLE_NAME" || { echo "Cannot cd into $BASE_DIR/$SAMPLE_NAME"; exit 1; }

# Input FASTQ files
READ3=$(ls *L8_1.fq.gz)
READ4=$(ls *L8_2.fq.gz)
# Output trimmed filenames
OUT_R1="${SAMPLE_NAME}_L3_trim_1.fq.gz"
OUT_R2="${SAMPLE_NAME}_L3_trim_2.fq.gz"

echo "[$(date)] Trimming adapters for $SAMPLE_NAME..."

# Run cutadapt
    cutadapt -m 20 -a "${ADAPTER}" -A "${ADAPTER1}" \
        -o "${BASE_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_L8_trim_1.fq.gz" -p "${BASE_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_L8_trim_2.fq.gz" \
        ${READ3} ${READ4}

echo "[$(date)] Running FastQC for $SAMPLE_NAME..."

# Run FastQC (your exact command)
fastqc "${BASE_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_L8_trim_1.fq.gz" "${BASE_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_L8_trim_2.fq.gz"

echo "[$(date)] Completed trimming and FastQC for $SAMPLE_NAME."
