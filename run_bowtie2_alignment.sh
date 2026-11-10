#!/bin/bash
#SBATCH --job-name=bowtie2_alignment
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=bowtie2_alignment_%A_%a.out
#SBATCH --error=bowtie_2_alignment_%A_%a.err
#SBATCH --array=0-7
#SBATCH --mail-type=ALL

# =============================================================================
# This script aligns paired-end FASTQ reads to a reference genome using
# Bowtie2. It automatically finds all read pairs and includes robust
# logic to skip samples that have already completed successfully.
# =============================================================================

# Exit immediately if a command fails or a variable is not set
set -eu

# --- Configuration ---
# The absolute path to the folder containing the raw data subdirectories
BASE_DIR="2026_3_30_CCL4_Nucleosome_Seq/01.RawData"

# The directory where you want to save the output SAM files
OUTPUT_DIR="2026_3_30_CCL4_Nucleosome_Seq/bowtie2_alignment"
mkdir -p "$OUTPUT_DIR"

# The directory to store the Bowtie2 summary log files
LOG_DIR="2026_3_30_CCL4_Nucleosome_Seq/bowtie2_logs"
mkdir -p "$LOG_DIR"

# --- !! CORRECTED PART !! ---
# The absolute path to the Bowtie2 index prefix.
# This points to the parent folder of 01.RawData, where the index is located.
INDEX_NAME="2026_2_5_H3_3_CR_NIH_H4A_KO/mm39_bowtie2_index/mm39_bowtie_build_index"

# Load bowtie2 module
module load bowtie2/2.5.2

# Load samtools
module load samtools/1.22.1

# Sample list and paths 
SAMPLES=($(cat 2026_3_30_CCL4_Nucleosome_Seq/sample_list.txt))
SAMPLE_NAME=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Go to sample folder
cd "$BASE_DIR/$SAMPLE_NAME" || { echo "Cannot cd into $BASE_DIR/$SAMPLE_NAME"; exit 1; }

# Input FASTQ files
READ3=$(ls *trim_1.fq.gz)
READ4=$(ls *trim_2.fq.gz)

# Output directory for this sample
SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE_NAME}"
mkdir -p "$SAMPLE_OUT"

# Output directory for logs
LOG_OUT="${LOG_DIR}/${SAMPLE_NAME}"
mkdir -p "$LOG_OUT"

# Construct the full paths for the output SAM and log files
OUTPUT_SAM="${SAMPLE_OUT}/${SAMPLE_NAME}_aligned.sam"
LOG_FILE_8="${LOG_OUT}/${SAMPLE_NAME}_bowtie2.L8.log"

# Run Bowtie2 and redirect stderr (2>) to the log file
bowtie2 --very-sensitive -X 1000 \
  -p $SLURM_CPUS_PER_TASK \
  -x "$INDEX_NAME" \
  -1 "$READ3" \
  -2 "$READ4" \
  2> "$LOG_FILE_8" | \
samtools view -@ $SLURM_CPUS_PER_TASK -b - | \
samtools sort -@ $SLURM_CPUS_PER_TASK \
  -o "${SAMPLE_OUT}/${SAMPLE_NAME}.L8.trimmed.sorted.bam" -
  
samtools index "${SAMPLE_OUT}/${SAMPLE_NAME}.L8.trimmed.sorted.bam"
