#!/bin/bash
#SBATCH --job-name=danpos
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=danpos_working_merged_ccl4%A_%a.out
#SBATCH --error=danpos_working_merged_ccl4%A_%a.err
#SBATCH --array=0
#SBATCH --mail-type=ALL

# Exit immediately if a command fails or a variable is not set
set -eu

# --- Configuration ---
# The absolute path to the folder containing the raw data subdirectories
BASE_DIR="2026_3_30_CCL4_Nucleosome_Seq/Merge_bam"

# The directory where you want to save the output SAM files
OUTPUT_DIR="2026_3_30_CCL4_Nucleosome_Seq/danpos_output_dpos_merged"

# Load bowtie2 module
module load r/4.4.1

# conda activate danpos2_env
source /cm/shared/apps/python/3.8.x-anaconda/etc/profile.d/conda.sh
conda activate danpos2_env

# Sample list and paths
SAMPLES=($(cat /2026_2_17_SNAP_H4C3_Nucleosome_Seq/sample_list.txt))
SAMPLE_NAME=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Output directory for this sample
SAMPLE_OUT="${OUTPUT_DIR}/Snap_Input_Input_merged_ccl4_oil"
mkdir -p "$SAMPLE_OUT"

# Construct the full paths for the output SAM and log files
Input_BAM="${BASE_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.merged.sorted.bam"

# Now, run DANPOS dpos
/danpos2_env/bin/python /2026_3_9_Danpos/DANPOS3/danpos.py dpos \
${BASE_DIR}/Snap_input_ccl4/SNAP.merged.ccl4.input.trimmed.sorted.bam:${BASE_DIR}/Snap_input_merged/SNAP.merged.input.trimmed.sorted.bam \
-m 1 \
-jd 120 \
--frsz 150 \
-n F \
-u 1e-10 \
-o ${SAMPLE_OUT}
