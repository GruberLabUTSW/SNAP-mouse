#!/bin/bash
#SBATCH --job-name=danpos
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=120:00:00
#SBATCH --output=danpos_profile_merged%A_%a.out
#SBATCH --error=danpos_profile_merged%A_%a.err
#SBATCH --array=0
#SBATCH --mail-user=Lavanya.Vumma@utsouthwestern.edu
#SBATCH --mail-type=ALL

# Exit immediately if a command fails or a variable is not set
set -eu

# --- Configuration ---
# The absolute path to the folder containing the raw data subdirectories
BASE_DIR="/work/InternalMedicine/s222632/2026_2_17_SNAP_H4C3_Nucleosome_Seq/danpos_output_dpos_merged"

# The directory where you want to save the output SAM files
OUTPUT_DIR="/work/InternalMedicine/s222632/2026_2_17_SNAP_H4C3_Nucleosome_Seq/danpos_output_dpos_merged/profile"

# module load R to prevent errors rpy2
module load r/4.4.1

# conda activate danpos2_env
source /cm/shared/apps/python/3.8.x-anaconda/etc/profile.d/conda.sh
conda activate danpos2_env

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Now, run DANPOS profile
/work/InternalMedicine/s222632/danpos2_env/bin/python /work/InternalMedicine/s222632/2026_3_9_Danpos/DANPOS3/danpos.py profile \
${BASE_DIR}/pooled/SNAP.merged.input.trimmed.sorted.Fnor.smooth.wig,${BASE_DIR}/pooled/SNAP.merged.IP.trimmed.sorted.Fnor.smooth.wig \
  --wigfile_aliases input,SNAP_H4C3_IP \
  --bed3file_paths /work/InternalMedicine/s222632/2026_2_17_SNAP_H4C3_Nucleosome_Seq/danpos_output_dpos_merged/2026_4_16_ref_positions.bed \
  --genomic_sites center \
  --flank_up 3000 \
  --flank_dn 3000 \
  --bin_size 20 \
  --name SNAP_H4_profile
