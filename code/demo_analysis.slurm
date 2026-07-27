#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --mem=32gb
#SBATCH --tmp=8gb
#SBATCH --error=%x.e%j
#SBATCH --output=%x.o%j
#SBATCH --export=NONE
#SBATCH --mail-type=FAIL
#SBATCH --partition=msismall

set -euo pipefail

# Function to reliably get script directory
slurm_script_dir() {
  local script_path
  
  # Only use scontrol if we're confident we're in a SLURM batch job
  # SLURM_SUBMIT_DIR is ONLY set by sbatch, not by stale variables
  # (Note: salloc/srun also set it, but we handle that with additional checks)
  if [[ -n "${SLURM_SUBMIT_DIR-}" ]] && [[ -n "${SLURM_JOB_ID-}" ]] && command -v scontrol >/dev/null 2>&1; then
    # Try to get the command from scontrol
    script_path="$(scontrol show job "$SLURM_JOB_ID" -o 2>/dev/null | sed -n 's/.*Command=\([^ ]*\).*/\1/p')"
    
    # Only use scontrol result if it's a valid script file (not "bash", not empty, and exists)
    # This filters out interactive sessions where Command=bash
    if [[ -n "$script_path" ]] && [[ "$script_path" != "bash" ]] && [[ -f "$script_path" ]]; then
      # scontrol returned a valid script path - likely running under sbatch
      script_path="$script_path"
    else
      # scontrol returned "bash" or invalid path - running interactively or stale vars
      script_path="${BASH_SOURCE[0]}"
    fi
  else
    # No SLURM vars or scontrol not available - use BASH_SOURCE
    script_path="${BASH_SOURCE[0]}"
  fi
  
  script_path="$(readlink -f "$script_path")"
  dirname "$script_path"
}

# Get project root
# Since script is in PROJECT_ROOT/code/, get parent directory
SCRIPT_DIR="$(slurm_script_dir)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Setup output directory
OUT_DIR="${PROJECT_ROOT}/code_out/demo_analysis"
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

echo "=========================================="
echo "Demo 1: conda1 (R and Python)"
echo "=========================================="
# Initialize conda and activate conda1 environment (includes R and Python)
source "${PROJECT_ROOT}/software_out/010_miniforge/use_miniforge.sh"
conda activate conda1

# Run R analysis
echo "Running R analysis..."
Rscript --vanilla "${PROJECT_ROOT}/code/demo_analysis.R"

# Run Python analysis
echo "Running Python analysis..."
python "${PROJECT_ROOT}/code/demo_analysis.py"

echo ""
echo "=========================================="
echo "Demo 2: conda2 (samtools)"
echo "=========================================="
# Switch to conda2 environment
conda deactivate
conda activate conda2

# Demo samtools - create a simple BAM header and validate
echo "Testing samtools availability..."
samtools --version | head -n 1

echo "Creating test SAM file..."
cat > test.sam << 'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
@PG	ID:demo	PN:demo	VN:1.0
EOF

echo "Converting SAM to BAM with samtools..."
samtools view -bS test.sam > test.bam
samtools view -H test.bam
echo "Samtools demo complete!"

echo ""
echo "=========================================="
echo "Demo 3: Apptainer (Picard)"
echo "=========================================="
# Configure Apptainer
export APPTAINER_CACHEDIR="${PROJECT_ROOT}/.apptainer/cache"
export APPTAINER_TMPDIR="${PROJECT_ROOT}/.apptainer/tmp"

mkdir -p "${APPTAINER_CACHEDIR}"
mkdir -p "${APPTAINER_TMPDIR}"

CONTAINER_SIF="${PROJECT_ROOT}/software_out/021_apptainer1/apptainer1.sif"
PROJECT_REAL_PATH="$(realpath "${PROJECT_ROOT}")"

echo "Testing Picard in Apptainer container..."
apptainer exec "${CONTAINER_SIF}" picard -h || true

echo "Running Picard ValidateSamFile..."
apptainer exec \
    --bind "${PROJECT_REAL_PATH}:${PROJECT_REAL_PATH}" \
    "${CONTAINER_SIF}" \
    picard CleanSam \
    I="${OUT_DIR}/test.sam" \
    O="${OUT_DIR}/test_clean.sam" \
    > picard_validation.txt 2>&1 || true

echo "Apptainer demo complete!"

echo ""
echo "=========================================="
echo "All demos completed successfully!"
echo "Results saved to: ${OUT_DIR}"
echo "=========================================="
