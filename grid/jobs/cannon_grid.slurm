#!/bin/bash
#SBATCH -J at12_grid
#SBATCH -n 1                       # one ATLAS12 process per task
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH -p conroy_priority,shared,sapphire
#SBATCH --mem-per-cpu=4000         # ATLAS12 is ~900 MB RSS; headroom for the
                                   # line lists it memory-maps during SELECTLINES
#SBATCH -o logs/at12_grid_%A_%a.log
#SBATCH -e logs/at12_grid_%A_%a.log
#
# Build the ATLAS12+SYNTHE model grid as a SLURM job array.
#
#   sbatch --array=0-99 cannon_grid.slurm
#
# Each task takes a strided slice of the node list, so tasks never collide and
# the union covers the grid exactly once. The build is idempotent: a node whose
# .spec is verified complete is skipped, so resubmitting the same array picks
# up whatever is missing. Re-run until `--dry-run` reports 0 to compute.
#
# Optional: --export=ALL,NTASKS=100,ONLY_MISSING=1
#   ONLY_MISSING=1 slices the MISSING nodes instead of all nodes, which
#   balances the work on a resubmission after a partial run.

date; hostname; pwd
module purge
module load python/3.10.9-fasrc01
source activate c3k

# --- paths ---------------------------------------------------------------
PROJECT_DIR=$SCRATCH/conroy_lab/Lab/$USER/ngsl
export ATLAS12=$SCRATCH/conroy_lab/Lab/$USER/atlas12
cd $PROJECT_DIR

NTASKS=${NTASKS:-$SLURM_ARRAY_TASK_COUNT}
ONLY_MISSING=${ONLY_MISSING:-0}
EXTRA=""
[ "$ONLY_MISSING" = "1" ] && EXTRA="--only-missing"

echo "ATLAS12=$ATLAS12"
echo "task $SLURM_ARRAY_TASK_ID of $NTASKS"

# One worker per task: ATLAS12 is single-threaded and memory-bandwidth bound,
# so packing several per task gains nothing. Let SLURM schedule the width.
python3 grid/build_grid.py --workers 1 \
        --task $SLURM_ARRAY_TASK_ID --ntasks $NTASKS $EXTRA

date
