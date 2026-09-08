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

# NTASKS is the WIDTH OF THE PARTITION, not the number of tasks in this
# submission. They are the same for a full array and differ when resubmitting a
# subset, which is the only time you must set it by hand:
#
#   sbatch --array=0-99                        NTASKS=100 automatically
#   sbatch --array=5,7,9 --export=ALL,NTASKS=100    resubmitting 3 of those 100
#
# Getting it wrong is not benign. Tasks take pool[task::NTASKS]; with
# --array=5,7,9 the default NTASKS=3 makes task 5 compute pool[5::3], which
# aliases onto a 3-way split. Each task would claim ~568 nodes instead of ~17
# and several would compute the SAME node concurrently, writing one .spec path
# from two processes.
NTASKS=${NTASKS:-$SLURM_ARRAY_TASK_COUNT}

if [ -z "$NTASKS" ] || [ "$SLURM_ARRAY_TASK_ID" -ge "$NTASKS" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID is not in [0,$NTASKS)."
    echo "       Resubmitting a subset of an array? Pass the ORIGINAL width:"
    echo "         sbatch --array=<ids> --export=ALL,NTASKS=<original> $0"
    exit 1
fi
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
