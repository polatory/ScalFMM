#!/usr/bin/env bash
## name of job
#SBATCH -J implicit-by-two
#SBATCH -p special
## Resources: (nodes, procs, tasks, walltime, ... etc)
#SBATCH -c 12
#SBATCH --time=0-00:30:00
#SBATCH --exclusive
# # output error message
#SBATCH -e implicit_%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT --mail-user=martin.khannouz@inria.fr
## modules to load for the job
source $HOME/env.sh

## variable for the job
export TREE_HEIGHT=8
export NB_NODE=$(($SLURM_JOB_NUM_NODES*2))
export STARPU_NCPU=12
export STARPU_FXT_PREFIX=$SLURM_JOB_ID
export FINAL_DIR="`pwd`/dir_$SLURM_JOB_ID"
export STARPU_COMM_STATS=1
export STARPU_WATCHDOG_TIMEOUT=30000000
export STARPU_WATCHDOG_CRASH=1
export STARPU_GENERATE_TRACE=0
export STARPU_SILENT=1
export STARPU_MPI_COMM=1
NUMACTL="numactl --interleave=all"
mkdir $FINAL_DIR
echo "my jobID: " $SLURM_JOB_ID > $FINAL_DIR/stdout
echo "Model: " $MODEL >> $FINAL_DIR/stdout
echo "Nb node: " $SLURM_JOB_NUM_NODES >> $FINAL_DIR/stdout
echo "Nb thread: " $STARPU_NCPU >> $FINAL_DIR/stdout
echo "Tree height: " $TREE_HEIGHT >> $FINAL_DIR/stdout
echo "Group size: " $GROUP_SIZE >> $FINAL_DIR/stdout
echo "Algorithm: implicit-two-by" >> $FINAL_DIR/stdout
echo "Max Task Submitted: " $STARPU_LIMIT_MAX_SUBMITTED_TASKS >> $FINAL_DIR/stdout
echo "Min Task Submitted: " $STARPU_LIMIT_MIN_SUBMITTED_TASKS >> $FINAL_DIR/stdout
echo "Total particles: " $NB_PARTICLE >> $FINAL_DIR/stdout
mpiexec -n $NB_NODE --map-by socket $NUMACTL ./Build/Tests/Release/testBlockedImplicitChebyshev -nb $NB_PARTICLE -bs $GROUP_SIZE -h $TREE_HEIGHT -no-validation -$MODEL >> $FINAL_DIR/stdout 2>&1

source $HOME/move_trace.sh
