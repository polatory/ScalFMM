#!/usr/bin/env bash
## name of job
#SBATCH -J starpu
#SBATCH -p defq
## Resources: (nodes, procs, tasks, walltime, ... etc)
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --time=02:00:00
#SBATCH --exclusive
# # output error message
#SBATCH -e starpu_%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT --mail-user=martin.khannouz@inria.fr
## modules to load for the job
source $HOME/env.sh

##Setting variable for the job
export TREE_HEIGHT=8
export NB_NODE=$SLURM_JOB_NUM_NODES
export STARPU_NCPU=24
export STARPU_FXT_PREFIX=$SLURM_JOB_ID
export FINAL_DIR="`pwd`/dir_$SLURM_JOB_ID"
export STARPU_WATCHDOG_TIMEOUT=20000000
export STARPU_WATCHDOG_CRASH=1
mkdir $FINAL_DIR
NUMACTL="numactl --interleave=all"

## Write data into an stdout file
echo "my jobID: " $SLURM_JOB_ID > $FINAL_DIR/stdout
echo "Model: " $MODEL >> $FINAL_DIR/stdout
echo "Nb node: " $NB_NODE >> $FINAL_DIR/stdout
echo "Nb thread: " $STARPU_NCPU >> $FINAL_DIR/stdout
echo "Tree height: " $TREE_HEIGHT >> $FINAL_DIR/stdout
echo "Group size: " $GROUP_SIZE >> $FINAL_DIR/stdout
echo "Algorithm: starpu" >> $FINAL_DIR/stdout
echo "Total particles: " $NB_PARTICLE >> $FINAL_DIR/stdout
$NUMACTL ./Build/Tests/Release/testBlockedChebyshev -nb $NB_PARTICLE -bs $GROUP_SIZE -h $TREE_HEIGHT -no-validation -$MODEL >> $FINAL_DIR/stdout

source $HOME/move_trace.sh
