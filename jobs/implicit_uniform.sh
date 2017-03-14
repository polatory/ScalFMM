#!/usr/bin/env bash
## name of job
#SBATCH -J implicit
#SBATCH -p special
## Resources: (nodes, procs, tasks, walltime, ... etc)
#SBATCH -c 24
#SBATCH --time=0-00:30:00
#SBATCH --exclusive
# # output error message
#SBATCH -e implicit_%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT --mail-user=martin.khannouz@inria.fr
## modules to load for the job
source $HOME/env.sh

## variable for the job
export TREE_HEIGHT=6
export NB_NODE=$SLURM_JOB_NUM_NODES
export STARPU_NCPU=24
export STARPU_FXT_PREFIX=$SLURM_JOB_ID
export FINAL_DIR="`pwd`/dir_$SLURM_JOB_ID"
export STARPU_COMM_STATS=1
export STARPU_WATCHDOG_TIMEOUT=20000000
export STARPU_WATCHDOG_CRASH=0
export STARPU_GENERATE_TRACE=0
export STARPU_SILENT=1
export STARPU_MPI_COMM=1
NUMACTL="numactl --interleave=all"
mkdir $FINAL_DIR
echo "my jobID: " $SLURM_JOB_ID > $FINAL_DIR/stdout
echo "Model: cube" >> $FINAL_DIR/stdout
echo "Nb node: " $NB_NODE >> $FINAL_DIR/stdout
echo "Nb thread: " $STARPU_NCPU >> $FINAL_DIR/stdout
echo "Tree height: " $TREE_HEIGHT >> $FINAL_DIR/stdout
echo "Group size: " $GROUP_SIZE >> $FINAL_DIR/stdout
echo "Algorithm: implicit_uniforme" >> $FINAL_DIR/stdout
echo "Total particles: " $NB_PARTICLE >> $FINAL_DIR/stdout
mpiexec -n $NB_NODE $NUMACTL ./Build/Tests/Release/testBlockedImplicitUniform -nb $NB_PARTICLE -bs $GROUP_SIZE -h $TREE_HEIGHT >> $FINAL_DIR/stdout 2>&1

source $HOME/move_trace.sh
