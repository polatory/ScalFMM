#!/usr/bin/env bash
## name of job
#SBATCH -J simple-mpi
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
export TREE_HEIGHT=8
export NB_NODE=$SLURM_JOB_NUM_NODES
export FINAL_DIR="`pwd`/dir_$SLURM_JOB_ID"
export STARPU_FXT_PREFIX=$SLURM_JOB_ID
NUMACTL="numactl --interleave=all"
mkdir $FINAL_DIR
echo "my jobID: " $SLURM_JOB_ID > $FINAL_DIR/stdout
echo "Model: " $MODEL >> $FINAL_DIR/stdout
echo "Nb node: " $NB_NODE >> $FINAL_DIR/stdout
echo "Nb thread: 24" >> $FINAL_DIR/stdout
echo "Tree height: " $TREE_HEIGHT >> $FINAL_DIR/stdout
echo "Group size: 1" >> $FINAL_DIR/stdout
echo "Algorithm: simple-mpi" >> $FINAL_DIR/stdout
echo "Total particles: " $NB_PARTICLE >> $FINAL_DIR/stdout
mpiexec -n $NB_NODE $NUMACTL ./Build/Tests/Release/testFmmAlgorithmThreadProc -nb $NB_PARTICLE -h $TREE_HEIGHT -$MODEL >> $FINAL_DIR/stdout 2>&1

source $HOME/move_trace.sh
