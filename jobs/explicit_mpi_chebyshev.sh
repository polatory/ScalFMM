#!/usr/bin/env bash
## name of job
#SBATCH -J explicit_chebyshev_50M_10_node
#SBATCH -p longq
## Resources: (nodes, procs, tasks, walltime, ... etc)
#SBATCH -N 10
#SBATCH -c 24
# #  standard output message
#SBATCH -o explicit_chebyshev_50M_10_node%j.out
# # output error message
#SBATCH -e explicit_chebyshev_50M_10_node%j.err
#SBATCH --mail-type=ALL --mail-user=martin.khannouz@inria.fr
## modules to load for the job
module purge
module load slurm
module add compiler/gcc/5.3.0 tools/module_cat/1.0.0 intel/mkl/64/11.2/2016.0.0
. /home/mkhannou/spack/share/spack/setup-env.sh
spack load fftw
spack load hwloc
spack load openmpi
spack load starpu@svn-trunk+fxt
## variable for the job
export GROUP_SIZE=500
export TREE_HEIGHT=8
export NB_NODE=$SLURM_JOB_NUM_NODES
export STARPU_NCPU=24
export NB_PARTICLE_PER_NODE=5000000
export STARPU_FXT_PREFIX=`pwd`/
echo "===== Explicit MPI ===="
echo "my jobID: " $SLURM_JOB_ID
echo "Model: cube"
echo "Nb node: " $NB_NODE
echo "Nb thread: " $STARPU_NCPU
echo "Tree height: " $TREE_HEIGHT
echo "Group size: " $GROUP_SIZE
echo "Algorithm: explicit"
echo "Particle per node: " $NB_PARTICLE_PER_NODE
echo "Total particles: " $(($NB_PARTICLE_PER_NODE*$NB_NODE))
mpiexec -n $NB_NODE ./Build/Tests/Release/testBlockedMpiChebyshev -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT -no-validation | grep Average

