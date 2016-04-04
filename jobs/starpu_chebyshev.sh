#!/usr/bin/env bash
## name of job
#SBATCH -J chebyshev_50M_1_node
#SBATCH -p longq
## Resources: (nodes, procs, tasks, walltime, ... etc)
#SBATCH -N 1
#SBATCH -c 24
# #  standard output message
#SBATCH -o chebyshev_50M_1_node%j.out
# # output error message
#SBATCH -e chebyshev_50M_1_node%j.err
#SBATCH --mail-type=ALL --mail-user=martin.khannouz@inria.fr
module purge
module load slurm
module add compiler/gcc/5.3.0 tools/module_cat/1.0.0 intel/mkl/64/11.2/2016.0.0
. /home/mkhannou/spack/share/spack/setup-env.sh
spack load fftw
spack load hwloc
spack load openmpi
spack load starpu@svn-trunk+fxt
## modules to load for the job
export GROUP_SIZE=500
export TREE_HEIGHT=8
export NB_NODE=$SLURM_JOB_NUM_NODES
export STARPU_NCPU=24
export NB_PARTICLE_PER_NODE=50000000
export STARPU_FXT_PREFIX=`pwd`/
echo "=====my job informations ===="
echo "Node List: " $SLURM_NODELIST
echo "my jobID: " $SLURM_JOB_ID
echo "Nb node: " $NB_NODE
echo "Particle per node: " $NB_PARTICLE_PER_NODE
echo "Total particles: " $(($NB_PARTICLE_PER_NODE*$NB_NODE))
echo "In the directory: `pwd`"
rm -f canard.fma > /dev/null 2>&1
mpiexec -n $NB_NODE ./Build/Tests/Release/testBlockedChebyshev -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT -no-validation | grep Kernel
