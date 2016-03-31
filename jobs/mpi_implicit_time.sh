#!/usr/bin/env bash
## name of job
#SBATCH -J Implicit_MPI_time
#SBATCH -p special
## Resources: (nodes, procs, tasks, walltime, ... etc)
#SBATCH -N 40
# #  standard output message
#SBATCH -o batch%j.out
# # output error message
#SBATCH -e batch%j.err
module purge
module load slurm
module add compiler/gcc/5.3.0 tools/module_cat/1.0.0 intel/mkl/64/11.2/2016.0.0
. /home/mkhannou/spack/share/spack/setup-env.sh
spack load fftw
spack load hwloc
spack load openmpi
spack load starpu@svn-trunk
## modules to load for the job
export GROUP_SIZE=500
export TREE_HEIGHT=5
export NB_NODE=$SLURM_JOB_NUM_NODES
export NB_PARTICLE_PER_NODE=100000
echo "=====my job informations ===="
echo "Node List: " $SLURM_NODELIST
echo "my jobID: " $SLURM_JOB_ID
echo "Nb node: " $NB_NODE
echo "Particle per node: " $NB_PARTICLE_PER_NODE
echo "In the directory: `pwd`"
rm -f canard.fma > /dev/null 2> /dev/null
mpiexec -n $NB_NODE ./Build/Tests/Release/testBlockedMpiAlgorithm -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT > loutre
cat loutre | grep Executing
cat loutre | grep Average
sleep 10
mpiexec -n $NB_NODE ./Build/Tests/Release/testBlockedImplicitAlgorithm -f canard.fma -bs $GROUP_SIZE -h $TREE_HEIGHT > loutre
cat loutre | grep Executing
cat loutre | grep Average
rm -f canard.fma > /dev/null 2> /dev/null
sleep 10
mpiexec -n $NB_NODE ./Build/Tests/Release/testBlockedMpiChebyshev -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT > loutre
cat loutre | grep Executing
cat loutre | grep Average
sleep 10
mpiexec -n $NB_NODE ./Build/Tests/Release/testBlockedImplicitChebyshev -f canard.fma -bs $GROUP_SIZE -h $TREE_HEIGHT > loutre
cat loutre | grep Executing
cat loutre | grep Average
