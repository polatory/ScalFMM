#!/usr/bin/env bash
## name of job
#SBATCH -J explicit_50M_4N
#SBATCH -p court
## Resources: (nodes, procs, tasks, walltime, ... etc)
#SBATCH -N 4
#SBATCH -c 24
#SBATCH --time=04:00:00
# # output error message
#SBATCH -e explicit_chebyshev_50M_10_node%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT --mail-user=martin.khannouz@inria.fr
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
export NB_PARTICLE_PER_NODE=12500000
export STARPU_FXT_PREFIX=$SLURM_JOB_ID
export FINAL_DIR="`pwd`/dir_$SLURM_JOB_ID"
export STARPU_COMM_STATS=1
$NUMACTL=numactl --interleave=all
mkdir $FINAL_DIR
echo "my jobID: " $SLURM_JOB_ID > $FINAL_DIR/stdout
echo "Model: cube" >> $FINAL_DIR/stdout
echo "Nb node: " $NB_NODE >> $FINAL_DIR/stdout
echo "Nb thread: " $STARPU_NCPU >> $FINAL_DIR/stdout
echo "Tree height: " $TREE_HEIGHT >> $FINAL_DIR/stdout
echo "Group size: " $GROUP_SIZE >> $FINAL_DIR/stdout
echo "Algorithm: explicit" >> $FINAL_DIR/stdout
echo "Particle per node: " $NB_PARTICLE_PER_NODE >> $FINAL_DIR/stdout
echo "Total particles: " $(($NB_PARTICLE_PER_NODE*$NB_NODE)) >> $FINAL_DIR/stdout
mpiexec -n $NB_NODE $NUMACTL ./Build/Tests/Release/testBlockedMpiChebyshev -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT -no-validation >> $FINAL_DIR/stdout

#Create argument list for starpu_fxt_tool
cd $FINAL_DIR
list_fxt_file=`ls ../$STARPU_FXT_PREFIX*`

#Clean to only keep trace.rec
mkdir fxt
for i in $list_fxt_file; do
	mv $i fxt
done
cd ..

##Move the result into a directory where all result goes
mv $FINAL_DIR jobs_result
