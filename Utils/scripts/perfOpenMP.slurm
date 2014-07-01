#!/usr/bin/env bash
#SBATCH --job-name=OpenMP_perf
#SBATCH  --mail-user=olivier.coulaud@inria.fr
#SBATCH -m a
#   
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#
project_dir=/projets/scalfmm/scalfmm-tests ;
project_dir=$HOME/Dev/src/ScalFMM/scalfmm ;

#
# PlaFRIM environment
#
module add build/cmake/2.8.12.2 compiler/intel/2013_sp1.2.144
#
EXEC="Examples/Release/ChebyshevInterpolationAlgorithm"
FILEPERF="RES_Chebyshev-openmp"
FILE="unifCube_2.fma"
#
cd $project_dir/Build/
#
#
PER_SIZE=`cat /proc/cpuinfo |grep processor |wc -l`
HOST=`hostname`
#
pwd
#export OMP_PROC_BIND=true
export KMP_AFFINITY=verbose,scatter
echo  "# Core TIME  ENERGY Pot_0 Pot_5000000 Pot_9999999"> $FILEPERF-${HOST}.out 

for l in `seq 1 $PER_SIZE `;
do
    OUTPUT=${FILEPERF}-${HOST}-${l}.out
    echo "Running  per = " ${l}
    $EXEC -f $FILE  -depth 7 -subdepth 4  -t $l >  $OUTPUT
    #
    TIME=`grep "@Algorithm"  $OUTPUT | awk '{print $4}'`
    Energy=`grep "Energy"  $OUTPUT | awk '{print $2}'`
    P0=`grep "Index 0  potential"  $OUTPUT | awk '{print $4}'`
    P1=`grep "Index 5000000  potential"  $OUTPUT | awk '{print $4}'`
    P2=`grep "Index 9999999  potential"  $OUTPUT | awk '{print $4}'`
    echo " " $l "   " $TIME "  " $Energy "  " $P0 "  " $P1  "  " $P2  >>  $FILEPERF-${HOST}.out
done
