#!/usr/bin/env bash
#SBATCH --job-name=OpenMPAlgo_perf
#SBATCH  --mail-user=olivier.coulaud@inria.fr
#SBATCH -m a
#
#SBATCH -p court
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#
data_dir=/projets/scalfmm/data/openmpImprovements ;
project_dir=$HOME/Dev/src/ScalFMM/scalfmm ;

#
# PlaFRIM environment
#
source $HOME/Config/bashrc.bash
module add compiler/gcc/5.1.0 compiler/intel/64/2016_beta; module li

#
EXEC="Utils/Release/ChebyshevInterpolationCmpAlgo"
FILEPERF="RES-Cheb-openmp"
OPTION="-show-info -show-compile"
#
cd $project_dir/BuildIntel2016

#
#
PER_SIZE=`cat /proc/cpuinfo |grep processor |wc -l`
PER_SIZE="4 8 12 16 20 24"
HOST=`hostname`
NBTEST=3  # 0..NBTEST
DISTRIB=(unitcube unitsphere prolate plummer)
DEPTH=(6           7           9            8)
ALGO=(basic task sectiontask)
#
pwd
#export OMP_PROC_BIND=true
export KMP_AFFINITY=scatter
NUM=`git rev-list HEAD --count`
echo $DISTRIB
for dist in `seq 0 $NBTEST`;
do
    echo "Distribution:  " $dist ${DISTRIB[$dist]}   "  depth octree "  ${DEPTH[$dist]}
    FILE=${data_dir}/${DISTRIB[$dist]}-1M-ref.bfma
    mkdir  ${DISTRIB[$dist]} 
    #
    #    Loop on algorithm 
    for a in `seq 0 2`;
    do 
echo  "# Core GlabalTIME P2PTIME M2LTIME  ENERGY PotentialError ForceError"> $FILEPERF-${DISTRIB[$dist]}-${ALGO[$a]}.out
	#  One computation
	for l in $PER_SIZE ;
	do
	    OUTPUT=${FILEPERF}-${DISTRIB[$dist]}-${ALGO[$a]}-${l}.out
	    echo "ScalFMM git version: " $NUM > $OUTPUT
	    echo "Running  on " ${l} " threads with algo "   ${ALGO[$a]}
	    $EXEC ${OPTION} -f $FILE  -depth  ${DEPTH[$dist]}  -subdepth 4  -t $l  -cmp  -algo ${ALGO[$a]}>>  $OUTPUT
	    #
	    TIME=`grep "@Algorithm"  $OUTPUT | awk '{print $4}'`
	    TIMEP2P=`grep "P2P"  $OUTPUT | awk '{print $4}'`
	    TIMEM2L=`grep "M2L"  $OUTPUT | grep "second" | awk '{print $2}'`
	    Energy=`grep "Energy"  $OUTPUT | awk '{print $3}'`
	    P1=`grep " Potential "  $OUTPUT | awk '{print $7}'`
	    P2=`grep " F "  $OUTPUT | awk '{print $7}'`
	    echo " " $l "   " $TIME "  " $TIMEP2P "  " $TIMEM2L "  " $Energy "     " $P1  "  " $P2
	    echo " " $l "   " $TIME "  " $TIMEP2P "  " $TIMEM2L "  " $Energy "     " $P1  "  " $P2  >>  $FILEPERF-${DISTRIB[$dist]}-${ALGO[$a]}.out
	    mv  $OUTPUT  ${DISTRIB[$dist]} 
	done
    echo "----------------------------------------------------------"
    done
    echo "-------------------------------------------------------------------"
done
echo " ----------------------   FINISH -------------------------------"
