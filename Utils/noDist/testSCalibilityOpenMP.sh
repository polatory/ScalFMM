#!/bin/bash -v
project_dir=$HOME/Dev/src/ScalFMM/scalfmmN/BuildOMP/run1
#
# PlaFRIM environment
#
#
#
#MaxCore=20
#  `cat /proc/cpuinfo |grep processor |wc -l`
# tester si linux alors
MaxCore=`getconf _NPROCESSORS_ONLN`
HOST=`hostname`
echo $HOST
EXEC="../Examples/Release/ChebyshevInterpolationFMM"
REP="Test"
FILEPERF="output-scab-Cheb-10M"
FILE="../../Data/unitCubeXYZQ20k.bfma"
FILE="casTest-20000.fma"
FILE="/projets/scalfmm/benchPlafrim/unifCube_2_10000000.bfma"
echo "=================================================================="
echo "Run "
echo "   pgm:   ${EXEC}"
echo "   file:  ${FILE}"
echo "   args  -depth 7 -subdepth 4 "
echo "   "
echo "   Host:    $HOST"
echo "   MaxCore: $MaxCore"
echo "   File:    $FILE"
echo "   "
echo "   Projectdir:    $project_dir"
echo "=================================================================="
#
mkdir -f $project_dir
cd $project_dir
#
pwd
export OMP_PROC_BIND=true
export KMP_AFFINITY=verbose,scatter
echo $FILEPERF-${HOST}.out 
echo  "# Core TIME  ENERGY Pot_0 Pot_5000000 Pot_9999999"> $FILEPERF-${HOST}.out 
P0=1
P1=1
P2=2
for l in `seq 1 $MaxCore `;
do
    OUTPUT=${FILEPERF}-${HOST}-${l}.out
    echo "Running  per = " ${l}
    $EXEC --datehost --flags -show-params -f $FILE  -depth 7 -subdepth 4  -t $l >  $OUTPUT
    #
    TIME=`grep "@Algorithm"  $OUTPUT | awk '{print $4}'`
    Energy=`grep "Energy"  $OUTPUT | awk '{print $2}'`
    P0=`grep "Index 0  potential"  $OUTPUT | awk '{print $4}'`
    P1=`grep "Index 5000000  potential"  $OUTPUT | awk '{print $4}'`
    P2=`grep "Index 9999999  potential"  $OUTPUT | awk '{print $4}'`
    echo " " $l "   " $TIME "  " $Energy "  " $P0 "  " $P1  "  " $P2  >>  $FILEPERF-${HOST}.out
done


