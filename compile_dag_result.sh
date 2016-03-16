#!/bin/sh
export SCALFMM_SIMGRIDOUT='scalfmm.out'
make testBlockedImplicitAlgorithm generateMapping testBlockedMpiAlgorithm compareDAGmapping  -j16
if [ $? -ne 0 ]; then
	exit
fi
mpiexec -n 8 ./Tests/Release/testBlockedMpiAlgorithm -nb 8 -bs 8 -h 3
if [ $? -ne 0 ]; then
	exit
fi

a=`ls $SCALFMM_SIMGRIDOUT\_*`
rm -f $SCALFMM_SIMGRIDOUT
echo $a
for i in $a; do
	echo $i
	cat $i >> $SCALFMM_SIMGRIDOUT
done

cp -f $SCALFMM_SIMGRIDOUT scalfmm_explicit.out
mpiexec -n 8 ./Tests/Release/generateMapping -nb 8 -bs 8 -h 3
mpiexec -n 8 ./Tests/Release/testBlockedImplicitAlgorithm -map mapping -f canard.fma -bs 8 -h 3
if [ $? -ne 0 ]; then
	exit
fi
cp -f scalfmm.out_0 scalfmm_implicit.out
./Tests/Release/compareDAGmapping -e scalfmm_explicit.out -i scalfmm_implicit.out -h 3 > output
