#!/bin/sh
export SCALFMM_SIMGRIDOUT='scalfmm.out'
export GROUP_SIZE=500
export TREE_HEIGHT=5
export NB_NODE=16
#export NB_PARTICLE_PER_NODE=$(( (`awk "BEGIN{print 8 ** ($TREE_HEIGHT-1)}"` / $NB_NODE) ))
export NB_PARTICLE_PER_NODE=15000

echo "GROUP_SIZE=$GROUP_SIZE"
echo "TREE_HEIGHT=$TREE_HEIGHT"
echo "NB_NODE=$NB_NODE"
echo "NB_PARTICLE_PER_NODE=$NB_PARTICLE_PER_NODE"

#Compile only what we need
make testBlockedImplicitChebyshev testBlockedMpiChebyshev testBlockedImplicitAlgorithm testBlockedMpiAlgorithm compareDAGmapping  -j $((`nproc`*2))
if [ $? -ne 0 ]; then
	exit
fi

#Execute explicit mpi version
sleep 10
mpiexec -n $NB_NODE ./Tests/Release/testBlockedMpiAlgorithm -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT 2>/dev/null
if [ $? -ne 0 ]; then
	echo
	echo " /!\\Error on explicit"
	echo
	exit
fi

#Aggregate task information from explicit execution
a=`ls $SCALFMM_SIMGRIDOUT\_*`
rm -f $SCALFMM_SIMGRIDOUT
for i in $a; do
	cat $i >> $SCALFMM_SIMGRIDOUT
done

#Get task information
cp -f $SCALFMM_SIMGRIDOUT scalfmm_explicit.out

#Generate mapping for implicite version
#mpiexec -n $NB_NODE ./Tests/Release/generateMapping -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT > /dev/null
#Execute implicit version
sleep 10
mpiexec -n $NB_NODE ./Tests/Release/testBlockedImplicitAlgorithm -f canard.fma -bs $GROUP_SIZE -h $TREE_HEIGHT 2>/dev/null
if [ $? -ne 0 ]; then
	echo
	echo " /!\\Error on implicit"
	echo
	exit
fi

#Get task information
cp -f scalfmm.out_0 scalfmm_implicit.out

#Compare DAGs
./Tests/Release/compareDAGmapping -e scalfmm_explicit.out -i scalfmm_implicit.out -h $TREE_HEIGHT > output

sleep 10
mpiexec -n $NB_NODE ./Tests/Release/testBlockedMpiChebyshev -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT 2>/dev/null
if [ $? -ne 0 ]; then
	echo
	echo " /!\\Error on explicit Chebyshev"
	echo
	exit
fi
sleep 10
mpiexec -n $NB_NODE ./Tests/Release/testBlockedImplicitChebyshev -f canard.fma -bs $GROUP_SIZE -h $TREE_HEIGHT 2>/dev/null
if [ $? -ne 0 ]; then
	echo
	echo " /!\\Error on implicit Chebyshev"
	echo
	exit
fi
