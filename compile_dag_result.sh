#!/bin/sh
export SCALFMM_SIMGRIDOUT='scalfmm.out'
export GROUP_SIZE=8
export TREE_HEIGHT=6
export NB_NODE=8
export NB_PARTICLE_PER_NODE=$((`awk "BEGIN{print 8 ** ($TREE_HEIGHT-1)}"` / $NB_NODE))
echo "GROUP_SIZE=$GROUP_SIZE"
echo "TREE_HEIGHT=$TREE_HEIGHT"
echo "NB_NODE=$NB_NODE"
echo "NB_PARTICLE_PER_NODE=$NB_PARTICLE_PER_NODE"

make testBlockedImplicitAlgorithm generateMapping testBlockedMpiAlgorithm compareDAGmapping  -j $((`nproc`*2))
if [ $? -ne 0 ]; then
	exit
fi
mpiexec -n $NB_NODE ./Tests/Release/testBlockedMpiAlgorithm -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT
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
mpiexec -n $NB_NODE ./Tests/Release/generateMapping -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT
mpiexec -n $NB_NODE ./Tests/Release/testBlockedImplicitAlgorithm -map mapping -f canard.fma -bs $GROUP_SIZE -h $TREE_HEIGHT
if [ $? -ne 0 ]; then
	exit
fi
cp -f scalfmm.out_0 scalfmm_implicit.out
./Tests/Release/compareDAGmapping -e scalfmm_explicit.out -i scalfmm_implicit.out -h $TREE_HEIGHT > output
