#!/bin/sh
export SCALFMM_SIMGRIDOUT='scalfmm.out'
export STARPU_STATS=1
export GROUP_SIZE=32
export TREE_HEIGHT=6
export NB_NODE=8
export NB_PARTICLE_PER_NODE=$((`awk "BEGIN{print 8 ** ($TREE_HEIGHT-1)}"` / $NB_NODE))
echo "GROUP_SIZE=$GROUP_SIZE"
echo "TREE_HEIGHT=$TREE_HEIGHT"
echo "NB_NODE=$NB_NODE"
echo "NB_PARTICLE_PER_NODE=$NB_PARTICLE_PER_NODE"

#Compile only what we need
make testBlockedImplicitAlgorithm generateMapping testBlockedMpiAlgorithm compareDAGmapping  -j $((`nproc`*2))
if [ $? -ne 0 ]; then
	exit
fi

#Execute explicit mpi version
mpiexec -n $NB_NODE ./Tests/Release/testBlockedMpiAlgorithm -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT
if [ $? -ne 0 ]; then
	exit
fi

#Aggregate task information from explicit execution
a=`ls $SCALFMM_SIMGRIDOUT\_*`
rm -f $SCALFMM_SIMGRIDOUT
echo $a
for i in $a; do
	cat $i >> $SCALFMM_SIMGRIDOUT
done

#Get task information
cp -f $SCALFMM_SIMGRIDOUT scalfmm_explicit.out
#Get task information from fxt
#a=`ls $STARPU_FXT_PREFIX/../installprof*`
#ARG_FXT_TOO=""
#for i in $a; do
	#ARG_FXT_TOO="$ARG_FXT_TOO -i $i"
#done
#$STARPU_DIR/bin/starpu_fxt_tool $ARG_FXT_TOO
#cp $STARPU_FXT_PREFIX/../trace.rec explicit.rec
#rm -f $STARPU_FXT_PREFIX/../installprof*

#Generate mapping for implicite version
mpiexec -n $NB_NODE ./Tests/Release/generateMapping -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT > /dev/null
#Execute implicit version
mpiexec -n $NB_NODE ./Tests/Release/testBlockedImplicitAlgorithm -map mapping -f canard.fma -bs $GROUP_SIZE -h $TREE_HEIGHT
if [ $? -ne 0 ]; then
	exit
fi

#Get task information
cp -f scalfmm.out_0 scalfmm_implicit.out
#Get task information from fxt
a=`ls $STARPU_FXT_PREFIX/../installprof*`
#ARG_FXT_TOO=""
#for i in $a; do
	#ARG_FXT_TOO="$ARG_FXT_TOO -i $i"
#done
#$STARPU_DIR/bin/starpu_fxt_tool $ARG_FXT_TOO
#cp $STARPU_FXT_PREFIX/../trace.rec implicit.rec
#rm -f $STARPU_FXT_PREFIX/../installprof*

#Compare DAGs
./Tests/Release/compareDAGmapping -e scalfmm_explicit.out -i scalfmm_implicit.out -h $TREE_HEIGHT > output


