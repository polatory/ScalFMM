#!/bin/sh
export SCALFMM_SIMGRIDOUT='scalfmm.out'
export GROUP_SIZE=50
export TREE_HEIGHT=5
export NB_NODE=4
#export NB_PARTICLE_PER_NODE=$(( (`awk "BEGIN{print 8 ** ($TREE_HEIGHT-1)}"` / $NB_NODE) ))
export NB_PARTICLE_PER_NODE=5000
export STARPU_NCPU=1
export STARPU_FXT_PREFIX=`pwd`/

echo "GROUP_SIZE=$GROUP_SIZE"
echo "TREE_HEIGHT=$TREE_HEIGHT"
echo "NB_NODE=$NB_NODE"
echo "NB_PARTICLE_PER_NODE=$NB_PARTICLE_PER_NODE"

#Compile only what we need
time make testBlockedImplicitChebyshev testBlockedMpiChebyshev testBlockedImplicitAlgorithm testBlockedMpiAlgorithm compareDAGmapping  -j $((`nproc`*2))
if [ $? -ne 0 ]; then
	exit
fi

test_kernel()
{
	#Execute explicit mpi version
	mpiexec -n $NB_NODE ./Tests/Release/testBlockedMpiAlgorithm -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT
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
		rm -f $i
	done

	#Get task information
	cp -f $SCALFMM_SIMGRIDOUT scalfmm_explicit.out

	#Execute implicit version
	mpiexec -n $NB_NODE ./Tests/Release/testBlockedImplicitAlgorithm -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT
	if [ $? -ne 0 ]; then
		echo
		echo " /!\\Error on implicit"
		echo
		exit
	fi

	#Get task information
	cp -f $SCALFMM_SIMGRIDOUT\_0 scalfmm_implicit.out
	rm -f $SCALFMM_SIMGRIDOUT\_*


	#Compare DAGs
	./Tests/Release/compareDAGmapping -e scalfmm_explicit.out -i scalfmm_implicit.out -h $TREE_HEIGHT > output
}
chebyshev_kernel()
{
	mpiexec -n $NB_NODE ./Tests/Release/testBlockedMpiChebyshev -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT
	if [ $? -ne 0 ]; then
		echo
		echo " /!\\Error on explicit Chebyshev"
		echo
		exit
	fi

	##Aggregate task information from explicit execution
	a=`ls $SCALFMM_SIMGRIDOUT\_*`
	rm -f $SCALFMM_SIMGRIDOUT
	for i in $a; do
		cat $i >> $SCALFMM_SIMGRIDOUT
		rm -f $i
	done

	#Get task information
	cp -f $SCALFMM_SIMGRIDOUT scalfmm_explicit.out

	mpiexec -n $NB_NODE ./Tests/Release/testBlockedImplicitChebyshev -nb $NB_PARTICLE_PER_NODE -bs $GROUP_SIZE -h $TREE_HEIGHT
	if [ $? -ne 0 ]; then
		echo
		echo " /!\\Error on implicit Chebyshev"
		echo
		exit
	fi

	#Get task information
	cp -f $SCALFMM_SIMGRIDOUT\_0 scalfmm_implicit.out
	rm -f $SCALFMM_SIMGRIDOUT\_*


	#Compare DAGs
	./Tests/Release/compareDAGmapping -e scalfmm_explicit.out -i scalfmm_implicit.out -h $TREE_HEIGHT > narval
}

test_kernel
chebyshev_kernel

