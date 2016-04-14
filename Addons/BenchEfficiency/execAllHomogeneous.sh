#!/bin/bash

echo "Perform the computation for :"
echo "$SCALFMM_NB particles"
echo "$SCALFMM_H tree height"
echo "Up to $SCALFMM_MAX_NB_CPU CPUs"
echo ""
echo "Using granularities:"
echo "$SCALFMM_BS_CPU_SEQ and $SCALFMM_BS_CPU_PAR"

for (( cpu=1 ; cpu<=$SCALFMM_MAX_NB_CPU ; cpu++)) ; do
    echo ">> CPU = $cpu"

    STARPU_NCPUS=$cpu
    STARPU_NCUDA=0
    ./Tests/Release/testBlockedUnifCudaBench -nb $SCALFMM_NB -h $SCALFMM_H -bs $SCALFMM_BS_CPU_SEQ
    rec_name="trace-nb_$SCALFMM_NB-h_$SCALFMM_H-bs_$SCALFMM_CPU_SEQ-CPU_$cpu.rec"
    mv trace.rec output/$rec_name
    python $SCALFMM_STARPU_DIR/bin/starpu_trace_state_stats.py -t output/$rec_name

    ./Tests/Release/testBlockedUnifCudaBench -nb $SCALFMM_NB -h $SCALFMM_H -bs $SCALFMM_BS_CPU_PAR
    rec_name="trace-nb_$SCALFMM_NB-h_$SCALFMM_H-bs_$SCALFMM_CPU_PAR-CPU_$cpu.rec"
    mv trace.rec output/$rec_name
    python $SCALFMM_STARPU_DIR/bin/starpu_trace_state_stats.py -t output/$rec_name
done

