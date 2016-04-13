#!/bin/bash

function GetExecTime() {
    >&2 echo "[LOG]     Exec $1 $2"
    local res_output=`$1 "$2"`
    >&2 echo "[LOG]     Try to find $3"
    local time_result=`echo "$res_output" | grep "$3" | cut -d'=' -f2 | cut -d's' -f1`
    # >&2 echo "[LOG] output : $res_output"
    >&2 echo "[LOG]     Done in $time_result"
    echo $time_result
}

function isLower() {
    awk -v n1=$1 -v n2=$2 'BEGIN{ if (n1<n2) exit 0; exit 1}'
}

function getMin(){
    if isLower $1 $2 ; then
        echo $1
    else
        echo $2
    fi
}

if [[ "$#" -ne "3" ]] ; then
    echo "Error you must provide 3 parameters:"
    echo "$0 test-commnand-line starting-bs ending-bs"
    return
fi

echo "You ask to find the best bs for:"
echo "Command: $1"
echo "From $2 to $3"

outputfile=./bs_bench.data

echo "# BS TIME" > $outputfile

left_bs="$2"
left_exec_time=`GetExecTime "$1" "$left_bs" "@EXEC TIME"`
echo "$left_bs $left_exec_time" >> $outputfile


best_bs=$left_bs
best_exec_time="$left_exec_time"

right_bs="$3"
right_exec_time=`GetExecTime "$1" "$right_bs" "@EXEC TIME"`
echo "$right_bs $right_exec_time" >> $outputfile

if isLower $right_exec_time $best_exec_time ; then
    best_bs=$right_bs
    best_exec_time="$right_exec_time"
fi

for (( depth=0; depth<10; depth++ ))
do
    if [[ $left_bs -ge $right_bs ]] ; then
        break;
    fi

    >&2 echo "[LOG] Depth $depth ..."

    sub_left_bs=$(( ($right_bs-$left_bs)/3 + $left_bs ))
    sub_left_exec_time=`GetExecTime "$1" "$sub_left_bs" "@EXEC TIME"`
    echo "$sub_left_bs $sub_left_exec_time" >> $outputfile

    if isLower $sub_left_exec_time $best_exec_time ; then
        best_bs=$sub_left_bs
        best_exec_time="$sub_left_exec_time"
    fi

    sub_right_bs=$(( 2*($right_bs-$left_bs)/3 + $left_bs ))
    sub_right_exec_time=`GetExecTime "$1" "$sub_right_bs" "@EXEC TIME"`
    echo "$sub_right_bs $sub_right_exec_time" >> $outputfile

    if isLower $sub_right_exec_time $best_exec_time ; then
        best_bs=$sub_right_bs
        best_exec_time="$sub_right_exec_time"
    fi
    
    min_from_left=`getMin $left_exec_time $sub_left_exec_time`
    min_from_right=`getMin $right_exec_time $sub_right_exec_time`

    if isLower $min_from_left $min_from_right ; then 
        right_bs=$sub_right_bs
        right_exec_time="$sub_right_exec_time"
    else
        left_bs=$sub_left_bs
        left_exec_time="$sub_left_exec_time"
    fi
done

echo "@BEST BS = $best_bs (in = $best_exec_time s)"

