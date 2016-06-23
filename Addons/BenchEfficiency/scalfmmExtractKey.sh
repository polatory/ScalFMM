#!/bin/bash

if [[ $# -ne 1 ]] ; then
    echo "You must pass a key as parameter"
    return
fi

input=$(cat)
res=`echo "$input" | grep "$1" | cut -d'=' -f2 | cut -d' ' -f2`
echo $res
