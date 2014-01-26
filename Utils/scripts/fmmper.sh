#!/bin/bash
#
#PBS -N MG_water
# 
# preciser le temps en heures, minutes, secondes
#PBS -l walltime=01:59:00
# nombre de noeuds et de coeurs
#PBS -l nodes=1:ppn=12
#PBS -l naccesspolicy=singlejob
#
# PlaFRIM environment
w
source $HOME/.bashrc
#
FMMPER_EXE="Tests/Release/testSphericalEwalAlgorithm"
FILE="../Data/forceNacl_2000_dlpolyPer.txt"
FILE="../Data/forceNacl_128_dlpolyPer.txt"
cd ~/Dev/src/ScalFMM/scalfmmT/BuildGCC
#
# RUN 1 influence de la taille de la boite sur la precision des calculs
#    Regarde size, Nx ny nz energie total, dipole et le temps de calcul
#
PER_SIZE="-1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18"
#
rm -f FMM_PER.out
 echo  "# PerSize DOMAIN_SIZE  ENERGY"> FMM_PER.out
pwd
for l in $PER_SIZE 
do
    echo "Running  per = " ${l}
    echo  "$FMMPER_EXE -f $FILE  -h 4 -sh 2 -per $l  >  $FILE.out "
    $FMMPER_EXE -f $FILE  -h 4 -sh 2 -per $l  >  $FILE.out 
    DOMAINSIZE=`grep "Simulated box:"  $FILE.out | awk '{print $3}'`
    ENERGY=`grep "Energy FMM="  $FILE.out |awk '{ print $3 }'`
    echo " " $l   $DOMAINSIZE "   " $ENERGY " "   >>  FMM_PER.out
done

