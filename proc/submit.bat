#!/bin/bash
set -e
#
# Define working directories
# --------------------------
#
dir=$HOME/LDAS_ISBA
path1=${dir}/data_out
path2=${dir}/src1
path3=${dir}/data_in
path4=${dir}/proc
if [[ $3 == "TB" ]] 
then
path2=${dir}/src2
fi
expid=$2
exptype=$1
#
#  Execute the model
#  ----------------
#
cp ${path4}/namelist fort.8 
#
# Reference and open loop runs
#
if [[ $exptype == "REF" || $exptype == "OL" ]]
then
${path2}/main-ref
fi
#
# OI - Extended Kalman filter and 2DVAR
#
if [[ $exptype == "OI_EC" || $exptype == "OI_MF" || $exptype == "EKF" || $exptype == "2DVAR" ]]
then
${path2}/main-ref
#
# Ensemble Kalman filter
#
fi
if [[ $exptype == "ENKF" ]]
then
${path2}/main-ref
fi
#
#  Store results in appropriate directory according to run type
#  ------------------------------------------------------------
#
#  a) Optimum interpolation
#
if [[ $exptype == "OI_EC" || $exptype == "OI_MF" ]]
then
for n in 22 23 24 25 26 71 72 73 55 56
do
mv fort.$n ${path1}/FIC${n}_${exptype}${expid}.dat
done
fi
#
#  b) Extended Kalman filter or simplified 2D-Var
#
if [[  $exptype == "EKF" || $exptype == "2DVAR" ]]
then
for n in 22 23 24 25 26 55 
do
mv fort.$n ${path1}/FIC${n}_${exptype}${expid}.dat
done
fi
#
#  c) Ensemble Kalman filter 
#
if [[ $exptype == "ENKF" ]]
then
for n in 22 23 24 25 26 27 
do
mv fort.$n ${path1}/FIC${n}_${exptype}${expid}.dat
done
fi
#
#  d) Open loop 
#
if [[ $exptype == "OL" ]]
then
for n in 22 23 24 25 26
do
mv fort.$n ${path1}/FIC${n}_${exptype}${expid}.dat
done
rm -f fort.30
fi
#
#  e) Reference run + storage of synthetic data
#
if [[ $exptype == "REF"  ]]
then
for n in 22 23 24 25 26 
do
mv fort.$n ${path1}/FIC${n}_${exptype}${expid}.dat
done
mv fort.30 ${path3}/OBS_SITE${expid}.dat
fi
#
#  Remove useless temporary files
#  ------------------------------
#
rm -f INPUT_OI.DAT *~ 
rm -f fort.400 fort.98 fort.97 fort.8
