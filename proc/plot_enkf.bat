#!/bin/bash
path=$HOME/LDAS_ISBA
path1=${path}/data_out
path2=${path}/graphs
lab=$1
exp=$2
set -evx
cat << EOF > fic1
set term post landscape noenhanced monochrome "Helvetica" 16
set output '${path2}/FC_ERRORS_${lab}_${exp}.ps'
set style data lines
set xrange [0:31]
#set yrange [0:0.014]
set xlabel "Days"
set title " DERIVED FORECAST ERRORS FROM ${lab} ${exp}"
set ylabel "[m3/m3]"
plot "${path1}/FIC27_${lab}${exp}.dat" using 1:2 title "WG" lw 4 , \
"${path1}/FIC27_${lab}${exp}.dat" using 1:3 title "W2" lw 4
EOF
gnuplot fic1
rm fic1 
