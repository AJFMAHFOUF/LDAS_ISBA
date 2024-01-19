#!/bin/bash
set -ex
# id of current experiment
lab=$2
# id of reference and open loop simulations
expid=2c
exptype=$1
echo $exptype
exp=${exptype}${lab}
dir=$HOME/LDAS_ISBA
path1=${dir}/data_out
path2=${dir}/graphs
cat << EOF > diff.f90
 program diff
 do 
  read (10,*,end=400) xj,x1,y1
  read (11,*) xj,x2,y2
  write (12,*) xj,x1-x2,y1-y2
 enddo
400 continue
 end program diff
EOF
gfortran diff.f90
for nvar in 22 23
do
cp  ${path1}/FIC${nvar}_${exp}.dat fort.10
cp  ${path1}/FIC${nvar}_REF${expid}.dat fort.11
a.out
cp fort.12 ${path1}/FIC${nvar}_dif1.dat
cp  ${path1}/FIC${nvar}_OL${expid}.dat fort.10
cp  ${path1}/FIC${nvar}_REF${expid}.dat fort.11
a.out
cp fort.12 ${path1}/FIC${nvar}_dif2.dat
done 
cat << EOF > fic1
set term post portrait noenhanced monochrome "Helvetica" 12
set output '${path2}/FIGURE_${exp}.ps'
set style data lines
set multiplot
set origin 0.00,0.75
set size 1.00,0.25
set xrange [0:31]
set xlabel "Days"
set title " SCREEN-LEVEL TEMPERATURE ($exptype $lab) "
set ylabel "[K]"
plot "${path1}/FIC22_dif1.dat" using 1:2 title 'ASSIM-TRUTH' lw 3, \
"${path1}/FIC22_dif2.dat" using 1:2 title 'OPEN LOOP-TRUTH' lw 2
set origin 0.00,0.50
set size 1.00,0.25
set title " SCREEN-LEVEL RELATIVE HUMIDITY ($exptype $lab) "
set ylabel "[]"
plot "${path1}/FIC22_dif1.dat" using 1:3 title 'ASSIM-TRUTH' lw 3, \
"${path1}/FIC22_dif2.dat" using 1:3 title 'OPEN LOOP-TRUTH' lw 2
set origin 0.00,0.25
set size 1.00,0.25
set title "SURFACE TEMPERATURE ($exptype $lab) "
set ylabel "[K]"
plot "${path1}/FIC23_dif1.dat" using 1:2 title 'ASSIM-TRUTH' lw 3, \
"${path1}/FIC23_dif2.dat" using 1:2 title 'OPEN LOOP-TRUTH' lw 2
set origin 0.00,0.00
set size 1.00,0.25
set title "MEAN SOIL TEMPERATURE (${exptype} ${lab}) "
set ylabel "[K]"
plot "${path1}/FIC23_dif1.dat" using 1:3 title 'ASSIM-TRUTH' lw 3, \
"${path1}/FIC23_dif2.dat" using 1:3 title 'OPEN LOOP-TRUTH' lw 2
set nomultiplot
set multiplot
set origin 0.00,0.66
set size 1.00,0.33
set title "SUPERFICIAL SOIL MOISTURE CONTENT ($exptype $lab) "
set ylabel "[m3/m3]"
set xrange [0:31]
plot "${path1}/FIC23_${exp}.dat" using 1:4 title 'ASSIM' lw 3, \
"${path1}/FIC23_REF${expid}.dat" using 1:4 title 'TRUTH' lw 2, \
"${path1}/FIC23_OL${expid}.dat" using 1:4 title 'OPEN LOOP' lw 2 
set origin 0.00,0.33
set size 1.00,0.33
set title "MEAN SOIL MOISTURE CONTENT ($exptype $lab) "
set ylabel "[m3/m3]"
plot "${path1}/FIC23_${exp}.dat" using 1:5 title 'ASSIM' lw 3 , \
"${path1}/FIC23_REF${expid}.dat" using 1:5 title 'TRUTH' lw 2 , \
"${path1}/FIC23_OL${expid}.dat" using 1:5 title 'OPEN LOOP' lw 2 
set origin 0.00,0.00
set size 1.00,0.33
set title "ACCUMULATED EVAPORATION ($exptype $lab)"
set ylabel "[mm]"
plot "${path1}/FIC25_${exp}.dat" using 1:2 title 'ASSIM' lw 3, \
"${path1}/FIC25_REF${expid}.dat" using 1:2 title 'TRUTH' lw 2, \
"${path1}/FIC25_OL${expid}.dat" using 1:2 title 'OPEN LOOP' lw 2 
set nomultiplot
EOF
gnuplot fic1
rm a.out fic1 diff.f90 ${path1}/FIC23_dif2.dat  ${path1}/FIC23_dif1.dat ${path1}/FIC22_dif1.dat ${path1}/FIC22_dif2.dat fort.10 fort.11 fort.12 
