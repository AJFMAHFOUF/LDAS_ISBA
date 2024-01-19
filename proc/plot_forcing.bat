#!/bin/bash
set -evx
dir1=$HOME/LDAS_ISBA
FC=gfortran
path1=${dir1}/data_in
path2=${dir1}/graphs
cat << EOF > conv.f90
program conv
do 
  read (10,*,end=100) i,rg,ra,pr,ta,ua,va,ps,qa
  day = real(i)/24.0
  umod = sqrt(ua**2 + va**2)
  pr = pr*86400.
  qa = qa*1000.0
  ta = ta - 273.15
  write (11,*) day,rg,ra,pr,ta,umod,qa
enddo
100 continue  
end program conv
EOF
$FC conv.f90
# number of chosen location
nb=2
# selected month
month=07
# label on title
lab="Site ${nb} - Month ${month}"
cp ${path1}/FORCINGb_SITE${nb}_${month}.dat fort.10
a.out
mv fort.11 ${path1}/FORCINGa_SITE${nb}_${month}.dat
set -evx
cat << EOF > fic1
set term post portrait monochrome "Helvetica" 12
set output '${path2}/FORCING_SITE${nb}_${month}.ps'
set style data lines
set multiplot
set origin 0.00,0.80
set size 1.00,0.20
set xrange [0:31]
set xlabel "Days"
set title "SOLAR AND THERMAL RADIATION (${lab})"
set ylabel "[W/m2]"
plot "${path1}/FORCINGa_SITE${nb}_${month}.dat" using 1:2 title  'SW'lw 3, \
"${path1}/FORCINGa_SITE${nb}_${month}.dat" using 1:3 title  'LW' lw 3 
set origin 0.00,0.60
set size 1.00,0.20
set style data boxes
set title "PRECIPITATION (${lab})"
set ylabel "[mm/day]"
plot "${path1}/FORCINGa_SITE${nb}_${month}.dat" using 1:4 notitle  lw 3 
set origin 0.00,0.40
set size 1.00,0.20
set style data lines
set title "AIR TEMPERATURE (${lab})"
set ylabel "[K]"
plot "${path1}/FORCINGa_SITE${nb}_${month}.dat" using 1:5 notitle  lw 3 
set origin 0.00,0.20
set size 1.00,0.20
set style data lines
set title "AIR SPECIFIC HUMIDITY (${lab})"
set ylabel "[g/kg]"
plot "${path1}/FORCINGa_SITE${nb}_${month}.dat" using 1:7 notitle  lw 3 
set origin 0.00,0.00
set size 1.00,0.20
set ylabel "[m/s]"
set title "WIND SPEED (${lab})"
set ylabel "[K]"
plot "${path1}/FORCINGa_SITE${nb}_${month}.dat" using 1:6 notitle  lw 3 
set nomultiplot
EOF
gnuplot fic1
rm fic1 fort.10 a.out ${path1}/FORCINGa* conv.f90
