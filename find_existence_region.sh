#!/bin/bash

####We find the existence region of periodic orbits in the plane eps-r.
###This is done by varying the ratio \tr/\te between 0 and and a given ratio<rho, and then
###calling function find_delta_max.sh. This should generate a file called delta_max.dat containing ratio ###epsilon (which is delta_max because \te=1) and 1-r.
###To initiate the method we allways choose bary for y_0, but t_0 depends on the actual ratio.

n=$(cat "system.dat" | head -1 |  awk '{ print $1 }')
m=$(cat "system.dat" | head -1 |  awk '{ print $2 }')
by0=$(cat "system.dat" | head -1 |  awk '{ print $3 }')
bt0=$(cat "system.dat" | head -1 |  awk '{ print $4 }')
deltaini=$(cat "system.dat" | head -1 |  awk '{ print $5 }')
#ratio=$(cat "system.dat" | head -1 |  awk '{ print $6 }')
omega=$(cat "system.dat" | head -1 |  awk '{ print $7 }')

ratiomin=0.0892
ratiomax=0.0917 ##This has to be smaller than rho, computed in, for instance, comparison_Hog89_4.mw.

cp system.dat system.dat.ini_existence_region

for ratio in `seq $ratiomin 0.0001 $ratiomax`; do
	echo "$n $m $by0 $bt0 $deltaini $ratio $omega" > system.dat.tmp
	./solve_wt0 "system.dat.tmp" "system.dat"
	#rm "system.dat.tmp"
	./find_delta_max.sh

done
