#!/bin/bash

##Here I'm reading r and epsilon from system.dat and by0 and bt0 from the impacts given by the last execution of find_nm_ini_cond.c. Then I run it again.


omega=$(cat "system.dat" | head -1 |  awk '{ print $7 }')
n=$(cat "system.dat" | head -1 |  awk '{ print $1 }')
m=$(cat "system.dat" | head -1 |  awk '{ print $2 }')
eps=$(cat "system.dat" | head -1 |  awk '{ print $5 }')
r=$(cat "system.dat" | head -1 |  awk '{ print $6 }')

by0=$(cat "impacts.dat" | head -1 |  awk '{ print $1 }')
bt0=$(cat "impacts.dat" | head -1 |  awk '{ print $2 }')


echo "$n $m $by0 $bt0 $eps $r $omega" > "system_tmp.dat"

./find_nm_ini_cond "system_tmp.dat" "output.dat" "new_impacts.dat" > "process"
tail -2 "process"
./simulate_impacts_linear "new_impacts.dat" "system_tmp.dat" "orbit.dat"
