#!/bin/bash


n=$(cat "system.dat" | head -1 |  awk '{ print $1 }')
m=$(cat "system.dat" | head -1 |  awk '{ print $2 }')
by0=$(cat "system.dat" | head -1 |  awk '{ print $3 }')
bt0=$(cat "system.dat" | head -1 |  awk '{ print $4 }')
deltaini=$(cat "system.dat" | head -1 |  awk '{ print $5 }')
deltafinal=0.6
ratio=$(cat "system.dat" | head -1 |  awk '{ print $6 }')
omega=$(cat "system.dat" | head -1 |  awk '{ print $7 }')
#h=$(echo "scale=10; ($deltafinal-$deltaini)/($npoints+1)" | bc)

hini=0.1
tol=0.01

delta=$deltaini

cp system.dat system.dat.copy

./find_nm_ini_cond > process
##We need the method to converge at this first execution.

prevh=$hini
prevdelta=$delta
delta=$(echo "scale=20; $prevdelta+$hini" | bc)
oldcondt=$(cat "impacts.dat" | head -1 |  awk '{ print $2 }')

while [ $(echo "scale=20; $delta < $deltafinal" | bc ) -eq 1 ]; do
#while [ 22 -eq 1 ]; do
	echo "$n $m $by0 $bt0 $delta $ratio $omega" > "system.dat"
	./use_last_initial_conditions.sh
	newcondt=$(cat "new_impacts.dat" | head -1 |  awk '{ print $2 }')
	echo "$newcondt"

	if [ "$newcondt" == "nan" ] || [ $(echo "scale=20; $newcondt - $oldcondt >$tol" | bc) -eq 1 ]; then
		echo "Failed: reducing the step"
		sleep 2s
		prevh=$(echo "scale=20; $prevh/2" | bc)
		delta=$(echo "scale=20; $prevdelta+$prevh" | bc)
	elif [ "$newcondt" == "nan" ] || [ $(echo "scale=20; $oldcondt - $newcondt >$tol" | bc) -eq 1 ]; then
		echo "Failed: reducing the step"
		sleep 2s
		prevh=$(echo "scale=20; $prevh/2" | bc)
		delta=$(echo "scale=20; $prevdelta+$prevh" | bc)
	else
		oldcondt=$(cat "new_impacts.dat" | head -1 |  awk '{ print $2 }')
		mv -f new_impacts.dat impacts.dat
		cp impacts.dat "impacts_$delta.dat"
		cp orbit.dat "orbit_$delta.dat"
		prevdelta=$delta
		delta=$(echo "scale=20; $prevdelta+$prevh" | bc)
		clear
		echo "prevh=$prevh"
	fi

done
#cp system.dat.backup system.dat

exit


eps=$epsini
echo "$n $m $by0 $bt0 $eps $omega" > "system_$n-$m-$eps.dat"
./find_nm_ini_cond "system_$n-$m-$eps.dat" "output_$n-$m-$eps.dat" "impacts_$n-$m-$eps.dat" > "tmp_$n-$m-$eps.dat"
tail "tmp_$n-$m-$eps.dat"
./simulate_impacts_linear "impacts_$n-$m-$eps.dat" "system_$n-$m-$eps.dat" "orbit_$n-$m-$eps.dat"
rm "system_$n-$m-$eps.dat" "output_$n-$m-$eps.dat" "impacts_$n-$m-$eps.dat" "tmp_$n-$m-$eps.dat"
preveps=$eps
eps=$(echo "scale=10; $preveps+$h" | bc)
epsini2=$eps

for eps in `seq $epsini2 $h $epsfinal`; do
#while [ $eps -le $epsfinal ]; do
	#by0=$(cat "impacts_$n-$m-$preveps.dat" | head -1 |  awk '{ print $1 }')
	#ty0=$(cat "impacts_$n-$m-$preveps.dat" | head -1 |  awk '{ print $2 }')
	echo "$n $m $by0 $bt0 $eps $omega" > "system_$n-$m-$eps.dat"
	./find_nm_ini_cond "system_$n-$m-$eps.dat" "output_$n-$m-$eps.dat" "impacts_$n-$m-$eps.dat" > "tmp_$n-$m-$eps.dat"
	tail "tmp_$n-$m-$eps.dat"
	./simulate_impacts_linear "impacts_$n-$m-$eps.dat" "system_$n-$m-$eps.dat" "orbit_$n-$m-$eps.dat"
	rm "system_$n-$m-$eps.dat" "output_$n-$m-$eps.dat" "impacts_$n-$m-$eps.dat" "tmp_$n-$m-$eps.dat"
	preveps=$eps
done


