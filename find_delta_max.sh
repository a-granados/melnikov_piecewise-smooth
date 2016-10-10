#!/bin/bash


# Function to evaluate the absolute value of a number
abs () {
# Check for numeric input
if expr $1 + 0 2>/dev/null 1>&2 ; then
# Is the number negative?
#if [ $1 -lt 0 ] ; then
#echo `expr 0 - $1`
if [$(echo "scale=20; $1 >0" | bc) -eq 1 ]; then
echo $(echo "scale=20; -$1" | bc)
else
echo $1
fi
return 0 # OK
else
return 1 # Not a number
fi
}



##This is the same function as follow_orbit.sh but it stops when the delta can not be increased
##The criteria that we are using here is to stop when the step is reduced fails_tol consecutive times.
##The result is written in a file called delta_max.dat and contains the following info:
####    $ratio $prevdelta $ratio*$prevdelta $by0 $inicondt $oldby0 $oldcond

n=$(cat "system.dat" | head -1 |  awk '{ print $1 }')
m=$(cat "system.dat" | head -1 |  awk '{ print $2 }')
by0=$(cat "system.dat" | head -1 |  awk '{ print $3 }')
bt0=$(cat "system.dat" | head -1 |  awk '{ print $4 }')
deltaini=$(cat "system.dat" | head -1 |  awk '{ print $5 }')
ratio=$(cat "system.dat" | head -1 |  awk '{ print $6 }')
omega=$(cat "system.dat" | head -1 |  awk '{ print $7 }')
#h=$(echo "scale=10; ($deltafinal-$deltaini)/($npoints+1)" | bc)

#hini=0.1
tol=0.001
fails_tol=4
deltafinal=10
delta=$deltaini
hini=$deltaini

cp system.dat system.dat.copy_delta_max

./find_nm_ini_cond > process
##We need the method to converge at this first execution.

prevh=$hini
prevdelta=$delta
delta=$(echo "scale=20; $prevdelta+$hini" | bc)
oldcondt=$(cat "impacts.dat" | head -1 |  awk '{ print $2 }')
inicondt=$(cat "impacts.dat" | head -1 |  awk '{ print $2 }')

fails_count=0




while [ $(echo "scale=20; $delta < $deltafinal" | bc ) -eq 1 ]; do
	echo "$n $m $by0 $bt0 $delta $ratio $omega" > "system.dat"
	./use_last_initial_conditions.sh
	newcondt=$(cat "new_impacts.dat" | head -1 |  awk '{ print $2 }')
	if [ "$newcondt" == "nan" ] || [ $(echo "scale=20; $newcondt - $oldcondt >$tol" | bc) -eq 1 ]; then
		if [ $fails_count -ge $fails_tol ]; then
			####We assume here that te=1, so eps=delta.
			echo "$ratio $prevdelta $(echo "scale=20; $ratio*$prevdelta" | bc) $by0 $inicondt $oldby0 $oldcondt" >> delta_max.dat
			exit	
		else
			echo "Failed: reducing the step"
#			sleep 2s
			prevh=$(echo "scale=20; $prevh/2" | bc)
			delta=$(echo "scale=20; $prevdelta+$prevh" | bc)
			fails_count=$(($fails_count+1))
		fi
	elif [ "$newcondt" == "nan" ] || [ $(echo "scale=20; $oldcondt - $newcondt>$tol" | bc) -eq 1 ]; then


		if [ $fails_count -ge $fails_tol ]; then
			####We assume here that te=1, so eps=delta.
			echo "$ratio $prevdelta $(echo "scale=20; $ratio*$prevdelta" | bc) $by0 $inicondt $oldby0 $oldcont" >> delta_max.dat
			exit	
		else
			echo "Failed: reducing the step"
#			sleep 2s
			prevh=$(echo "scale=20; $prevh/2" | bc)
			delta=$(echo "scale=20; $prevdelta+$prevh" | bc)
			fails_count=$(($fails_count+1))
		fi
	else
		fails_count=0
		oldby0=$(cat "new_impacts.dat" | head -1 |  awk '{ print $1 }')
		oldcondt=$(cat "new_impacts.dat" | head -1 |  awk '{ print $2 }')
		mv -f new_impacts.dat impacts.dat
		#cp impacts.dat "impacts_$delta.dat"
		#cp orbit.dat "orbit_$delta.dat"
		prevdelta=$delta
		delta=$(echo "scale=20; $prevdelta+$prevh" | bc)
		clear
		echo "ratio=$ratio"
		echo "prevh=$prevh"
		echo "delta=$delta"
	fi

done
#cp system.dat.backup system.dat


