#!/usr/bin/env bash

MY_PTH=$PWD

function plot_freq {
    new_name=${1%%.*}
    new_name=${new_name##*_}
    echo "
	FILES='$@'
	FILES1='$1'
	FILES2='$2'
	print FILES
	ST=900
	EV=100
	set title \"FREQ\"
#	set terminal wxt size 1800,900; unset key;
#	set multiplot layout 1,4
	plot [:] for [ data in FILES ] data u 2:((\$21+\$22)/(\$10)) w lp t 'eTr',\
	for [ data in FILES ] data u 2:((\$23+\$24)/(\$14)) w lp t 'eA',\
	for [ data in FILES ] data u 2:((\$25+\$26)/(\$18)) w lp t 'eB'
#	for [ data in FILES ] data u 2:((\$19+\$20)/(\$6) ) w lp t 'eV',\

	pause -1
#	unset multiplot
    " > to_plot
    gnuplot to_plot
}



files=`ls total_dR2.dat`
echo $files

#plot_N $files

for i in $files; do
    j=${i%.avg}.step
#    echo $i $j
#    while true; do
    plot_freq $i
#    end=="q"; read end;
#    if [ "$end" == "q" ]; then break; fi;
#    done
done
