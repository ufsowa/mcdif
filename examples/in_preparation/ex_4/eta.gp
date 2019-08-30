#!/usr/bin/env bash

MY_PTH=$PWD
rm *.jpeg
#cd ..

function plot_d {
    new_name=${1%%.*}
    new_name=${new_name##*_}
#    if [ $new_name -le 100 ]; then
    echo "
#	plot [0.1:0.9][:] '$1' u 4:((\$5+\$6)) w l
#	plot '$1' u 4:5 w lp,'$1' u 4:6 w lp,'$1' u 4:(1.0/\$8) w lp

#	plot [0.1:0.9] '$1' u (1.0-\$4):(-(16.0*\$6/\$5)/\$2) w lp
	plot [:] '$1' u (\$3):(\$7/\$2) w lp
#	plot [0.1:0.9] '$1' u 4:((\$7 - 1.0/\$8)) w lp

	pause -1
    " > to_plot
    gnuplot to_plot

#$MY_PTH/$new_name.jpeg
#    fi
}

function plot_Ns {
    new_name=${1%%.*}
    new_name=${new_name##*_}
    echo "
	FILES='$@'
	FILES1='$1'
	FILES2='$2'
	print FILES
	ST=900
	EV=100
	set title \"Concentration\"
#	set terminal wxt size 1800,900; unset key;
#	set multiplot layout 1,4
	plot [:] for [ data in FILES ] data every :EV::ST u 5:(\$6/(\$6 + \$7 + \$8 + \$9)) w lp t 'V',\
	for [ data in FILES ] data every :EV::ST u 5:(\$7/(\$6 + \$7 + \$8+ \$9)) w lp t 'Tr',\
	for [ data in FILES ] data every :EV::ST u 5:(\$8/(\$6 + \$7 + \$8+ \$9)) w lp t 'A',\
	for [ data in FILES ] data every :EV::ST u 5:(\$9/(\$6 + \$7 + \$8+ \$9)) w lp t 'B'
	pause -1
#	unset multiplot

    " > to_plot
    gnuplot to_plot
}


function plot_Js {
    new_name=${1%%.*}
    new_name=${new_name##*_}
    echo "
	FILES='$@'
	FILES1='$1'
	FILES2='$2'
	print FILES
	ST=900
	EV=100
	set title \"Flux\"
#	set terminal wxt size 1800,900; unset key;
#	set multiplot layout 1,4
	plot [:] for [ data in FILES ] data every :EV::ST u 5:(\$10/\$1) w lp t 'JV',\
	for [ data in FILES ] data every :EV::ST u 5:(\$11/\$1) w lp t 'JTr',\
	for [ data in FILES ] data every :EV::ST u 5:(\$12/\$1) w lp t 'JA',\
	for [ data in FILES ] data every :EV::ST u 5:(\$13/\$1) w lp t 'JB'
	pause -1
#	unset multiplot
    " > to_plot
    gnuplot to_plot
}


function plot_eta {
    new_name=${1%%.*}
    new_name=${new_name##*_}
    echo "
	FILES='$@'
	FILES1='$1'
	FILES2='$2'
	print FILES
	ST=900
	EV=100
	set title \"LRO\"
#	set terminal wxt size 1800,900; unset key;
#	set multiplot layout 1,4
	plot [:] for [ data in FILES ] data u 2:( abs(\$3-\$4)/(\$3+\$4) ) w lp t 'eV',\
	for [ data in FILES ] data u 2:(abs(\$5-\$6)/(\$5+\$6)) w lp t 'eTr',\
	for [ data in FILES ] data u 2:(abs(\$7-\$8)/(\$7+\$8)) w lp t 'eA',\
	for [ data in FILES ] data u 2:(abs(\$9-\$10)/(\$9+\$10)) w lp t 'eB'
	pause -1
#	unset multiplot
    " > to_plot
    gnuplot to_plot
}



files=`ls *N.dat`
echo $files

#plot_N $files

for i in $files; do
    j=${i%.avg}.step
#    echo $i $j
#    while true; do
    plot_eta $i
#    end=="q"; read end;
#    if [ "$end" == "q" ]; then break; fi;
#    done
done
