#!/bin/bash


plot(){
    echo -e "
	file='$1'
	equi='$2'
	set hidden3d
	set ticslevel 0
	set zrange [:]
#        splot file every :5 u 2:((\$4+\$5)/2.):(\$6/(\$6+\$7+\$8)) w l
        splot file every :100:20:100:80 u 2:((\$4+\$5)/2.):((\$9 - \$12)/\$2) w l,\
		equi every :100:2:100:37 u 2:((\$4+\$5)/2.):((\$12)/\$2) w l

#        splot file every :100:20:100:80 u 2:((\$4+\$5)/2.):((\$9 - \$12)/\$2) w l,\
#		file every :100:20:100:80 u 2:((\$4+\$5)/2.):((\$9 + \$12)/\$2) w l

#		file every :100:20:100:80 u 2:((\$4+\$5)/2.):(\$9/\$2) w l
#		file every :10:20:100:80 u 2:((\$4+\$5)/2.):(\$10/\$2) w l
#		file every :10:20:100:80 u 2:((\$4+\$5)/2.):(\$11/\$2) w l

	pause -1
    " > to_plot
    gnuplot to_plot
}


hist=`ls *hist.dat`
biny=`ls *bin.dat`

flux=($hist)
equi=($biny)


for (( i=0; i<=${#hist[*]}; i++ )); do
    echo ${flux[i]} ${equi[i]}
    plot ${flux[i]} ${equi[i]}
done;

