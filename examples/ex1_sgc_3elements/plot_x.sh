#!/bin/bash

module load gnuplot
counter=0

while [ ! -f u*N.dat ] ;
do
      sleep 2
      ((counter++));
      if [ $counter -gt 10 ];
      then
        echo $counter
        echo "No file to plot"
        break;
      fi;
done

filename=`ls u*N.dat`

echo "
set xlabel \"MC steps\"
set ylabel \"x_{Ni} \"
set key bottom

plot '$filename' u 1:((\$5+\$6)/(\$5+\$6+\$7+\$8)) w p pt 6 ps 2 t 'x_{Ni}'
pause 5
reread
" > plot_x.gp

gnuplot plot_x.gp &