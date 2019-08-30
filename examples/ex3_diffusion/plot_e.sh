#!/bin/bash

module load gnuplot
counter=0

while [ ! -f 1E.dat ] ;
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

filename=`ls 1E.dat`

echo "
set xlabel \"MC steps\"
set ylabel \"U [eV]\"
set key bottom

plot '$filename' u 1:3 w p pt 6 ps 2

pause 5
reread
" > plot_e.gp

gnuplot plot_e.gp &