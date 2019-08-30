#!/bin/bash

module load gnuplot
counter=0

while [ ! -f 1N.dat ] ;
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

filename=`ls 1N.dat`

echo "
set xlabel \"MC steps\"
set ylabel \"No. of atoms \"
set key bottom

plot '$filename' u 1:(\$5+\$6) w p pt 6 ps 2 t 'Ni',\\
 '$filename' u 1:(\$7+\$8) w p pt 6 ps 2 t 'Al',\\
 '$filename' u 1:(\$3+\$4) w p pt 6 ps 2 t 'Vac'
pause 5
reread
" > plot_n.gp

gnuplot plot_n.gp &