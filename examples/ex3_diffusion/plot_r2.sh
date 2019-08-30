#!/bin/bash

module load gnuplot
counter=0

while [ ! -f 1dR2.dat ] ;
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

filename=`ls 1dR2.dat`

echo "
set xlabel \"MC steps\"
set ylabel \"R^{2}/N [A^2/N]\"
set key left

timescale=10**-8   #s
plot '$filename' u (\$2*timescale):((\$7+\$8+\$9)/\$10)  w p pt 6 ps 2 t 'R^2_{Ni}',\\
'$filename' u (\$2*timescale):((\$11+\$12+\$13)/\$14)  w p pt 6 ps 2 t 'R^2_{Al}'

pause 5
reread
" > plot_r2.gp

gnuplot plot_r2.gp &