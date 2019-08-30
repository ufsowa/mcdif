
set xlabel "MC steps"
set ylabel "x_{Ni} "
set key bottom

plot 'u-0.4153u-0.3735N.dat' u 1:(($5+$6)/($5+$6+$7+$8)) w p pt 6 ps 2 t 'x_{Ni}'
pause 5
reread

