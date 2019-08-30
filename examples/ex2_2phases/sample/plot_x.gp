
set xlabel "MC steps"
set ylabel "x_{Ni} "
set key bottom

plot 'u-0.3956u-0.393N.dat' u 1:(($5+$6)/($5+$6+$7+$8)) w p pt 6 ps 2 t 'x_{Ni}'
pause 5
reread

