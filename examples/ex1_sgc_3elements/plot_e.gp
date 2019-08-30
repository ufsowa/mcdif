
set xlabel "MC steps"
set ylabel "U [eV]"
set key bottom

plot 'u-0.4153u-0.3735E.dat' u 1:3 w p pt 6 ps 2

pause 5
reread

