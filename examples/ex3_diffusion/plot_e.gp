
set xlabel "MC steps"
set ylabel "U [eV]"
set key bottom

plot '1E.dat' u 1:3 w p pt 6 ps 2

pause 5
reread

