
set xlabel "MC steps"
set ylabel "No. of atoms "
set key bottom

plot '1N.dat' u 1:($5+$6) w p pt 6 ps 2 t 'Ni',\
 '1N.dat' u 1:($7+$8) w p pt 6 ps 2 t 'Al',\
 '1N.dat' u 1:($3+$4) w p pt 6 ps 2 t 'Vac'
pause 5
reread

