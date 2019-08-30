
set xlabel "MC steps"
set ylabel "No. of atoms "
set key bottom

plot 'u-0.395243u-0.392446N.dat' u 1:($5+$6) w p pt 6 ps 2 t 'Ni',\
 'u-0.395243u-0.392446N.dat' u 1:($7+$8) w p pt 6 ps 2 t 'Al',\
 'u-0.395243u-0.392446N.dat' u 1:($3+$4) w p pt 6 ps 2 t 'Vac'
pause 5
reread

