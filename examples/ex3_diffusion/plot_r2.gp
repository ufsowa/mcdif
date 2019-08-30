
set xlabel "MC steps"
set ylabel "R^{2}/N [A^2/N]"
set key left

timescale=10**-8   #s
plot '1dR2.dat' u ($2*timescale):(($7+$8+$9)/$10)  w p pt 6 ps 2 t 'R^2_{Ni}',\
'1dR2.dat' u ($2*timescale):(($11+$12+$13)/$14)  w p pt 6 ps 2 t 'R^2_{Al}'

pause 5
reread

