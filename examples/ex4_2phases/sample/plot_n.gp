
name=system("ls u*N.dat")
plot name u 1:($5+$6) t 'Ni', name u 1:($7+$8) t 'Al', name u 1:($3+$4) t 'Vac'
pause 5
reread