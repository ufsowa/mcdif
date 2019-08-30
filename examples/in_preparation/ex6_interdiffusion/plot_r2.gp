
name=system("ls 1dR2.dat")
timescale=10**-8   #s
plot name u ($2*timescale):(($7+$8+$9)/$10) t 'R^2_{Ni}', name u ($2*timescale):(($11+$12+$13)/$14) t 'R^2_{Al}'
pause 5
reread