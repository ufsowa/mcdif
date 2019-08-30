
name=system("ls 1dR2.dat")
timescale=10**-8   #s
plot name u ($2*timescale):(timescale*($7+$8+$9)/(6.0*$10*$2)) t 'D_{Ni}', name u ($2*timescale):(timescale*($11+$12+$13)/(6.0*$14*$2)) t 'D_{Al}'
pause 5
reread