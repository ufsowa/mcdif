name=system("ls u*N.dat")
plot name u 1:(($5-$6)/($5+$6)) t 'Etha_{Ni}', name u 1:(($8-$7)/($7+$8)) t 'Etha_{Al}'
pause 5
reread