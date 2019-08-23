
name=system("ls *R2.dat")
plot name u 1:($7+$8+$9) t 'R^2_{Ni}', name u 1:($11+$12+$13) t 'R^2_{Al}'
pause 5
reread