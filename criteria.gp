#!/usr/bin/gnuplot -persist

unset key
#set xtics 0.1
#set ytics 0.1
set grid
f(x,y)=1.0 - abs(y-x)
n(x,y,z)= abs(y-z)*0.2

L=1.0
P=0.0

set title sprintf("%f/%f", L,P)
plot [0:1] f(x,L) lc -1, f(x,P) lc 1,  f(x,L)+n(x,L,P) lc -1, f(x,L)-n(x,L,P) lc -1,f(x,P)+n(x,L,P) lc 1, f(x,P)-n(x,L,P) lc 1
