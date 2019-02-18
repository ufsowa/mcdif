data := $(shell date +%d%m%Y)	
NAME=mcdif_${data}
CC=g++
CFLAGS=-std=c++11 -g -Wall -c
#mcdif_${data} 

${NAME}: generators.o site.o lattice.o mc.o potential.o mymath.o opcja.o plaster.o pairjump.o
	$(CC) -fopenmp mymath.o generators.o site.o pairjump.o potential.o lattice.o mc.o opcja.o plaster.o -o ${NAME}
generators.o: generators.cpp generators.h
	$(CC) $(CFLAGS) generators.cpp
site.o: site.cpp site.h pairjump.h wektor.h
	$(CC) $(CFLAGS) site.cpp
lattice.o: lattice.cpp lattice.h site.h generators.h potential.h plaster.h
	$(CC) -fopenmp $(CFLAGS) lattice.cpp
mc.o: mc.cpp mc.h lattice.h site.h generators.h pairjump.h opcja.h
	$(CC) -fopenmp $(CFLAGS) mc.cpp
potential.o: potential.cpp potential.h
	$(CC) $(CFLAGS) potential.cpp
opcja.o: opcja.cpp opcja.h lattice.h plaster.h
	$(CC) $(CFLAGS) opcja.cpp
plaster.o: plaster.cpp plaster.h site.h 
	$(CC) $(CFLAGS) plaster.cpp
mymath.o: mymath.cpp mymath.h
	$(CC) $(CFLAGS) mymath.cpp
pairjump.o: pairjump.cpp pairjump.h site.h
	$(CC) $(CFLAGS) pairjump.cpp
	
clean:
	rm -rf *o
	rm mcdif_*


#debug mode: https://stackoverflow.com/questions/3718998/fixing-segmentation-faults-in-c
# http://valgrind.org/docs/manual/quick-start.html#quick-start.intro
# http://www.gnu.org/software/gdb/documentation/
# 