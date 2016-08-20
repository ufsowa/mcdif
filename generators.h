#include <math.h>

#ifndef TIME
#define TIME
#include <time.h>
#endif

#ifndef IOSTREAM
#define IOSTREAM
#include <iostream>
#endif

#ifndef FSTREAM
#define FSTREAM
#include <fstream>
#endif

#ifndef STDLIB
#define STDLIB
#include <stdlib.h>
#endif

#ifndef SYS
#define SYS
#include <sys/time.h>
#endif

using namespace std;

double ran01();
long ranZ(long n);
double rownomierny(double a,double b);
void initialize();
void initialize_seed();
static double rnd(){return ((double)rand())/((double)RAND_MAX+1.0);}
