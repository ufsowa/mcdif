#include "generators.h"

double start,start2;

void initialize()
{
//extern ofstream history;

start=time(0);
start2=time(0);


long int seed1;
time(&seed1);
timeval seed2;
struct timezone zone1;
gettimeofday(&seed2,&zone1);
srand(seed1+seed2.tv_usec);
srandom (seed1+seed2.tv_usec);
//history<<"seed: "<<seed1+seed2.tv_usec<<" start: "<<start<<" "<<start2<<endl;
for(long j=0;j<10000;j++)
ran01();

}

void initialize_seed()
{
	
	long int seed1;
	time(&seed1);
	timeval seed2;
	struct timezone zone1;
	gettimeofday(&seed2,&zone1);
	srand(seed1+seed2.tv_usec);

}


double rnd(){
	return ((double)rand())/((double)RAND_MAX+1.0);
	}

double ran01()
{
double rnda,IA=65539,IM=2147483648.0,x;
rnda=fmod(IA*start,IM);
start=rnda;
x=rnda/IM;
return x;
}

/*--------------------------------------------------------------*/

long ranZ(long n)
{
int rnda;
double IM=6075,IA=1366,IC=1283,ISEED=1;

ISEED=fmod(((double)start2)*IA+IC,IM);
start2=ISEED;
rnda=(long)(((n)*ISEED)/IM);

return rnda;
}

/*--------------------------------------------------------*/

double rownomierny(double a,double b) 
{	
		
	double x1;
	double y1;
	
	
	x1 = double(random())/double(RAND_MAX);//ładujemy losowa liczbe
	y1 = (double)((b-a)*x1+a);//sprowadzamy do rozkładu rownomiernego od <a,b>
			
	return y1;	//zwracamy jedna z wartosci
};


