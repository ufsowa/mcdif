#ifndef IOSTREAM
#define IOSTREAM
#include <iostream>
#endif

using namespace std;
/**********************class site start****************************/
class mark
{
public:
int x,y,z;
int atom;
int x0,y0,z0,nr_jump; 			//wspolrzedne wakancji
long double R2;
double aktualR2;

mark ()
	{
	x=0;
	y=0;
	z=0;
	atom=-3;
	};
	
mark(int _x,int _y,int _z,int _atom)
	{
	x=_x;
	y=_y;
	z=_z;
	atom=_atom;
	};
	
mark(int _x,int _y,int _z)
	{
	x=_x;
	y=_y;
	z=_z;
	};	
	
void show_mark()
{
cout<<"x = "<<x<<" y = "<<y<<" z = "<<z<<" atom = "<<atom<<endl;
}

void write_mark_cord(long war)
{
	ofstream output("cord.xyz", ios :: app);
	
	if(pow(-1,x)>0 && pow(-1,y)>0 && pow(-1,z)>0)
	output<<"Ni "<<x<<" "<<y<<" "<<z<<endl;
	else
	output<<"Al "<<x<<" "<<y<<" "<<z<<endl;
	
};

~mark()
	{
	};

bool operator ==(mark &A)
{
if(x==A.x && y==A.y && z==A.z && atom==A.atom)
return true;
else
return false;
};

};
/*************************class site end*************************/
