#ifndef MYMATH_H
#define MYMATH_H
#include "mymath.h"
#endif

//http://4programmers.net/C/Artyku%C5%82y/Operacje_na_wektorach
//http://matematyka.pisz.pl/strona/1630.html

using namespace std;
/**********************class site start****************************/
class wektor
{
public:
double x,y,z;


wektor(){
	x=0.0;
	y=0.0;
	z=0.0;
	};
	
wektor(double _x,double _y,double _z){
	x=_x;
	y=_y;
	z=_z;
	};
	
wektor(const wektor &A){
	x=A.x;
	y=A.y;
	z=A.z;
};


~wektor(){
	};
	
void show()
	{
		control_output<<x<<" "<<y<<" "<<z<<endl;
		//cout<<x<<" "<<y<<" "<<z<<endl;
	};
	
double lenght() const {							//WARUNKI BRZEGOWE
	return sqrt( x*x + y*y + z*z );
	};

double IloczynSkalarny( const wektor &B) const
{
	double x1 = B.x * this->x;
	double y1 = B.y * this->y;
	double z1 = B.z * this->z;
	double C = x1+y1+z1;
	return C;
};

wektor IloczynWektorowy(const wektor &Wektor2) const
{
        wektor vWynik;

        vWynik.x = this->y * Wektor2.z - this->z * Wektor2.y;
        vWynik.y = this->z * Wektor2.x - this->x * Wektor2.z;
        vWynik.z = this->x * Wektor2.y - this->y * Wektor2.x;

        return vWynik;
}


double cos_AB(const wektor &B) const
{
	double AL = this->lenght();
	double BL = B.lenght();
	double IS = this->IloczynSkalarny(B);
	double C = IS/(AL*BL);
	return C;
};

double sin_AB(const wektor &B) const
{
	wektor C = this->IloczynWektorowy(B);
	double AL = this->lenght();
	double BL = B.lenght();
	double CL = C.lenght();
	double sinL = CL/(AL*BL);
	return sinL;
};

double kat(const wektor &B) const
{
	double cosW = this->cos_AB(B);
	double L = acos(cosW);		// * 180.0 / PI
	return L;		//in radians
};

double kat_from_sin(const wektor &B) const
{
	double sinW = this->sin_AB(B);
	double L = asin(sinW);		// * 180.0 / PI
	return L;		//in radians
};

wektor norma(){
	double k = this->lenght();            
	if(k==0){
		control_output<<"ERROR in wektor::norma(). Wektor has lenght: "<<k<<endl; exit(1);
	}
	wektor C;
	C.x = this->x / k ;
	C.y = this->y / k ;
	C.z = this->z / k ;
	return C;
};

wektor wersor(){
	wektor C;
	if(x == 0){	C.x = 0.0 ;}else{ 	C.x = this->x / ( fabs(this->x) ) ;}
	if(y == 0){	C.y = 0.0 ;}else{ 	C.y = this->y / ( fabs(this->y) ) ;}
	if(z == 0){	C.z = 0.0 ;}else{ 	C.z = this->z / ( fabs(this->z) ) ;}
	return C;
};


bool operator ==(const wektor &A){
	if(x==A.x && y==A.y && z==A.z)
	return true;
	else
	return false;
};

bool lower(const wektor &B){

	bool x =( (B.x) > (this->x) );
	bool y =( (B.y) > (this->y) );
	bool z =( (B.z) > (this->z) );
	//control_output<<x<<" "<<y<<" "<<z<<" "<<(x && y && z)<<endl;
	return (x && y && z);
};


wektor operator >(const wektor &A){
	wektor wynik;
	wynik.x = int ( fabs(this->x) > fabs(A.x) );
	wynik.y = int ( fabs(this->y) > fabs(A.y) );
	wynik.z = int ( fabs(this->z) > fabs(A.z) );
	return wynik;
};

wektor operator <(const wektor &A){
	wektor wynik;
	wynik.x = int ( fabs(this->x) < fabs(A.x) );
	wynik.y = int ( fabs(this->y) < fabs(A.y) );
	wynik.z = int ( fabs(this->z) < fabs(A.z) );
	return wynik;
};

wektor& operator =(const wektor &A){
	this->x = A.x;
	this->y = A.y;
	this->z = A.z;
	return *this;
};

wektor& operator ()(double a, double b, double c){
	this->x = a;
	this->y = b;
	this->z = c;
	return *this;
};

wektor operator *(const wektor &A){
	wektor C;
	C.x = this->x * A.x;
	C.y = this->y * A.y;
	C.z = this->z * A.z;
	return C;
};

wektor operator *(double k){	//skalowanie
	wektor C;
	C.x = this->x * k ;
	C.y = this->y * k ;
	C.z = this->z * k ;
	return C;
};


wektor operator -(const wektor &A){
	wektor C;
	C.x = this->x - A.x;
	C.y = this->y - A.y;
	C.z = this->z - A.z;
	return C;
};


wektor operator +(const wektor &A){
	wektor C;
	C.x = this->x + A.x;
	C.y = this->y + A.y;
	C.z = this->z + A.z;
	return C;

};


wektor& operator +=(const wektor &A){
	this->x = this->x + A.x;
	this->y = this->y + A.y;
	this->z = this->z + A.z;
	return *this;
};

wektor& operator -=(const wektor &A){
	this->x = this->x - A.x;
	this->y = this->y - A.y;
	this->z = this->z - A.z;
	return *this;
};


double &operator[](int i){
	if(i==1){
		return x;
	}
	else if(i==2){
		return y;
	}
	else if(i==3){
		return z;
	}
	else{
		cout << "Index out of bounds for wektor: "<<i<<endl; 
		show();
		exit(1);
	}

};

};
/*************************class site end*************************/
