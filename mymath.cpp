
#include "mymath.h"

int myRound( double fValue )
{
 //   return fValue < 0 ? ceil( fValue - 0.5 )
  //      : floor( fValue + 0.5 );
	const double sd = 10; //for accuracy to 3 decimal places
	return int(fValue*sd + (fValue<0? -0.5*sd : 0.5*sd))/sd;     
	//return int(fValue*sd + (fValue<0? -0.5 : 0.5))/sd;     
        
}


double set_prec(double x, double prec)
{
	//zaokraglanie do 4 miejsca po przecinku
	
	int scale = pow(10,(prec+2));
	
  long int y = x * scale; // przesuwamy przecinek o 6 miejsca i ucinamy reszte za przecinkiem - y jest calkowite
 // control_output.precision(7);
//	control_output<<"Zaookraglam: "<<std::fixed<<x<<" "<<std::fixed<<y<<" modulo "<<(y % 10);
	
	
	
	if (x<0)		// zalezy czy y jest ujemne czy dodatnie. Zaokraglam do 5 
	{
	   if (y % 10 <= -5){ y = y - 10;}
	}
	else
	{
		if (y % 10 >= 5){ y = y + 10;}
	} // jezeli cyfra jednosci >= 5
  
  y=y/10;		// i ucinam reszte do 5 miejsca po przecinku
 // control_output<<" ypo1 "<<y;
  
  if (x<0)  //zaokraglam do 4 miejsca
	{
	   if (y % 10 <= -5){ y = y -10;}
	}
	else
	{
		if (y % 10 >= 5){ y =y + 10;}
	}
  // i ucinam reszte za 4
  y=(y / 10);
  double ilo = pow(10,prec);
  double new_value = (y/ilo);
  
//  control_output<<" ypo2 "<<y<<" "<<new_value<<std::endl;
//	int o; std::cin >> o;
  return new_value; // usuwamy ostatnia cyfre i zamieniamy na liczbe zmiennoprzecinkowa
}

double integral_data(std::vector<double> &X, std::vector<double> &Y, double a){
	double SUM=0;
	double SUM_F=0;
	int sign=1;
	for(int i=1;i<X.size();i++){
		double diff = (X[i]) - (X[i-1]);
		if( (X[i]-a) < 0  and (X[i-1]-a) < 0 ){diff = diff * (-1); sign=-1;}
		else if( (X[i]-a) > 0  and (X[i-1]-a) > 0 ){diff = diff * (1);sign=1;}
		else{diff=0;sign=0;}
		
		SUM += CalculateTrapezoidArea(Y[i],Y[i-1], diff);
		SUM_F += sign*Y[i];
//		std::cout<<a<<" "<<X[i-1]<<" "<<X[i]<<" "<<(X[i-1] - a)<<" "<<(X[i]-a)<<" "<<diff<<" "<<Y[i-1]<<" "<<Y[i]<<" "<<SUM<<" "<<SUM_F<<std::endl;
	}
	return SUM;
	
}

float CalculateTrapezoidArea(float sideA, float sideB, float height){
	return ( ( (sideA+sideB) * height  ) / 2.0f );
}

double rad2degree(double kat_in){return kat_in * 180.0 / PI;}
double degree2rad(double kat_in){return kat_in * PI / 180.0;}

