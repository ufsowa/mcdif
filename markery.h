/* KLASA MARKER */

#ifndef MARK_H
#define MARK_H
#include "mark.h"
#endif

#ifndef FSTREAM
#define FSTREAM
#include <fstream>
#endif

#include <sstream>
#include <stdio.h>
#include <stdlib.h>

class markery
{
	private:
	int sizex;
	int sizey;
	int sizez;
	
	public:
	mark ***marker;
	mark cord_mark;
	int rozkladA[10000];
	int get_sizex();
	int get_sizey();
	int get_sizez();
	void show_R2byatom(double *stat);
	void pritn_rozklad(int step);
	double r2(int x,int y,int z);
	void zlicz(double r);
	markery(int _sizex,int _sizey,int _sizez);
	~markery();
};


void markery :: zlicz(double r)
{
	for(int i=0;i<r;i++)
	{
		rozkladA[i]=0;
		}
	
	}

int markery :: get_sizex()
{
int _size=sizex;

return _size;

}

int markery :: get_sizey()
{
int _size=sizey;

return _size;

}

int markery :: get_sizez()
{
int _size=sizez;

return _size;

}

markery :: markery(int _sizex,int _sizey,int _sizez)
	{
		sizex=_sizex;
		sizey=_sizey;
		sizez=_sizez;
		cout<<sizex<<sizey<<sizez<<endl;
		marker = new mark**[(sizex+3)];
	
		for(int i=0;i<(sizex+3);i++)
		{
			marker[i]=new mark*[(sizey+3)];
			for(int j=0;j<(sizey+3);j++)
			marker[i][j]=new mark[(sizez+3)];
		}
	for(int i=0;i<sizex;i++)
		{
			for(int j=0;j<sizey;j++)         
			{
				for(int k=0;k<sizez;k++)
				{
				marker[i][j][k]=mark();
				}
			}  
		}	
	}

markery :: ~markery()
{
	int _sizex=get_sizex();
	int _sizey=get_sizey();
	int _sizez=get_sizez();
	cout<<"destruktor"<<_sizex<<_sizey<<_sizez<<endl;
	for (int i = 0; i < (_sizex+3); i++) 
	{
		for (int j = 0; j < (_sizey+3); j++)
			{
			delete [] marker[i][j];
			//cout<<" i_j:"<<i<<" "<<j;
			}	
		delete [] marker[i];
	}
	delete [] marker;
	cout<<"skasowane"<<endl;
}


double markery :: r2(int x,int y,int z)
{
	int r2x=x,r2y=y,r2z=z;
	return (r2x*r2x + r2y*r2y + r2z*r2z);
}

void markery :: pritn_rozklad(int step)
{
	int Nx=get_sizex();	
	int Ny=get_sizey();
	int Nz=get_sizez();
	long double R2=0,x=0,y=0,z=0;
	string nrstep;
	string namefile1,namefile2,namefile3;
	stringstream s;
	s<<step;
	nrstep=s.str();
	namefile1="R2A"+nrstep+".dat";
	namefile2="R2B"+nrstep+".dat";
	namefile3="R2V"+nrstep+".dat";
	ofstream outputA(namefile1.c_str(), ios :: app);
	ofstream outputB(namefile2.c_str(), ios :: app);
	ofstream outputV(namefile3.c_str(), ios :: app);
	s.clear();
	s.str("");
	
	
	////////////////////////////////////////////////    

    for(int i=0;i<Nx;i++)
	{
	for(int j=0;j<Ny;j++)
		{
		for(int k=0;k<Nz;k++)
			{
				if(marker[i][j][k].atom==1)
					{	
						x=marker[i][j][k].x;
						y=marker[i][j][k].y;
						z=marker[i][j][k].z;
						R2=r2(x,y,z);
						outputA<<R2<<endl;
					}
				else if(marker[i][j][k].atom==2)
					{	
						x=marker[i][j][k].x;
						y=marker[i][j][k].y;
						z=marker[i][j][k].z;
						R2=r2(x,y,z);
						outputB<<R2<<endl;
					}
				else if(marker[i][j][k].atom==0)
					{	
						x=marker[i][j][k].x;
						y=marker[i][j][k].y;
						z=marker[i][j][k].z;
						R2=r2(x,y,z);
						outputV<<R2<<endl;
					}
					     
			}
		}
    }
	

	}

void markery :: show_R2byatom(double *stat)
{
	
int Nx=get_sizex(),x=0,y=0,z=0;
int Ny=get_sizey();
int Nz=get_sizez();
int nr_markA=0,nr_markB=0,nr_markV=0;
long double R2=0,R2x=0,R2y=0,R2z=0;
long double sumR2A=0,sumR2Ax=0,sumR2Ay=0,sumR2Az=0;
long double sumR2B=0,sumR2Bx=0,sumR2By=0,sumR2Bz=0;
long double sumR2V=0,sumR2Vx=0,sumR2Vy=0,sumR2Vz=0;
long double mR2A=0,mR2Ax=0,mR2Ay=0,mR2Az=0;
long double mR2B=0,mR2Bx=0,mR2By=0,mR2Bz=0;
long double mR2V=0,mR2Vx=0,mR2Vy=0,mR2Vz=0;
////////////////////////////////////////////////////    

    for(int i=0;i<Nx;i++)
	{
	for(int j=0;j<Ny;j++)
		{
		for(int k=0;k<Nz;k++)
			{
				if(marker[i][j][k].atom==1)
					{	
						x=marker[i][j][k].x;
						y=marker[i][j][k].y;
						z=marker[i][j][k].z;
						R2=r2(x,y,z);
						R2x=x*x;
						R2y=y*y;
						R2z=z*z;
						sumR2A=sumR2A+R2;						
						sumR2Ax=sumR2Ax+R2x;
						sumR2Ay=sumR2Ay+R2y;
						sumR2Az=sumR2Az+R2z;
					
						nr_markA++;
						//cout<<i<<j<<k<<" "<<x<<" "<<y<<" "<<z<<" "<<marker[i][j][k].atom<<" "<<R2<<endl;
						
					}
				else if(marker[i][j][k].atom==2)
					{	
						x=marker[i][j][k].x;
						y=marker[i][j][k].y;
						z=marker[i][j][k].z;
						R2=r2(x,y,z);
						R2x=x*x;
						R2y=y*y;
						R2z=z*z;
						sumR2B=sumR2B+R2;						
						sumR2Bx=sumR2Bx+R2x;
						sumR2By=sumR2By+R2y;
						sumR2Bz=sumR2Bz+R2z;
						//cout<<i<<j<<k<<" "<<x<<" "<<y<<" "<<z<<" "<<marker[i][j][k].atom<<" "<<R2<<endl;
						
						nr_markB++;
					}
				else if(marker[i][j][k].atom==0)
					{	
						x=marker[i][j][k].x;
						y=marker[i][j][k].y;
						z=marker[i][j][k].z;
						R2=r2(x,y,z);
						R2x=x*x;
						R2y=y*y;
						R2z=z*z;
						sumR2Vx=sumR2Vx+R2x;
						sumR2Vy=sumR2Vy+R2y;
						sumR2Vz=sumR2Vz+R2z;
						//cout<<i<<j<<k<<" "<<x<<" "<<y<<" "<<z<<" "<<marker[i][j][k].atom<<" "<<R2<<endl;
						sumR2V=sumR2V+R2;
						nr_markV++;
					}
					     
			}
		}
    }
	
	if(nr_markA>0)
	{
	mR2A=sumR2A/nr_markA;
	mR2Ax=sumR2Ax/nr_markA;
	mR2Ay=sumR2Ay/nr_markA;
	mR2Az=sumR2Az/nr_markA;
	}
	else
	{
	mR2A=0;
	mR2Ax=0;
	mR2Ay=0;
	mR2Az=0;
	}
	
	if(nr_markB>0)
	{
	mR2B=sumR2B/nr_markB;
	mR2Bx=sumR2Bx/nr_markB;
	mR2By=sumR2By/nr_markB;
	mR2Bz=sumR2Bz/nr_markB;
	}
	else
	{
	mR2B=0;
	mR2Bx=0;
	mR2By=0;
	mR2Bz=0;
	}
	
	if(nr_markV>0)
	{
	mR2V=sumR2V/nr_markV;
	mR2Vx=sumR2Vx/nr_markV;
	mR2Vy=sumR2Vy/nr_markV;
	mR2Vz=sumR2Vz/nr_markV;
	}
	else
	{
	mR2V=0;
	mR2Vx=0;
	mR2Vy=0;
	mR2Vz=0;
	}
    //cout<<mR2A<<endl;
	
  stat[0]=mR2A;
  stat[1]=mR2Ax;
  stat[2]=mR2Ay;
  stat[3]=mR2Az;
  stat[4]=mR2B;
  stat[5]=mR2Bx;
  stat[6]=mR2By;
  stat[7]=mR2Bz;
  stat[8]=mR2V;
  stat[9]=mR2Vx;
  stat[10]=mR2Vy;
  stat[11]=mR2Vz;
	
}
