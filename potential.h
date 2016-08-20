#ifndef SITE_H
#define SITE_H
#include "site.h"
#endif

class potential {

unsigned int coordination_zones;
std::vector<std::vector<std::vector<double> > > V;
//tablica 3D v[k][i][j] zawiera potencjaly Vij dla kolejnych stref k 
std::vector <double> rmin;
std::vector <double> rmax;
//zawieraja zasieg oddzialywan


public:
potential();
~potential();
unsigned int get_coordination_number();
void get_interaction_zone(double &_rmin,double &_rmax, int i);

int check_coordination_zone(site *A, site *B);

void init(int atom_type_size);
void set_boundary_condAt_to_pot(wektor &boundary);
void set_boundary_condEn_to_pot(wektor &boundary);
void show();
//double getV1(int _at1,int _at2);
//double getV2(int _at1,int _at2);

//double get_energy(int atom, vector <site*> &jumper_neighbour); 
double get_energy(site *atom1, site* atom2); 
double get_energy(site *atom1, int typ1, site* atom2, int typ2);
double get_energy(site *atom);
double get_energy(site *atom, int typ);

};
