#ifndef SITE_H
#define SITE_H
#include "site.h"
#endif

#ifndef SSTREAM
#define SSTREAM
#include <sstream>
#endif

class potential {

bool SAVE;
unsigned int coordination_zones, atoms_type, sublattices;
std::vector<std::vector<std::vector<double> > > V;
//tablica 3D v[k][i][j] zawiera potencjaly Vij dla kolejnych stref k 
std::vector <double> rmin;
std::vector <double> rmax;
//zawieraja zasieg oddzialywan
vector <vector <vector<double> > > bars;								//keeps energy of the jumps barrier

vector <vector <vector <wektor > > > save_bar;

public:
potential();
~potential();
unsigned int get_coordination_number();
void get_interaction_zone(double &_rmin,double &_rmax, int i);

int check_coordination_zone(site *A, site *B);
unsigned int get_zone(double r);

void init(unsigned int atom_type_size, unsigned int lattices);
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

void read_bars(string file_input);
double get_barier(site* node1, site* node2);

void set_save(){SAVE = true;};
void unset_save(){SAVE = false;};
bool save_status(){return SAVE;};
void add_barrier(const pairjump &jump);
void save_barriers(double time, double step, string name);

};

//http://stackoverflow.com/questions/19062874/how-to-represent-a-3d-array-with-map
