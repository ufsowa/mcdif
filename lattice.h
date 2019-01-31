#ifndef MYOMP_H
	#define MYOMP_H
	#ifdef _OPENMP
		#include <omp.h>
	#else
		#define omp_get_thread_num() 0
	#endif
#endif

#ifndef POTENTIAL_H
#define POTENTIAL_H
#include "potential.h"
#endif

#ifndef GENERATORS_H
#define GENERATORS_H
#include "generators.h"
#endif

#ifndef BOX
#define BOX
#include "box.h"
#endif

#ifndef PLASTER_H
#define PLASTER_H
#include "plaster.h"
#endif

#ifndef SSTREAM
#define SSTREAM
#include <sstream>
#endif

#ifndef SET_H
#define SET_H
#include <set>
#endif

using namespace std;

class lattice
{
private:
site sizex[10]; 	
 		
potential *POTENCIALY;
vector < vector <double> > *BARIERY;
list <pairjump> *EVENTY;

vector <vector<double> > V;
vector <vector<vector<box> > > matrix;		//dziala, ale jest nie uzywana UWAGA
vector <site*> atom_list;
vector <site*> sim_atom_list;
vector <wektor> sublatt_typ;		//przechowuje kombinacje sublattice i atomow na sublattice
vector <vector<site> > cells;

vector <int> atoms_type;			//zawiera typy atomow 0 is reserved for vacancy
vector <string> atoms_name;			//ZAMIENIC NA KEY MAP


unsigned int x_size;			// obszar wykorzystany na zbudowanie siatki w pamieci
unsigned int y_size;
unsigned int z_size;
						//region probki - calosc
wektor st_region;
wektor end_region;
						//obszar symulacji
wektor st_sim_area;
wektor end_sim_area;  
bool TRANSPARENT,MOVE_SIM_REGION;										//true if boundary_plane is transparent for atoms/sim_area is allow to move
						//warunki brzegowe
wektor boundary_con_at;
wektor boundary_con_en;

double latt_const;
vector <double> x_trans;		//rozmiar komorki elementarnej. Wykorzytywany przy liczeniu dyfuzji r2. 
vector <double> y_trans;
vector <double> z_trans;

double Rmin;			//zasieg skokow atomow - definiowany w conf.dat main()
double Rmax;

wektor interaction_zone;		//zawiera maksymalny rozmiar zakresu sasiadow
int max_coordination_number;	//maksymalna liczba sasiadow atomu w probce

public:
unsigned int sublattice;
int local_control_atom;
lattice(int _sizex,int _sizey,int _sizez);
~lattice();

//initialization of lattice
void simulation_set(double _r_min, double _r_max,wektor a,wektor b,wektor c,wektor d, wektor _max_zone, wektor e, wektor f);
void simulation_initialize();
void init_sim_boundary(vector <double> &parameters);
void atoms_list_init();
void sim_atoms_list_init();
void jumps_shell_init();
void interaction_shell_init();
void sites_zone_init(int typ,double rmin, double rmax,site *atom_list, vector <site*> &tmp_atom_list);
void set_atoms_list(vector <site *> &kontener, int typ);
void set_atoms_list(set <site *> &kontener, int typ);
void set_atoms_map(vector <site *> &kontener);
void read_structure(string file_name,wektor start,wektor end,wektor set_st, int lattice_nr);
void set_alg_objects(list <pairjump> &evt, vector < vector <double> > &bar, potential &pot);

//opearation of lattice
void sort_atoms(vector <site> &atoms);
site* search_site(site *A);
void refresh_structure(string file_name);
bool reinit_sim_area(wektor a, wektor b, set<site*> &vatoms);
void update_vac_list( set<site*> &ADD,  set<site*> &OLD);
double calc_energy();
double calc_energy_global();

//Saving results
string get_file_name(string name,string format);
void save_energy(double Time, double Step, string name_of_file, int setON=0);
void save_Natoms(double Time, double step, string name_of_file, int setON=0);
void save_NandE(double Time, double Step, string name, int setON=0);
void save_SRO(double Time, double Step, string name);
void save_SRO_deep(double Time, double Step, string name);
void clear_dR();
void save_dR(double Time, long step, string name_of_file, int setON=0);
void save_hist_dR(string file_name,int typ, double Time, double st_bin, double size_bin, double end_bin);
void makepic(long step,long step_break,wektor st, wektor end, string name_of_file);
void pic_stech(long step,double stech, wektor make_pic_vec_st, wektor make_pic_vec_ed, string name_of_file);
void pic_diff(long step,long step_break,wektor st, wektor end, string name_of_file);

//initialization of sites
void add_atom_type(int _atom, string _name);
void add_sublatt_typ(int sublatt, int atomtyp);
void put_site(int x,int y,int z,site* Site);
void put_atom(int x,int y,int z,site* Site);

//operation on sites
void exchange_sites( site* A, site* B);
void get_sites(vector <plaster> &tmp);
void get_sites(plaster &tmp);
site* get_site(long pozition);
int get_atom_name(int typ);
int get_atom_spin(int typ);
int get_atom_type(string name);
int get_atom_type(int spin);
double move(double x2, double x1, int dir);

//reading lattice parameters
int get_size(int typ);
int get_vec_lattice_typ_size();
unsigned int get_atom_typ_numbers();
long get_atoms_number();
int get_sub_lattice_type(double X, double Y, double Z, int lattice_nr);
long get_sim_atom_number();
double get_latice_const(int direction, int i);
double get_latice_transition(int direction);
wektor get_st_sim_wektor();
wektor get_end_sim_wektor();
wektor get_PB();
int get_max_coordination_number();
potential& get_potentials();
double get_stech(int typ1,int typ2);
void get_window(site* A, site* V, vector <site> &tab);
void get_sity_from_nnbox(int x,int y, int z,int latt_num,vector <site*> &tmp_atom_list);

//controls
bool check_site_belonging_to_region(site *A);
bool check_site_belonging_to_sim_area(site *A);
bool check_site_mobile(site* node);
bool check_boundary_conditions(wektor *wsk);
void check_neighbours(int typ);
void check_atoms();

//alghoritms
void init_events_list(set <site *> &kontener);
void update_events(site* sajt);

private:
void update_site_events(site* sajt);	
void clear_events_index(site* sajt);
void create_events_index(site* siteA, vector <pairjump> &events);
void create_events_trans(site* siteA, vector <pairjump> &events);

};

