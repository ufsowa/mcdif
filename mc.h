#ifndef MYOMP_H
	#define MYOMP_H
	#ifdef _OPENMP
		#include <omp.h>
	#else
		#define omp_get_thread_num() 0
	#endif
#endif

#ifndef LATTICE_H
#define LATTICE_H
#include "lattice.h"
#endif

#ifndef TASK_H
#define TASK_H
#include "task.h"
#endif

#ifndef OPCJA_H
#define OPCJA_H
#include "opcja.h"
#endif

#include <typeinfo>
#include <utility>



std::string name_of_control_file="control_file.dat";
std::ofstream control_output(name_of_control_file.c_str(),std::ios :: app);


//Schedule
int execute_task(task &comenda, vector <task> &savings, lattice *sample);
int save_results(lattice *sample, vector <task> &savings, string name, double a=0, double b=0);

//Methods
void widom(lattice *sample, long step, long number_of_steps, double T);
void widom_rnd(lattice *sample, long step, long number_of_steps, double T);
void sgcmc(lattice *sample,long number_of_steps,double T,vector <double> &chem);
void exchange_mechanism(lattice *sample,long steps,double T);
void direct_exchange(lattice *sample,long number_of_steps,double T);
double residence_time(lattice *sample,long number_of_steps,double T, double _barr1, double _barr2);
double residence_time_energy(lattice *sample,long number_of_steps,double T, int numer_plik);
double RTA_random_alloy(lattice *sample,long number_of_steps,double T, double barr1, double barr2);
double vac_mechanism(lattice *sample,long number_of_steps,long direct_step,double T, double _barr1, double _barr2);
//Operations
void insert_atoms(int number, int from_typ, int to_typ,lattice *sample);
void exchange_atoms(long atoms,int from,int to,lattice *sample);

//alghoritms
double time_increment(double norma);
void init_events( lattice* sample, list <pairjump> &tablica_skokow , double beta);
void update_site_events(lattice *sample, list <pairjump> &EVENTS,site* sajt);	
void update_events(lattice *sample, list <pairjump> &EVENTS, site* sajt);
void create_events_index(lattice *sample, site* siteA, vector <pairjump> &events);
int try_jump(lattice* sample,double T,site* vac_to_jump, site* atom_to_jump);
void make_jump(lattice* sample, site* vac_to_jump, site* atom_to_jump);

int main(int arg,char *argc[]);

