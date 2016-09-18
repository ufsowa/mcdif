#ifndef WEKTOR_H
#define WEKTOR_H
#include "wektor.h"
#endif

#ifndef IOSTREAM
#define IOSTREAM
#include <iostream>
#endif

#ifndef VECTOR
#define VECTOR
#include <vector>
#endif

#ifndef STDLIB
#define STDLIB
#include <stdlib.h>
#endif

#ifndef PAIRJUMP_H
#define PAIRJUMP_H
#include "pairjump.h"
#endif

#ifndef LIST_H
#define LIST_H
#include <list>
#endif

class site
{
private:
double x,y,z; 			//aktualne wspolrzedne atomu
int atom;				//typ atomu
int sub_latt_name;		//nazwa podsieci
int latt_number;		//numer sieci krystalicznej
double dx;
double dy;		//liczba skokow w danym kerunku
double dz;		
std :: vector <long int> nr_jump; 			// licza wszystkich skokow situ
std :: vector <site*> site_at_neigh; 		//sasiedzie dla bcc przy liczeniu skokow
std :: vector <std::vector <site*> > site_en_neigh;		//sasiedzie dla bcc przy liczeniu energii w kolejnych strefach
std :: vector < list <pairjump>::iterator > site_events;
long int Vindex;
int hist_index;
int block_index;
int rez_index;
 
public:
site ();
site (double _x,double _y,double _z, int _atom, int sub_lat_name);
site (double _x,double _y,double _z, int _atom, int sub_lat_name, int _latt_num);
site(site* Site);
site(const site &Site);
~site();
//site events
void add_events_index(list<pairjump>::iterator skok);
void get_events_index(vector < list <pairjump>::iterator > &skoki);
void clear_events_index();
//site parameters
void clear_neighbours(int _typ);
void put_neigh( site *Site, int typ, int zone);
void put_neighbours(std :: vector <site*> &neig_vect, int typ);
void set_atom(int _atom);
void set_x(double _x);
void set_y(double _y);
void set_z(double _z);
void set_drx(double _x);
void set_dry(double _y);
void set_drz(double _z);
void set_jumps(long int _z, int zon=0);
void set_jumps( std :: vector <long int> &jumps);
void set_z_coordination(int z_cor, int typ);
void set_vindex(long int i);
void set_hist_index(int i);
void set_block_index(int i);
void set_rez_index(int i);
//get site parameters
int get_latt_number();
int get_z_coordination(int typ);
double get_x();
double get_y();
double get_z();
double get_position(int dir);
wektor get_position();
int get_atom();
int get_sub_latt();
double get_drx();
double get_dry();
double get_drz();
long int get_jumps(int zon=0);
void get_jumps( std::vector <long int > &jumps);
int get_no_zones();
long get_vindex();
int get_hist_index();
int get_block_index();
int get_rez_index();
void reset_index(string type);
void reset_vindex();
int events_size();

//site operations
void reset_site ();
void refresh_site();
void read_site_neighbours(std :: vector <site*> &vac_neighbour, int typ, int zone = 0);
void show_neigh(int typ);
void show_site();
void show_site_jumps();
//logical site operators 
bool operator == (const site &A);
site& operator =(const site &A);
void change_to(const site &source);

};

