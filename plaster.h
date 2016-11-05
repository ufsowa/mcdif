#ifndef SITE_H
#define SITE_H
#include "site.h"
#endif

#ifndef GENERATORS_H
#define GENERATORS_H
#include "generators.h"
#endif

#ifndef LIST_H
#define LIST_H
#include <list>
#endif

//http://stackoverflow.com/questions/20854068/removing-object-pointer-from-stl-list
//http://www.cplusplus.com/reference/list/list/remove_if/

class plaster {
	
	private:
	vector <site*> PL_REF_TO_SITES;			//przechowuje wszystkie sity w plastrze
	vector < list < site* > > PL_SITES_TYP;
	
	string PL_NAME;
	unsigned int PL_TYPES,PL_DIRECTION,PL_INDEX;
	double PL_P0,PL_P1;		//P - pozycja poczatkowa i koncowa plastra |....| size bin
	vector <long> PL_EQ_FLUX;
	vector <long> PL_NET_FLUX;
	vector <long> PL_PROB_ADD;
	vector <long> PL_PROB_DEL;
	vector<double> PL_ATOMS;
	vector<double> PL_AVG_PARS;	//p0,p1,A,B,V,flux_A,...,eq_fluxA,..., counter;
	unsigned long PL_JUMPS;
	unsigned long PL_JUMPS_EQ;
	
	struct is_equal{
		is_equal(site* to_find) : to_find(to_find) {}
		bool operator()(site* const in_list){
//			control_output<<to_find<<" ?? "<<in_list<<" "; in_list->show_site();
			if(in_list == to_find){
//				control_output<<"found and removed: "<<in_list<<" ";in_list->show_site();
				unsigned int typ = to_find->get_atom();
				unsigned int atyp = in_list->get_atom();
				if(atyp != typ){
					control_output<<"ERROR: plaster::is_equal():46: "<<typ<<" "<<atyp<<endl;exit(1);
				}
				return true;
			}else{ 
//				control_output<<"Site: "<<to_find<<" not found"<<endl;
				return false;}
		}
		site* to_find;
	};

	public:
	
	
//	plaster();
//	plaster(const plaster &obj);
//	~plaster(){control_output<<"plaster destructed"<<endl;};
	

	plaster(int i, int atoms_type, int direction, int id, double x0 , double x1 , string name );	

	void cumulate();
	void call_avg(vector<double>& results);
	void set_atoms_list(vector <site *> &kontener, int typ);
	void set_atoms_list(list <site *> &kontener, int typ);
	void init_calc(int FLAG=0);
	bool check(unsigned int typ_min, unsigned int typ_max, int delta);
	site* choose_atom(unsigned int typ);
	int choose_typ( const vector <int> &exclude= vector<int>());
	int	get_direction(){return PL_DIRECTION;};
	double	get_st(){return PL_P0;};
	double	get_end(){return PL_P1;};
	int	get_index(){return PL_INDEX;};
	void reset_indexes();
	string get_name(){return PL_NAME;};
	void update_plaster(site* node, bool status);
	void update_hist(site* node, bool status);
	void copy_fluxes(plaster &source);
	void push_back( site* item ){PL_REF_TO_SITES.push_back(item);};
	unsigned int size(){return PL_REF_TO_SITES.size();};
	unsigned int size(int typ){return PL_SITES_TYP[typ].size();};
	unsigned int get_size_types(){return PL_SITES_TYP.size();};
	site* get_site(long pozition){return PL_REF_TO_SITES[pozition];};
	void jump_occured(){PL_JUMPS++;};
	void jump_occured_dislocation(){PL_JUMPS_EQ++;};
	void show();
	bool check_types();
	
	site* get_site(int typ,int nr){
		list<site*>::iterator it= PL_SITES_TYP[typ].begin();
		if( !(PL_SITES_TYP[typ].empty()) ){
			advance(it, nr);
		}else{
			control_output<<"ERROR: plaster::get_site(). Type list empty: "<<typ<<endl;
			exit(1);
		}
		int atyp = (*it)->get_atom();
		if(atyp != typ){
			control_output<<"ERROR: plaster::get_site():97: "<<typ<<" "<<atyp<<endl;
			exit(1);
		}
		return *it;
	};
	
	void	delete_site(unsigned int typ, unsigned long numer){
		unsigned int st_size=PL_SITES_TYP[typ].size();
		list<site*>::iterator item = PL_SITES_TYP[typ].begin();
		if( !(PL_SITES_TYP[typ].empty()) ){
			advance(item, numer);
			unsigned int atyp = (*item)->get_atom();
			if(atyp != typ){
				control_output<<"ERROR: plaster::delete_site():109: "<<typ<<" "<<atyp<<endl;show();exit(1);
				
			}
			PL_SITES_TYP[typ].erase(item);
			eq_flux_delta(typ,0);
		}else{
			control_output<<"ERROR: plaster::delete_site(). Type list empty: "<<typ<<endl;
			show();
			exit(1);
		}
		if( (st_size - PL_SITES_TYP[typ].size()) != 1){
			control_output<<"ERROR: plaster::delete_site(). site not removed "<<typ<<endl;
			control_output<<st_size<<" "<<PL_SITES_TYP[typ].size()<<" "<<(st_size - PL_SITES_TYP[typ].size())<<" "<<abs(st_size - PL_SITES_TYP[typ].size())<<typ<<endl;
			show();
			exit(1);		
		}

	};
	
	void	add_site(int typ, site* new_site){
		unsigned int st_size=PL_SITES_TYP[typ].size();
		int atyp=new_site->get_atom();
		if(typ != atyp){
			control_output<<"ERROR: plaster::add_site(). wrong site type  "<<typ<<" "<<atyp<<endl;
			new_site->show_site();
			exit(1);					
		}
		PL_SITES_TYP[typ].push_back(new_site);
		eq_flux_delta(typ,1);
		
		if( ( PL_SITES_TYP[typ].size() - st_size ) != 1){
			control_output<<"ERROR: plaster::add_site(). site not added "<<typ<<endl;
			control_output<<st_size<<" "<<PL_SITES_TYP[typ].size()<<" "<<( PL_SITES_TYP[typ].size() - st_size )<<" "<<abs(st_size - PL_SITES_TYP[typ].size())<<typ<<endl;
			show();
			exit(1);		
		}

	};


	void	plaster_delete_site(site* node){
		unsigned int typ = node->get_atom();
		unsigned int st_size=PL_SITES_TYP[typ].size();
		if( PL_SITES_TYP[typ].empty() ){
			control_output<<"ERROR: plaster::plaster_delete_site(). Type list empty: "<<typ<<endl;
			show();
			exit(1);
		}																//		control_output<<" del site in plaster "<<node<<" ";
																		//		control_output<<typ<<" | "<<size()<<" | "<<size(typ)<<" ";node->show_site();
		PL_SITES_TYP[typ].remove_if(is_equal(node));
		eq_flux_delta(typ,0);
		prob_update(typ,0);
		if( (st_size - PL_SITES_TYP[typ].size()) != 1){
			control_output<<"ERROR: plaster::plaster_delete_site(). site not removed "<<typ<<endl;
			control_output<<st_size<<" "<<PL_SITES_TYP[typ].size()<<" "<<(st_size - PL_SITES_TYP[typ].size())<<" "<<abs(st_size - PL_SITES_TYP[typ].size())<<typ<<endl;
			show();
			exit(1);		
		}
																		//		control_output<<typ<<" | "<<size()<<" | "<<size(typ)<<endl;
	};
	
	void	plaster_add_site(site* node){
//		control_output<<" add site in plaster "<<node<<" ";
		unsigned int typ = node->get_atom();
		unsigned int st_size=PL_SITES_TYP[typ].size();

//		control_output<<typ<<" | "<<size()<<" | "<<size(typ)<<" ";node->show_site();
		PL_SITES_TYP[typ].push_back(node);
		eq_flux_delta(typ,1);
		prob_update(typ,1);
//		control_output<<typ<<" | "<<size()<<" | "<<size(typ)<<endl;
		if( (PL_SITES_TYP[typ].size() - st_size ) != 1){
			control_output<<"ERROR: plaster::plaster_add_site(). site not added ";
			control_output<<st_size<<" "<<PL_SITES_TYP[typ].size()<<" "<<(PL_SITES_TYP[typ].size() - st_size)<<" "<<abs(st_size - PL_SITES_TYP[typ].size())<<typ<<endl;
			show();
			exit(1);		
		}
	};
	
	void eq_flux_delta(unsigned int typ, int flaga){
		//flaga =0 -> site remved; flaga=1 -> site created in plaster
		if(typ>=PL_EQ_FLUX.size()){cout<<"ERROR in plaster::eq_flux: "<<typ<<endl;exit(1);}
		if(flaga==1){
			PL_EQ_FLUX[typ]++;
		}else if(flaga==0){
			PL_EQ_FLUX[typ]--;		
		}
		else{cout<<"ERROR in plaster::eq_flux"<<endl;exit(1);}
	};

	void flux_net_delta(unsigned int typ, int flaga){
		//flaga =0 -> site remved; flaga=1 -> site created in plaster
		if(typ>=PL_NET_FLUX.size()){cout<<"ERROR in plaster::net_flux: "<<typ<<endl;exit(1);}
		if(flaga==1){
			PL_NET_FLUX[typ]++;
		}else if(flaga==0){
			PL_NET_FLUX[typ]--;		
		}
		else{cout<<"ERROR in plaster::net_flux"<<endl;exit(1);}
	};

	void prob_hist_l(unsigned int typ, int flaga){
		//flaga =0 -> site remved; flaga=1 -> site created in plaster
		if(typ>=PL_PROB_ADD.size()){cout<<"ERROR in plaster::flux: "<<typ<<endl;exit(1);}
		if(flaga==1){
			PL_PROB_ADD[typ]++;
		}else if(flaga==0){
			PL_PROB_ADD[typ]--;		
		}
		else{cout<<"ERROR in plaster::flux"<<endl;exit(1);}
	};

	void prob_hist_r(unsigned int typ, int flaga){
		//flaga =0 -> site remved; flaga=1 -> site created in plaster
		if(typ>=PL_PROB_DEL.size()){cout<<"ERROR in plaster::flux: "<<typ<<endl;exit(1);}
		if(flaga==1){
			PL_PROB_DEL[typ]++;
		}else if(flaga==0){
			PL_PROB_DEL[typ]--;		
		}
		else{cout<<"ERROR in plaster::flux"<<endl;exit(1);}
	};


	
	void prob_update(unsigned int typ, int flaga){
		//flaga =0 -> site remved; flaga=1 -> site created in plaster
		if(flaga==1){
			if(typ>=PL_PROB_ADD.size()){cout<<"ERROR in plaster::flux: "<<typ<<endl;exit(1);}
			PL_PROB_ADD[typ]++;
		}else if(flaga==0){
			if(typ>=PL_PROB_DEL.size()){cout<<"ERROR in plaster::flux: "<<typ<<endl;exit(1);}
			PL_PROB_DEL[typ]++;		
		}
		else{cout<<"ERROR in plaster::flux"<<endl;exit(1);}
	};

	long eq_flux_get(unsigned int typ){
		if(typ>=PL_EQ_FLUX.size()){cout<<"ERROR in plaster::eq_flux: "<<typ<<endl;exit(1);}
		return PL_EQ_FLUX[typ];
	};

	long net_flux_get(unsigned int typ){
		if(typ>=PL_NET_FLUX.size()){cout<<"ERROR in plaster::eq_flux: "<<typ<<endl;exit(1);}
		return PL_NET_FLUX[typ];
	};
	
	long flux_net_get(unsigned int typ){
		if(typ>=PL_NET_FLUX.size()){cout<<"ERROR in plaster::net_flux: "<<typ<<endl;exit(1);}
		return PL_NET_FLUX[typ];
	};

	site* &operator[](unsigned int i){
		if( i > PL_REF_TO_SITES.size() ){
			cout << "Index out of bounds" <<endl; 
            // return first element.
            return PL_REF_TO_SITES[0];
		}
		return PL_REF_TO_SITES[i];
	};

	double get_stech(){
		double stech = static_cast<double> (PL_ATOMS[1])/(static_cast<double> (PL_ATOMS[2]) + static_cast<double> (PL_ATOMS[1]));	
		return stech;
	};

	double get_vac(){
		long int size = PL_REF_TO_SITES.size();	
		double vac = static_cast<double> (PL_ATOMS[0])/static_cast<double> (size);
		return vac;
	};


};

//http://stackoverflow.com/questions/8192570/use-of-list-inside-map-c
//http://stackoverflow.com/questions/10450477/iterate-through-a-map-of-lists
//http://stackoverflow.com/questions/10720359/iterating-over-c-mapint-listint
//http://stackoverflow.com/questions/34526710/how-to-insert-value-in-c-map-stdmapstring-liststring-in-c
