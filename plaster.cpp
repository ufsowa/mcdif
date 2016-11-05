#include "plaster.h"

plaster ::	plaster(int i, int atoms_type, int dir, int id, double x0, double x1, string name){

		PL_NAME=name;
		PL_REF_TO_SITES.reserve(i);
		PL_SITES_TYP.reserve(10);
		PL_TYPES=atoms_type;
		PL_INDEX=id;
		PL_DIRECTION=dir;
		PL_P0=x0;
		PL_P1=x1;
		
		for(unsigned int i=0;i<PL_TYPES;i++){
			PL_ATOMS.push_back(0);
			PL_EQ_FLUX.push_back(0);
			PL_NET_FLUX.push_back(0);
			PL_PROB_ADD.push_back(0);
			PL_PROB_DEL.push_back(0);
		}
		PL_JUMPS=0;
		PL_JUMPS_EQ=0;
		PL_AVG_PARS.reserve(20);
		PL_AVG_PARS.push_back(0);
		PL_AVG_PARS.push_back(0);
		PL_AVG_PARS.push_back(0);
		for(unsigned int i=0;i<PL_ATOMS.size();i++){PL_AVG_PARS.push_back(PL_ATOMS[i]);};
		for(unsigned int i=0;i<PL_NET_FLUX.size();i++){PL_AVG_PARS.push_back(PL_NET_FLUX[i]);};
		for(unsigned int i=0;i<PL_EQ_FLUX.size();i++){PL_AVG_PARS.push_back(PL_EQ_FLUX[i]);};
		for(unsigned int i=0;i<PL_PROB_ADD.size();i++){PL_AVG_PARS.push_back(PL_PROB_ADD[i]);};
		for(unsigned int i=0;i<PL_PROB_DEL.size();i++){PL_AVG_PARS.push_back(PL_PROB_DEL[i]);};
		PL_AVG_PARS.push_back(0);
		PL_AVG_PARS.push_back(0);
//		control_output<<"plaster initied"<<endl;

}

/*
plaster ::	plaster(){

		PL_NAME="null";
		PL_REF_TO_SITES.reserve(5000);
		PL_SITES_TYP.reserve(10);
		PL_TYPES=0;
		PL_INDEX=0;
		PL_DIRECTION=0;
		PL_P0=-1;
		PL_P1=-1;
		
		for(unsigned int i=0;i<PL_TYPES;i++){
			PL_ATOMS.push_back(0);
			PL_EQ_FLUX.push_back(0);
			PL_NET_FLUX.push_back(0);
			PL_PROB_ADD.push_back(0);
			PL_PROB_DEL.push_back(0);
		}
		PL_JUMPS=0;
		PL_JUMPS_EQ=0;
		PL_AVG_PARS.reserve(20);
		PL_AVG_PARS.push_back(0);
		PL_AVG_PARS.push_back(0);
		PL_AVG_PARS.push_back(0);
		for(unsigned int i=0;i<PL_ATOMS.size();i++){PL_AVG_PARS.push_back(PL_ATOMS[i]);};
		for(unsigned int i=0;i<PL_NET_FLUX.size();i++){PL_AVG_PARS.push_back(PL_NET_FLUX[i]);};
		for(unsigned int i=0;i<PL_EQ_FLUX.size();i++){PL_AVG_PARS.push_back(PL_EQ_FLUX[i]);};
		for(unsigned int i=0;i<PL_PROB_ADD.size();i++){PL_AVG_PARS.push_back(PL_PROB_ADD[i]);};
		for(unsigned int i=0;i<PL_PROB_DEL.size();i++){PL_AVG_PARS.push_back(PL_PROB_DEL[i]);};
		PL_AVG_PARS.push_back(0);
		PL_AVG_PARS.push_back(0);
		control_output<<"plaster() initied"<<endl;

}

plaster :: plaster(const plaster &obj){

		PL_NAME=obj.PL_NAME;
		PL_REF_TO_SITES=obj.PL_REF_TO_SITES;
		PL_SITES_TYP=obj.PL_SITES_TYP;
		PL_TYPES=obj.PL_TYPES;
		PL_INDEX=obj.PL_INDEX;
		PL_DIRECTION=obj.PL_DIRECTION;
		PL_P0=obj.PL_P0;
		PL_P1=obj.PL_P1;
		
		PL_ATOMS=obj.PL_ATOMS;
		PL_EQ_FLUX=obj.PL_EQ_FLUX;
		PL_NET_FLUX=obj.PL_NET_FLUX;
		PL_PROB_ADD=obj.PL_PROB_ADD;
		PL_PROB_DEL=obj.PL_PROB_DEL;

		PL_JUMPS=obj.PL_JUMPS;
		PL_JUMPS_EQ=obj.PL_JUMPS_EQ;
		PL_AVG_PARS=obj.PL_AVG_PARS;

		control_output<<"plaster copied"<<endl;
}

*/
 
void plaster :: cumulate(){
	PL_AVG_PARS[0]+=(PL_INDEX);
	PL_AVG_PARS[1]+=(PL_P0);
	PL_AVG_PARS[2]+=(PL_P1);
	int j = 3;
	for(unsigned int i=0;i<PL_ATOMS.size();i++,j++){PL_AVG_PARS[j]+=(PL_ATOMS[i]);};
	for(unsigned int i=0;i<PL_NET_FLUX.size();i++,j++){PL_AVG_PARS[j]+=(PL_NET_FLUX[i]);};
	for(unsigned int i=0;i<PL_EQ_FLUX.size();i++,j++){PL_AVG_PARS[j]+=(PL_EQ_FLUX[i]);};
	for(unsigned int i=0;i<PL_PROB_ADD.size();i++,j++){PL_AVG_PARS[j]+=(PL_PROB_ADD[i]);};
	for(unsigned int i=0;i<PL_PROB_DEL.size();i++,j++){PL_AVG_PARS[j]+=(PL_PROB_DEL[i]);};
	PL_AVG_PARS[j]=PL_JUMPS;j++;
	PL_AVG_PARS[j]=PL_JUMPS_EQ;
}

void plaster :: call_avg(vector<double>& results){
	vector<double> tmp;

	for(unsigned int i=0;i<PL_AVG_PARS.size();i++){
		tmp.push_back(PL_AVG_PARS[i]);
		};
	results=tmp;
	
//	for(unsigned int i=0;i<PL_ATOMS.size();i++){PL_ATOMS[i]=0;};
//	for(unsigned int i=0;i<PL_NET_FLUX.size();i++){PL_NET_FLUX[i]=0;};
//	for(unsigned int i=0;i<PL_EQ_FLUX.size();i++){PL_EQ_FLUX[i]=0;};
//	for(unsigned int i=0;i<PL_PROB_ADD.size();i++){PL_PROB_ADD[i]=0;};
//	for(unsigned int i=0;i<PL_PROB_DEL.size();i++){PL_PROB_DEL[i]=0;};
//	for(unsigned int i=0;i< (3 + PL_ATOMS.size());i++){PL_AVG_PARS[i]=0;};

	for(unsigned int i=0;i< (PL_AVG_PARS.size());i++){PL_AVG_PARS[i]=0;};
	
}

bool plaster :: check(unsigned int typ_min, unsigned int typ_max, int delta){
	
	//delta musi byc wieksze od 0
	bool status = false;
	if ( typ_max < PL_SITES_TYP.size()){
		for (unsigned int i=typ_min;i <= typ_max;i++){
			int rozmiar = PL_SITES_TYP[i].size();
			if(rozmiar < delta){
				status = true;
				break;
			}
		}
	}
	else{
		control_output<<"Over types in plaster :: check"<<endl;
		exit(1);
	}
	return status;
}

void plaster :: copy_fluxes(plaster &source){
	this->PL_EQ_FLUX=source.PL_EQ_FLUX;
	this->PL_NET_FLUX=source.PL_NET_FLUX;
}

site* plaster :: choose_atom(unsigned int typ){
	unsigned long N1=(long)(rnd()*(size(typ)));
	site* node = get_site(typ,N1);
	int atyp = node->get_atom();
	if(typ != typ){
		control_output<<"ERROR:plaster::choose_atom:problem with types: "<<typ<<" "<<atyp<<endl;exit(1);
	}
	return node;
}

int plaster :: choose_typ(const vector <int> &exclude){	
	
	unsigned int SUM = 0; unsigned int i=1; 	int TYP = -1;
	typedef vector <pair <unsigned int,int> > mykey;
	mykey target;
	for(; i<PL_SITES_TYP.size(); i++ ){
		if( ! inlist(exclude.begin(),exclude.end(),i) ){
			unsigned int s = size(i);
			target.push_back(make_pair(SUM, i));
			SUM = SUM + s;	
		}
	}
	target.push_back(make_pair( SUM, i ));
	
	mykey::iterator event = target.begin();
	mykey::iterator next_event = target.begin();

	if(SUM==0){
		control_output<<"Print target: "<<target.size()<<endl;
		int counter=0;
		for(event=target.begin(); event != target.end(); ++event){
			control_output<<counter<<" p: ";
			control_output<<(*event).first<<" "<<(*event).second<<endl;
			counter++;
		}
		control_output<<"WARRNING:plaster::choose_typ:no type available"<<endl;
	}else{
		unsigned int R = rnd()*SUM; 
		event = target.begin();
		next_event = target.begin();
		for( ++next_event ; next_event != target.end(); ++event, ++next_event){	
			unsigned int Lvalue = (*event).first;
			unsigned int Rvalue = (*next_event).first;	
			//		control_output<<Lvalue<<" "<<R<<" "<<Rvalue<<endl;
			if( R>=Lvalue and R < Rvalue){
				TYP = (*event).second;
			//			control_output<<"Find event: "<<Lvalue<<" "<<R<<" "<<Rvalue<<" "<<(*event).second<<" "<<TYP<<endl;
				}
		}
	}

	return TYP;			//TYP < 0 means that there is no atoms in the bin
}


void plaster :: set_atoms_list(vector <site *> &kontener, int typ)
{
	site *pointer=0;
	
	kontener.clear();
	
	for(unsigned int i=0;i<PL_REF_TO_SITES.size();i++)
	{
		pointer=PL_REF_TO_SITES[i];
		
		
		int atom = PL_REF_TO_SITES[i]->get_atom();
		
		if(atom==typ)
		{
			kontener.push_back(pointer);
		}
		
	}
	
//	control_output<<"set atom list typ/size: "<<typ<<" / "<<kontener.size()<<endl;
	
	//for(int i=0;i<kontener.size();i++)
	//{
	//	 control_output<<"nr "<<i<<" "<<kontener[i]<<" "<<endl;
	//	kontener[i]->show_site();
	//	kontener[i]->show_neigh(1);
	//}	
	
}	

void plaster :: set_atoms_list(list <site*> &kontener, int typ)
{
	site *pointer=0;
	
	kontener.clear();
	
	for(unsigned int i=0;i<PL_REF_TO_SITES.size();i++)
	{
		pointer=PL_REF_TO_SITES[i];
		
		
		int atom = PL_REF_TO_SITES[i]->get_atom();
		
		if(atom==typ)
		{
			kontener.push_back(pointer);
		}
		
	}

	
//	control_output<<"set atom list typ/size: "<<typ<<" / "<<kontener.size()<<endl;
	
//	list<site*>::iterator iter;
//	int i=0;
//	for(iter=kontener.begin(); iter!=kontener.end();iter++,i++)
//	{
//		 control_output<<"nr "<<i<<" "<< *iter <<" "<<endl;
//(*iter)->show_site();
//(*iter)->show_neigh(1);
//	}	
	
}	

void plaster :: show(){
	
	control_output<<PL_NAME<<" "<<PL_TYPES<<" "<<PL_DIRECTION<<" "<<PL_INDEX<<" "<<PL_P0<<" "<<PL_P1<<" ";
	control_output<<PL_REF_TO_SITES.size()<<" "<<PL_SITES_TYP.size()<<" "<<PL_ATOMS.size()<<endl;
	
	for(unsigned int i=0; i<PL_SITES_TYP.size(); i++){
		control_output<<i<<" "<<PL_SITES_TYP[i].size()<<" ";
	}
	control_output<<endl;
	
	for(unsigned int i=0; i<PL_ATOMS.size(); i++){
		control_output<<i<<" "<<PL_ATOMS[i]<<" ";
	}
	control_output<<endl;
	
	check_types();
//	vector <long> PL_EQ_FLUX;
//	vector <long> PL_NET_FLUX;
//	vector <long> PL_PROB_ADD;
//	vector <long> PL_PROB_DEL;
//	vector<double> PL_AVG_PARS;	//p0,p1,A,B,V,flux_A,...,eq_fluxA,..., counter;
//	unsigned long PL_JUMPS;
//	unsigned long PL_JUMPS_EQ;
	
}

bool plaster :: check_types(){

	for(unsigned int i=0; i<PL_SITES_TYP.size(); i++){
		list<site*>::iterator item = PL_SITES_TYP[i].begin();
		int counter=0;
		for(; item != PL_SITES_TYP[i].end(); ++item,counter++){
			unsigned int typ = (*item)->get_atom();
			if(typ != i){
				control_output<<"ERROR:plaster::check_types: "<<i<<" "<<counter<<" "<<typ<<endl;
				(*item)->show_site();
				//exit(1);
			}
		}
	}
	
	return true;
}

//UWAGA: sprawdzic -> zaokraglanie, dzielenie przez 0 -> dodac warunki
void plaster :: init_calc(int FLAG){
	
		PL_SITES_TYP.clear();
		PL_SITES_TYP.reserve(4);
		
		for (unsigned int i=0; i < PL_TYPES; i++)
		{ 	
			list <site*> tmp;
			tmp.clear();
	//		control_output<<"do set_atom_list in init_calc for type: "<<i<<endl; 
			set_atoms_list(tmp,i);
			PL_SITES_TYP.push_back(tmp);
		}
	//	control_output<<"in init_calc calculate stech,vac: "; 
		for (unsigned int i=0; i < PL_TYPES; i++){
			PL_ATOMS[i] = PL_SITES_TYP[i].size();
		}

		if(FLAG){	

		control_output<<fixed;
		control_output<<PL_INDEX<<" "<<PL_P0<<" "<<PL_P1<<" ";
		control_output<<"V: "<<PL_ATOMS[0]<<" A: "<<PL_ATOMS[1]<<" B: "<<PL_ATOMS[2]<<" ";
		for(unsigned int i=0;i<PL_NET_FLUX.size();i++){
			control_output<<"||"<<PL_NET_FLUX[i];
			}
		for(unsigned int i=0;i<PL_EQ_FLUX.size();i++){
			control_output<<"|"<<PL_EQ_FLUX[i];
			}
		for(unsigned int i=0;i<PL_PROB_ADD.size();i++){
			control_output<<"||"<<PL_PROB_ADD[i];
			}
		for(unsigned int i=0;i<PL_PROB_DEL.size();i++){
			control_output<<"||"<<PL_PROB_DEL[i];
			}
		control_output<<"|"<<endl; 
		}
	}

void plaster :: reset_indexes(){
	for (vector <site*>::iterator item = PL_REF_TO_SITES.begin(); item != PL_REF_TO_SITES.end(); ++item){
		(*item)->reset_index(PL_NAME);
	}
}

void plaster :: update_plaster(site* node, bool status){

	if(status){
//		control_output<<"stat: "<<status;;
		plaster_add_site(node);
	}else{
//		control_output<<"stat: "<<status;
		plaster_delete_site(node);
	}
}
