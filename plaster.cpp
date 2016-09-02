#include "plaster.h"

plaster ::	plaster(int i, int atoms_type, int dir, int id, double x0, double x1, string name){

		PL_NAME=name;
		PL_REF_TO_SITES.reserve(i);
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
}

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
	return node;
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
