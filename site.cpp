#include "site.h"


using namespace std;
site :: ~site()
{
	//cout<<"destruktor sita"<<endl;
	//delete [] site_neigh;
}

site :: site ()
{
	x=0.0;
	y=0.0;
	z=0.0;
	atom=-3;
	dx=0;
	dy=0;
	dz=0;
	Vindex=-1;
	hist_index=-1;
	block_index=-1;
	rez_index=-1;
	sub_latt_name=-1;
	
	site_events.reserve(500);
	site_at_neigh.reserve(500);
	site_en_neigh.reserve(20);
	nr_jump.reserve(5);
	for (int i=0;i<3;i++){nr_jump.push_back(0);}

};

site :: site(double _x,double _y,double _z,int _atom, int sub_lat_name)
{
	x=_x;
	y=_y;
	z=_z;
	atom=_atom;
	sub_latt_name=sub_lat_name;
	latt_number=0;
	dx=0;
	dy=0;
	dz=0;
	Vindex=-1;
	hist_index=-1;
	block_index=-1;
	rez_index=-1;

	site_events.reserve(500);
	site_at_neigh.reserve(500);
	site_en_neigh.reserve(20);
	nr_jump.reserve(5);
	for (int i=0;i<3;i++){nr_jump.push_back(0);}

}

site :: site(double _x,double _y,double _z,int _atom, int sub_lat_name, int latt_num)
{
	x=_x;
	y=_y;
	z=_z;
	atom=_atom;
	sub_latt_name=sub_lat_name;
	latt_number=latt_num;
	dx=0;
	dy=0;
	dz=0;
	Vindex=-1;
	hist_index=-1;
	block_index=-1;
	rez_index=-1;

	site_events.reserve(500);
	site_at_neigh.reserve(500);
	site_en_neigh.reserve(20);
	nr_jump.reserve(5);
	for (int i=0;i<3;i++){nr_jump.push_back(0);}

}  

site :: site(site* Site)	//przesylam adres komorki
{
	dx=0;
	dy=0;
	dz=0;
	site_events.reserve(500);
	site_at_neigh.reserve(500);
	site_en_neigh.reserve(20);
	nr_jump.reserve(5);
	x=Site->get_x();
	y=Site->get_y();
	z=Site->get_z();
	atom=Site->get_atom();
	sub_latt_name=Site->sub_latt_name;
	latt_number=Site->latt_number;
	dx=Site->dx;
	dy=Site->dy;
	dz=Site->dz;
	nr_jump=Site->nr_jump; 	
	Vindex=Site->Vindex;
	hist_index=Site->hist_index;
	block_index=Site->block_index;
	rez_index=Site->rez_index;
}

site :: site(const site &Site)
{
	dx=0;
	dy=0;
	dz=0;
	site_events.reserve(500);
	site_at_neigh.reserve(500);
	site_en_neigh.reserve(20);
	nr_jump.reserve(5);
	x=Site.x;
	y=Site.y;
	z=Site.z;
	atom=Site.atom;
	sub_latt_name=Site.sub_latt_name;
	latt_number=Site.latt_number;
	dx=Site.dx;
	dy=Site.dy;
	dz=Site.dz;
	nr_jump=Site.nr_jump; 	
	Vindex=Site.Vindex;
	hist_index=Site.hist_index;
	block_index=Site.block_index;
	rez_index=Site.rez_index;
}

/*----------------------------------------------*/

void site :: refresh_site ()		
{
	for(unsigned int i=0; i< nr_jump.size();i++){nr_jump[i]=0;}
	dx=0;dy=0; dz=0;
}

void site :: reset_site ()		
{
	for(unsigned int i=0; i< nr_jump.size();i++){nr_jump[i]=0;}
	Vindex=-10;
	dx=0; dy=0; dz=0;
}

void site :: read_site_neighbours(std :: vector <site*> &vac_neighbour, int typ, int zone)
{
//	control_output<<"jestem site1 "<<endl;
	if(typ){
		vac_neighbour=site_at_neigh;
//		control_output<<"jestem site at "<<endl;
	}else{
		int zones = get_no_zones();
		if(zone >= zones){
			cout<<"ERROR: site :: read_site_neighbours -> coordination szhell failed"<<endl;
			cout<<" ZONE in site neigh.. : "<<zones<<" ZONE called in function: "<<zone<<endl;
			exit(0);
		}
		else{
	//	control_output<<"jestem site en 1 zone: "<<zone<<" corr: "<<zones<<endl;
//		this->show_site();
//		this->show_neigh(0);
	//	control_output<<"size neigh: "<<site_en_neigh[zone].size()<<endl;
		
		vac_neighbour=site_en_neigh[zone];
		//vac_neighbour.clear();
		//for(int i=0;i<site_en_neigh[zone].size();i++)
		//{
		//vac_neighbour.push_back(site_en_neigh[zone][i]);
		//}
	//	control_output<<"jestem site en 2 "<<endl;
		}
	}
	
	//cout<<site_en_neigh.size()<<" n";
	
}

//void site :: show_site()
//{
//control_output<<"x/dx = "<<x<<"/";
//for(int i=0; i< dx.size();i++){control_output<<dx[i]<<"/";}
//control_output<<" y/dy = "<<y<<"/";
//for(int i=0; i< dy.size();i++){control_output<<dy[i]<<"/";}
//control_output<<" z/dz = "<<z<<"/";
//for(int i=0; i< dy.size();i++){control_output<<dz[i]<<"/";}
//control_output<<" atom/subl = "<<atom<<"/"<<sub_latt_name<<endl;
//}

void site :: show_site()
{
control_output<<"x/dx = "<<x<<"/";
control_output<<dx<<"/";
control_output<<" y/dy = "<<y<<"/";
control_output<<dy<<"/";
control_output<<" z/dz = "<<z<<"/";
control_output<<dz<<"/";
control_output<<" atom/subl = "<<atom<<"/"<<sub_latt_name<<" ";
control_output<<" V/H/B/R/E: "<<Vindex<<"/"<<hist_index<<"/"<<block_index<<"/"<<rez_index<<"/"<<events_size()<<endl;

}


void site :: show_site_jumps()
{

control_output<<"jumps = ";
for(unsigned int i=0; i< nr_jump.size();i++){control_output<<nr_jump[i]<<"/";}
control_output<<endl;

}

int site :: get_no_zones()
{
	return site_en_neigh.size();
}

void site :: show_neigh(int typ)
{
		
			if(typ)
			{	
				control_output<<"Sasiedzi atomowi: "<<site_at_neigh.size()<<endl;
				for(unsigned int i=0; i<site_at_neigh.size();i++)
				{
				control_output<<i<<" adress: "<<site_at_neigh[i]<<" ";
				site_at_neigh[i]->show_site();
				}
			}
			else
			{
				int zones=get_no_zones();
				control_output<<"Sasiedzi energii: Liczba stref: "<<zones<<endl;
				for(int j=0;j<zones;j++){
					control_output<<"Zone: "<<j<<" liczba sasiadow: "<<site_en_neigh[j].size()<<endl;
				for(unsigned int i=0; i<site_en_neigh[j].size();i++)
				{
				control_output<<i<<" adress "<<site_en_neigh[j][i]<<" ";
				site_en_neigh[j][i]->show_site();
				}
			}
			}
		
}

double site :: cal_stech(int typ1){
	double x=-1.0;
	
	map <int,int> TYPY;

	for(unsigned int i=0;i<site_at_neigh.size();i++){
		int typ =(site_at_neigh[i])->get_atom();
		
		if(TYPY.count(typ)){
			TYPY[typ]++;
		}else{
			pair< map<int,int>::iterator,bool> ret;
			ret = TYPY.insert(pair<int,int>(typ,1));
			if (ret.second==false) {
				control_output << "ERROR: site::cal_stech. Element "<<typ<<" already existed";
				control_output << " with a value of " << ret.first->second << '\n';exit(1);
			}
		}
	}	
	
	int licz = TYPY[typ1];
	int SUM = 0;
	
	for(auto it = TYPY.begin(); it != TYPY.end(); ++it){
		if( (it->first > 0) ){
			SUM += it->second;
		}
	}
	
	if(SUM > 0){
		x = static_cast<double> (licz) / static_cast<double> (SUM);
		return x;
	}else{
		return 0;
	}
}

void site :: change_to(const site &A){
	
	atom=A.atom;
	dx=A.dx;
	dy=A.dy;
	dz=A.dz;
	nr_jump=A.nr_jump; 		
	Vindex=A.Vindex;	

}



void site :: clear_neighbours(int neigh_typ)
{
	if(neigh_typ)
	{
		site_at_neigh.clear();
	}
	else
	{
		int zones = get_no_zones();
		for (int i=0;i<zones;i++)
		{
			site_en_neigh[i].clear();
		}
		site_en_neigh.clear();
	}
}


void site :: add_events_index( list<pairjump>::iterator skok){
	
		site_events.push_back(skok);
}

void site :: get_events_index(vector < list <pairjump>::iterator > &skoki){
		skoki=site_events;
}

void site :: clear_events_index(){	
		site_events.clear();
}

int site :: events_size(){	
	return site_events.size();
}

void site :: put_neighbours( vector <site*> &neig_vect, int neigh_typ)
{
	//cout<<neig_vect.size()<<" ve";
	if(neigh_typ)
	{
		site_at_neigh.clear();
		site_at_neigh=neig_vect;
	}
	else
	{
		//site_en_neigh.clear();				//TUTAJ WARUNEK NA CZYSZCZENIE W ZALEZNOSCI OD ROZMIARU UWAGA
		site_en_neigh.push_back(neig_vect);
		//=neig_vect;
	//	cout<<site_en_neigh.size()<<" ve"<<endl;
	//	this->show_site();	
	//	int o;
	//	cin>>o;
	}
}

void site :: put_neigh(site *Site, int typ, int zone)
{	
		if(typ)			// jesli 1 to sasiedzi atomu 
		{
			site_at_neigh.push_back(Site);
		}
		else 			// jesli 0 to sasiedzi do liczenia energii + sity z poza obszaru symulacji
		{
			site_en_neigh[zone].push_back(Site);
		}
	}

double site :: get_position(int direction){
		double d=0.0;
		if(direction==1){
			d=get_x();
		}else if(direction==2){
			d=get_y();
		}else if(direction==3){
			d=get_z();
		}else{
			cout<<"Wrong direction number in lattice::get_sites parameters| x-1|y-2|z-3"<<endl;
			exit(1);
		}
		return d;
	}

wektor site :: get_position(){
	wektor r(x,y,z);
	return r;
}

/*--------------------------------------------------------------*/

void site :: set_atom(int _atom){atom=_atom;}
void site :: set_x(double _x){x=_x;}
void site :: set_y(double _y){y=_y;}
void site :: set_z(double _z){z=_z;}
void site :: set_drx(double _x){dx=_x;}
void site :: set_dry(double _y){dy=_y;}
void site :: set_drz(double _z){dz=_z;}
void site :: set_jumps(long int _z, int zon){nr_jump[zon]=_z;}
void site :: set_jumps( std :: vector <long int> &input){nr_jump=input;}

void site :: set_vindex(long _i){if( (_i >=0) ){
	Vindex=_i;
	}else{
		control_output<<"ERROR in site::set_Vindex()"<<endl;exit(1);}
	}

void site :: set_hist_index(int _i){if( (_i >=0) ){	
//	control_output<<"Setting hist index "<<_i<<" to "; show_site(); 
	hist_index=_i;
	}else{control_output<<"ERROR in site::set_hist_index()"<<endl;exit(1);}}

void site :: set_block_index(int _i){if( (_i >=0) ){
//	control_output<<"Setting block index "<<_i<<" to "; show_site(); 
	block_index=_i;}else{control_output<<"ERROR in site::set_bloks_index()"<<endl;exit(1);}}

void site :: set_rez_index(int _i){if( (_i >=0) ){
//	control_output<<"Setting rez index "<<_i<<" to "; show_site(); 
	rez_index=_i;}else{control_output<<"ERROR in site::set_rez_index()"<<endl;exit(1);}}
/*--------------------------------------------------------------*/
double site :: get_x(){return x;}
double site :: get_y(){return y;}
double site :: get_z(){return z;}
double site :: get_drx(){return dx;}
double site :: get_dry(){return dy;}
double site :: get_drz(){return dz;}
long int site :: get_jumps(int zon){return nr_jump[zon];}
void site :: get_jumps( std :: vector <long int> &output){output=nr_jump;}
int site :: get_latt_number(){return latt_number;}
unsigned int site :: get_atom(){if(atom<0){ cout<<"ERROR: atom type < 0. Bad lattice introduction."<<endl;exit(0);}return atom;}
int site :: get_sub_latt(){return sub_latt_name;}
long site :: get_vindex(){return Vindex;}
int site :: get_hist_index(){return hist_index;}
int site :: get_block_index(){return block_index;}
int site :: get_rez_index(){return rez_index;}

void site :: reset_vindex(){ Vindex = -1;}

void site :: reset_index(string type){
	if(type=="hist"){
		hist_index=-1;
	}else if(type=="block"){
		block_index=-1;
	}else if(type=="rez"){
		rez_index=-1;
	}else{
		control_output<<"ERROR in site::reset_index(): "<<type<<endl;exit(1);
	}


}

/*--------------------------------------------------------------*/
bool site :: operator ==(const site &A)
{
if(set_prec(x,1)==set_prec(A.x,1) && set_prec(y,1)==set_prec(A.y,1) && set_prec(z,1)==set_prec(A.z,1))
return true;
else
return false;
}

site& site :: operator= (const site &A)
{
	nr_jump.reserve(3);
	x=A.x;
	y=A.y;
	z=A.z;
	atom=A.atom;
	dx=A.dx;
	dy=A.dy;
	dz=A.dz;
	nr_jump=A.nr_jump; 		
	site_at_neigh=A.site_at_neigh; 
	site_en_neigh=A.site_en_neigh; 
	sub_latt_name=A.sub_latt_name;
	Vindex=A.Vindex;	
	hist_index=A.hist_index;
	block_index=A.block_index;

	return *this;
}

