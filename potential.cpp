#include "potential.h"

potential :: potential(){};
potential :: ~potential(){};

using namespace std;

void potential :: init(int atom_type_size){
	
	control_output<<"Initialize potentials..."<<endl;
	string text;
	ifstream energy_file("energy.in");
	if (energy_file)
	{	control_output<<"energy_file reading from energy.in..."<<endl;}	
	else{
	control_output<<"nie ma pliku energy.in"<<endl;	
	exit(0);
	}  
	
	vector< vector<int> >::iterator iter_ii;
    vector<int>::iterator iter_jj;
    int liczba_stref=0;

	//tyle ile bylo typow atomow uzytych w strukturze, tyle par energii
	int size_V = atom_type_size;
	//vector Vrow(size_V);
	//V.push_back;
	//cout<<"Ile typow: "<<size_V<<endl;
	//int o;
	//cin>>o;

		
	energy_file>>text;
	energy_file>>liczba_stref;
	coordination_zones=liczba_stref;
	//vector<vector<double> > V_tmp(size_V, vector<double> (size_V));
	//for(int p=0;p<size_V;p++)
	//{	
	//	for(int d=0;d<size_V;d++)
	//		V_tmp[p][d]=1000.0;
	//}
	int ls;
	for (ls=0;ls<liczba_stref;ls++)
	{
		V.push_back(vector <vector <double> >(size_V, vector <double>(size_V,1000.0) ) );
		
		for(int k=0;k<(size_V*size_V+1);k++)  // +1 ze wzgledu na rmin_rmax
		{	
		int i=0,j=0;
		double pot=0.0,r1=0.0,r2=0.0;
		energy_file>>text;
		if(text=="rmin_rmax"){
		energy_file>>r1>>r2;
		rmin.push_back(r1);
		rmax.push_back(r2);
		//cout<<text<<endl;
		//cout<<"Promien oddzialywania: "<<r1<<" "<<r2<<endl;
		//int o;
		//cin>>o;
		}
		else
		{energy_file>>i>>j>>pot;
			if(V[ls][i][j]!=1000.0)
			{
			control_output<<"You try overwritte energy in energy.in file!"<<endl;
		 control_output<<i<<" "<<j<<" "<<pot<<endl;
			exit(0);
			}
			if( (i >= size_V) or (j >= size_V))
			{
			control_output<<"You have more types in energy.in than declared in structure.in"<<endl;
		 control_output<<i<<" "<<j<<" "<<pot<<endl;
			exit(0);
			}
			
		//cout<<text<<endl;
		//cout<<"Oddzialywanie: i j V "<<i<<" "<<j<<" "<<pot<<endl;
		//int o;
		//cin>>o;
			
			V[ls][i][j]=pot;
			}
		}	
	}
	
	show();
	
	control_output<<"Potentials initialized."<<endl;
		
	
	};
	
	
void potential :: show(){
	
	for(unsigned int jj=0; jj<rmin.size(); jj++)
		{ 
			control_output<<"Zapisane odziaływania: rmin "<<rmin[jj]<<" rmax "<<rmax[jj]<<endl;
			
		}	
		
	for(unsigned int kk=0;kk<V.size();kk++){
	for(unsigned int jj=0; jj<V[kk].size(); jj++)
	      {
	         for(unsigned int ii=0; ii<V[kk][jj].size(); ii++)
	         {
	            control_output<<"V^(k)_ij "<<kk<<" "<<ii<<" "<<jj<<" "<< V[kk][ii][jj] << endl;
	         }
	      }}
}	

unsigned int potential :: get_coordination_number(){
	return coordination_zones;
	};	
	
void potential :: get_interaction_zone(double &_rmin,double &_rmax, int i){
	_rmin=rmin[i];
	_rmax=rmax[i];
}

	
double potential :: get_energy(site *atom)	
{
//	cout<<"get_energy() for "<<endl;
	double Vsite=0.0;
//	int zone=0;
	int A=atom->get_atom(); //typ atomu
//	atom->show_site();
//	atom->show_neigh(0);
//	int o;
//	cin>>o;
//	control_output<<"Rozmiar V.size: "<<V.size()<<endl;
	for(unsigned int j=0;j<V.size();j++)
	{
//	control_output<<"jestem pot1 "<<endl;

		vector <site*> neigh;
//		control_output<<"jestem pot2 "<<endl;
		atom->read_site_neighbours(neigh,0,j);//jego sasiedzi energetyczni

//	control_output<<"jestem pot3 "<<endl;

	//	cout<<"Strefa numer: "<<j<<endl;
	//	for (int t=0;t<neigh.size();t++)
	//	{
	//		neigh[t]->show_site();
	//	}
		//licz energie konfiguracji 
		for(unsigned int i=0;i < neigh.size();i++)
		{
			int B = neigh[i]->get_atom();
			//atom->show_site();
			//neigh[i]->show_site();
		//zone=check_coordination_zone(atom, neigh[i]);
			Vsite=Vsite+V[j][A][B];
//		cout<<"V[j][A][B]: "<<V[j][A][B]<<endl;
		}
	}
//	int o;
	//cin>>o;
	
//	V=get_energy(A,neigh);	
	return Vsite;
}
	
double potential :: get_energy(site *atom, int typ)	
{
//	cout<<"get_energy() for "<<endl;
	double Vsite=0.0;
	//int zone=0;
	int A=typ; //typ atomu
//	atom->show_site();
//	atom->show_neigh(0);
//	int o;
//	cin>>o;
//	control_output<<"Rozmiar V.size: "<<V.size()<<endl;
	for(unsigned int j=0;j<V.size();j++)
	{
//	control_output<<"jestem pot1 "<<endl;

		vector <site*> neigh;
//		control_output<<"jestem pot2 "<<endl;
		atom->read_site_neighbours(neigh,0,j);//jego sasiedzi energetyczni

//	control_output<<"jestem pot3 "<<endl;

	//	cout<<"Strefa numer: "<<j<<endl;
	//	for (int t=0;t<neigh.size();t++)
	//	{
	//		neigh[t]->show_site();
	//	}
		//licz energie konfiguracji 
		for(unsigned int i=0;i < neigh.size();i++)
		{
			int B = neigh[i]->get_atom();
			//atom->show_site();
			//neigh[i]->show_site();
		//zone=check_coordination_zone(atom, neigh[i]);
			Vsite=Vsite+V[j][A][B];
//		cout<<"V[j][A][B]: "<<V[j][A][B]<<endl;
		}
	}
//	int o;
	//cin>>o;
	
//	V=get_energy(A,neigh);	
	return Vsite;
}
	


int potential :: check_coordination_zone(site *A, site *B){
	
	static int control_zone=0;
	int zone=-1,zones=0;
	
	zones=A->get_no_zones();
	for(int i=0;i<zones;i++)
	{
		vector <site*> neigh;
		A->read_site_neighbours(neigh,0,i);
		for(unsigned int j=0;j<neigh.size();j++)
		{
			if(B == neigh[j])
			{
				zone=i;
	//			cout<<"Adres atomu z funkcji: "<<B<<" Adres sasiada: "<<neigh[j]<<endl;
				
			}
		}
	}
	
	if(zone<0)
	{
		if(control_zone==0){
		control_output<<"WARNNING: potencial::check_coordination_zone, zone<0, atoms do not interract "<<endl;
		control_zone=1;
		}
	}

//	A->show_site();
//	B->show_site();
//	cout<<"Zone: "<<zone<<endl;
//	A->show_neigh(0);

	//int o;
	//cin>>o;
	return zone;
}

//double potential :: get_energy(int atom, vector <site*> &_neighbour)
//{
//	double V=0.0;
//	
//	for(int i=0;i < _neighbour.size();i++)
//	{
//		V=V+getV1(atom,_neighbour[i]->get_atom());	
//	}
	//cout<<_neighbour.size()<<" N";
	//cout<<V<<" V";
//	return V;
//}

double potential :: get_energy(site *atom1, site* atom2){
//	cout<<"get_energy(site site)"<<endl;
	double Vinter=0.0;
	int zone=-1;
	int A=atom1->get_atom();
	int B=atom2->get_atom();
	vector <site*> neigh;
	//int zones = atom1->get_no_zones();
	zone=check_coordination_zone(atom1, atom2);
	if(zone<0){Vinter=0.0;}
	else{Vinter=V[zone][A][B];}
	//atom1->show_site();	
	//atom2->show_site();	
	//cout<<V[zone][A][B]<<endl;
	
	return Vinter;
}

double potential :: get_energy(site *atom1, int typ1, site* atom2, int typ2){
//	cout<<"get_energy(site site)"<<endl;
	double Vinter=0.0;
	int zone=-1;
	int A=typ1;
	int B=typ2;
	vector <site*> neigh;
	//int zones = atom1->get_no_zones();
	zone=check_coordination_zone(atom1, atom2);
	if(zone<0){Vinter=0.0;}
	else{Vinter=V[zone][A][B];}

	//atom1->show_site();	
	//atom2->show_site();	
	//cout<<V[zone][A][B]<<endl;
	
	return Vinter;
}

//double potential :: getV1(int _at1,int _at2)
//{
	//cout<<V[_at1][_at2]<<" VV";
//	return V[0][_at1][_at2];
//	}
	
//double potential :: getV2(int _at1,int _at2)
//{
	//cout<<V[_at1][_at2]<<" VV";
//	return V[0][_at1][_at2];
//	}
