#include "potential.h"

potential :: potential(){
	SAVE=false;
	model="Hamiltonian model not assigned!";
	coordination_zones=0;
	atoms_type=0;
	sublattices=0;
	};
potential :: ~potential(){};

using namespace std;

void potential :: init(unsigned int atom_type_size, unsigned int lattices){

	control_output<<"Initialize potentials..."<<endl;
	string text, used_model;
	ifstream energy_file("energy.in");
	if (energy_file)
	{	control_output<<"energy_file reading from energy.in..."<<endl;}
	else{
	control_output<<"nie ma pliku energy.in"<<endl;
	exit(0);
	}

	energy_file>>text;
	energy_file>>used_model;
	if(used_model=="CVM" or used_model=="Ising"){    //set variable model for cmv
		model=used_model;
		control_output<<"Used energy model: "<<model<<endl;
		}else{
		control_output<<"Bad model used in energy.in: "<<used_model<<endl;
		control_output<<model<<endl;
		exit(0);
	}
	//from this point must distinguise for CVM or Ising becasue of different files format
	if(model=="CVM"){
		init_cvm(energy_file, atom_type_size, lattices);
	}else if(model=="Ising"){
		init_ising(energy_file, atom_type_size, lattices);
	}else{
		control_output<<"Bad model used in energy.in: "<<used_model<<endl;
		control_output<<model<<endl;
		exit(0);
	}
}

void potential :: init_cvm(std::ifstream &energy_file, unsigned int size_V, unsigned int lattices){
	string text;
	vector< vector<int> >::iterator iter_ii;
    vector<int>::iterator iter_jj;
    unsigned int cluster_size=0, liczba_stref=0;
    unsigned int functions_size=0;

	double r1=0.0,r2=0.0;
	energy_file>>text>>liczba_stref;
	coordination_zones=liczba_stref;
	if(text=="rmin_rmax:"){
		for(unsigned int i=0; i<coordination_zones; i++){
			energy_file>>r1>>r2;
			rmin.push_back(r1);
			rmax.push_back(r2);
		}
	}else{
		control_output<<"Bad format in energy.in: "<<text<<endl;}

	energy_file>>text>>cluster_size;
	if(text!="Clusters:"){
		control_output<<"Bad format in energy.in. 'Clusters:'/="<<text<<"'"<<endl;
		control_output<<"You have declared more coordination zones than regions in energy.in"<<endl;
		control_output<<rmin.size()<<" "<<coordination_zones<<" "<<rmin.back()<<" "<<rmax.back()<<endl;exit(0);
	}

	atoms_type = size_V;
	sublattices = lattices;

	for (unsigned int ls=0;ls<cluster_size;ls++){
		ECI.push_back(vector <double>(size_V,0.0) );		//[cluster_size][eci1,eci2,eci3]

		for(unsigned int k=0;k<(size_V);k++){
			double pot=0.0;
			energy_file>>pot;
			//cout<<"Oddzialywanie: i j V "<<ls<<" "<<k<<" "<<pot<<endl;int o;cin>>o;
			if(ECI[ls][k]!=0.0){
				control_output<<"WARRNING:You try overwritte energy in energy.in file!"<<endl;
				control_output<<ls<<" "<<k<<" "<<pot<<endl;
			}
			ECI[ls][k]=pot;
		}
	}
	//here need also read in othogonal functions for every cluster and atom type
	//fun(a,b,x) = a*(x+d)^c + b all function has the same shape with one variable
	//so need just read all parameters a and b for evry cluster/type
	//functions[cluster_size][atom type][a,b]
	energy_file>>text>>functions_size;
	if(text!="FunParameters:"){control_output<<"Bad format in energy.in. 'FunParameters:'/="<<text<<"'"<<endl;
		control_output<<"You have more types in energy.in than declared in structure.in"<<endl;exit(1);}
	for (unsigned int ls=0;ls<functions_size;ls++){
		ORT_FUN.push_back(vector <double>(4,0.0) );		//[cluster_size][eci1,eci2,eci3]
		double a=0.0, b=0.0, c=0.0, d=0.0;
		energy_file>>a>>b>>c>>d;
		//cout<<"fun(s):a b c d "<<ls<<": "<<a<<" "<<b<<" "<<c<<" "<<d<<endl;int o;cin>>o;
		if( ORT_FUN[ls][0]!=0.0 or ORT_FUN[ls][1]!=0.0 or ORT_FUN[ls][2]!=0.0 or ORT_FUN[ls][3]!=0.0){
			control_output<<"WARRNING:You try overwritte function parameters in energy.in file!"<<endl;
			control_output<<ls<<" "<<ORT_FUN[ls][0]<<ORT_FUN[ls][1]<<ORT_FUN[ls][2]<<ORT_FUN[ls][3]<<endl;
		}
		ORT_FUN[ls][0]=a;
		ORT_FUN[ls][1]=b;
		ORT_FUN[ls][2]=c;
		ORT_FUN[ls][3]=d;
	}

	energy_file>>text;
	if(text!="END"){control_output<<"Bad format in energy.in. 'END'/="<<text<<"'"<<endl;
		control_output<<"You have more types in energy.in than declared in structure.in"<<endl;exit(1);}

	for(unsigned int jj=0; jj<rmin.size(); jj++){
		control_output<<"Zapisane odziaływania: rmin "<<rmin[jj]<<" rmax "<<rmax[jj]<<endl;
	}

	for(unsigned int kk=0;kk<ECI.size();kk++){
	for(unsigned int jj=0; jj<ECI[kk].size(); jj++){
	        control_output<<"ECI(s)_i "<<kk<<" "<<jj<<" "<< ECI[kk][jj] << endl;
	}}

	for(unsigned int kk=0;kk<ORT_FUN.size();kk++){
	for(unsigned int jj=0; jj<ORT_FUN[kk].size(); jj++){
	        control_output<<"fun(s)_i "<<kk<<" "<<jj<<" "<< ORT_FUN[kk][jj] << endl;
	}}

}

void potential :: init_ising(std::ifstream &energy_file, unsigned int size_V, unsigned int lattices){
	string text;
	vector< vector<int> >::iterator iter_ii;
    vector<int>::iterator iter_jj;
    unsigned int liczba_stref=0;

	energy_file>>text;
	energy_file>>liczba_stref;
	coordination_zones=liczba_stref;
	atoms_type = size_V;
	sublattices = lattices;
	//vector<vector<double> > V_tmp(size_V, vector<double> (size_V));
	//for(int p=0;p<size_V;p++)
	//{
	//	for(int d=0;d<size_V;d++)
	//		V_tmp[p][d]=1000.0;
	//}
	unsigned int ls;
	for (ls=0;ls<liczba_stref;ls++){
		V.push_back(vector <vector <double> >(size_V, vector <double>(size_V,0.0) ) );

		for(unsigned int k=0;k<(size_V*size_V+1);k++){							// +1 ze wzgledu na rmin_rmax
			unsigned int i=0,j=0;
			double pot=0.0,r1=0.0,r2=0.0;
			energy_file>>text;
			if(text=="rmin_rmax"){
				energy_file>>r1>>r2;
				rmin.push_back(r1);
				rmax.push_back(r2);
								//cout<<text<<endl;
								//cout<<"Promien oddzialywania: "<<r1<<" "<<r2<<endl																	//int o;//cin>>o;
			}else{
				energy_file>>i>>j>>pot;
					if( (i >= size_V) or (j >= size_V)){
						control_output<<"You have more interactions in energy.in than declared types in structure.in"<<endl;
						control_output<<i<<" "<<j<<" "<<pot<<endl;exit(0);}
					if(V[ls][i][j]!=0.0){
						control_output<<"WARRNING:You try overwritte energy in energy.in file!"<<endl;
						control_output<<i<<" "<<j<<" "<<pot<<endl;exit(1);}
					if(rmin.size() != (ls+1) ) {
						control_output<<"You have declared more coordination zones than interactions in energy.in"<<endl;
						control_output<<"rmin.size()"<<" "<<"zone no."<<" "<<"rmin.back()"<<" "<<"rmax.back()"<<endl;
						control_output<<rmin.size()<<" "<<ls<<" "<<rmin.back()<<" "<<rmax.back()<<endl;exit(1);}

				V[ls][i][j]=pot;}
		}
	}
	energy_file>>text;
	if(text != "END"){
		control_output<<"You have bad format in energy.in"<<endl;exit(1);}

	show();
	control_output<<"Ising potentials initialized-> "<<atoms_type<<" "<<sublattices;

	for (unsigned int i=0;i<atoms_type;i++){
		save_bar.push_back( vector<vector<wektor>>(sublattices, vector<wektor>(sublattices, wektor() )));
	}
	control_output<<" save_bar set: "<<save_bar.size()<<endl;

	control_output<<save_bar.size()<<" ";
	for(unsigned int ib = 0; ib < save_bar.size(); ib++){
		control_output<<save_bar[ib].size()<<" ";
		for(unsigned int jb=0; jb < save_bar[ib].size(); jb++){
		control_output<<save_bar[ib][jb].size()<<endl;
		for(unsigned int kb=0; kb < save_bar[ib][jb].size(); kb++){
			control_output<<ib<<" "<<jb<<" "<<kb<<" "; (save_bar[ib][jb][kb]).show();
	}}}

}

void potential :: add_barrier(const pairjump &jump){

	if(SAVE){
	double bariera = jump.get_barier();
	double Em = ( jump.get_e1() + jump.get_e2() )/2.0 ;

	site* vac_to_jump = jump.get_vac_to_jump();
	site* atom_to_jump = jump.get_atom_to_jump();

	int ATOM = atom_to_jump->get_atom("add_barier");
	int NR_LATa = atom_to_jump->get_sub_latt();
	int NR_LATv = vac_to_jump->get_sub_latt();

	save_bar[ATOM][NR_LATa][NR_LATv].x = save_bar[ATOM][NR_LATa][NR_LATv].x + Em;
	save_bar[ATOM][NR_LATa][NR_LATv].y = save_bar[ATOM][NR_LATa][NR_LATv].y + bariera;
	save_bar[ATOM][NR_LATa][NR_LATv].z = save_bar[ATOM][NR_LATa][NR_LATv].z + 1;
	}
}

void potential :: save_barriers(double time, double step, string name){

	string file="";
	stringstream total(name);

	int word_count=0 ;
    string word;
    while( total >> word ) ++word_count;

	if(word_count == 1){
		file=name+"bar";
	}else if(word_count == 2){
		stringstream ss(name);
		int log=0;
		while(ss>>word){
			if(log==0){
			file=word+"bar";log++;}
			else if(log>=1){log++;}
			else {control_output<<"ERROR1 in potential:save "<<log<<endl;exit(1);}
			}
	}else if(word_count>2){control_output<<"ERROR2 in potential::save->file_name: "<<word_count<<endl;exit(1);}
	else{control_output<<"ERROR3 in potential::save->file_name: "<<word_count<<endl;exit(1);}
	string name_of_file= file + ".dat";
	ofstream file_out(name_of_file.c_str(),ios :: app);

	file_out<<step<<" "<<time;
	for(unsigned int ib = 0; ib < save_bar.size(); ib++){
		for(unsigned int jb=0; jb < save_bar[ib].size(); jb++){
		for(unsigned int kb=0; kb < save_bar[ib][jb].size(); kb++){
			 file_out<<" "<<save_bar[ib][jb][kb].x<<" "<<save_bar[ib][jb][kb].y<<" "<<save_bar[ib][jb][kb].z;
			 save_bar[ib][jb][kb].x=0.0; save_bar[ib][jb][kb].y=0.0; save_bar[ib][jb][kb].z=0.0;
	}}}
	file_out<<endl;
	file_out.close();
}

void potential :: read_bars(string file_input){

	control_output<<"Initialize "<<file_input;
	if(atoms_type > 0 and coordination_zones > 0){
		for (unsigned int ls=0;ls<coordination_zones;ls++){
			bars.push_back(vector <vector <double> >(atoms_type, vector <double>(atoms_type,1000.0) ) );
		}
		if(bars.size()==0){
			control_output<<"ERROR:potential::read_bar():95->Cooridination zones not defined"<<endl;exit(0);
		}
	}else{
		control_output<<"ERROR:potential::read_bar():93->Atoms type//zone not defined: "<<atoms_type<<"/"<<coordination_zones<<endl;exit(0);
	}
	string text;
	ifstream file_in(file_input,ios :: in);

	if( file_in.good() ){
		std::string napis;												//   std::cout << "\nZawartosc pliku:" << std::endl;
        while( !file_in.eof() ){
            getline( file_in, napis );									//    control_output<<"String Line: "<< napis << std::endl;
            double data; vector <double> line;
			istringstream string_line(napis);							//			int line_count=0;
			while(string_line>>data){									//		cout<<typeid(data).name()<<" "<<data<<endl;
				line.push_back(data);									//		cout<<line.size()<<endl;
			}
			unsigned int zone = line[0], typ1 = line[1], typ2=line[2];
			double bar = line[3];
			if(line.size() != 4){control_output<<"ERROR:potential::read_bar()->wrong format: zone typ1 typ2 bar"<<endl;exit(0);}
			bars[zone][typ1][typ2]=bar;
			bars[zone][typ2][typ1]=bar;
        }
        file_in.close();
	}else{
		control_output<< "Error! Missing "<<file_input<<" file!" << std::endl;exit(1);
	}

	control_output<<"Barriers initialized ..."<<endl;
	for(unsigned int i=0; i<bars.size();i++){
		for(unsigned int j=0; j < (bars[i]).size();j++){
			for(unsigned int k=0; k < (bars[i][j]).size();k++){

				control_output<<"z "<<i<<" t1 "<<j<<" t2 "<<k<<": "<<bars[i][j][k]<<endl;
	}}}

}

double potential::get_barier(site* node1, site* node2){
	unsigned int zone = check_coordination_zone(node1,node2);
	int atom1 = node1->get_atom("get_barrier");
	int atom2 = node2->get_atom("get_barrier");
	if(bars.size() > (unsigned) zone){
		return bars[zone][atom1][atom2];}
	else{
		control_output<<"ERROR: V size is lower then en_neigh zone: "<<bars.size()<<"<="<<zone<<endl;exit(1);
	}
}

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

unsigned int potential::get_zone(double r){

	int zon = -1;
	for(unsigned int i=0; i<rmin.size();i++){
		if( rmin[i] < r and r < rmax[i] ){
			zon = i;
		}
	}

	if(zon<0){
		control_output<<"ERROR:potential::get_zone-> no such zone: "<<zon<<endl;exit(0);
	}

	return zon;
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

double potential :: get_energy(site *atom){
//	control_output<<"get_energy(site)"<<endl;
	double Vsite=0.0;
//	int zone=0;
	int A=atom->get_atom(); //typ atomu
//	atom->show_site();
//	atom->show_neigh(0);
//	int o;
//	cin>>o;
//	control_output<<"get_energy(site) V.size: "<<V.size()<<" "<<A<<endl;

	for(unsigned int j=0;j<V.size();j++){
		vector <site*> neigh;
//		control_output<<"get_energy(site) for zone: "<<j<<endl;
		atom->read_site_neighbours(neigh,0,j);

	//	cout<<"Strefa numer: "<<j<<endl;
	//	for (int t=0;t<neigh.size();t++)
	//	{
	//		neigh[t]->show_site();
	//	}
		//licz energie konfiguracji
		for(unsigned int i=0;i < neigh.size();i++){
			int B = neigh[i]->get_atom();
			//atom->show_site();
			//neigh[i]->show_site();
			//zone=check_coordination_zone(atom, neigh[i]);
			if(B<0 or A < 0){
				control_output<<"get_energy(site) for neigh: "<<i<<" "<<A<<" "<<B<<" "<<j<<endl;
				control_output<<"Initialization of structure went wrong. You have atom with type < 0"<<endl;
				exit(1);
			}
			Vsite=Vsite+V[j][A][B];
			//cout<<"V[j][A][B]: "<<V[j][A][B]<<endl;
		}
	}
	return Vsite;
}

double potential :: get_energy(site *atom, int typ){
//	cout<<"get_energy() for "<<endl;
	double Vsite=0.0;
	//int zone=0;
	int A=typ; //typ atomu
//	atom->show_site();
//	atom->show_neigh(0);
//	int o;
//	cin>>o;
//	control_output<<"Rozmiar V.size: "<<V.size()<<endl;
	for(unsigned int j=0;j<V.size();j++){
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
			int B = neigh[i]->get_atom("get_energy(site typ)");
			//atom->show_site();
			//neigh[i]->show_site();
			if(B<0 or A < 0){
				control_output<<"get_energy(site) for neigh: "<<i<<" "<<A<<" "<<B<<" "<<j<<endl;
				control_output<<"Initialization of structure went wrong. You have atom with type < 0"<<endl;
				exit(1);
			}
			Vsite=Vsite+V[j][A][B];
//		cout<<"V[j][A][B]: "<<V[j][A][B]<<endl;
		}
	}
//	int o;
	//cin>>o;

//	V=get_energy(A,neigh);
	return Vsite;
}

double potential :: get_energy(site *atom1, site* atom2){
//	control_output<<"get_energy(site site)"<<endl;
	double Vinter=0.0;
	int zone=-1;
	int A=atom1->get_atom("energy for siteA");
	int B=atom2->get_atom("energy for siteB");
	if(B<0 or A < 0){
		control_output<<"get_energy(site) for neigh: "<<A<<" "<<B<<endl;
		control_output<<"Initialization of structure went wrong. You have atom with type < 0"<<endl;
		exit(1);
	}

	vector <site*> neigh;
	//int zones = atom1->get_no_zones();
	zone=check_coordination_zone(atom1, atom2);
//	control_output<<"get_energy(site site): "<<A<<" "<<B<<" "<<zone<<endl;
	if(zone<0){Vinter=0.0;}
	else if(V.size() > (unsigned) zone){Vinter=V[zone][A][B];}
	else{
		control_output<<"ERROR: V size is lower then en_neigh zone: "<<V.size()<<"<="<<zone<<endl;
		control_output<<"Probably energy.in is corrupted"<<endl;exit(1);
	}
	//atom1->show_site();
	//atom2->show_site();
	//cout<<V[zone][A][B]<<endl;
//	control_output<<"get_energy(site site) END"<<endl;

	return Vinter;
}

double potential :: get_energy(site *atom1, int typ1, site* atom2, int typ2){
//	cout<<"get_energy(site site)"<<endl;
	double Vinter=0.0;
	int zone=-1;
	int A=typ1;
	int B=typ2;
	if(B<0 or A < 0){
		control_output<<"get_energy(site) for neigh: "<<A<<" "<<B<<endl;
		control_output<<"Initialization of structure went wrong. You have atom with type < 0"<<endl;
		exit(1);
	}

	vector <site*> neigh;
	//int zones = atom1->get_no_zones();
	zone=check_coordination_zone(atom1, atom2);
	if(zone<0){Vinter=0.0;}
	else if(V.size() > (unsigned) zone){Vinter=V[zone][A][B];}
	else{
		control_output<<"ERROR: V size is lower then en_neigh zone: "<<V.size()<<"<="<<zone<<endl;exit(1);
	}

	//atom1->show_site();
	//atom2->show_site();
	//cout<<V[zone][A][B]<<endl;

	return Vinter;
}