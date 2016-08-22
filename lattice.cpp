#include "lattice.h"

using namespace std;

//definicja obszaru probki do symulacji



lattice :: lattice(int _xsize,int _ysize,int _zsize )
{ 
	//rezerwuje pamiec na symulacje. Buduje siec boxow.
	x_size=_xsize;
	y_size=_ysize;
	z_size=_zsize;
	sublattice=0;
	//cout<<"Creating lattice in adress: "<<&matrix<<endl;
	long iter=0;
	matrix.clear();
	matrix.reserve(50000);
	for(int i=0;i < int(x_size);i++)
    {
		//cout<<i<<" ";
		matrix.push_back( vector<vector<box> >(y_size,vector <box> (z_size, box(iter))));
	//	matrix.push_back( vector<vector<box> > ());
		
	}

/*	for(int i=0;i < int(x_size);i++)
	{
		matrix[i].clear();
		matrix[i].reserve(100);
		for(int j=0;j < int(y_size);j++)
		{
			matrix[i].push_back(vector<box> ());
		}
	}
	
	for(int i=0;i < int(x_size);i++)
	{
	for(int j=0;j < int(y_size);j++)
	{
		matrix[i][j].clear();
		matrix[i][j].reserve(100);
		for(int k=0;k < int(z_size);k++)
		{box item(iter);
		cout<<i<<" "<<j<<" "<<k<<" ";
		matrix[i][j].push_back(item);
		cout<<&matrix[i][j][k]<<endl;
		}
	}
	}

*/	
	for(int i=0;i < int(x_size);i++)
    {
		for(int j=0;j < int(y_size);j++)
    {
			for(int k=0;k < int(z_size);k++)
    {
		iter++;
	matrix[i][j][k].set_boxid(iter);
//	cout<<i<<" "<<j<<" "<<k<<" ";
//	cout<<&matrix[i][j][k]<<endl;

	}}}
	
	
		
//	int o;
	string text;
	//wczytuje atomy bazowe struktury
	ifstream unit_cell("structure.in");
	if (unit_cell)
	{	control_output<<"unit_cell reading from structure.in..."<<endl;
	control_output<<"All conditions are taken with '/<' not '<=' "<<endl;}	
	else{
	control_output<<"nie ma pliku structure.in"<<endl;	
	exit(0);
	}
	
	
	int nr_struct;
	int nr_atoms;
	//rozmiary komorki elementarnej i wektor translacji tej komorki
	double x_left_border;
	double x_right_border;
	double x_translation;
	double y_left_border;
	double y_right_border;
	double y_translation;
	double z_left_border;
	double z_right_border;
	double z_translation;
	
	cells.reserve(5);	
	cells.clear();
	atoms_type.reserve(8);
	atoms_type.clear();
	//liczba struktur do wczytania
	unit_cell>>text>>nr_struct;
	
	//iteruje dla kazdej struktury
	for(int l=0;l<nr_struct;l++)
	{
		unit_cell>>text>>x_left_border>>y_left_border>>z_left_border;
		unit_cell>>x_right_border>>y_right_border>>z_right_border;
		unit_cell>>x_translation>>y_translation>>z_translation;
	
		control_output<<"structure size: Lbord/Rbord/trans"<<endl;
		control_output<<x_left_border<<" "<<x_right_border<<" "<<x_translation<<endl;
		control_output<<y_left_border<<" "<<y_right_border<<" "<<y_translation<<endl;
		control_output<<z_left_border<<" "<<z_right_border<<" "<<z_translation<<endl;
		
		//zapisanie rozmiaru komorki elementarnej
		x_trans.push_back(x_translation);		
		y_trans.push_back(y_translation);
		z_trans.push_back(z_translation);
		
		
		
		
		// tablica vecorow do przechowywania sitow. Pelni role sitow/boxow w ktorych siedza atomy o double cord: 1<=1.44 <2
		//wypelnione -3. atom <0 means that site is not used in simulation
		
		
		
		//cin>>o;
		//liczba atom w tej komorce do wczytania
		unit_cell>>text>>nr_atoms;
		control_output<<"atoms in unit cell: "<<nr_atoms<<endl;
		//cin>>o;
		
		
	//	cout<<"LATTICE"<<endl;
	//	matrix[2][2][2].show_site();
		
		vector <site> atoms_incell;
		atoms_incell.reserve(10);	
		atoms_incell.clear();
		unit_cell>>text;
		
		//wczytuje uzyte atomy komorki do tablicy 
		for (int i=0;i<nr_atoms;i++)
			{
				double x=0.0,y=0.0,z=0.0;
				int atom=1;
				string name=""; 
				int sub_lat_name;
				unit_cell>>x>>y>>z>>atom>>name>>sub_lat_name;
				control_output<<x<<" "<<y<<" "<<z<<" "<<atom<<" "<<name<<" "<<sub_lat_name<<" || ";
				control_output<<endl;
				  
				if(x<x_translation and y<y_translation and z<z_translation)
				{
					add_sublatt_typ(sub_lat_name, atom);
					add_atom_type(atom,name);	//dodaje typy uzytych atomow				
					atoms_incell.push_back(site(x,y,z,atom,sub_lat_name,l));
				} 
				else 
			 	{
					control_output<<"Erorr: Atoms in cell behind cell. Site cord must be < trans"<<endl;
					exit(0);
				}
			  
			}
		
		control_output<<"atoms in cell: "<<atoms_incell.size()<<endl;
		
		cells.push_back(atoms_incell);   
		
		//cin>>o;
		for(double x=x_left_border;x<x_right_border;x=x+x_translation)
		{
		for(double y=y_left_border;y<y_right_border;y=y+y_translation)
		{
			for(double z=z_left_border;z<z_right_border;z=z+z_translation)
			{
				vector <site> :: iterator K;
		//wpiasnie atom komorki elementarnej do vectora atomow
				for(K=atoms_incell.begin();K!=atoms_incell.end();K++)
				{
					double X=K->get_x();
					double Y=K->get_y();
					double Z=K->get_z();  
					int ATOM=K->get_atom();
					int name=K->get_sub_latt();
					  
					//control_output<<" "<<name;
							
					if(set_prec(x+X)<set_prec(x_right_border) and set_prec(y+Y)<set_prec(y_right_border) and set_prec(z+Z)<set_prec(z_right_border))
					{
					//	control_output<<set_prec(x+X)<<" "<<set_prec(y+Y)<<" "<<set_prec(z+Z)<<" "<<ATOM<<" / ";
					//	control_output<<x_trans[l]<<" "<<y_trans[l]<<" "<<z_trans[l]<<" || ";
						site new_site(set_prec(x+X),set_prec(y+Y),set_prec(z+Z),ATOM,name,l);
					//	cout<<"Adres new_site w funkcji: "<<&new_site<<endl;
//						int a = myRound((x+X)/x_trans[l]);
	//					int b = myRound((y+Y)/y_trans[l]);
		//				int c = myRound((z+Z)/z_trans[l]);
						int a = int((x+X)/x_trans[l]);
						int b = int((y+Y)/y_trans[l]);
						int c = int((z+Z)/z_trans[l]);
		//				control_output<<a<<" "<<b<<" "<<c<<endl;  
						matrix[a][b][c].put_site(new_site);
					//	cout<<"Jeszcze raz sprawdzam sity w boxie z lattice"<<endl;
			//			matrix[a][b][c].show_sity();  
					//	int o;
					//	cin>>o; 
					  
						
					}
					//site *_new_site=&new_site;
					//site *_old_site=&new_site;
					
					//cout<<"adres sita przed: "<<_old_site<<endl;
					//_old_site=search_site(_new_site);
					//cout<<"adres sita po: "<<_old_site<<endl;
					
					//search_site(new_site);  skanuje aktualna lista sitow i porownuje new_site z sitami
					//nadpisac zawartosc komorki-zeby adresu nie zmieniac
					
					
							 
					 //wpisuje atom do boxa		
					
					
					
					/*
					if(_old_site)
					{	cout<<endl;
						cout<<"You are traing overwrite existing site! \n To continue press 1 To end program press 0"<<endl;
						cin>>o;
						if(o==0)
						{exit(0);}
						//_old_site->show_site();
						_old_site->set_atom(ATOM);
						//cout<<"Nadpisano sita"<<endl;
						//_old_site->show_site();
						
						}
					else
						{
					//		cout<<"Dodaj sita"<<endl;
							atom_list.push_back(new_site); 
							
																//klasa vector tworzy obiekt typu vector dla typow X
					//		cout<<"Sit dodano"<<endl;			//biorac adres danego obiektu typu X 
																//i tworzy nowy obiekt typu X korzystajac z adresu ob. typu X.
						}										//i ten nowy utworzony obiekt dodaje do vectora
					*/											//innymi slowy, korzysta z konstruktora site(cons site &obiekt)
					
				}

				
		//		control_output<<atom_list.size()<<endl;
				
			}
		}
		}
		
		
		
		control_output<<"structure " <<l<<" readed ok"<<endl;
	
	//	for(int i=0;i<int(x_size);i++)
	//	{
	//		for(int j=0;j<int(y_size);j++)
	//		{
	//			for(int k=0;k<int(z_size);k++)
	//			{
	//				control_output<<"i j k "<<i<<" "<<j<<" "<<k<<endl;
	//				matrix[i][j][k].show_sity();
	//			}
	//		}
	//	}
		
//cin>>o;
	
	
	
	
	}
	
	//wczytuje energie do tablicy potencjalow V znajduje sie w klasie potencial
	
		
}

long lattice :: get_atoms_number()
{
	return atom_list.size();
	
	}

void lattice :: set_atoms_map(vector <site *> &kontener)
{
	site *pointer=0;
	
	kontener.clear();
	
	for(unsigned int i=0;i<atom_list.size();i++)
	{
		pointer=atom_list[i];
		
		if(check_site_belonging_to_sim_area(pointer))
		{
		int atom = atom_list[i]->get_atom();
		
		if(atom >= 0)
		{
			kontener.push_back(pointer);
			}
		}
	}
//		control_output<<"set atom map: "<<endl;
//	for(int i=0;i<kontener.size();i++)
//	{
//		 control_output<<"nr "<<i<<" "<<kontener[i]<<" "<<endl;
//		kontener[i]->show_site();
//		}	
	
}	

site* lattice :: get_site(long pozition)
{
	return sim_atom_list[pozition];
}
	
void lattice :: set_atoms_list(vector <site *> &kontener, int typ)
{
	site *pointer=0;
	long int counter=0;
	kontener.clear();
	
	for(unsigned int i=0;i<atom_list.size();i++)
	{
		pointer=atom_list[i];
		
		if(check_site_belonging_to_sim_area(pointer))
		{
		int atom = atom_list[i]->get_atom();
		
		if(atom==typ){
			if(typ==0){pointer->set_vindex(counter); counter++;}	// jak wakancja to ustaw vindex
			kontener.push_back(pointer);
			}
		}
	}
	
	if(typ==0)
	{control_output<<"set atom list typ/size: "<<typ<<" / "<<kontener.size()<<endl;
	
//	for(int i=0;i<kontener.size();i++)
//	{
//		control_output<<"nr "<<i<<" "<<kontener[i]<<" "<<kontener[i]->get_vindex()<<endl;
//		kontener[i]->show_site();
//		kontener[i]->show_neigh(1);
//	}	
	}
}	

void lattice :: add_sublatt_typ(int sublatt, int atom_typ)
	{
		wektor typ_sublat(atom_typ,sublatt,0);
		
		vector <wektor> :: iterator I;
		wektor new_sublatt= typ_sublat;
		int add_new_type = 1;
		int new_sublattice = 1;
		//if(sublatt_typ.size() <1)
		//{
		//	sublatt_typ.push_back(typ_sublat);
		//}
		//else
		{
			for(I=sublatt_typ.begin();I!=sublatt_typ.end();I++)
			{
				wektor old_sublatt = wektor(*I);
				
				if(sublatt==old_sublatt.y)
					{
						new_sublattice=0;
					}
				
				if(old_sublatt == new_sublatt)	// sprawdza czy typ jest juz w lisice
				{
					add_new_type = 0;	
				}
			}
		}
		
	if(new_sublattice)
	{	
		sublattice++;
	}
	
	if(add_new_type)
	{
		sublatt_typ.push_back(typ_sublat);
	}
	
	control_output<<"latt/typ:"<< sublatt<<"/"<<atom_typ<<" |sublatt_typ size "<<sublatt_typ.size()<<" nr_sublattice "<<sublattice<<endl;
}
	
	
void lattice :: add_atom_type(int _atom, string name)
{
	int iter=0;
	int new_atom = _atom;
	int add_new_type = 1;
	vector <int> :: iterator I;
	//control_output<<"type size: "<<atoms_type.size()<<endl;
	
	if(atoms_type.size() < 1)	//zawsze na poczatek dodaje wakancje do pustej listy
	{
		atoms_type.push_back(0);
		atoms_name.push_back("Fe");
	//	control_output<<"in add at. new typ: "<< 0 <<" "<<"Fe"<<endl;	

	}
				// jesli juz cos jest dodane to dodaje kolejny
	for(I=atoms_type.begin();I!=atoms_type.end();I++)
		{
			int old_atom = int(*I);		// wczytuje typ bedacy w liscie
		//	control_output<<"in add at. old typ: "<< old_atom<<" "<<atoms_name[iter]<<endl;
			iter++;
			if(old_atom == new_atom)	// sprawdza czy typ jest juz w lisice
				{
				add_new_type = 0;	
				}
		}
	
	
	if(add_new_type)
	{atoms_type.push_back(new_atom);
	atoms_name.push_back(name);
	//control_output<<"in add at. new typ: "<< new_atom<<" "<<name<<endl;	
	}
	
	control_output<<"size: |typ "<<atoms_type.size()<<" |name: "<<atoms_name.size()<<endl;
}

int lattice :: get_size(int typ){
	if(typ ==1){return x_size;}
	else if(typ ==2){return y_size;}
	else if(typ ==3){return z_size;}
	else{control_output<<"ERROR in lattice::get_size() "<<typ<<endl;exit(1);}
}

unsigned int lattice :: get_atom_typ_numbers()
{
	unsigned int vec_size=atoms_type.size();
	return vec_size;
}

int lattice :: get_atom_type(int typ)
{
	return atoms_type[typ];
}

int lattice :: get_atom_type(string name)
{
	int pozycja=0;
	for(unsigned int i=0;i<atoms_name.size();i++)
	{
		if(name == atoms_name[i])
			break;
			
		pozycja++;
	}
	
	return atoms_type[pozycja];
}

int lattice :: get_vec_lattice_typ_size()
{
	int vec_size=sublattice;
	return vec_size;	
}

void lattice :: get_window(site* A, site* V, vector <site> &tab){
	
	vector <site*> Aneigh;
	vector <site*> Vneigh;
	vector <site> tmp_window;
	tmp_window.reserve(10);
	tmp_window.clear();
	V->read_site_neighbours(Vneigh,1,0);  //- atomy sasiedzi	
	A->read_site_neighbours(Aneigh,1,0); // - 0 enegia
	
	for(unsigned int vn=0; vn<Vneigh.size();vn++){
		for(unsigned int jn=0; jn<Aneigh.size();jn++){
			if(*(Vneigh[vn])==*(Aneigh[jn])){
				//int typ_win = A[jn]->get_atom();
				tmp_window.push_back(*(Vneigh[vn]));	//typ_window
			}
		}
	}
	
	double boundary_conx = boundary_con_at.x;
	double boundary_cony = boundary_con_at.y;
	double boundary_conz = boundary_con_at.z;

	double latt_constx = get_latice_const(1,0);
	double latt_consty = get_latice_const(2,0);
	double latt_constz = get_latice_const(3,0);

//	double st_areax=st_sim_area.x;
//	double st_areay=st_sim_area.y;
//	double st_areaz=st_sim_area.z;
	
//	double ed_areax=end_sim_area.x;
//	double ed_areay=end_sim_area.y;
//	double ed_areaz=end_sim_area.z;

	double N0X=A->get_x();
	double N0Y=A->get_y();
	double N0Z=A->get_z();	

	double N1X=V->get_x();
	double N1Y=V->get_y();
	double N1Z=V->get_z();

//	control_output<<"Row:"<<endl;
//	control_output<<"ATOM0: "<<N0X<<" "<<N0Y<<" "<<N0Z<<endl;
//	control_output<<"ATOM1: "<<N1X<<" "<<N1Y<<" "<<N1Z<<endl;
//	for (int ii=0;ii<tmp_window.size();ii++){
//		tmp_window[ii].show_site();
//	}

													//sprawdz warunki brzegowe i przesun

	//		(set_prec(x0-r_maxx) < set_prec(st_areax))						
	//		double r2tr = (x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 


	
//	double RN0X=N0X-N0X;
//	double RN0Y=N0Y-N0Y;
//	double RN0Z=N0Z-N0Z;

	double RN1X=N1X-N0X;
	double RN1Y=N1Y-N0Y;
	double RN1Z=N1Z-N0Z;
	
	double leftx=(-1.5*latt_constx);
	double lefty=(-1.5*latt_consty);
	double leftz=(-1.5*latt_constz);

	double rightx=(1.5*latt_constx);
	double righty=(1.5*latt_consty);
	double rightz=(1.5*latt_constz);

//	control_output<<leftx<<" "<<rightx<<" "<<latt_constx<<endl;
//	control_output<<lefty<<" "<<righty<<" "<<latt_consty<<endl;
//	control_output<<leftz<<" "<<righty<<" "<<latt_constz<<endl;

	
	if(RN1X < leftx){RN1X += boundary_conx;}
	if(RN1Y < lefty){RN1Y += boundary_cony;}
	if(RN1Z < leftz){RN1Z += boundary_conz;}

	if(RN1X > rightx){RN1X -= boundary_conx;}
	if(RN1Y > righty){RN1Y -= boundary_cony;}
	if(RN1Z > rightz){RN1Z -= boundary_conz;}
	
	//wektory
	double Rxy=sqrt(pow(RN1X,2) + pow(RN1Y,2));
	double Rxyz=sqrt(pow(RN1X,2) + pow(RN1Y,2) + pow(RN1Z,2));
	//policz katy o jakie przekrecic uklad L i B

	double sinL=(RN1Z)/(Rxyz);
	double cosL=(Rxy)/(Rxyz);
	double sinB=(RN1Y)/(Rxy);
	double cosB=(RN1X)/(Rxy);

	if(Rxy==0){sinB=0;cosB=1;}

	//przeliczyc window na wzgledne pozycje
	for (unsigned int i=0;i<tmp_window.size();i++){
													//sprawdz warunki brzegowe i przesun

		double WX=(tmp_window[i].get_x())-N0X;
		double WY=(tmp_window[i].get_y())-N0Y;
		double WZ=(tmp_window[i].get_z())-N0Z;

		if(WX < leftx){WX += boundary_conx;}
		if(WY < lefty){WY += boundary_cony;}
		if(WZ < leftz){WZ += boundary_conz;}

		if(WX > rightx){WX -= boundary_conx;}
		if(WY > righty){WY -= boundary_cony;}
		if(WZ > rightz){WZ -= boundary_conz;}

		double r=((WX)*(cosB) + (WY)*(sinB));

		double new_x=((r)*(cosL) + (WZ)*sinL);
		double new_y=((-WX)*(sinB) + (WY)*(cosB));
		double new_z=((-r)*(sinL) + (WZ)*(cosL));
		tmp_window[i].set_x(new_x);
		tmp_window[i].set_y(new_y);
		tmp_window[i].set_z(new_z);		


//		tmp_window[i].set_x(WX);
//		tmp_window[i].set_y(WY);
//		tmp_window[i].set_z(WZ);		
		}

//	control_output<<"Absloute:"<<endl;
//	control_output<<"ATOM0: "<<RN0X<<" "<<RN0Y<<" "<<RN0Z<<endl;
//	control_output<<"ATOM1: "<<RN1X<<" "<<RN1Y<<" "<<RN1Z<<endl;
//	for (int ii=0;ii<tmp_window.size();ii++){
//		tmp_window[ii].show_site();
//	}

//	for (int i=0;i<tmp_window.size();i++){

//		double WX=(tmp_window[i].get_x());
//		double WY=(tmp_window[i].get_y());
//		double WZ=(tmp_window[i].get_z());
//		double r=((WX)*(cosB) + (WY)*(sinB));
		
//		double new_x=((r)*(cosL) + (WZ)*sinL);
//		double new_y=((-WX)*(sinB) + (WY)*(cosB));
//		double new_z=((-r)*(sinL) + (WZ)*(cosL));
//		tmp_window[i].set_x(new_x);
//		tmp_window[i].set_y(new_y);
//		tmp_window[i].set_z(new_z);		
//	}

//		control_output<<"Rotated:"<<endl;
//		control_output<<"Rxy/Rxyz: "<<Rxy<<" "<<Rxyz<<endl;
//		control_output<<"sinL/cosL/sinB/cosB: "<<sinL<<" "<<cosL<<" "<<sinB<<" "<<cosB<<endl;
//		for (int ii=0;ii<tmp_window.size();ii++){
//			tmp_window[ii].show_site();
//		}
	

	
	tab=tmp_window;
		
}

void lattice :: sort_atoms(vector <site> &atoms){

//	control_output<<"Sorting..."<<endl;
	vector <double> Z;
	Z.reserve(10);
	Z.clear();
	vector <double> Y;
	Y.reserve(10);
	Y.clear();
	vector <int> I;
	I.reserve(10);
	I.clear();
	//mapowanie sitow na z i y
	for (unsigned int i=0;i<atoms.size();i++){
		double z = atoms[i].get_z();
		double y = atoms[i].get_y();
		Z.push_back(z);
		Y.push_back(y);
		I.push_back(i);
	}
//	control_output<<"Sizes:"<<Z.size()<<" "<<Y.size()<<" "<<I.size()<<endl;
	
//	for(int i=0;i<Z.size();i++){
//		control_output<<Z[i]<<" ";
//	}
//	control_output<<endl;
//	for(int i=0;i<Y.size();i++){
//		control_output<<Y[i]<<" ";
//	}
//	control_output<<endl;
//	for(int i=0;i<I.size();i++){
//		control_output<<I[i]<<" ";
//	}
//	control_output<<endl;
	
	unsigned int licznik=0;
	while(licznik < I.size()){
		double minZ=Z[licznik];
		double minY=Y[licznik];
		int minI=I[licznik];
		unsigned int index=licznik;
//		control_output<<"przed: "<<licznik<<" "<<index<<" "<<minZ<<" "<<minY<<" "<<minI<<endl;

	
		for(unsigned int i=licznik;i<Z.size();i++){
			if(minZ > Z[i]){
				minZ=Z[i];
				minY=Y[i];
				minI=I[i];
				index=i;
			}
			else if(minZ==Z[i]){
				if(minY > Y[i]){
					minY=Y[i];
					minI=I[i];
					index=i;
				}
			}	
		}
//porownalem element[licznik] z wszystkimi w tablicy, mam minimum teraz podmien elementy z licznik
		if(index != licznik){
//			control_output<<"podmieniam..."<<endl;

			Z[index]=Z[licznik];
			Z[licznik]=minZ;
			Y[index]=Y[licznik];
			Y[licznik]=minY;
			I[index]=I[licznik];
			I[licznik]=minI;
		}
		licznik++;
//		control_output<<"po: "<<licznik<<" "<<index<<" "<<minZ<<" "<<minY<<" "<<minI<<endl;

	}
//zapisz atoms do tmp_atoms wedlug I, ktore przechowuje index wedlug posortowanych z i y
//	control_output<<"Sorted:"<<endl;
//	for(int i=0;i<Z.size();i++){
//		control_output<<Z[i]<<" ";
//	}
//	control_output<<endl;
//	for(int i=0;i<Y.size();i++){
//		control_output<<Y[i]<<" ";
//	}
//	control_output<<endl;
//	for(int i=0;i<I.size();i++){
//		control_output<<I[i]<<" ";
//	}
//	control_output<<endl;


	vector <site> tmp;
	tmp.reserve(10);
	tmp.clear();
//	control_output<<"Saving..."<<endl;
	for(unsigned int i=0;i<I.size();i++){
		int index = I[i];
//		control_output<<index<<endl;
		site atom(atoms[index]);
//		atom.show_site();
		tmp.push_back(atom);
	}

//	control_output<<tmp.size()<<endl;
//	for (int ii=0;ii<tmp.size();ii++){
//		tmp[ii].show_site();
//	}

	atoms.clear();
	atoms=tmp;
	
//	for (int ii=0;ii<atoms.size();ii++){
//		atoms[ii].show_site();
//	}

	
}

void lattice::sim_atoms_list_init()
{	//int o;
	control_output<<"sim atom list init.."<<endl;
	//cin>>o;
	vector <site*> tmp_vector;
	sim_atom_list.reserve(50000);
	sim_atom_list.clear();
	
	

	for(unsigned int i=0;i<x_size;i++)
	{
	for(unsigned int j=0;j<y_size;j++)
	{
		
		for(unsigned int k=0;k<z_size;k++)
		{
			tmp_vector.clear();
			tmp_vector.reserve(50);
		//	control_output<<i<<" "<<j<<" "<<k<<endl;  
	//		matrix[i][j][k].show_sity();
			matrix[i][j][k].get_sity_in_box(tmp_vector);
			
			
			for(unsigned int l=0;l<tmp_vector.size();l++) 
			{
				//int ATOM=tmp_vector[l]->get_atom();
	//			control_output<<i<<" "<<j<<" "<<k;
		
				
				
					site *wsk_site=0;
					wsk_site=tmp_vector[l];
					
				//	control_output<<"wskaznik "<<wsk_site<<endl;
					if(check_site_belonging_to_sim_area(wsk_site))
					{
						int atom = wsk_site->get_atom();
						
					//	control_output<<"Atom in sim: "<<atom<<endl;
						
						if(atom >= 0)
						{
							sim_atom_list.push_back(wsk_site);
							//control_output<<" "<<atom<<" "<<tmp_vector[l]<<endl;
						}
					}			
			}
}}}

//	for(int i=0;i<sim_atom_list.size();i++)
//		{
//		control_output<<i<<" "<<sim_atom_list[i]<<" ";
//		sim_atom_list[i]->show_site();
//		}
control_output<<"Atoms list - ok: "<<sim_atom_list.size()<<endl;
	//cin>>o;	
}	



void lattice::atoms_list_init()
{	//int o;
	control_output<<"atom list init.."<<endl;
	//cin>>o;
	vector <site*> tmp_vector;
	atom_list.reserve(500000);
	atom_list.clear();
	
	

	for(unsigned int i=0;i<x_size;i++)
	{
	for(unsigned int j=0;j<y_size;j++)
	{
		
		for(unsigned int k=0;k<z_size;k++)
		{
			tmp_vector.clear();
	//		control_output<<i<<" "<<j<<" "<<k<<endl;
	//		matrix[i][j][k].show_sity();
			matrix[i][j][k].get_sity_in_box(tmp_vector);
			
			
			for(unsigned int l=0;l<tmp_vector.size();l++) 
			{
	//			int ATOM=tmp_vector[l]->get_atom();
	//			control_output<<i<<" "<<j<<" "<<k;
	//			control_output<<" "<<ATOM<<" "<<tmp_vector[l]<<endl;
				
				
					site * wsk_site=tmp_vector[l];
					atom_list.push_back(wsk_site);
				
			}
}}}

	//for(int i=0;i<atom_list.size();i++)
//		{
//	control_output<<i<<" "<<atom_list[i]<<" ";
	//	atom_list[i]->show_site();
	//	}
	control_output<<"Atoms list - ok: "<<atom_list.size()<<endl;
	//cin>>o;	
}	


void lattice :: set_alg_objects(list <pairjump> &evt, vector < vector <double> > &bar, potential &pot){
	
	POTENCIALY=&pot;
	BARIERY=&bar;
	EVENTY=&evt;

	cout<<"Adres obiekt potencjal w lattice"<<POTENCIALY<<endl;
	cout<<"Adres obiekt bar w lattice"<<BARIERY<<endl;
	cout<<"Adres obiekt event w lattice"<<EVENTY<<endl;

}

potential& lattice :: get_potentials(){
	return *POTENCIALY;
}

long lattice :: get_sim_atom_number()
{
	return sim_atom_list.size();
}

site * lattice :: search_site(site *A)
{	
	//cout<<"w search jestem"<<endl;
	site *_old_site=0;
	
	//vector <site> :: iterator K;
	int iter=0;
/*
	for(K=atom_list.begin();K!=atom_list.end();K++)
	{
		if( (site(*K)) == (site(*A)) )
			{	cout<<"Znalazlem taki site "<<endl;;
				cout<<"adres "<<&K<<endl;
				_old_site=&(*K);
				iter++;
				cout<<iter<<endl;
			}
	}
	*/
	
	
	for(unsigned int i=0;i<atom_list.size();i++)
		{
	
			if( site(atom_list[i]) == (site(A)) )	//wykorzystany tu jest domyslny konstruktor site(const site &A) jesli site(*A)
			{	
				_old_site = atom_list[i];
				iter++;
				//cout<<"Szukanie Pierwsze|Drugie: "<<_old_site<<"| "<<&atom_list[i]<<endl;
			}
	}
	
	
	if(iter>1){
		control_output<<"In function search_site in lattice"<<endl;
		control_output<<"ERROR: Many identical sites in sample. Check structure definition!"<<endl;
		exit(0);
	}
	
	return _old_site;
	}

double lattice :: get_latice_transition(int direction)
{
	double a=0.0;
	
	
	if(direction==1)
	a=boundary_con_at.x;
	
	if(direction==2)
	a=boundary_con_at.y;
	
	if(direction==3)
	a=boundary_con_at.z;
	
	
	return a;
	}

double lattice :: get_latice_const(int direction, int i)
{	
	// i okresla ktora stala sieci z jakiej struktury, kolejnosc jak w structure.in
	double a=0.0;
	
	if(direction==1)
	a=x_trans[i];	//stala sieci w x pierwsza dodana z pliku structure.in
	
	if(direction==2)
	a=y_trans[i];
	
	if(direction==3)
	a=z_trans[i];
	
	return a;
	
	
}

bool lattice :: check_site_belonging_to_region(site *A)
		{
			
			   
			double x = (A->get_x());
			double y = (A->get_y());
			double z = (A->get_z());  
			
			//control_output<<" "<<x<<" "<<y<<" "<<z<<endl;
		//	control_output<<st_sim_area.x<<" "<<end_sim_area.x<<endl;
		//	control_output<<st_sim_area.y<<" "<<end_sim_area.y<<endl;
		//	control_output<<st_sim_area.z<<" "<<end_sim_area.z<<endl;
				if((set_prec(x)>=set_prec(st_region.x)) and (set_prec(x)<set_prec(end_region.x)))
				{
				if((set_prec(y)>=set_prec(st_region.y)) and (set_prec(y)<set_prec(end_region.y)))
				{
				if((set_prec(z)>=set_prec(st_region.z)) and (set_prec(z)<set_prec(end_region.z)))	
				{
					
					return true;
				}}}
			
			return false;
		}

bool lattice :: check_site_belonging_to_sim_area(site *A)
		{
			
			   
			double x = (A->get_x());
			double y = (A->get_y());
			double z = (A->get_z());  
			
			//control_output<<" "<<x<<" "<<y<<" "<<z<<endl;
		//	control_output<<st_sim_area.x<<" "<<end_sim_area.x<<endl;
		//	control_output<<st_sim_area.y<<" "<<end_sim_area.y<<endl;
		//	control_output<<st_sim_area.z<<" "<<end_sim_area.z<<endl;
				if((set_prec(x)>=set_prec(st_sim_area.x)) and (set_prec(x)<set_prec(end_sim_area.x)))
				{
				if((set_prec(y)>=set_prec(st_sim_area.y)) and (set_prec(y)<set_prec(end_sim_area.y)))
				{
				if((set_prec(z)>=set_prec(st_sim_area.z)) and (set_prec(z)<set_prec(end_sim_area.z)))	
				{
					
					return true;
				}}}
			
			return false;
		}

void lattice::get_sity_from_nnbox(int x,int y, int z,int latt_num,vector <site*> &tmp_atom_list)
	{
		if(local_control_atom==control_atom)
		control_output<<"Getting sity from nnboxes... "<<endl;
		
		tmp_atom_list.clear();
		vector <long> used_boxes(27);
		
		int X=0;
		int Y=0;
		int Z=0;
		int border=0;
		//obcinam do jednosci bo komorki pamieci wskazuje
		int r_maxx=int(interaction_zone.x);
		int r_maxy=int(interaction_zone.y);
		int r_maxz=int(interaction_zone.z);
		
		int xsize=x_size;	//int((end_region.x-st_region.x)/x_trans[latt_num]);
		int ysize=y_size;	//int((end_region.y-st_region.y)/y_trans[latt_num]);
		int zsize=z_size;	//int((end_region.z-st_region.z)/z_trans[latt_num]);
		
		int left_borderx = 0;	//int(st_region.x/x_trans[latt_num]);
		int left_bordery = 0;	//int(st_region.y/y_trans[latt_num]);
		int left_borderz = 0;	//int(st_region.z/z_trans[latt_num]);
		
		
		int right_borderx = x_size;	//int(end_region.x/x_trans[latt_num]);
		int right_bordery = y_size;	//int(end_region.y/y_trans[latt_num]);
		int right_borderz = z_size;	//int(end_region.z/z_trans[latt_num]);
		
	if(local_control_atom==control_atom)
{	
	control_output<<" "<<x<<" "<<y<<" "<<z<<endl;
	control_output<<"bef "<<tmp_atom_list.size()<<endl;
	control_output<<" "<<left_borderx<<" "<<left_bordery<<" "<<left_borderz;
	control_output<<" "<<right_borderx<<" "<<right_bordery<<" "<<right_borderz;
	control_output<<" "<<xsize<<" "<<ysize<<" "<<zsize<<endl;
}	
		for(int i=-r_maxx;i<=r_maxx;i++)
		{
			for(int j=-r_maxy;j<=r_maxy;j++)
			{
				for(int k=-r_maxz;k<=r_maxz;k++) 
				{	
					border=0;
					X=x+i;
					Y=y+j;
					Z=z+k;
					
					if(local_control_atom==control_atom)					
					control_output<<"XYZ: "<<(X)<<" "<<(Y)<<" "<<(Z)<<endl;
 
					if((X)>=0 and (Y)>=0 and (Z)>=0)
						{
							if((unsigned)X<x_size and (unsigned)Y<y_size and unsigned(Z)<z_size)
							{
								if(local_control_atom==control_atom)
								control_output<<"hole"<<endl;
								
								int wasnt_used=1;
								long boxid=matrix[X][Y][Z].get_box_id();
								vector <long> :: iterator K;
								
								if(local_control_atom==control_atom)
								control_output<<" box id: "<<boxid<<endl;
							
								
								for(K=used_boxes.begin();K!=used_boxes.end();K++)							
									{
									if(long(*K)==boxid)
										{
											wasnt_used=0;
										}
									}
								
								if(wasnt_used)
								{matrix[X][Y][Z].get_sity_in_box(tmp_atom_list);
								used_boxes.push_back(boxid);
								if(local_control_atom==control_atom)
									matrix[X][Y][Z].show_sity();

								}
							}
						}

					
					if(X<left_borderx)
						{X=X+xsize;
						border=1;}
					 
					if(Y<left_bordery)
						{Y=Y+ysize;
						border=1;}
					
					if(Z<left_borderz)
						{Z=Z+zsize;
					border=1;}
					
					if(X>=right_borderx)
						{X=X-xsize;
					 border=1;}
					 
					if(Y>=right_bordery)
						{Y=Y-ysize;
					border=1;}
					
					if(Z>=right_borderz)
						{Z=Z-zsize;
					border=1;}
					
					if(local_control_atom==control_atom)				
					{control_output<<"XYZ aft borders: "<<(X)<<" "<<(Y)<<" "<<(Z)<<endl;
				//	if((x+i)>=0 and (y+j)>=0 and (z+k)>=0)
					}
						if(border)
						{	
							if(local_control_atom==control_atom)
							control_output<<"border: "<<endl;
							long boxid=matrix[X][Y][Z].get_box_id();
							if(local_control_atom==control_atom)
							control_output<<" box id: "<<boxid<<endl;
							
							int wasnt_used=1;
							vector <long> :: iterator K;
								
							for(K=used_boxes.begin();K!=used_boxes.end();K++)							
							{
								if(long(*K)==boxid)
								{
									wasnt_used=0;
									}
								
								
								}
							if(wasnt_used)
							{matrix[X][Y][Z].get_sity_in_box(tmp_atom_list);
							used_boxes.push_back(boxid);
							if(local_control_atom==control_atom)
								matrix[X][Y][Z].show_sity();

							}
							}
		
		if(local_control_atom==control_atom)
		control_output<<"aft "<<tmp_atom_list.size()<<endl;	
					
		}	}	}
	if(local_control_atom==control_atom)
	control_output<<"aft "<<tmp_atom_list.size()<<endl;	
	} 
	
void lattice :: jumps_shell_init(){
	control_output<<"Initialized jump zone: ";
//	cout<<"Init jump zone "<<endl;
	vector <int> cor_zone;
	cor_zone.reserve(4);
	int zone=-1;
	double rmin=Rmin,rmax=Rmax;
	for(unsigned int K=0;K<atom_list.size();K++)
	{
		vector <site*> atom_neigh;
		sites_zone_init(1,rmin,rmax,atom_list[K],atom_neigh);
		atom_list[K]->put_neighbours(atom_neigh,1);

		int z =atom_neigh.size();
		if(z != zone)
		{
			if(cor_zone.size() == 0)
			{
				cor_zone.push_back(z);
				zone=z;
			}
			else
			{
				int r=1;
				for(unsigned int i=0; i<cor_zone.size();i++)
				{
					if(cor_zone[i] == z)
					{
						r=0;
					}
				}
				if(r)
				{
					cor_zone.push_back(z);
					zone=z;
				//	cout<<"Dodaje nowy zone"<<endl;
				//	control_atom=atom_list[K]->get_atom();
				//	atom_list[K]->show_site();
				//	sites_zone_init(1,rmin,rmax,atom_list[K],atom_neigh);
				//	int o;
				//	cin>>o;
				}
			}
		}
		
	//Sprawdza
	//atom_list[K]->show_neigh(typ);
		double ATOM=atom_list[K]->get_atom();
		if(ATOM==control_atom)
		{
			atom_list[K]->show_site();
			for(unsigned int i=0; i<atom_neigh.size();i++)
			{
				control_output<<endl;
				control_output<<i<<" ";
				atom_neigh[i]->show_site();
			}		
		}	
	}
	
	control_output<<"with coordination zones: ";
	
	max_coordination_number=0;
	for(unsigned int i=0; i<cor_zone.size();i++)
			{	
				int cor_nr=cor_zone[i];
				control_output<<cor_nr<<" ";
				if(max_coordination_number<cor_nr)
				{
					max_coordination_number=cor_nr;
				}
			}
			control_output<<endl;
}	

int lattice :: get_max_coordination_number(){
	return max_coordination_number;
}

void lattice :: interaction_shell_init(){
	control_output<<"Initiated interaction zone ";	
	//TUTAJ PETLA OD pot->rmin[i]
	unsigned  int zones_number= POTENCIALY->get_coordination_number();
	double rmin=0.0,rmax=0.0;
	vector <int> cor_zone;
	cor_zone.reserve(4);
	int zone=-1;
	
	
	for(unsigned int K=0;K<atom_list.size();K++)
	{
		for(unsigned int i=0;i<zones_number;i++)
		{
	//	cout<<"iter: "<<i<<endl;
		POTENCIALY->get_interaction_zone(rmin,rmax,i);
		vector <site*> atom_neigh;
		sites_zone_init(0,rmin,rmax,atom_list[K],atom_neigh);
		atom_list[K]->put_neighbours(atom_neigh,0);
	//	atom_list[K]->show_neigh(0);
		int z =atom_neigh.size();
		if(z != zone)
		{
			if(cor_zone.size() == 0)
			{
				cor_zone.push_back(z);
				zone=z;
			}
			else
			{
				int r=1;
				for(unsigned int i=0; i<cor_zone.size();i++)
				{
					if(cor_zone[i] == z)
					{
						r=0;
					}
				}
				if(r)
				{
					cor_zone.push_back(z);
					zone=z;
				}
			}
		}

	//Sprawdza
	//atom_list[K]->show_neigh(typ);
		double ATOM=atom_list[K]->get_atom();
		if(ATOM==control_atom)
		{
			atom_list[K]->show_site();
			atom_list[K]->show_neigh(0);
			//for(int i=0; i<atom_neigh.size();i++)
			//{
			//	cout<<i<<" "<<atom_neigh[i]<<" ";
			//	atom_neigh[i]->show_site();
			//}		
		}
			
	}
	}
	control_output<<"with coordination zones: ";
	for(unsigned int i=0; i<cor_zone.size();i++)
	{
		control_output<<cor_zone[i]<<" ";
	}
	control_output<<endl;
	
}	

void lattice :: sites_zone_init(int typ, double rmin, double rmax, site *atom_list, vector <site*> &atom_neigh)
{
	//wczytuje sasiadow atomu atom_list do tablicy sasiadow w klasie site
	
	//control_output<<"Init interaction zone ";
	//cout<<"Init interaction zone "<<endl;

	double r_maxx=rmax;	//interaction_zone.x;
	double r_maxy=rmax;	//interaction_zone.y;
	double r_maxz=rmax;	//interaction_zone.z;
		
	double boundary_conx=0.0;
	double boundary_cony=0.0;
	double boundary_conz=0.0;	
		
	if(typ)  // typ= 1 oznacza atomy
	{	
	//	control_output<<"for atom exchange"<<endl;
		boundary_conx = boundary_con_at.x;
		boundary_cony = boundary_con_at.y;
		boundary_conz = boundary_con_at.z;
	}
	else
	{
		
	//	control_output<<"for energy calculation"<<endl;
		boundary_conx = boundary_con_en.x;
		boundary_cony = boundary_con_en.y;
		boundary_conz = boundary_con_en.z;
	}
	

	double st_areax=st_region.x;
	double st_areay=st_region.y;
	double st_areaz=st_region.z;
	
	double ed_areax=end_region.x;
	double ed_areay=end_region.y;
	double ed_areaz=end_region.z;
	
	//control_output<<st_sim_area.x<<" "<<end_sim_area.x<<endl;
	//control_output<<st_sim_area.y<<" "<<end_sim_area.y<<endl;
	//control_output<<st_sim_area.z<<" "<<end_sim_area.z<<endl;
		
	int set_x = 1; 
	int set_y = 1;
	int set_z = 1;
	
//	int set_writex = 1;
//	int set_writey = 1;
//	int set_writez = 1;
	

//	for(int K=0;K<atom_list.size();K++)
	//{
		
		vector <site*> tmp_neigbours;
		tmp_neigbours.clear();
		tmp_neigbours.reserve(50);
	//	atom_list->clear_neighbours(typ);  //atom
		double x0=atom_list->get_x();
		double y0=atom_list->get_y();
		double z0=atom_list->get_z();
		double ATOM=atom_list->get_atom();
		int latt_num=atom_list->get_latt_number();

		//sprawdzam czy atom jest w obszarze symulacji   simul_atoms_init
	if(check_site_belonging_to_region(atom_list)){
	
	if(ATOM==control_atom)
	{	control_output<<"ATOM: "<<ATOM<<endl;
		local_control_atom=control_atom;
		atom_list->show_site();
	}	
	else
	{
		local_control_atom=-5;	
	}
			//warunki brzegowe dla atomu w obszarze symulacji, dla ktorego bedziemy szukac sasiadow!!!
			//okreslam jego pozycje wzgledem warunkow brzegowych
							
		if((set_prec(x0-r_maxx) < set_prec(st_areax)) or (set_prec(x0+r_maxx) >= set_prec(ed_areax)))  // do funkcji wektor lattice:: find_atom_BC(x0,y0,z0)
		{	
				if(ATOM==control_atom)
				control_output<<"atom na granicy x"<<endl;
			
			if(set_prec(boundary_conx) < 0)
			{set_x = 2;
			
				if(ATOM==control_atom)
				control_output<<"set x take nn site for all: 2"<<endl;
			
			}	// do 1 sasiada
			else if (set_prec(boundary_conx) == 0)
			{set_x = 0;
			
			if(ATOM==control_atom)
			control_output<<"set x no boundary condition: 0"<<endl;
			
			}	// anuluj wpisywanie
			else if (set_prec(boundary_conx) > 0)
			{
				int control=0;
				
				if(set_prec(x0-r_maxx) < set_prec(st_areax)) 		//przesun x0 o wektor +translate
				{set_x = -1;
				
				if(ATOM==control_atom)
				control_output<<"set x coppy cell on right: -1"<<endl;							//UWAGA na blad oba if prawdziwe tylo ze nie ma else to 2 if jest robiony
				
				control++;
				}
				if (set_prec(x0+r_maxx) >= set_prec(ed_areax))	//przesun o wektor -translate
				{set_x = 1;
			
				if(ATOM==control_atom)
				control_output<<"set x coppy cell on left: 1"<<endl;
			
				control++;
				}
				if(control>=2)
				{
					control_output<<"Bad interractions radius -> cell_vector to big  "<<endl;
					exit(0);}
			}
					
			}
		else
		{set_x = 2;
		
		if(ATOM==control_atom)
		control_output<<"x atom w srodku"<<endl;
	
		}	// do 1 sasiada
		
			//to samo dla y i z
		if((set_prec(y0-r_maxy) < set_prec(st_areay)) or (set_prec(y0+r_maxy) >= set_prec(ed_areay)))
		{
			if(ATOM==control_atom)
			control_output<<"atom na granicy y"<<endl;
			if(set_prec(boundary_cony) < 0)
			{set_y = 2;
			if(ATOM==control_atom)
			control_output<<"set y take nn site for all: 2"<<endl;
			}	// do 1 sasiada
			else if (set_prec(boundary_cony) == 0)
			{set_y = 0;
			if(ATOM==control_atom)
			control_output<<"set y no boundary condition: 0"<<endl;
			}	// anuluj wpisywanie
			else if (set_prec(boundary_cony) > 0)
			{
				int control=0;
				if(set_prec(y0-r_maxy) < set_prec(st_areay)) 		//przesun o wektor +translate
				{set_y = -1;
				if(ATOM==control_atom)
				control_output<<"set y coppy cell on right: -1"<<endl;
				control++;
				}
				if (set_prec(y0+r_maxy) >= set_prec(ed_areay))	//przesun o wektor -translate
				{set_y = 1;
				if(ATOM==control_atom)
				control_output<<"set y coppy cell on left: 1"<<endl;
				control++;
				}
				if(control>=2)
				{
				control_output<<"Bad interractions radius -> cell_vector to big  "<<endl;
					exit(0);}
			}
					
			}
		else
		{set_y = 2;
		if(ATOM==control_atom)
		control_output<<"y atom w srodku"<<endl;
		}
		
		if((set_prec(z0-r_maxz) < set_prec(st_areaz)) or (set_prec(z0+r_maxz) >= set_prec(ed_areaz)))
		{
			if(ATOM==control_atom)
			control_output<<"atom na granicy z"<<endl;
			if(set_prec(boundary_conz) < 0)
			{set_z = 2;
			if(ATOM==control_atom)
			control_output<<"set z take nn site for all: 2"<<endl;
			}	// do 1 sasiada
			else if (set_prec(boundary_conz) == 0)
			{set_z = 0;
			if(ATOM==control_atom)
			control_output<<"set z no boundary condition: 0"<<endl;
			}	// anuluj wpisywanie
			else if (set_prec(boundary_conz) > 0)
			{
				int control=0;
				if(set_prec(z0-r_maxz) < set_prec(st_areaz)) 		//przesun z0 o wektor +translate
				{set_z = -1;
				if(ATOM==control_atom)
				control_output<<"set z coppy cell on right: -1"<<endl;
				control++;
				}
				if (set_prec(z0+r_maxz) >= set_prec(ed_areaz))	//przesun z0 o wektor -translate
				{set_z = 1;
				if(ATOM==control_atom)
				control_output<<"set z coppy cell on left: 1"<<endl;
				control++;
				}
			if(control>=2)
				{
					control_output<<"Bad interractions radius -> cell_vector to big  "<<endl;
					exit(0);}
			}
					
			}
		else
		{set_z = 2;

		if(ATOM==control_atom)
		control_output<<"z atom w srodku "<<endl;
		}
		//okreslilem pozycje wybranego atomu i teraz
		//dla wybranego atomu skanuje cala liste atomow w poszukiwaniu jego sasiadow
		
		//nowy obiekt box pelniacy role situ, bedzie zawieral sity mieszczace sie w danym boxie
		// vector sitow
		// put_site i get_site
		
		vector <site*> tmp_atom_list;
		tmp_atom_list.reserve(200);  //lista z adresami do bliskich boxow bedacych w sferze odidzialywania.
		//int o;
	//	cin>>o;
		
		if(control_atom==ATOM)
		{
			control_output<<"Checking myRound for get_sity_from_nnbox "<<endl;
			control_output<<int(x0/x_trans[latt_num])<<" "<<int(y0/y_trans[latt_num])<<" "<<int(z0/z_trans[latt_num])<<" "<<latt_num<<endl;
		}
		
		get_sity_from_nnbox(int(x0/x_trans[latt_num]),int(y0/y_trans[latt_num]),int(z0/z_trans[latt_num]),latt_num,tmp_atom_list);

	//	if(ATOM==control_atom)
//		control_output<<"tmp_atom_list size "<<tmp_atom_list.size()<<endl;

		// matrix[][][].get_box(tmp_atom_list);
		// matrix +1 razy 26 kombinacje +/-1 | +/-1 | +/-1 bo tyle jest najblizszych boxow
		//site *neighbour=0; 
//				neighbour = &matrix[i];
		//if (r2<=(rmax*rmax) and r2 >= (rmin*rmin))	
		// 		tmp_atom_list.pusch_back();	
		
		for(unsigned int i=0;i<tmp_atom_list.size();i++)
			{
				//wczytuje wspolrzedne potencjalnego sasiada
				double x=tmp_atom_list[i]->get_x();//I->get_x();
				double y=tmp_atom_list[i]->get_y();//I->get_y();
				double z=tmp_atom_list[i]->get_z();//I->get_z();
				site *neighbour=0;
				neighbour = tmp_atom_list[i];		//zapisz adres komorki ktora moze okazac sie sasiadem
		
				if(ATOM==control_atom)
				control_output<< set_x<<" "<<set_y<<" "<<set_z<<endl;
				
	//			int set_writex = 0;
	//			int set_writey = 0;
	//			int set_writez = 0;	
				
				if (check_site_belonging_to_region(neighbour))	//sprawdza czy site jest w obszerze symulacji
					{
					double r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
					
					// warunek dla wylaczonych warunkow brzegowych boundary =0
					// atom jest we wnetrzu obszaru symulacji
					if (set_prec(r2)<=set_prec(rmax*rmax) and set_prec(r2) >= set_prec(rmin*rmin))	//jesli sasiad znajduje sie w sferze odzidzialywania
					{			
							tmp_neigbours.push_back(neighbour);						//dodaj go do listy sasiadow danego sita glownego
					}
					
					////////////////atom jet na brzegu, warunki brzegowe aktywne, boundary >0
					if(boundary_conx > 0 )	// przesun komorke symulacji Ox o wektro boundary_conx w lewo lub prawo/ ste_x
					{

						double r2tr = (x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))	//sprawdz sfere odidzialywan
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0)	// y - sciany
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if( boundary_conz > 0 )
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x-x0)*(x-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);			
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0)
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(z-z0)*(z-z0); 
				
						//cout<<"r2 "<<r2tr<<endl;
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_conz > 0)	//przesun komorke Ox i Oz
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0 and boundary_conz > 0)	//y i z		- krawedzie
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0 and boundary_conz > 0) // z,y,z - rogi
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
				}	
				// to byli sasiezie jesli warungi brzegowe byly aktywne
				//teraz jesli sa nieaktywne, boundary<0
				///////////////// warunki dla sasiedniego atomu   SPRAWDZIC
				//~(check_site_belonging_to_sim_area(neighbour))
				
				if(boundary_conx < 0)	// bierze sasiadow spoza strefy symulacji do nn
				{
				if((x<st_sim_area.x) or (x>=end_sim_area.x))	// w kierunku x - sciana
				{
				if((y>=st_sim_area.y) and (y<end_sim_area.y))
				{
				if((z>=st_sim_area.z) and (z<end_sim_area.z))	
				{	
				
				double r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
					// wnetrze
					if (set_prec(r2)<=set_prec(rmax*rmax) and set_prec(r2) >= set_prec(rmin*rmin))
					{
						tmp_neigbours.push_back(neighbour);	
					}
					
							if(boundary_conx > 0 )	// przesun komorke symulacji Ox o wektro boundary_conx w lewo lub prawo/ ste_x
					{

						double r2tr = (x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))	//sprawdz sfere odidzialywan
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0)	// y - sciany
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if( boundary_conz > 0 )
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x-x0)*(x-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);			
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0)
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(z-z0)*(z-z0); 
				
						//cout<<"r2 "<<r2tr<<endl;
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_conz > 0)	//przesun komorke Ox i Oz
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0 and boundary_conz > 0)	//y i z		- krawedzie
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0 and boundary_conz > 0) // z,y,z - rogi
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
				}}}
				}
				
				if(boundary_cony < 0)
				{
				if((y<st_sim_area.y) or (y>=end_sim_area.y))		//sciana y +/-
				{
				if((x>=st_sim_area.x) and (x<end_sim_area.x))
				{
				if((z>=st_sim_area.z) and (z<end_sim_area.z))	
				{	
				double r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
					// wnetrze
					if (set_prec(r2)<=set_prec(rmax*rmax) and set_prec(r2) >= set_prec(rmin*rmin))
					{
						tmp_neigbours.push_back(neighbour);	
					}
					
										if(boundary_conx > 0 )	// przesun komorke symulacji Ox o wektro boundary_conx w lewo lub prawo/ ste_x
					{

						double r2tr = (x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))	//sprawdz sfere odidzialywan
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0)	// y - sciany
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if( boundary_conz > 0 )
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x-x0)*(x-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);			
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0)
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(z-z0)*(z-z0); 
				
						//cout<<"r2 "<<r2tr<<endl;
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_conz > 0)	//przesun komorke Ox i Oz
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0 and boundary_conz > 0)	//y i z		- krawedzie
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0 and boundary_conz > 0) // z,y,z - rogi
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
				}}}
				}
				
				
				if(boundary_conz < 0)
				{
				if((z<st_sim_area.z) or (z>=end_sim_area.z))
				{
				if((x>=st_sim_area.x) and (x<end_sim_area.x))
				{
				if((y>=st_sim_area.y) and (y<end_sim_area.y))
				{
				double r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
					// wnetrze
					if (set_prec(r2)<=set_prec(rmax*rmax) and set_prec(r2) >= set_prec(rmin*rmin))
					{
						tmp_neigbours.push_back(neighbour);	
					}
					
										if(boundary_conx > 0 )	// przesun komorke symulacji Ox o wektro boundary_conx w lewo lub prawo/ ste_x
					{

						double r2tr = (x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))	//sprawdz sfere odidzialywan
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0)	// y - sciany
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if( boundary_conz > 0 )
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x-x0)*(x-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);			
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0)
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(z-z0)*(z-z0); 
				
						//cout<<"r2 "<<r2tr<<endl;
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_conz > 0)	//przesun komorke Ox i Oz
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0 and boundary_conz > 0)	//y i z		- krawedzie
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0 and boundary_conz > 0) // z,y,z - rogi
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
				}}}
				}
				
				
				if(boundary_conx < 0 and boundary_cony < 0)		// krawedz z translate +/-x,+/-y
				{
				if((x<st_sim_area.x) or (x>=end_sim_area.x))
				{
				if((y<st_sim_area.y) or (y>=end_sim_area.y))
				{
				if((z>=st_sim_area.z) and (z<end_sim_area.z))	
				{	
				double r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
					// wnetrze
					if (set_prec(r2)<=set_prec(rmax*rmax) and set_prec(r2) >= set_prec(rmin*rmin))
					{
						tmp_neigbours.push_back(neighbour);	
					}
					
										if(boundary_conx > 0 )	// przesun komorke symulacji Ox o wektro boundary_conx w lewo lub prawo/ ste_x
					{

						double r2tr = (x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))	//sprawdz sfere odidzialywan
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0)	// y - sciany
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if( boundary_conz > 0 )
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x-x0)*(x-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);			
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0)
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(z-z0)*(z-z0); 
				
						//cout<<"r2 "<<r2tr<<endl;
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_conz > 0)	//przesun komorke Ox i Oz
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0 and boundary_conz > 0)	//y i z		- krawedzie
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0 and boundary_conz > 0) // z,y,z - rogi
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
				}}}
				}
			
				if(boundary_conx < 0 and boundary_conz < 0) // y
				{
				if((x<st_sim_area.x) or (x>=end_sim_area.x))
				{
				if((z<st_sim_area.z) or (z>=end_sim_area.z))
				{
				if((y>=st_sim_area.y) and (y<end_sim_area.y))
				{
				double r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
					// wnetrze
					if (set_prec(r2)<=set_prec(rmax*rmax) and set_prec(r2) >= set_prec(rmin*rmin))
					{
						tmp_neigbours.push_back(neighbour);	
					}
					
										if(boundary_conx > 0 )	// przesun komorke symulacji Ox o wektro boundary_conx w lewo lub prawo/ ste_x
					{

						double r2tr = (x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))	//sprawdz sfere odidzialywan
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0)	// y - sciany
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if( boundary_conz > 0 )
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x-x0)*(x-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);			
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0)
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(z-z0)*(z-z0); 
				
						//cout<<"r2 "<<r2tr<<endl;
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_conz > 0)	//przesun komorke Ox i Oz
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0 and boundary_conz > 0)	//y i z		- krawedzie
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0 and boundary_conz > 0) // z,y,z - rogi
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
				}}}
				}
				
				if(boundary_cony < 0 and boundary_conz < 0) //krawedzie x
				{
				if((y<st_sim_area.y) or (y>=end_sim_area.y))
				{
				if((z<st_sim_area.z) or (z>=end_sim_area.z))
				{
				if((x>=st_sim_area.x) and (x<end_sim_area.x))
				{
				
				double r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
					// wnetrze
					if (set_prec(r2)<=set_prec(rmax*rmax) and set_prec(r2) >= set_prec(rmin*rmin))
					{
						tmp_neigbours.push_back(neighbour);	
					}
					
										if(boundary_conx > 0 )	// przesun komorke symulacji Ox o wektro boundary_conx w lewo lub prawo/ ste_x
					{

						double r2tr = (x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))	//sprawdz sfere odidzialywan
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0)	// y - sciany
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if( boundary_conz > 0 )
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x-x0)*(x-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);			
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0)
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(z-z0)*(z-z0); 
				
						//cout<<"r2 "<<r2tr<<endl;
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_conz > 0)	//przesun komorke Ox i Oz
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0 and boundary_conz > 0)	//y i z		- krawedzie
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0 and boundary_conz > 0) // z,y,z - rogi
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
				}}}
				}
				
					
				if(boundary_conx < 0 and boundary_cony < 0 and boundary_conz < 0) // wierzcholki 
				{
				if((y<st_sim_area.y) or (y>=end_sim_area.y))
				{
				if((z<st_sim_area.z) or (z>=end_sim_area.z))
				{
				if((x<st_sim_area.x) or (x>=end_sim_area.x))
				{
				
				double r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
					// wnetrze
					if (set_prec(r2)<=set_prec(rmax*rmax) and set_prec(r2) >= set_prec(rmin*rmin))
					{
						tmp_neigbours.push_back(neighbour);	
					}
					
										if(boundary_conx > 0 )	// przesun komorke symulacji Ox o wektro boundary_conx w lewo lub prawo/ ste_x
					{

						double r2tr = (x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))	//sprawdz sfere odidzialywan
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0)	// y - sciany
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0)+(z-z0)*(z-z0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if( boundary_conz > 0 )
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x-x0)*(x-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);			
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0)
					{
						double r2tr = (y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(z-z0)*(z-z0); 
				
						//cout<<"r2 "<<r2tr<<endl;
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_conz > 0)	//przesun komorke Ox i Oz
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y-y0)*(y-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_cony > 0 and boundary_conz > 0)	//y i z		- krawedzie
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0)+(x-x0)*(x-x0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
					if(boundary_conx > 0 and boundary_cony > 0 and boundary_conz > 0) // z,y,z - rogi
					{
						double r2tr = (z+(set_z*boundary_conz)-z0)*(z+(set_z*boundary_conz)-z0)+(x+(set_x*boundary_conx)-x0)*(x+(set_x*boundary_conx)-x0)+(y+(set_y*boundary_cony)-y0)*(y+(set_y*boundary_cony)-y0); 
				
						if (set_prec(r2tr)<=set_prec(rmax*rmax) and set_prec(r2tr) >= set_prec(rmin*rmin))
						{
							tmp_neigbours.push_back(neighbour);	
						}
					}
					
				}}}
				}
			
			/////////////
							
				
			}
	}
			
	//		control_output<<"NEIGH at: "<<ATOM<<endl;
	//		atom_list->show_neigh(1);
	//		control_output<<"NEIGH en: "<<ATOM<<endl;
	//		atom_list->show_neigh(0);
			//cout<<"koniec SITA"<<endl;
	atom_neigh=tmp_neigbours;
//	atom_list->put_neighbours(tmp_neigbours, typ);
	//atom_list->show_neigh(typ);
		

//	if(ATOM==control_atom)
	//{
		//atom_list->show_site();
	//	for(int i=0; i<tmp_neigbours.size();i++)
		//{
	//		control_output<<i<<" ";
		//	tmp_neigbours[i]->show_site();
	//	}	
//	}
		
		
	
		
	/*
	 * zeby zmienic sita w boxie potrzeba funkcji w box
	 * matrix[1][1][1].sity[1].show_site();
		
	 * zeby odwolac sie do konkretnego sita w boxie trzeba podac koordynaty
	 * site* site_in_box_at(double x,double  y,double  z)
	 * {
		 * site Site(x,y,z);
		 * site * wsk=0;
		 * vector <site*> tmp;
		 * tmp.clear();
		 * tmp.reserve(2);
		 * 
		 * matrix[int(x)][][].get_sity_in_box(tmp);
		 *
		 *  for(int i=0;i<tmp.size();i++)
		 * 		{
			 * 		if ( site(tmp[i]) == Site )
		 * 			{
			 * 			wsk=tmp[i]	
			 * 		}		
			 * 		else
			 * 		{
				 * 		exit(0);
				 * 		cout<<"No atoms at: "<<x<<" "<<y<<" "<<z<<endl;
				 * 	}
		 * 		}
		 * return wsk;
		 * }
	 * */
	
//	}
	//koniec petli for K po atom_list
	//control_output<<endl;
}



void lattice :: check_neighbours(int typ)
{
	control_output<<"Printings neighbours... "<<endl;
	vector <site*> atom_inbox;
	for(int i=0;i<int(x_size);i++)
		{
			for(int j=0;j<int(y_size);j++)
			{
				for(int k=0;k<int(z_size);k++)
				{
					control_output<<"i j k "<<i<<" "<<j<<" "<<k<<endl;
					atom_inbox.clear();
					matrix[i][j][k].get_sity_in_box(atom_inbox);
					
					for(unsigned int l=0;l<atom_inbox.size();l++)
					{
						atom_inbox[l]->show_site();
						atom_inbox[l]->show_neigh(typ);
					}
					
				}
			}
		}
	
}

void lattice :: check_atoms()
{
	control_output<<"Printing atoms... "<<endl;
	vector <site*> atom_inbox;
	for(int i=0;i<int(x_size);i++)
		{
			for(int j=0;j<int(y_size);j++)
			{
				for(int k=0;k<int(z_size);k++)
				{
					control_output<<"i j k "<<i<<" "<<j<<" "<<k<<endl;
					atom_inbox.clear();
					matrix[i][j][k].get_sity_in_box(atom_inbox);
					
					for(unsigned int l=0;l<atom_inbox.size();l++)
					{
						atom_inbox[l]->show_site();
					}
					
				}
			}
		}
	
}


lattice :: ~lattice()
{
	}
	
	
/*------------------------------------------------------*/

void lattice::simulation_initialize(double _r_min, double _r_max, wektor a,wektor b,wektor c,wektor d, wektor _max_zone, wektor e, wektor f)
{
	control_output<<"Initialize boundary conditions... "<<endl;
	st_sim_area=a;
	end_sim_area=b;
	boundary_con_at=c;
	boundary_con_en=d;
	st_region=e;
	end_region=f;
	interaction_zone=_max_zone; 

	atoms_list_init();
	sim_atoms_list_init();

	Rmin=_r_min;
	Rmax=_r_max;

	control_output<<"Boundary conditions initialized "<<endl;

	wektor* boundary_pointer=0;

	boundary_pointer=&boundary_con_at;
	check_boundary_conditions(boundary_pointer);
	boundary_pointer=&boundary_con_en;
	check_boundary_conditions(boundary_pointer);

//	cout<<"simulation area: "<<endl;
//	st_sim_area.show_wektor();
//	end_sim_area.show_wektor();
//	cout<<"boundary translations |e and a|: "<<endl;
//	boundary_con_en.show_wektor();
//	boundary_con_at.show_wektor();

	
}


bool lattice :: check_boundary_conditions(wektor *wsk)
{
	control_output<<"Checking boundary conditions... "<<endl;
	int status=0;
	
//	double r_maxx=interaction_zone.x;
//	double r_maxy=interaction_zone.y;
//	double r_maxz=interaction_zone.z;
		
	double boundary_conx=0.0;
	double boundary_cony=0.0;
	double boundary_conz=0.0;	
		
	//	control_output<<"for atom exchange"<<endl;
	boundary_conx = wsk->x;
	boundary_cony = wsk->y;
	boundary_conz = wsk->z;
	
	double st_areax=st_region.x;
	double st_areay=st_region.y;
	double st_areaz=st_region.z;
	
	double ed_areax=end_region.x;
	double ed_areay=end_region.y;
	double ed_areaz=end_region.z;

	if(boundary_conx>0 and boundary_conx != ((ed_areax-st_areax)))
	{
		control_output<< "Wrong x translation vector in boundary condition. "<<endl;
		control_output<< "Or too large region area. "<<endl;
		control_output<<"Must be equal trans in structure.in file. if want PBC! "<<endl;
		exit(0);
		}
	if (boundary_cony>0 and boundary_cony != (ed_areay-st_areay))
	{
		control_output<< "Wrong y translation vector in boundary condition. "<<endl;
				control_output<< "Or too large region area. "<<endl;
		control_output<<"Must be equal trans in structure.in file. if want PBC! "<<endl;
		exit(0);
		}
	if (boundary_conz>0 and boundary_conz != (ed_areaz-st_areaz))
	{	
		control_output<< "Wrong z translation vector in boundary condition. "<<endl;
		control_output<< "Or too large region area. "<<endl;
		control_output<<"Must be equal trans in structure.in file. if want PBC! "<<endl;
		exit(0);
		}
	else
	{control_output<<"Boundary condition <=0 or PBC are fine."<<endl;
	status = 1;
	}
	control_output<<endl;
	return status;
}

/*------------------------------------------------------*/

void lattice :: refresh_structure(string file_name)
{
	ifstream file_structure(file_name.c_str(),ios :: in);
	string atom;
	double X;
	double Y;
	double Z;
	
//	int i=0;
//	int j=0;
//	int k=0; 
	long n =0;
	file_structure>>n;
	
	for(long l=0;l<n;l++)
	{
	file_structure>>atom>>X>>Y>>Z;
	
	
	int typ=get_atom_type(atom);
	int sublatt=get_sub_lattice_type(X,Y,Z,0);
	site ATOM(X,Y,Z,typ,sublatt,0);
	
	site* wsk=0;
	wsk=&ATOM;
	put_atom(int(X/x_trans[0]),int(Y/y_trans[0]),int(Z/z_trans[0]),wsk);
		
	
	}
	
	file_structure.close();
} 

/*-----------------------------------------------------------*/
	
void lattice :: reinit_sim_area(wektor a, wektor b){
	control_output<<"reinit_sim_area:"<<endl;
	st_sim_area.show();
	end_sim_area.show();
	st_sim_area += a;
	end_sim_area += b;
	st_sim_area.show();
	end_sim_area.show();

	sim_atoms_list_init();  	
}
	
/*-----------------------------------------------------------*/

void lattice::read_structure(string file_name,wektor start,wektor end,wektor set_st, int lattice_nr)
{	
		
//	int o;
	wektor set_vec=set_st;
	wektor st_vec=start;
	wektor end_vec=end;
	//set_vec.show_wektor();
	//st_vec.show_wektor();
	//end_vec.show_wektor();
	wektor delta= (end_vec-st_vec);
	//delta.show_wektor();
	wektor size = (set_vec+delta);
	        
	//size.show_wektor();
	//set_vec.show_wektor();
	//st_vec.show_wektor();
	//end_vec.show_wektor();
	//delta.show_wektor();
	
	long n =0;
	ifstream file_structure(file_name.c_str(),ios :: in);
	string atom;
	double X,dx=0.0;
	double Y,dy=0.0;
	double Z,dz=0.0;
	vector <long int> jumps(3,0);

	control_output<<"czytam sobie "<<file_name<<endl;
	file_structure>>n;
	control_output<<n<<endl;
	int i=int(st_vec.x/x_trans[0]);
	int j=int(st_vec.y/y_trans[0]);
	int k=int(st_vec.z/z_trans[0]);
	control_output<<"st "<<i<<" "<<j<<" "<<k<<endl;
	i=int(end_vec.x/x_trans[0]);
	j=int(end_vec.y/y_trans[0]);   
	k=int(end_vec.z/z_trans[0]);
	control_output<<"end "<<i<<" "<<j<<" "<<k<<endl; 
	int kounter=0;
//	int kounter2=0;
	
	for(long i=0;i<n;i++)
	{
	string line;
	while (getline(file_structure,line)){
	if(!line.empty()){
	double buf;
	vector <double> values;	
	stringstream ss(line);
	ss >> atom;
	while ( ss >> buf ){values.push_back(buf);}

//	control_output<<"wczytano: "<<values.size()<<" "<<atom;
//	for (int vi=0;vi<values.size();vi++){ control_output<<" "<<values[vi];}
//	control_output<<endl;
//	kounter++;	
	X=values[0];Y=values[1];Z=values[2];

	if(values.size() > 3){
		dx=values[3];dy=values[4];dz=values[5]; 
		for(unsigned int ni=6;ni<values.size();ni++){jumps[(ni-6)]=values[ni];}
	}
	
	if((X>=set_vec.x)&&(Y>=set_vec.y)&&(Z>=set_vec.z))
		{	//cout<<"wczytano1 "<<atom<<" "<<X<<" "<<Y<<" "<<Z<<endl;
			//cin>>o;
		if((X<size.x)&&(Y<size.y)&&(Z<size.z)){
			X=X-set_vec.x+st_vec.x;
			Y=Y-set_vec.y+st_vec.y;
			Z=Z-set_vec.z+st_vec.z;
			int typ=get_atom_type(atom);
			int sublatt=get_sub_lattice_type(X,Y,Z,lattice_nr);
			site ATOM(X,Y,Z,typ,sublatt,lattice_nr);
			
			if(values.size() > 3){		
				ATOM.set_drx(dx);ATOM.set_dry(dy);ATOM.set_drz(dz);
				ATOM.set_jumps(jumps);
			}
			site* wsk=0;
			wsk=&ATOM;
//				control_output<<"wczytano2 "<<X<<" "<<Y<<" "<<Z<<" "<<atom<<" "<<typ;
//				for(int ni=0;ni<jumps.size();ni++){control_output<<" "<<jumps[ni];}
//				control_output<<endl;
				kounter++;
			put_atom(int(X/x_trans[0]),int(Y/y_trans[0]),int(Z/z_trans[0]),wsk);
		//	cin>>o;
			
		//	if(atom=="Al")	
		//	matrix[X][Y][Z].set_atom(2);
		//	if(atom=="Fe")
		//	matrix[X][Y][Z].set_atom(0);
			//cout<<"wczytaj "<<atom<<endl;
			}
	
		}
	}}}
	file_structure.close();
	control_output<<"Total readed: "<<kounter<<endl;
	}	
  
/*----------------------------------------------------------*/


int lattice :: get_sub_lattice_type(double X, double Y, double Z, int lattice_nr)
{
	int sublatice_type=0;
//	int atom_type=0;
	double x = fmod(X,x_trans[lattice_nr]);
	double y = fmod(Y,y_trans[lattice_nr]);
	double z = fmod(Z,z_trans[lattice_nr]);
	
//	cout<<"For "<<x<<" "<<y<<" "<<z<<" find sublatice type"<<endl;
	
//	cells[lattice_nr] vector zawiera sity z komurki elementarnej dla sieci lattice_nr
	
	for(unsigned int i=0;i<cells[lattice_nr].size();i++)	//iteruje po sitach w cells
	{											//dla kazdego situ sprawdz czy x,y,z rownaja sie tym z bazy x,y,z
		
		double cx=cells[lattice_nr][i].get_x();	//pobierz wspolrzedne w cell
		double cy=cells[lattice_nr][i].get_y();
		double cz=cells[lattice_nr][i].get_z();
		
		//porownaj z dokladnoscia
		
		if(set_prec(x)==set_prec(cx) and set_prec(y)==set_prec(cy) and set_prec(z)==set_prec(cz) )	//jesli sie rowna to pobierz sub_lattice_type
		{
			sublatice_type=cells[lattice_nr][i].get_sub_latt();
//			atom_type=cells[lattice_nr][i].get_atom();
		}
	}
//	cout<<"sub lat: "<<sublatice_type<<" atom typ: "<<atom_type<<endl;
	return sublatice_type;
}

void lattice :: put_site(int x, int y, int z, site *Site)
{
//	control_output<<"put: "<<x<<" "<<y<<" "<<z<<endl;
	matrix[x][y][z].put_site(*Site);
	}
	
void lattice :: put_atom(int x, int y, int z, site *Site)
{
	//cout<<"put: "<<x<<" "<<y<<" "<<z<<endl;
//	matrix[x][y][z].show_sity();
	matrix[x][y][z].put_atom(*Site);
//	matrix[x][y][z].show_sity();
	}
  

void lattice :: makepic(long step,long step_break, wektor make_pic_vec_st, wektor make_pic_vec_ed, string name_of_file)
{
	stringstream total(name_of_file);
	int all=0,sum=1;
	
	int word_count=0 ;
    string word;
    while( total >> word ) ++word_count;
    
	if(word_count == 1)
	{
		sum=0;
	}
	else if(word_count == 2)
	{
		stringstream ss(name_of_file);
		int log=0;
		while(ss>>all){
			if(log==0){all-=1;log++;}
			else if(log>=1){log++;}
			else {control_output<<"ERROR in lattice::pic_diff: "<<all<<endl;exit(1);}
			sum *= all;
			}
	}else if(word_count>2){control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	else{control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	
		sum+=step;
		
	
    
	if( ((sum != 0) and (sum % step_break == 0)) or ((step == 0) and (step_break == 0)) )
	{
	double xs=make_pic_vec_st.x;
	double ys=make_pic_vec_st.y;
	double zs=make_pic_vec_st.z;

	double xe=make_pic_vec_ed.x;
	double ye=make_pic_vec_ed.y;
	double ze=make_pic_vec_ed.z;

//	int Nx=xe-xs;
//	int Ny=ye-ys;
//	int Nz=ze-zs;

	stringstream s;

	string name;
	s<<(sum);
	name=s.str()+"pic.xyz";
	ofstream file(name.c_str());

	file<<endl;	//UWAGA: dzielic przez vektor translacji x_trans , ... 
	file<<endl;
	file<<endl;
	file<<endl;
	file<<endl;
	file<<endl;
	file<<endl;
	//Nx=Nx-1;
	//Ny=Ny-1;
	//Nz=Nz-1;
	long atoms=0;
	//vector <site> :: iterator K;

	for(unsigned int K=0;K<atom_list.size();K++)
	{
		//K->show_site();
		double i=0.0;
		double j=0.0;
		double k=0.0;
		int atom = -3;
		
		i=atom_list[K]->get_x();
		j=atom_list[K]->get_y();
		k=atom_list[K]->get_z();
		atom=atom_list[K]->get_atom();

		if((i>=xs and i<=xe) and (j>=ys and j<=ye) and (k>=zs and k<=ze))
		{
			if(atom>-1)
			{	
			file<<atoms_name[atom]<<" "<<i<<" "<<j<<" "<<k<<endl;	
			atoms++;
			}
		}	
	
	}
file.seekp(0);	
file<<atoms<<endl;
file.close();	
}		

}

/*-------------------------------------------------------------------*/


	
wektor lattice::get_end_sim_wektor()
{
	wektor a=end_sim_area;
	
	return a; 	
	}
	
wektor lattice::get_PB(){
	return boundary_con_at;	
}

/*-------------------------------------------------------------------*/

double lattice :: get_stech(int typ1,int typ2){
	
	double stech=0;
	vector <site* > vtyp1;
	vtyp1.reserve(20000);
	vector <site* > vtyp2;
	vtyp2.reserve(20000);

	set_atoms_list(vtyp1, typ1);
	set_atoms_list(vtyp2, typ2);
	
	double N1=vtyp1.size();
	double N2=vtyp2.size();
	
	stech=N1/(N1+N2);
	return stech;

}
  

void lattice :: pic_stech(long step,double stech, wektor make_pic_vec_st, wektor make_pic_vec_ed, string name_of_file)
{
	
	double stechinsample=get_stech(1,2);
//	cout<<step<<" "<<stech<<" "<<stechinsample<<" "<<set_prec(stech,3)<<" "<<set_prec(stechinsample,3)<<endl;
	if(set_prec(stech,3) == set_prec(stechinsample,3) )
	{
	stringstream total(name_of_file);
	int all=0,sum=1;
	
	int word_count=0 ;
    string word;
    while( total >> word ) ++word_count;
    
	if(word_count == 1)
	{
		sum=0;
	}
	else if(word_count == 2)
	{
		stringstream ss(name_of_file);
		int log=0;
		while(ss>>all){
			if(log==0){all-=1;log++;}
			else if(log>=1){log++;}
			else {control_output<<"ERROR in lattice::pic_diff: "<<all<<endl;exit(1);}
			sum *= all;
			}
	}else if(word_count>2){control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	else{control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	
		sum+=step;
		
	double xs=make_pic_vec_st.x;
	double ys=make_pic_vec_st.y;
	double zs=make_pic_vec_st.z;

	double xe=make_pic_vec_ed.x;
	double ye=make_pic_vec_ed.y;
	double ze=make_pic_vec_ed.z;

//	int Nx=xe-xs;
//	int Ny=ye-ys;
//	int Nz=ze-zs;

	stringstream s;

	string name;
	s<<(sum);
	name=s.str()+"stech.xyz";
	ofstream file(name.c_str());

	file<<endl;	//UWAGA: dzielic przez vektor translacji x_trans , ... 
	file<<endl;
	file<<endl;
	file<<endl;
	file<<endl;
	file<<endl;
	file<<endl;
	//Nx=Nx-1;
	//Ny=Ny-1;
	//Nz=Nz-1;
	long atoms=0;
	//vector <site> :: iterator K;

	for(unsigned int K=0;K<atom_list.size();K++)
	{
		//K->show_site();
		double i=0.0;
		double j=0.0;
		double k=0.0;
		int atom = -3;
		
		i=atom_list[K]->get_x();
		j=atom_list[K]->get_y();
		k=atom_list[K]->get_z();
		atom=atom_list[K]->get_atom();

		if((i>=xs and i<=xe) and (j>=ys and j<=ye) and (k>=zs and k<=ze))
		{
			if(atom>-1)
			{	
			file<<atoms_name[atom]<<" "<<i<<" "<<j<<" "<<k<<endl;	
			atoms++;
			}
		}	
	
	}
file.seekp(0);	
file<<atoms;
file.close();	
}		

}

/*-------------------------------------------------------------------*/


void lattice :: pic_diff(long step,long step_break, wektor make_pic_vec_st, wektor make_pic_vec_ed, string name_of_file)
{
	stringstream total(name_of_file);
	int all=0,sum=1;
	
	int word_count=0 ;
    string word;
    while( total >> word ) ++word_count;
    
	if(word_count == 1)
	{
		sum=0;
	}
	else if(word_count == 2)
	{
		stringstream ss(name_of_file);
		int log=0;
		while(ss>>all){
			if(log==0){all-=1;log++;}
			else if(log>=1){log++;}
			else {control_output<<"ERROR in lattice::pic_diff: "<<all<<endl;exit(1);}
			sum *= all;
			}
	}else if(word_count>2){control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	else{control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	
		sum+=step;
	
	if( ((sum != 0) and (sum % step_break == 0)) or ((step == 0) and (step_break == 0)) )
	{
	double xs=make_pic_vec_st.x;
	double ys=make_pic_vec_st.y;
	double zs=make_pic_vec_st.z;

	double xe=make_pic_vec_ed.x;
	double ye=make_pic_vec_ed.y;
	double ze=make_pic_vec_ed.z;

//	int Nx=xe-xs;
//	int Ny=ye-ys;
//	int Nz=ze-zs;

	stringstream s;

	string name;
	s<<(sum);
	name=s.str()+"diff.xyz";
	ofstream file(name.c_str());

	file<<endl;	//UWAGA: dzielic przez vektor translacji x_trans , ... 
	file<<endl;
	file<<endl;
	file<<endl;
	file<<endl;
	file<<endl;
	file<<endl;
	//Nx=Nx-1;
	//Ny=Ny-1;
	//Nz=Nz-1;
	long atoms=0;
	//vector <site> :: iterator K;

	for(unsigned int K=0;K<atom_list.size();K++)
	{
		//K->show_site();
		double i=0.0;
		double j=0.0;
		double k=0.0;
		int atom = -3;
		
		i=atom_list[K]->get_x();
		j=atom_list[K]->get_y();
		k=atom_list[K]->get_z();
		atom=atom_list[K]->get_atom();
		
		double dx = atom_list[K]->get_drx();
		double dy = atom_list[K]->get_dry();
		double dz = atom_list[K]->get_drz();
		vector <long int> jumps;
		atom_list[K]->get_jumps(jumps);
	//	cout<<"jumps size: "<<jumps.size()<<endl;
	//	int o; cin >>o;
	//	for (int nj=0; nj<jumps.size();nj++){cout<<" "<<jumps[nj]; cin>>o;}
	//	cout<<endl;
		if((i>=xs and i<=xe) and (j>=ys and j<=ye) and (k>=zs and k<=ze))
		{
			if(atom>-1)
			{	
			file<<atoms_name[atom]<<" "<<i<<" "<<j<<" "<<k<<" "<<dx<<" "<<dy<<" "<<dz;
			for (unsigned int nj=0; nj<jumps.size();nj++){file<<" "<<jumps[nj];}
			file<<endl;	
			atoms++;
			}
		}	
	
	}
	file.seekp(0);	
	file<<atoms<<endl;
	file.close();	
	}		

}

/*-------------------------------------------------------------------*/



wektor lattice::get_st_sim_wektor()
{
	wektor a=st_sim_area;
	
	return a; 
	}

/*-------------------------------------------------------------------*//*------------------------------------------------------------------*/


void lattice :: get_sites(plaster &tmp){

	double L_B=0.0;
	double R_B=0.0;
//	double PBC=0.0;
	int direction = tmp.get_direction();
	double st = tmp.get_st();
	double end = tmp.get_end();

	if(direction == 1){L_B=0;R_B=x_size*2;}
	else if (direction == 2){L_B=0;R_B=y_size*2;}
	else if (direction == 3){L_B=0;R_B=z_size*2;}
	else{control_output<<"ERROR in lattice::get_sites -> wrong direction > 3!! "<<endl; exit(0);}

//	PBC=R_B-L_B;
//	if(st < L_B){st += PBC;}
//	if(st > R_B){st -= PBC;}
//	if(end < L_B){end += PBC;}
//	if(end > R_B){end -= PBC;}

	if(st < L_B){control_output<<"ERROR in lattice::get_sites -> reservuars reached border "<<st<<endl; exit(0);}
	if(st > R_B){control_output<<"ERROR in lattice::get_sites -> reservuars reached border "<<st<<endl; exit(0);}
	if(end < L_B){control_output<<"ERROR in lattice::get_sites -> reservuars reached border "<<end<<endl; exit(0);}
	if(end > R_B){control_output<<"ERROR in lattice::get_sites -> reservuars reached border "<<end<<endl; exit(0);}

	
	for(unsigned int i=0;i<atom_list.size();i++)
	{

		site *wsk_to_site = atom_list[i];
		double d=0;		//atom displacement
		
		if(direction==1)
		{
			d=atom_list[i]->get_x();
		}
		else if(direction==2)
		{
			d=atom_list[i]->get_y();
		}
		else if(direction==3)
		{
			d=atom_list[i]->get_z();
		}
		else
		{
			cout<<"Wrong direction number in lattice::get_sites parameters| x-1|y-2|z-3"<<endl;
			exit(1);
		}
	//		UWAGA zaokraglanie 	
//		control_output<<"spr: "<<st<<" "<<d<<" "<<end<<endl;
		if((st <= d )and (d < end) ){
			unsigned int id=tmp.get_index();
			if(tmp.get_name() == "block"){wsk_to_site->set_block_index(id);}
			if(tmp.get_name() == "hist"){wsk_to_site->set_hist_index(id);}	
			if(tmp.get_name() == "rez"){wsk_to_site->set_rez_index(id);}
			tmp.push_back(wsk_to_site);
		}
	}
}

void lattice :: get_sites(vector <plaster> &tmp){

	int direction = tmp[0].get_direction();

	
	for(unsigned int i=0;i<atom_list.size();i++)
	{

		site *wsk_to_site = atom_list[i];
		double d=0;		//atom placement
		
		if(direction==1)
		{
			d=wsk_to_site->get_x();
		}
		else if(direction==2)
		{
			d=wsk_to_site->get_y();
		}
		else if(direction==3)
		{
			d=wsk_to_site->get_z();
		}
		else
		{
			cout<<"Wrong direction number in lattice::get_sites parameters| x-1|y-2|z-3"<<endl;
			exit(1);
		}
		
		for(unsigned int ik=0; ik< tmp.size(); ik++){
	
			if(tmp[ik].get_st() <= d  and d < tmp[ik].get_end()){
		//		control_output<<"dla "<<d<<" zakres "<<miarka[ik]<<" "<<miarka[ik+1]<<" ";
		//		control_output<<"dodalem";
		//		control_output<<endl;
				unsigned int id=tmp[ik].get_index();
				if(id==ik){
					if(tmp[ik].get_name() == "block"){wsk_to_site->set_block_index(id);}
					if(tmp[ik].get_name() == "hist"){wsk_to_site->set_hist_index(id);}
					if(tmp[ik].get_name() == "rez"){wsk_to_site->set_rez_index(id);}
				}
				tmp[ik].push_back(wsk_to_site);
			}
		}


	}
}


void lattice :: save_hist_dR(string file_name, int direction, double Time, double st_bin, double size_bin, double end_bin)
{
	int bins_nr = int((end_bin-st_bin)/size_bin)+1;		//zero dodac
	
	
		vector < vector<long> > hist(atoms_type.size(),vector <long> (bins_nr,0)) ;
	//2D : typ atomu i histogram
	
	
	for(unsigned int i=0;i<sim_atom_list.size();i++)
	{
		int atom=sim_atom_list[i]->get_atom();
		double d=0;		//atom displacement
		
		if(direction==1)
		{
			d=sim_atom_list[i]->get_drx();
		}
		else if(direction==2)
		{
			d=sim_atom_list[i]->get_dry();
		}
		else if(direction==3)
		{
			d=sim_atom_list[i]->get_drz();
		}
		else
		{
			cout<<"Wrong direction number in HIST parameters| x-1|y-2|z-3"<<endl;
			exit(1);
		}
		
	if(atom>0)	//wylaczam wakancje ze statystki
	{	
		int d_bin=int(d/size_bin);		//sprawdzic zaokraglanie
		
		if(d_bin>st_bin and d_bin<end_bin)		// zliczam tylko to co jest w zakresie
		{
			
		//	cout<<"Bin range in HIST too small-> max_d: "<<d_bin<<endl;
		//	cout<<sim_atom_list[i];
		//	sim_atom_list[i]->show_site();
		//	exit(1);
			 
		
		if(d>=0)
		{
			d_bin=d_bin-int(st_bin/size_bin);
			hist[atom][d_bin]++;
		}
		
		if(d<0)
		{
			d_bin=d_bin-int(st_bin/size_bin);
			hist[atom][d_bin]++;
		}
		
	}
	}
	}
	
	
	
	
	for(unsigned int i=1;i<atoms_type.size();i++)		//wylaczylem wakancje. W lattice w vectore atoms_type wakancje maja pozycje 0 i nazwe Fe
	{
	stringstream dir;
	dir<<direction;
	
	string name_of_file=file_name+atoms_name[i]+dir.str()+".dat";
	ofstream out_data(name_of_file.c_str(),ios :: app);
	for(unsigned int j=0;j<hist[i].size();j++)
	{
		out_data<<Time<<" "<<(j*size_bin+st_bin)<<" "<<hist[i][j]<<endl;
	}
	out_data<<endl;	
	}
}

double lattice :: calc_energy()
{
	double totE=0.0;
//	omp_set_dynamic(1);	
	int MAX_THREADS = omp_get_max_threads();
	int LOCAL_THREADS = int(MAX_THREADS/2);

	#pragma omp parallel shared(totE) num_threads(LOCAL_THREADS)	
	{	
//		if(omp_get_thread_num()==0){cout<<"Threat numbers in E: "<<omp_get_num_threads()<<endl;}
		double E=0.0;
		#pragma omp for nowait schedule(runtime)
		for(unsigned int i=0;i<sim_atom_list.size();i++)
		{
			E += POTENCIALY->get_energy(sim_atom_list[i]);
		}
		#pragma omp critical(collectE)
		{
//		control_output<<"threadE: "<<omp_get_thread_num()<<" "<<E<<endl;
			totE += E;
		}
	}
//	omp_set_dynamic(0);
//	control_output<<"tot thread "<<totE<<endl;

	return (totE/2.0);
}

double lattice :: calc_energy_global()
{
	double totE=0.0;
//	omp_set_dynamic(1);	
	int MAX_THREADS = omp_get_max_threads();
	int LOCAL_THREADS = int(MAX_THREADS/2);

	#pragma omp parallel shared(totE) num_threads(LOCAL_THREADS)	
	{	
//		if(omp_get_thread_num()==0){cout<<"Threat numbers in E: "<<omp_get_num_threads()<<endl;}
		double E=0.0;
		#pragma omp for nowait schedule(runtime)
		for(unsigned int i=0;i<atom_list.size();i++)
		{
			E += POTENCIALY->get_energy(atom_list[i]);
		}
		#pragma omp critical(collectE)
		{
//		control_output<<"threadE: "<<omp_get_thread_num()<<" "<<E<<endl;
			totE += E;
		}
	}
//	omp_set_dynamic(0);
//	control_output<<"tot thread "<<totE<<endl;

	return (totE/2.0);
}


void lattice :: save_energy(double Time, double Step, string name, int setON)
{
	string file="";	
	stringstream total(name);
	
	int word_count=0 ;
    string word;
    while( total >> word ) ++word_count;
    
	if(word_count == 1)
	{
		file=name+"E";
	}
	else if(word_count == 2)
	{
		stringstream ss(name);
		int log=0;
		while(ss>>word){
			if(log==0){
			file=word+"E";log++;}
			else if(log>=1){log++;}
			else {control_output<<"ERROR in lattice::save_Natoms: "<<log<<endl;exit(1);}
			}
	}else if(word_count>2){control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	else{control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}

	
	double Etotal=calc_energy();
	string name_of_file= file + ".dat";
	ofstream out_data(name_of_file.c_str(),ios :: app);
	out_data<<" "<<Step<<" "<<Time<<" "<<Etotal<<endl;
	
	if(setON==1){
		Etotal=calc_energy_global();
		name_of_file= file + "global.dat";
		ofstream out_data(name_of_file.c_str(),ios :: app);
		out_data<<" "<<Step<<" "<<Time<<" "<<Etotal<<endl;
	}
	
}

void lattice :: save_Natoms(double Time, double Step, string name, int setON)
{		

	string file="";	
	stringstream total(name);
	
	int word_count=0 ;
    string word;
    while( total >> word ) ++word_count;
    
	if(word_count == 1)
	{
		file=name+"N";
	}
	else if(word_count == 2)
	{
		stringstream ss(name);
		int log=0;
		while(ss>>word){
			if(log==0){
			file=word+"N";log++;}
			else if(log>=1){log++;}
			else {control_output<<"ERROR in lattice::save_Natoms: "<<log<<endl;exit(1);}
			}
	}else if(word_count>2){control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	else{control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	
	string name_of_file = file + ".dat";
	ofstream out_data(name_of_file.c_str(),ios :: app);
	int size = atoms_type.size();
	vector <vector<long > > results(size,vector <long> (sublattice,0));
	vector <vector<long > > results_global(size,vector <long> (sublattice,0));
	int MAX_THREADS = omp_get_max_threads();
	int LOCAL_THREADS = int(MAX_THREADS/2);
//	omp_set_dynamic(1);
//	omp_set_num_threads(LOCAL_THREADS);
//	omp_set_num_threads(MAX_THREADS);

	#pragma omp parallel shared(results) num_threads(LOCAL_THREADS) 	
	{	
//		if(omp_get_thread_num()==0){control_output<<"Threat numbers in N: "<<omp_get_num_threads()<<endl;}
		#pragma omp for schedule(runtime)
		for(unsigned int i=0;i<atom_list.size();i++){
			int atom=atom_list[i]->get_atom();			
			int podsiec=atom_list[i]->get_sub_latt();
			#pragma omp critical(collectN)
			{
			//	cout<<"threadN: "<<omp_get_thread_num()<<" "<<results[0][0]<<endl;
				if(check_site_belonging_to_sim_area(atom_list[i])){
					results[atom][podsiec]++;
					}
				if(setON){
					results_global[atom][podsiec]++;
				}
			}
		}//koniec for
	}//koniec parallel
//	cout<<"tot_thread: "<<tot_results[0][0]<<endl;

	out_data<<" "<<Step<<" "<<Time;
	for(int i=0;i<size;i++){
		for(int j=0;j<sublattice;j++){
					out_data<<" "<<results[i][j];
	}}
	out_data<<endl;	
	out_data.close();
	
	if(setON==1){
		name_of_file= file + "global.dat";
		ofstream out_data(name_of_file.c_str(),ios :: app);
		out_data<<" "<<Step<<" "<<Time;
		for(int i=0;i<size;i++){
			for(int j=0;j<sublattice;j++){
				out_data<<" "<<results_global[i][j];
		}}
		out_data<<endl;	
		out_data.close();
	}

}

void lattice :: save_NandE(double Time, double Step, string name, int setON)
{		
	omp_set_nested(1);
//	omp_set_dynamic(0);
//	int MAX_THREADS = omp_get_max_threads();

	#pragma omp parallel sections shared(Time,Step,name) num_threads(2)
	{
//		int LOCAL_THREADS = int(MAX_THREADS/omp_get_num_threads());
//		omp_set_num_threads(LOCAL_THREADS);
//		if(omp_get_thread_num()==0){cout<<"MAIN threat numbers: "<<omp_get_num_threads()<<endl;}

		#pragma omp section
		{
			save_Natoms(Time,Step,name,setON);
		}//koniec one section

		#pragma omp section
		{
			save_energy(Time,Step,name,setON);
		}//koniec second section

		}//koniec omp parallel sections
		#pragma omp barrier
	omp_set_nested(0);
//	omp_set_dynamic(1);
//	omp_set_num_threads(MAX_THREADS);	 
}



/*------------------------------------------------------------------*/

string lattice :: get_file_name(string name, string format){

	string file="";	
	stringstream total(name);
	
	int word_count=0 ;
    string word;
    while( total >> word ) ++word_count;
    
	if(word_count == 1)
	{
		file=name+format;
	}
	else if(word_count == 2)
	{
		stringstream ss(name);
		int log=0;
		while(ss>>word){
			if(log==0){
			file=word+format;log++;}
			else if(log>=1){log++;}
			else {control_output<<"ERROR in lattice::save_Natoms: "<<log<<endl;exit(1);}
			}
	}else if(word_count>2){control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	else{control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	
	return file;
}

void lattice :: save_SRO(double Time, double Step, string name){
	
	
	string name_of_file = get_file_name(name,"sro") + ".dat";
	ofstream out_data(name_of_file.c_str(),ios :: app);

	unsigned int size = get_atom_typ_numbers();
	unsigned  int zones= POTENCIALY->get_coordination_number();	
	vector < vector <vector<long > > > results(zones, vector <vector<long > > (size,vector <long> (size,0) ) );
	
	//cout<<results.size()<<" "<<results[0].size()<<" "<<results[0][0].size()<<endl;
	
	int MAX_THREADS = omp_get_max_threads();
	
	#pragma omp parallel shared(results) num_threads(MAX_THREADS) 	
	{	
//		if(omp_get_thread_num()==0){control_output<<"Threat numbers in N: "<<omp_get_num_threads()<<endl;}
		#pragma omp for schedule(runtime)
		for(unsigned int i=0;i<atom_list.size();i++){
			if(check_site_belonging_to_sim_area(atom_list[i])){
				vector <site*> neighbours;
				atom_list[i]->read_site_neighbours(neighbours,1);
				int typ1 = atom_list[i]->get_atom();			
				
				for(unsigned int k =0;k<neighbours.size();k++){
					if(check_site_belonging_to_sim_area(neighbours[k])){	
						int typ2 = neighbours[k]->get_atom();
						unsigned int zone = POTENCIALY->check_coordination_zone(atom_list[i],neighbours[k]);		
						#pragma omp critical(collectSRO)
						{
//							cout<<"threadN: "<<omp_get_thread_num()<<endl;
						results[zone][typ1][typ2]++;
						}		
					}
				}
			}
		}
			 
	}//koniec parallel

	out_data<<" "<<Step<<" "<<Time;
	for(unsigned int i=0;i<zones;i++){
		for(unsigned int j=0;j<size;j++){
			for(unsigned int k=0;k<size;k++){
				out_data<<" "<<results[i][j][k];
	}}}
	out_data<<endl;	
	out_data.close();

}


void lattice :: save_dR(double Time, long Step, string name)
{
	string name_of_fileR="";	
	string name_of_fileR2="";	
	stringstream total(name);
	
	int word_count=0 ;
    string word;
    while( total >> word ) ++word_count;
    
	if(word_count == 1)
	{
			name_of_fileR=name+"dR.dat";
			name_of_fileR2=name+"dR2.dat";
	}
	else if(word_count == 2)
	{
		stringstream ss(name);
		int log=0;
		while(ss>>word){
			if(log==0){
			name_of_fileR=word+"dR.dat";
			name_of_fileR2=word+"dR2.dat";log++;}
			else if(log>=1){log++;}
			else {control_output<<"ERROR in lattice::save_Natoms: "<<log<<endl;exit(1);}
			}
	}else if(word_count>2){control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	else{control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
	
	

	ofstream out_data(name_of_fileR.c_str(),ios :: app);
	out_data.precision(10);

	ofstream out_data1(name_of_fileR2.c_str(),ios :: app);
	out_data1.precision(10);
	
	int size=atoms_type.size();
	
	vector <wektor> results1(size,wektor());
	vector <wektor> results3(size,wektor());
	vector <double> results2_0(size,0.0);
	vector <double> results2_1(size,0.0);
	vector <double> results2_2(size,0.0);
	vector <double> results4(size,0.0);

#pragma omp parallel shared(results1,results3,results4,results2_0,results2_1,results2_2)	
	{
//	if(omp_get_thread_num()==0){cout<<"Threat numbers R: "<<omp_get_num_threads()<<endl;}

	#pragma omp for nowait schedule(runtime)

	for(unsigned int i=0;i<sim_atom_list.size();i++){	
		int atom=sim_atom_list[i]->get_atom();
		if(atom<0){cout<<"Atom type <0 in save_dR "<<endl;exit(0);}

		double dx=sim_atom_list[i]->get_drx();
		double dy=sim_atom_list[i]->get_dry();
		double dz=sim_atom_list[i]->get_drz();
		long int j0=sim_atom_list[i]->get_jumps(0);
		long int j1=sim_atom_list[i]->get_jumps(1);
		long int j2=sim_atom_list[i]->get_jumps(2);

		#pragma omp critical(collectR)
		{
//		cout<<"thread: "<<omp_get_thread_num()<<" "<<results1[0].x<<endl;
			results1[atom].x=results1[atom].x+dx;
			results1[atom].y=results1[atom].y+dy;
			results1[atom].z=results1[atom].z+dz;
			results2_0[atom]=results2_0[atom]+j0;
			results2_1[atom]=results2_1[atom]+j1;
			results2_2[atom]=results2_2[atom]+j2;
			results3[atom].x=results3[atom].x+dx*dx;
			results3[atom].y=results3[atom].y+dy*dy;
			results3[atom].z=results3[atom].z+dz*dz;
			results4[atom]=results4[atom]+1;
		}
	}//koniec for
	}//koniec omp parallel


	#pragma omp parallel sections num_threads(2)
	{

	#pragma omp section
	{
	out_data<<Step<<" "<<Time;
		for(unsigned int i=0;i<results2_0.size();i++){
			out_data<<" "<<results1[i].x<<" "<<results1[i].y<<" "<<results1[i].z<<" "<<results2_0[i]<<" "<<results2_1[i]<<" "<<results2_2[i];
		}
		out_data<<endl;	
		out_data.close();
	}
	#pragma omp section
	{	
		out_data1<<Step<<" "<<Time;
		for(unsigned int i=0;i<results4.size();i++){
			out_data1<<" "<<results3[i].x<<" "<<results3[i].y<<" "<<results3[i].z<<" "<<results4[i];
		}
		out_data1<<endl;
		out_data1.close();
	}
	}//koniec omp sections	
}
	
void lattice :: clear_dR()
{
	//cout<<"Clear displacements"<<endl;
	for(unsigned int i=0;i<sim_atom_list.size();i++)
	{
		sim_atom_list[i]->refresh_site();
	}
}


void lattice :: update_events(site* sajt){
	
	sajt->show_site();
	control_output<<"prze events: "<<EVENTY->size()<<endl;
	if( check_site_belonging_to_sim_area(sajt) ){	
	update_site_events(sajt);
	//search for vacancy in neigh for sajt and update
	vector <site*> neighs;
	typedef vector <site*>::iterator iters;
	sajt->read_site_neighbours(neighs,1,0);
	for( iters it=neighs.begin(); it != neighs.end();++it){
		if( check_site_belonging_to_sim_area((*it)) ){	
			update_site_events( (*it) );			
		}
	}
	}
	control_output<<"po events: "<<EVENTY->size()<<endl;
}

void lattice :: update_site_events(site* sajt){
	
//	control_output<<"\nupdating site events: "<<EVENTY->size()<<endl;
//	sajt->show_site();

	typedef list <pairjump>::iterator it2list;
	vector < it2list > to_del;
	sajt->get_events_index(to_del);
	

	for( unsigned int i=0; i<to_del.size();i++){
//		control_output<<"del event: "<<&(*to_del[i])<<endl;
//		(to_del[i])->show();

		(*EVENTY).erase(to_del[i]);
		}
	sajt->clear_events_index();

//	control_output<<"deleted site events: "<<EVENTY->size()<<endl;
	

	typedef vector <pairjump>::iterator iterev;
	vector <pairjump> tmp_events; 
	int typ = sajt->get_atom();
	if(typ==0){
		create_events_index(sajt,tmp_events);
	}

	list <pairjump>::iterator point2l; 
	for( iterev it=tmp_events.begin(); it != tmp_events.end();++it){
		site* tmp=(*it).get_vac_to_jump();
		if(sajt == tmp){
		point2l = EVENTY->insert( EVENTY->end(),(*it));
		tmp->add_events_index(point2l);
		
//		control_output<<"added events: "<<&(*point2l)<<endl;

//		it->show();
//		point2l->show();

		}else{
		cout<<"Error in mc::update_site_events. Atom tries jump to atom."<<endl;exit(1);
		}
	}
//	control_output<<"added site events: "<<EVENTY->size()<<endl;
	
}

void lattice :: create_events_index(site* siteV, vector <pairjump> &tmp_events){

//	control_output<<"Create events"<<endl;
	if(check_site_belonging_to_sim_area(siteV)){
	vector <pairjump> skoki;
	skoki.reserve(50);
	vector <site*> neighbour;
	siteV->read_site_neighbours(neighbour,1,0);
	
	for(unsigned int k =0;k<neighbour.size();k++){
		if(check_site_belonging_to_sim_area(neighbour[k])){	
			int atom = neighbour[k]->get_atom();		//wczytaj typ atomu sasiada atomowego wakancji
			unsigned int zone = ( POTENCIALY->check_coordination_zone(siteV,neighbour[k]) );	
			double E1= ( POTENCIALY->get_energy(siteV) ) + ( POTENCIALY->get_energy(neighbour[k]) ) - ( POTENCIALY->get_energy(siteV,neighbour[k]) );
			double E2= ( POTENCIALY->get_energy(siteV,atom) ) + ( POTENCIALY->get_energy(neighbour[k],0) ) 
			- ( POTENCIALY->get_energy(siteV,0,neighbour[k],0) ) - ( POTENCIALY->get_energy(siteV,atom,neighbour[k],atom) ) 
			+ ( POTENCIALY->get_energy(siteV,atom,neighbour[k],0) );
			double barrier=1000.0;
			barrier=(*BARIERY)[atom][zone];				
			double bariera=(E1+E2)/2+barrier-E1;		// 	--policz bariere
			pairjump tmp(siteV,neighbour[k],E1,E2,barrier,bariera);
			skoki.push_back(tmp);
		}
	}
	tmp_events=skoki;
}
//	control_output<<"After create events: "<<tmp_events.size()<<endl;
}



