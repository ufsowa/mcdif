
//Program do wstestowania modelu symulacji MC  da ukladu otwartego dla sieci kwadratowej*/
    
    /*
    02_11_2015: increase tab variable for input structure from 5 to 15
    05_11_2015: fix widom_rnd()/execute()/
    13_11_2015: change int to long step in direct,widom, remove sub_steps
    */
    
   
/****************************************************
 
  - dodac wczytywanie podsieci do site(,,,,sublatice_typ) w lattice - ok
  - sprawdzic warunki brzegowe i inna konfiguracje = ok
  - uniezaleznic liczenie dyfuzji od stalej sieci -ok
  - przerobic wszystko zeby pracowalo na matrix a nie na atom_list. To znaczy adresy maja pokazywac na komorki matrixa. - ok
  - dodac wczytywanie struktur = ok
  - dodac exchange_atoms = ok
  - dodac resident time zmodyfikowac zeby dowolne bariery mogl wczytywac do obiektu task = ok  // przeniesc do structure z conf.in
  * ZMIENIC RESIDENT TIME, skoki tylko licz = ok // czestotliwosci przeskokow
  - dodac obliczanie parametrow = ok
  - dodac sgcmc
  - dodac direc_exchange
  - dodac schedule
  * sprawdzic adresy = ok
  * sprawdzic czy dobrze przepisuje = ok
* svn
 ************************************************/   
int debug_mode =0;   
   
#include "mc.h"
    
//using namespace std;
  


						// kontenery na ktorych idzie cala symulacja
						
vector <site *> Vatoms;
set <site *> V_LIST;

list <pairjump> EVENTS;
unsigned long EVENTS_size = EVENTS.size();
//list <pairjump> *wsk_EVENTS=&EVENTS;

potential pot;
vector <vector <double> > simple_bars;
//vector <vector <double> > *wsk_simple_bars=&simple_bars;
opcja opt_equi(pot, simple_bars,V_LIST);
//3D
double bar_list[3][2][2][3][3][3][3][3][3];
static long RTA_energy_executions=-1;

/*--------------------------------------------------------------*/
void exchange_atoms(long atoms,int from,int to,lattice *sample)
{
control_output<<"Exchanging atoms..."<<endl;
vector <site *> Aatoms;
vector <site *> Batoms;


sample->set_atoms_list(Aatoms,from);
sample->set_atoms_list(Batoms,to);

long Asize=Aatoms.size();
long Bsize=Batoms.size();

if(atoms<=Asize)
{
    for(long i=0;i<=atoms;i++)
    {
    long R=rnd()*Asize;
    Aatoms[R]->set_atom(to);
    }

}

if(atoms<=Bsize)
{
    for(long i=0;i<=atoms;i++)
    {
    long R=rnd()*Bsize;
    Batoms[R]->set_atom(from);

    }

}
control_output<<"...removed..."<<endl;

//setA(sample);
//setB(sample);
//setV(sample);

}

void insert_atoms(int number,int from_typ, int to_typ,lattice *sample)
{
	control_output<<"Removing atoms..."<<endl;
	vector <site *> Aatoms;
	//vector <site *> Batoms;
	int iter=0;
	sample->set_atoms_list(Aatoms,from_typ);
	//sample->set_atoms_list(Batoms,to_typ);

	long Asize=Aatoms.size();
	//long Bsize=Batoms.size();
	
	while(iter < number)
	{
		long R=rnd()*Asize;
		int atom=Aatoms[R]->get_atom();
		if(atom != to_typ)
		{
			Aatoms[R]->set_atom(to_typ);
			iter++;
			sample->set_atoms_list(Aatoms,from_typ);
			Asize=Aatoms.size();
		}
	}
	
	
	//cout<<number<<" "<<iter<<" atoms "<<from_typ<<" exchange to "<<to_typ<<endl;

	control_output<<iter<<" atoms "<<from_typ<<" exchange to "<<to_typ<<endl;

}

/*

void setmap(lattice *sample)
{
	
Aatoms.clear();
Batoms.clear();
Vatoms.clear();
map.clear();	
	
wektor Nstart = sample->get_st_sim_wektor();
wektor Nend = sample->get_end_sim_wektor();
int Nxs=Nstart.x;
int Nys=Nstart.y;
int Nzs=Nstart.z;

int Nxe=Nend.x;
int Nye=Nend.y;
int Nze=Nend.z;


long size_atom_list = sample->get_atoms_number();

for (long i=0;i<size_atom_list;i++)
{
	
	}


for(int i=Nxs;i<Nxe;i++)
{
for(int j=Nys;j<Nye;j++)
    {
    for(int k=Nzs;k<Nze;k++)
        {
        int atom=sample->get_atom(i,j,k);
        if(atom==1)
        {
        Aatoms.push_back(site(i,j,k,atom));
	}
        if(atom==2)
        {
        Batoms.push_back(site(i,j,k,atom));
	}

        if(atom==0)
        {
        Vatoms.push_back(site(i,j,k,atom));
	}
	if(atom>-1)
	map.push_back(site(i,j,k,atom));
	}    
    }

}

cout<< "map"<<endl;
vector <site> :: iterator K;

for(K=map.begin();K!=map.end();K++)
{

int X=K->get_x();
int Y=K->get_y();
int Z=K->get_z();
int ATOM=K->get_atom();

cout<<X<<" "<<Y<<" "<<Z<<" "<<ATOM<<" || ";

}
cout<<endl;

}
*/


/*--------------------------------------------------------------

void set_atoms_list(lattice *sample, int type)
{
if(type==1)
Aatoms.clear();
if(type==2)
Batoms.clear();
if(type==0)
Vatoms.clear();



wektor Nstart = sample->get_st_sim_wektor();
wektor Nend = sample->get_end_sim_wektor();
int Nxs=Nstart.x;
int Nys=Nstart.y;
int Nzs=Nstart.z;

int Nxe=Nend.x;
int Nye=Nend.y;
int Nze=Nend.z;


for(int i=Nxs;i<Nxe;i++)
{
for(int j=Nys;j<Nye;j++)
    {
    for(int k=Nzs;k<Nze;k++)
        if(sample->get_atom(i,j,k)==type)
        
        {
		if(type==1)
		Aatoms.push_back(site(i,j,k,1));
		if(type==2)
		Batoms.push_back(site(i,j,k,2));
		if(type==0)
		Vatoms.push_back(site(i,j,k,0));
	}
    }

}

cout<< "V map"<<endl;
vector <site> :: iterator K;
for(K=Vatoms.begin();K!=Vatoms.end();K++)
{

int X=K->get_x();
int Y=K->get_y();
int Z=K->get_z();
int ATOM=K->get_atom();

cout<<X<<" "<<Y<<" "<<Z<<" "<<ATOM<<" || ";

}

cout<<"Ustawil A"<<endl;
}



-----------------------------------------------------------*/
void sgcmc(lattice *sample,long number_of_steps,double T, vector <double> &chem)
{
	double beta=1.0/(kB*T);	
	//cout<<"Robie sgcmc: pot chem: "<<chem.size()<<" ";
	//for( int i=0;i<chem.size();i++)
	//{
	//	cout<<chem[i]<<" ";
	//}
	long Nsize=sample->get_sim_atom_number();//sample->get_size();
	int Ntyp = sample->get_atom_typ_numbers();
	//cout<<"atoms: "<<Nsize<<" types: "<<Ntyp<<endl;
	
	for(long m=0;m<number_of_steps;m++)
	{
//losuj atom z listy atomow symulacji
		long N=(long)(rnd()*Nsize);
		//cout<<<<" ";
//pobierz adres wylosowanego atomu
		site *rnd_site=0;
		rnd_site=sample->get_site(N);
		int old_typ = rnd_site->get_atom();
		//cout<<"Adres rnd_sita: "<<&rnd_site<<" Adres trzymany w rnd_site: "<<rnd_site<<endl;
		//rnd_site->show_site();
		
//losuj typ atomu z listy typow atomow
		int typ=(int)(rnd()*Ntyp);
		int new_typ=0;
		new_typ=sample->get_atom_type(typ);
		//cout<<typ<<" "<<rnd_typ<<" ";
		if(new_typ != old_typ)					//warunek wyklucza podstawienie identyczne
		{
//policz energie przed podmiana
			double old_E=pot.get_energy(rnd_site);
//zamien typ
			rnd_site->set_atom(new_typ);
//policz energie po podmianie
			double new_E=pot.get_energy(rnd_site);
//ustal potencial chemiczny

										//vector chem /vector pot_list istnieje lokalnie w execute_task/
										//po zakonczeniu SGCMC jest kasowany 
										//chem_list.size() - ile wierszy jest w pliku (punktow do symulacji)
										//chem_list[i].size() - ile typow atomow (pot. chem) jest zadeklarowanych
			//cout<<endl;
			//pozycja pot.chemicznego w vectorze chem odpowiada typowi atomu
			//bierze sie to z tad jak w kolumnach zapisane sa pot.chem w pliku chem.in
			// uAV uBV ...

			double old_chem=chem[old_typ];
			double new_chem=chem[new_typ];
			double mi=new_chem-old_chem;
			
			//cout<<"Zamiana: "<<old_typ<<" "<<new_typ<<endl;
			//cout<<old_chem<<" - "<<new_chem<<" = "<<mi<<endl;

//licz prawdopodobienstwo: exp(-dE/kT) -> dE=E2-E1 -> E = U-u
		    double P1=exp(beta*(mi-(new_E-old_E)));

		//zmien z powrotem na old_typ jesli zdarzenie to nie zostalo trafione rnd()
			if(P1<rnd())
			{
			rnd_site->set_atom(old_typ);
			}
		
		} 
//		else
//		{
		//	cout<<"te same typy"<<endl;
			//m--;
			
//		}		
	}
//	cout<<"Sub-run completed..."<<endl;
}

void widom(lattice *sample,long next, long steps, double T)
{
	double beta=1.0/(kB*T);
	int Ntyp = sample->get_atom_typ_numbers();
	string file_name="widom.dat";
	ofstream output(file_name.c_str(),ios :: app);

	vector <double> p((Ntyp),0.0);
	vector <double> e((Ntyp),0.0);


	direct_exchange(sample,steps,T);

	vector <site *> Aatoms;
	vector <site *> Batoms;
	vector <site *> Vatoms;
	sample->set_atoms_list(Vatoms,0);
	sample->set_atoms_list(Aatoms,1);
	sample->set_atoms_list(Batoms,2);

	long a=0;	
	long b=0;	
	long v=0;	

	for (unsigned int i=0;i<Aatoms.size();i++){
		
		if(Aatoms[i]->get_atom() <= 0 ){
		Aatoms[i]->show_site();
		Aatoms[i]->show_neigh(0);
		Aatoms[i]->show_neigh(1);
		control_output<<"vacancy in the atom list"<<endl;
		exit(1);}
		
		double E1=0.5*(pot.get_energy(Aatoms[i]));
		double E2=0.5*(pot.get_energy(Aatoms[i],0));
		double dE=E2-E1;
		e[1]+=dE;
		p[1]+=exp(-beta*dE);
		a++;
	}

	for (unsigned int i=0;i<Batoms.size();i++){
		
		if(Batoms[i]->get_atom() == 0 ){
		Batoms[i]->show_site();
		Batoms[i]->show_neigh(0);
		Batoms[i]->show_neigh(1);
		control_output<<"vacancy in the atom list"<<endl;
		exit(1);}
		
		double E1=0.5*(pot.get_energy(Batoms[i]));
		double E2=0.5*(pot.get_energy(Batoms[i],0));
		double dE=E2-E1;
		e[1]+=dE;
		p[1]+=exp(-beta*dE);
		b++;
	}

	for (unsigned int i=0;i<Vatoms.size();i++){

		if(Vatoms[i]->get_atom() != 0 ){
		Vatoms[i]->show_site();
		Vatoms[i]->show_neigh(0);
		Vatoms[i]->show_neigh(1);
		control_output<<"atom in the vacancy list"<<endl;
		exit(1);}

		double E=0.5*(pot.get_energy(Vatoms[i]));
		e[0]+=E;
		p[0]+=exp(-beta*E);
		v++;
	}
		

	output<<next<<" "<<T<<" "<<v<<" "<<a<<" "<<b;
	for (unsigned int j=0;j<e.size();j++){ output<<" "<<e[j];}
	for (unsigned int j=0;j<p.size();j++){ output<<" "<<p[j];}
	output<<endl;
	
}

void widom_rnd(lattice *sample,long next, long steps, double T)
{
	double beta=1.0/(kB*T);
	int Ntyp = sample->get_atom_typ_numbers();
	string file_name="widom.dat";
	ofstream output(file_name.c_str(),ios :: app);

	vector <double> Eadd((Ntyp-1),0.0);
	vector <double> Ermv((Ntyp-1),0.0);
	


	direct_exchange(sample,steps,T);

	vector <site *> Aatoms;
	vector <site *> Batoms;
	vector <site *> Vatoms;
	sample->set_atoms_list(Vatoms,0);
	sample->set_atoms_list(Aatoms,1);
	sample->set_atoms_list(Batoms,2);

    

	
	{
	long N=(long)(rnd()*Aatoms.size());
	double E1=pot.get_energy(Aatoms[N]);
	double E2=pot.get_energy(Aatoms[N],0);
	double dE=E2-E1;
	double p=exp(-beta*dE);
		if(Aatoms[N]->get_atom() == 0 ){
		control_output<<N<<" "<<E1<<" "<<E2<<" "<<dE<<endl;
		Vatoms[N]->show_site();
		Vatoms[N]->show_neigh(0);
		Vatoms[N]->show_neigh(1);
		control_output<<"vacancy in the atom list"<<endl;
		exit(1);
//		cout<<"stop"<<endl;
//		cin>>o;
		}

	Ermv[0]=p;
	}
	
	
	{
	    long N=(long)(rnd()*Batoms.size());
		double E1=pot.get_energy(Batoms[N]);
		double E2=pot.get_energy(Batoms[N],0);
		double dE=E2-E1;
		double p=exp(-beta*dE);
		if(Batoms[N]->get_atom() == 0 ){
		control_output<<N<<" "<<E1<<" "<<E2<<" "<<dE<<endl;
		Vatoms[N]->show_site();
		Vatoms[N]->show_neigh(0);
		Vatoms[N]->show_neigh(1);
		control_output<<"vacancy in the atom list"<<endl;
		exit(1);
//		cout<<"stop"<<endl;
//		cin>>o;
		}

		Ermv[1]=p;
		}

	{
	    long N=(long)(rnd()*Vatoms.size());
		double E1=pot.get_energy(Vatoms[N]);
		double E2=pot.get_energy(Vatoms[N],1);
		double dE=E2-E1;
		double p=exp(-beta*dE);
		if(Vatoms[N]->get_atom() != 0 ){
		control_output<<N<<" "<<E1<<" "<<E2<<" "<<dE<<endl;
		Vatoms[N]->show_site();
		Vatoms[N]->show_neigh(0);
		Vatoms[N]->show_neigh(1);
		control_output<<"atom in the vacancy list"<<endl;
		exit(1);
//		cout<<"stop"<<endl;
//		cin>>o;
		}
		Eadd[0]=p;}
		
	if(Batoms.size()>0){
	{
	    long N=(long)(rnd()*Vatoms.size());
		double E1=pot.get_energy(Vatoms[N]);
		double E2=pot.get_energy(Vatoms[N],2);
		double dE=E2-E1;
		double p=exp(-beta*dE);
		if(Vatoms[N]->get_atom() != 0 ){
		control_output<<N<<" "<<E1<<" "<<E2<<" "<<dE<<endl;
		Vatoms[N]->show_site();
		Vatoms[N]->show_neigh(0);
		Vatoms[N]->show_neigh(1);
		control_output<<"atom in the vacancy list"<<endl;
		exit(1);
//		cout<<"stop"<<endl;
//		cin>>o;
		}

		Eadd[1]=p;}
	}
		output<<next<<" "<<T<<" "<<Vatoms.size()<<" "<<Aatoms.size()<<" "<<Batoms.size();
		for (unsigned int j=0;j<Eadd.size();j++){ output<<" "<<Eadd[j];}
		for (unsigned int j=0;j<Ermv.size();j++){ output<<" "<<Ermv[j];}
		output<<endl;
	
}



void direct_exchange(lattice *sample,long steps,double T)
{
	double beta=1.0/(kB*T);	
	long Nsize=sample->get_sim_atom_number();
	for (long i=0; i<steps; i++){
		long N1=(long)(rnd()*Nsize);
		long N2=(long)(rnd()*Nsize);		
		site *site1=0;	
		site *site2=0;
		site1=sample->get_site(N1);
		site2=sample->get_site(N2);
		int typ1 = site1->get_atom();
		int typ2 = site2->get_atom();
		double E1=pot.get_energy(site1)+pot.get_energy(site2)-pot.get_energy(site1,site2);
		site1->set_atom(typ2);
		site2->set_atom(typ1);
		double E2=pot.get_energy(site1)+pot.get_energy(site2)-pot.get_energy(site1,site2);
		double dE=E2-E1;
		if(dE > 0 ){
			double p = exp(-beta*dE);
			double R = rnd();
			if(R >= p){
				site1->set_atom(typ1);
				site2->set_atom(typ2);
			}
		}
	}

}	

double vac_mechanism(lattice *sample,long number_of_steps,long direct_step, double T, int file_nr)
{	
	RTA_energy_executions++; double time=0.0;
	int equi_step=opt_equi.get_steps();			//return <-1 or ==0 or >+1
	if (equi_step > 0){ opt_equi.set_temperature(T);}
	
	//poczatek petli rta for
	for(long n=0;n<number_of_steps;n++)
	{
		if (equi_step > 0 and !((RTA_energy_executions*number_of_steps+n) % equi_step)){
			opt_equi.equilibrate();
		}
		
		//try make a jump for all vacancies in the system
		for(unsigned int tmpV=0;tmpV<Vatoms.size();tmpV++){
			site* vac_to_jump=0;
			site* atom_to_jump=0;
			vac_to_jump=Vatoms[tmpV];
			vector <site*> vac_neighbour;
			vac_to_jump->read_site_neighbours(vac_neighbour,1,0);  //- atomy sasiedzi
			long Vsize = vac_neighbour.size(); 
			long N2=(long)(rnd()*Vsize);	atom_to_jump=vac_neighbour[N2]; 
			
		//	bool NO_AT=false;int counter=0; 
			if(sample->check_site_belonging_to_sim_area(atom_to_jump) ){

				int atom_typ = atom_to_jump->get_atom();
		//		int counter=0; 
				int jump_occured=0;
		//		while( atom_typ==0){
		//			counter++;
		//			N2=(long)(rnd()*Vsize);		
		//			atom_to_jump=vac_neighbour[N2];
		//			atom_typ = atom_to_jump->get_atom();
		//			if(counter>Vsize){break;}
		//		}

			if(atom_typ>0){
				jump_occured=try_jump(sample,T,vac_to_jump,atom_to_jump);
			if(equi_step != 0 and jump_occured){
				opt_equi.call_flux(atom_to_jump,vac_to_jump);
			}}
			}
		
		}
		time++;
	}	//end of for n loop step

	return time;
}

//koniec petli vac_mechanism

int try_jump(lattice* sample, double T, site* vac_to_jump, site* atom_to_jump){

	double beta=1.0/(kB*T);
	int atom = atom_to_jump->get_atom();		//wczytaj typ atomu sasiada atomowego wakancji
	int jump_occured=0;
	unsigned int zone = pot.check_coordination_zone(vac_to_jump,atom_to_jump);	
	double E1=pot.get_energy(vac_to_jump)+pot.get_energy(atom_to_jump)-pot.get_energy(vac_to_jump,atom_to_jump);
	double E2=pot.get_energy(vac_to_jump,atom)+pot.get_energy(atom_to_jump,0)-pot.get_energy(vac_to_jump,0,atom_to_jump,0)-pot.get_energy(vac_to_jump,atom,atom_to_jump,atom) + pot.get_energy(vac_to_jump,atom,atom_to_jump,0);
	double barrier=1000.0;
	barrier=simple_bars[atom][zone];				
	double bariera=(E1+E2)/2+barrier-E1;		// 	--policz bariere

	double P = exp(-beta*(bariera));
	if(rnd() < P){
		make_jump(sample,vac_to_jump, atom_to_jump);
		jump_occured=1;
	}
	return jump_occured;
}

void make_jump(lattice* sample, site* vac_to_jump, site* atom_to_jump){
	
	double latt_constx = sample->get_latice_const(1,0);
	double latt_consty = sample->get_latice_const(2,0);
	double latt_constz = sample->get_latice_const(3,0);
	double boundary_con_x = sample->get_latice_transition(1);
	double boundary_con_y = sample->get_latice_transition(2);
	double boundary_con_z = sample->get_latice_transition(3);	
	int atom = atom_to_jump->get_atom();		//wczytaj typ atomu sasiada atomowego wakancji
	
//	site *adressbufor=NULL;	
//	adressbufor=vac_to_jump;
//wczytaj pozycje poczatkowe atomu i wakancji
	double xjumper=vac_to_jump->get_x();			//xjumper - pozycja na osi x wakancji przed skokiem
	double yjumper=vac_to_jump->get_y();
	double zjumper=vac_to_jump->get_z();
	long int Vindex = vac_to_jump->get_vindex();
//	control_output<<"Vindex: "<<Vindex<<endl;
	double xvac=atom_to_jump->get_x();		//xvac - pozycja na osi x wakancji po przeskoku
	double yvac=atom_to_jump->get_y();
	double zvac=atom_to_jump->get_z();
//policz kierunki skoku
	double xjump = move(xjumper,xvac,latt_constx,boundary_con_x);	//WARUNKI BRZEGOWE 
	double yjump = move(yjumper,yvac,latt_consty,boundary_con_y);	//UWZGLEDNIONE
	double zjump = move(zjumper,zvac,latt_constz,boundary_con_z);	//W move {return (xjumper-xvac)}

	//wykonaj monitor

//buforowanie
	vector <long int> Vjp;
	double Vdx=vac_to_jump->get_drx();
	double Vdy=vac_to_jump->get_dry();
	double Vdz=vac_to_jump->get_drz();
	vac_to_jump->get_jumps(Vjp);//liczba skokow tej wakancji			
	vector <long int> Ajp;
	double Adx=atom_to_jump->get_drx();
	double Ady=atom_to_jump->get_dry();
	double Adz=atom_to_jump->get_drz();				
	atom_to_jump->get_jumps(Ajp);//liczba skokow tej wakancji			
//identyfikacja rodzaju skoku na podstawie r^2
	double r2=sqrt(xjump*xjump + yjump*yjump + zjump*zjump);
//	int strefa=-1;			
	if(r2 < 1.8){		//skoki NN typu 1,1,1 -> sqrt((latt_constx/2.0)^2 + ...) UWAGA - > ustawione na sztywno dla tego przypadku
		//strefa=1;
		//vac		//total sum
		Vjp[0] += 1;
		Vdx += xjump;
		Vdy += yjump;
		Vdz += zjump;
		//NN shel
		Vjp[1] += 1;
		//atom
		Ajp[0] += 1;
		Adx -= xjump;
		Ady -= yjump;
		Adz -= zjump;
		//NN shel
		Ajp[1] += 1;
	}
	else if (r2 > 1.8 and r2 < 2.1 ){	//skoki NNN typu 2,0,0
		//strefa=2;
		//vac		//total sum
		Vjp[0] += 1;
		Vdx += xjump;
		Vdy += yjump;
		Vdz += zjump;
		//NNN shel
		Vjp[2] += 1;
		//atoms
		Ajp[0] += 1;
		Adx -= xjump;
		Ady -= yjump;
		Adz -= zjump;
		//NNN shel
		Ajp[2] += 1;
	}
	else{control_output<<"ERROR: mc::resident_time_energy -> no possible jumps occurs -> problem with coordinates "<<endl; exit(0);
	}
//mam cale wektory skokow vac i atomu zmienione, teraz trzeba je podmienic				//zamieniamy wakancje na atom
	vac_to_jump->set_atom(atom);//do adresu gdzie jest wakancja nadpisac typ na atom 
	vac_to_jump->set_drx(Adx);	//nadpisac liczbe skokow w kierunku x
	vac_to_jump->set_dry(Ady);
	vac_to_jump->set_drz(Adz);
	vac_to_jump->set_jumps(Ajp);
	vac_to_jump->reset_vindex();

	//w tym momencie kontener Vatoms przechowuje adres do atomu (symulacja gubi wakancje)
	//zamieniamy atom na wakacje
	atom_to_jump->set_atom(0);
	atom_to_jump->set_drx(Vdx);
	atom_to_jump->set_dry(Vdy);
	atom_to_jump->set_drz(Vdz);
	atom_to_jump->set_jumps(Vjp);
	atom_to_jump->set_vindex(Vindex);
	
	V_LIST.erase(vac_to_jump);
	V_LIST.insert(atom_to_jump);		
	
	//collect flux data

	opt_equi.call_flux(atom_to_jump,vac_to_jump);	
	
}

double residence_time_energy(lattice *sample,long number_of_steps,double T, int file_nr)
{	
	static long RTA_energy_executions;
	RTA_energy_executions++;
	
//	int o;
	string file_name1;
	string file_name2;  
	string file_name3;  
	stringstream s;
	s<<"";
	s<<file_nr;
	//file_name1=s.str()+"bar_all.dat";
	file_name2=s.str()+"bar_rel.dat";
	file_name3=s.str()+"bar_rel_ene.dat";
			
			
			
	s<<"";
	//ofstream barriers_all(file_name1.c_str(),ios :: app);
	ofstream barriers_jump(file_name2.c_str(),ios :: app);
	ofstream energy_jump(file_name3.c_str(),ios :: app);
	//wczytaj parametry sieci krystalicznej
	double latt_constx = sample->get_latice_const(1,0);
	double latt_consty = sample->get_latice_const(2,0);
	double latt_constz = sample->get_latice_const(3,0);
	
	double boundary_con_x = sample->get_latice_transition(1);
	double boundary_con_y = sample->get_latice_transition(2);
	double boundary_con_z = sample->get_latice_transition(3);	
	
	//definicja zmiennych
	
	double time=0.0;
	int nr_atoms_type=sample->get_atom_typ_numbers();
	int nr_lattice_type=sample->get_vec_lattice_typ_size();
	int Z = sample->get_max_coordination_number();
	//cout<<nr_atoms_type<<" "<<nr_lattice_type<<endl;

	//deklaracja tablicy 8-Z liczba koordynacyjna w I strefie (przyblzenie najblizszychsasiadow)
	double *target = new double[Z*Vatoms.size()+10];
	//vector <vector <wektor> > atoms_barriers(nr_atoms_type,vector <wektor> (nr_lattice_type,wektor()) );
	vector < vector <vector <wektor> > > atoms_barriers_jump(nr_atoms_type, vector < vector <wektor > > (nr_lattice_type, vector <wektor> (nr_lattice_type, wektor()) ) );
	vector <vector <vector <wektor> > > atoms_energy_jump(nr_atoms_type, vector < vector <wektor > > (nr_lattice_type, vector <wektor> (nr_lattice_type, wektor()) ) );
	//vector <> atoms_barriers_counts((nr_atoms_type*nr_lattice_type),0.0);
	
	int equi_step=0;
	if (opt_equi.get_steps()){ 
		equi_step=opt_equi.get_steps();
		opt_equi.set_temperature(T);
		}
	
	//poczatek petli rta for
	for(long n=0;n<number_of_steps;n++)
	{
		//rownowazenie wakancji
		//cout<<(RTA_energy_executions*number_of_steps+n)<<" "<<((RTA_energy_executions*number_of_steps+n) % equi_step)<<" "<<(!((RTA_energy_executions*number_of_steps+n) % equi_step))<<endl;
		if (equi_step){
		if(!((RTA_energy_executions*number_of_steps+n) % equi_step))
		{
			control_output<<"rownowaze... "<<(RTA_energy_executions*number_of_steps+n)<<endl;
			opt_equi.equilibrate();
			
		}}
		//UWAGA NA WAKANCJE TRZEBA ODSWIEZYC po rownowazeniu liste Vatoms
	
		//deklaracje kontenerow vector
		vector <long> vac_to_jump;
		vector <site*> atom_to_jump;
		
		//definicja zmiennych
		long itarget=0;
//		long Vsize=0;
		
		//inicjalizacja zmiennych
		for(unsigned long i=0;i<=(Z*Vatoms.size()+9);i++)
		target[i]=0.0;
		
		unsigned int i,k;
		double beta=1.0/(kB*T);
		vector <pairjump> tablica_skokow;
		tablica_skokow.clear();
		tablica_skokow.reserve(Z*Vatoms.size());
		
	

#pragma omp parallel default(none) private(i) shared(sample,cout,control_output,Vatoms,tablica_skokow,beta,simple_bars,bar_list,pot)
		{
	//	numit=0;
		i=0;
	//	id = omp_get_thread_num();
	//	printf ("Thread no %d starting... \n", id);
		

		//poczatek petli kroku rta (dla kazdej wakancji)
#pragma omp for schedule(runtime) private(i,k)							//shared(vac_to_jump,atom_to_jump,itarget,target)  private(i)
		for(i=0;i<Vatoms.size();i++)
		{	
			
		//	cout<<id<<" ";
		//	numit++;
		//	id = omp_get_thread_num();
	//		tablica.push_back(id);
			//deklaracje i inicjalizacja kontenerow do przechowywania sasiadow wakancji
			
			vector <site*> vac_neighbour;
			//vector <site*> vac_neighbour_en;
			Vatoms[i]->read_site_neighbours(vac_neighbour,1,0);  //- atomy sasiedzi
			//Vatoms[i]->read_site_neighbours(vac_neighbour_en,0); // - energia sasiedzi
			
			//poczatek petli dla kazdego atomu sasiada wakancji
			for(k =0;k<vac_neighbour.size();k++)			// dla kazdego sasiada atomowego
			{
			//	control_output<<"jestem1 "<<i<<" "<<k<<endl;
			//	vac_neighbour[k]->show_site();				
			//	int nr_lattice = vac_neighbour[k]->get_sub_latt();
			//	control_output<<"jestem2 "<<i<<" "<<k<<endl;
				int atom;
				double E1,E2;	
				atom = vac_neighbour[k]->get_atom();		//wczytaj typ atomu sasiada atomowego wakancji
//				int latA=vac_neighbour[k]->get_sub_latt();
//				int latV=Vatoms[i]->get_sub_latt();

		
			//	control_output<<"jestem3 "<<i<<" "<<k<<endl;
//#pragma omp critical(switching_atoms)
//{
			//	Vatoms[i]->set_atom(atom);
			//	vac_neighbour[k]->set_atom(0);
				E1=pot.get_energy(Vatoms[i])+pot.get_energy(vac_neighbour[k])-pot.get_energy(Vatoms[i],vac_neighbour[k]);
				E2=pot.get_energy(Vatoms[i],atom)+pot.get_energy(vac_neighbour[k],0)-pot.get_energy(Vatoms[i],0,vac_neighbour[k],0)-pot.get_energy(Vatoms[i],atom,vac_neighbour[k],atom) + pot.get_energy(Vatoms[i],atom,vac_neighbour[k],0);

			//	Vatoms[i]->set_atom(0);
			//	vac_neighbour[k]->set_atom(atom);
//}				
				//licz bariere przeskoku
			//	int o;
			//	cout<<"energie:"<<E1<<" "<<E2<<" dla i,k: "<<i<<k<<endl;
				//cin >> o;

				
				double barrier=20000.0;//sample->get_barrier(*I,*K);
				
			//2D simple barriers
				unsigned int zone = pot.check_coordination_zone(Vatoms[i],vac_neighbour[k]);	
				barrier=simple_bars[atom][zone];
		
		//3D full barriers
//			vector <int> okno(6,0);
	
			//for(int o=0; o<6;o++){
			//   okno[o]=0;
			//}
			
//			vector <site> window;
//			window.reserve(10);
//			window.clear();

//			sample->get_window(vac_neighbour[k],Vatoms[i],window);
//			sample->sort_atoms(window);
			
//			if(okno.size() < window.size()){ cout<<"ERROR: in RTA"; exit(0);}
				
	//		control_output <<"writting window..."<<endl;		
			
//			for(int o=0; o<window.size();o++){
//				okno[o]=(window[o].get_atom());
//			}

//			control_output <<"OKNO:"<<endl;		
//			for(int o=0; o<okno.size();o++){
//				control_output<<o<<" "<<okno[o]<<endl;
//			}

			
			
//			barrier=bar_list[atom][latA][latV][(okno[0])][(okno[1])][(okno[2])][(okno[3])][(okno[4])][(okno[5])];
		
							
				
				double bariera=(E1+E2)/2+barrier-E1;		// 	--policz bariere

				//vac_neighbour[k]->add_atom barrier(EBS)	--dodaj barriere do sita
		//		int miejsce_w_tablicy=atom;
		//		atoms_barriers[atom][nr_lattice].x=atoms_barriers[atom][nr_lattice].x+bariera;
		//		atoms_barriers[atom][nr_lattice].y=atoms_barriers[atom][nr_lattice].y+1;
				
	//			itarget++;		//krok zdarzenia +1
				pairjump tmp(Vatoms[i],vac_neighbour[k],E1,E2,barrier,bariera);
		//		tmp.show();
#pragma omp critical(adding_jump)
				tablica_skokow.push_back(tmp);
			//	target[itarget]=target[itarget-1]+exp(-beta*((E1+E2)/2.0+barrier-E1));//prawdopodobienstwo kolejnego zdarzenia zapisz w tablicy
			//	vac_to_jump.push_back(i);//dla kolejnego zdarzenia numer wakancji (ktory odpowiada nr z listy wakancji zapisz do kontenera
			//	atom_to_jump.push_back(vac_neighbour[k]);//adres atomu kolejnego zdarzenia zapisz do kontenera	
			
			
			}//koniec petli dla kazdego atomu sasiada wakancji 
		}//koniec petli dla kazdej wakancji
		
	//	printf("Thread %d performed %d iterations\n",id,numit);
	}
//	control_output<<"Rozmiar tab_skokow: "<<tablica_skokow.size()<<endl;
	for(unsigned long i=0;i<tablica_skokow.size();i++)
	{
		//control_output<<i<<" ";
		//tablica_skokow[i].show();
		itarget++;
		double bar=tablica_skokow[i].get_barier();
		target[itarget]=target[itarget-1]+exp(-beta*(bar));//prawdopodobienstwo kolejnego zdarzenia zapisz w tablicy
    
//		control_output<<itarget<<" "<<bar<<" "<<target[itarget]<<endl;

	}
//	control_output<<endl;
		
	
	
	//	control_output<<"jestem4 "<<endl;
	//	cin>>o;
		// sprawdzilem wszystkie przeskoki i zapisalem typ przeskoku i jego prawdopodobienstwo
		time=time-log(rnd())/target[itarget];//policz krok czasowy zdarzenia
		double R=rnd()*target[itarget];//wylosuj zdarzenie z tablicy zdarzen target)
		for(long m=0;m<=itarget-1;m++)			// szukam ktory przeskok bedzie zrealizowany skanujac wszystkie taget
		{//poczatek petli target
			if(R>target[m] && R<=target[m+1])
			{//poczatek if target
				//cin>>o; 
//				control_output<<"Wyknuje skok:"<<endl;
//				tablica_skokow[m].show();
				site* vac_to_jump=tablica_skokow[m].get_vac_to_jump();
				site* atom_to_jump=tablica_skokow[m].get_atom_to_jump();
				//long vac_nr=vac_to_jump[m];//wczytaj nr wakancji do przeskoku
				int ATOM=atom_to_jump->get_atom();
				int NR_LATa = atom_to_jump->get_sub_latt();
				int NR_LATv = vac_to_jump->get_sub_latt();
				//vector <site*> jumper_neighbour_en;				//kontener do przechowywania sasiadow energii atomu
				//vector <site*> vac_neighbour_en;
				//atom_to_jump[m]->read_site_neighbours(jumper_neighbour_en,0); // - 0 enegia
				//Vatoms[vac_nr]->read_site_neighbours(vac_neighbour_en,0); 
				
			//	double E1=sample->get_energy(ATOM,jumper_neighbour_en)+sample->get_energy(0,vac_neighbour_en)-sample->getV(0,ATOM);
				//double E2=sample->get_energy(0,jumper_neighbour_en)+sample->get_energy(ATOM,vac_neighbour_en)+sample->getV(0,ATOM)-sample->getV(ATOM,ATOM)-sample->getV(0,0);
		//omp section		
				double E1=0;
				double E2=0;			
				E1=pot.get_energy(vac_to_jump)+pot.get_energy(atom_to_jump)-pot.get_energy(vac_to_jump,atom_to_jump);
				E2=pot.get_energy(vac_to_jump,ATOM)+pot.get_energy(atom_to_jump,0)-pot.get_energy(vac_to_jump,0,atom_to_jump,0)-pot.get_energy(vac_to_jump,ATOM,atom_to_jump,ATOM) + pot.get_energy(vac_to_jump,ATOM,atom_to_jump,0);

								
			//	double     barrier=E1;//sample->get_barrier(*I,*K);
				
	//			unsigned int zone = pot.check_coordination_zone(vac_to_jump,atom_to_jump);
				//pobierz bariere 2 1 0 1 1 1 2 2 2
			//	barrier=bar_list[2][1][0][1][1][1][2][2][2];
//				barrier=bar_list[1][0][0][1][1][1][1][1][2];
			//	barrier=barr[ATOM][zone];
			//	switch(ATOM)
			//	{
			//		case 1 : barrier=barr1;break;
			//		case 2 : barrier=barr2;break;
			//		case 0 : barrier=1000000;break;
			//	} 
				
				double bariera = log(target[m+1]-target[m])/(-beta);
				double barrier = bariera + E1 - (E1+E2)/2 ;
				atoms_barriers_jump[ATOM][NR_LATa][NR_LATv].x=atoms_barriers_jump[ATOM][NR_LATa][NR_LATv].x+bariera;
				atoms_barriers_jump[ATOM][NR_LATa][NR_LATv].y=atoms_barriers_jump[ATOM][NR_LATa][NR_LATv].y+1;
				atoms_energy_jump[ATOM][NR_LATa][NR_LATv].x=atoms_energy_jump[ATOM][NR_LATa][NR_LATv].x+E1;
				atoms_energy_jump[ATOM][NR_LATa][NR_LATv].y=atoms_energy_jump[ATOM][NR_LATa][NR_LATv].y+E2;			
				atoms_energy_jump[ATOM][NR_LATa][NR_LATv].z=atoms_energy_jump[ATOM][NR_LATa][NR_LATv].z+barrier;			
				
		//		cout<<"e1:"<< E1<<" e2:"<<E2<<" w:"<<barrier<<" b:"<<bariera<<endl;
				//energie dzialaja
//				site *adressbufor=NULL;	//definiuje wskaznik do przechowywania adresu tymczasowo
	//			adressbufor=vac_to_jump;//zapisz adres komorki wakancji z listy wakancji, zostanie on zastapiony niwym adresem atomu
				
//wczytaj pozycje poczatkowe atomu i wakancji
				double xjumper=vac_to_jump->get_x();			//xjumper - pozycja na osi x wakancji przed skokiem
				double yjumper=vac_to_jump->get_y();
				double zjumper=vac_to_jump->get_z();
				long int Vindex = vac_to_jump->get_vindex();
							
				double xvac=atom_to_jump->get_x();		//xvac - pozycja na osi x wakancji po przeskoku
				double yvac=atom_to_jump->get_y();
				double zvac=atom_to_jump->get_z();

//policz kierunki skoku
				double xjump = move(xjumper,xvac,latt_constx,boundary_con_x);	//WARUNKI BRZEGOWE 
				double yjump = move(yjumper,yvac,latt_consty,boundary_con_y);	//UWZGLEDNIONE
				double zjump = move(zjumper,zvac,latt_constz,boundary_con_z);	//W move {return (xjumper-xvac)}
//buforowanie
				vector <long int> Vjp;
				double Vdx=vac_to_jump->get_drx();
				double Vdy=vac_to_jump->get_dry();
				double Vdz=vac_to_jump->get_drz();

				vac_to_jump->get_jumps(Vjp);//liczba skokow tej wakancji
				
				vector <long int> Ajp;
				double Adx=atom_to_jump->get_drx();
				double Ady=atom_to_jump->get_dry();
				double Adz=atom_to_jump->get_drz();				
				atom_to_jump->get_jumps(Ajp);//liczba skokow tej wakancji
				
//identyfikacja rodzaju skoku na podstawie r^2
				double r2=sqrt(xjump*xjump + yjump*yjump + zjump*zjump);
//				int strefa=-1;
				
				if(r2 < 1.8){		//skoki NN typu 1,1,1 -> sqrt((latt_constx/2.0)^2 + ...) UWAGA - > ustawione na sztywno dla tego przypadku
//					strefa=1;
			//vac		//total sum
					Vjp[0] += 1;
					Vdx += xjump;
					Vdy += yjump;
					Vdz += zjump;
					//NN shel
					Vjp[1] += 1;
			//atom
					Ajp[0] += 1;
					Adx -= xjump;
					Ady -= yjump;
					Adz -= zjump;
					//NN shel
					Ajp[1] += 1;
				}
				else if (r2 > 1.8 and r2 < 2.1 ){	//skoki NNN typu 2,0,0
//					strefa=2;
			//vac		//total sum
					Vjp[0] += 1;
					Vdx += xjump;
					Vdy += yjump;
					Vdz += zjump;
					//NNN shel
					Vjp[2] += 1;
			//atoms
					Ajp[0] += 1;
					Adx -= xjump;
					Ady -= yjump;
					Adz -= zjump;
					//NNN shel
					Ajp[2] += 1;
				}
				else{control_output<<"ERROR: mc::resident_time_energy -> no possible jumps occurs -> problem with coordinates "<<endl; exit(0);}

//mam cale wektory skokow vac i atomu zmienione, teraz trzeba je podmienic				//zamieniamy wakancje na atom
				vac_to_jump->set_atom(ATOM);//do adresu gdzie jest wakancja nadpisac typ na atom 
				vac_to_jump->set_drx(Adx);	//nadpisac liczbe skokow w kierunku x
				vac_to_jump->set_dry(Ady);
				vac_to_jump->set_drz(Adz);
				vac_to_jump->set_jumps(Ajp);
				vac_to_jump->set_vindex(0);

				//w tym momencie kontener Vatoms przechowuje adres do atomu (symulacja gubi wakancje)
				//zamieniamy atom na wakacje
				atom_to_jump->set_atom(0);
				atom_to_jump->set_drx(Vdx);
				atom_to_jump->set_dry(Vdy);
				atom_to_jump->set_drz(Vdz);
				atom_to_jump->set_jumps(Vjp);
				atom_to_jump->set_vindex(Vindex);
				
		//		control_output<<"vac przed "<<Vatoms[Vindex]<<" "<<Vatoms[Vindex]->get_vindex()<<endl;
			//	Vatoms[Vindex]->show_site();
//aktualny adres wakancji w liscie wakancji
				Vatoms[Vindex]=atom_to_jump;

		//		control_output<<"vac po "<<Vatoms[Vindex]<<" "<<Vatoms[Vindex]->get_vindex()<<endl;
		//		Vatoms[Vindex]->show_site();

//				for(long v=0;v<Vatoms.size();v++){
	//				if(Vatoms[v]==vac_to_jump){	Vatoms[v]=atom_to_jump;	//nadpisuje adres w Vatoms na aders wakancji}}



				//w tym momencie symulacja odzyskuje wakancje
			//nadpisalem tylko to co sie zmienilo, podmienilem adresy i reszte parametrow sita pozostawiam bez zmian
			}//koniec if target
		}//koniec petli target		
	}// koniec petli rta
	
//	barriers_all<<time;
	barriers_jump<<time;
	energy_jump<<time;
	for(unsigned int ib=0;ib<atoms_barriers_jump.size();ib++)
	{
		for(unsigned int jb=0;jb<atoms_barriers_jump[ib].size();jb++)
		{
		for(unsigned int kb=0;kb<atoms_barriers_jump[ib][jb].size();kb++)
		{

			
			// barriers_all<<" "<<atoms_barriers[ib][jb].x<<" "<<atoms_barriers[ib][jb].y;
			 barriers_jump<<" "<<atoms_barriers_jump[ib][jb][kb].x<<" "<<atoms_barriers_jump[ib][jb][kb].y;
			 energy_jump<<" "<<atoms_energy_jump[ib][jb][kb].x<<" "<<atoms_energy_jump[ib][jb][kb].y<<" "<<atoms_energy_jump[ib][jb][kb].z;
			 
		}}
	}
	
//	barriers_all<<endl;
	barriers_jump<<endl;
	energy_jump<<endl;
//	control_output<<"koniec RTA"<<endl;
	delete(target);//zwalniamy pamiec 
	
	return time/100.0;//zwracamy wartosc
}//koniec petli rta_energy	


/******************************koniec residence_time_energy**********************************/	

double time_increment(double norma){
	double sid = rnd();
	double sidL = log(sid); 
	double dt = sidL/norma; 
	while (sid == 0.0){
		control_output<<"TIME ERROR in mc:resident_time(); "<<endl;
		control_output<<sid<<" "<<sidL<<" "<<norma<<" "<<dt<<endl;
	
		sid = rnd();
		sidL = log(sid); 
		dt = sidL/norma;
	}	
	return -dt;
}	

/*****************************residence_time*************************************************/	
	
double residence_time(lattice *sample,long number_of_steps, double T, int file_nr)
{	
	RTA_energy_executions++;
	static int warrinig_jump_event=0;
	double time=0.0;
//	int Z = sample->get_max_coordination_number();
	double beta=1.0/(kB*T);
	int equi_step=opt_equi.get_steps();			//return <-1 or ==0 or >+1
	if (equi_step > 0){ 
		opt_equi.set_temperature(T);
		}
//	control_output<<"RTA start: "<<RTA_energy_executions<<endl;
//	control_output<<"EVENTS size: "<<EVENTS.size()<<endl;

	if(EVENTS.empty()){
		sample->init_events_list(V_LIST);
		if(EVENTS_size!=EVENTS.size()){control_output<<"Warrning! EVENTS changed after init from:"<<EVENTS_size<<" to: "<<EVENTS.size();
		control_output<<" in step: "<<RTA_energy_executions<<endl;
					EVENTS_size=EVENTS.size();}
	}
	
	for(long n=0;n<number_of_steps;n++){

		if (equi_step > 0){
		if(!((RTA_energy_executions*number_of_steps+n) % equi_step))
		{
		//	control_output<<"rownowaze... "<<(RTA_energy_executions*number_of_steps+n)<<endl;
			opt_equi.equilibrate();
			if(EVENTS_size!=EVENTS.size()){control_output<<"Warrning! EVENTS changed after equilibrate from:"<<EVENTS_size<<" to: "<<EVENTS.size();
				control_output<<" in step: "<<RTA_energy_executions<<endl;
				EVENTS_size=EVENTS.size();
			}
		}}

	double SUM=0;	
	typedef vector <pair <double,pairjump> > mykey;
	mykey target;
	if(!EVENTS.empty()){
		target.push_back(make_pair(SUM, EVENTS.front()));
		double bar= EVENTS.front().get_barier();
		SUM = SUM + exp(-beta*(bar));
	}
	
	if(!EVENTS.empty()){
		list<pairjump>::iterator it = EVENTS.begin(); ++it;
		for(; it != EVENTS.end(); ++it){	
			target.push_back(make_pair( SUM, (*it) ));
			double bar= (*it).get_barier();
			SUM = SUM + exp(-beta*(bar));
		}
		target.push_back(make_pair( SUM, EVENTS.back() ));
	}

	mykey::iterator event;
//	control_output<<"Print target: "<<target.size()<<endl;
//	int counter=0;
//	for(event=target.begin(); event != target.end(); ++event){
//		control_output<<counter<<" p: ";
//		control_output<<(*event).first<<endl;
//		(*event).second.show();
//		counter++;
//	}	
	
	site* vac_to_jump=0;
	site* atom_to_jump=0;
	double R=rnd()*SUM; 
	mykey::iterator next_event=target.begin();
	event = target.begin();

//	control_output<<"Schoot target at: "<<R<<endl;
//	control_output<<"Searching event: "<<target.size()<<endl;
	
	if(!target.empty()){	
		for( ++next_event ; next_event != target.end(); ++event, ++next_event){	
			double Lvalue = (*event).first;
			double Rvalue = (*next_event).first;	
//	control_output<<Lvalue<<" "<<R<<" "<<Rvalue<<endl;
			if( R>=Lvalue and R < Rvalue){
//	control_output<<"Find event: "<<Lvalue<<" "<<R<<" "<<Rvalue<<endl;
//	(*event).second.show();
				vac_to_jump=(*event).second.get_vac_to_jump();
				atom_to_jump=(*event).second.get_atom_to_jump();
				make_jump(sample,vac_to_jump,atom_to_jump);		//zaminia miejscami typy
				sample->update_events(vac_to_jump);	
				sample->update_events(atom_to_jump);
				if( (warrinig_jump_event <= 4) and (EVENTS_size!=EVENTS.size()) ){control_output<<"Warrning! EVENTS changed after jump from:"<<EVENTS_size<<" to: "<<EVENTS.size();
					control_output<<" in step: "<<RTA_energy_executions<<" "<<warrinig_jump_event<<endl;
					EVENTS_size=EVENTS.size();
					(*event).second.show();
					warrinig_jump_event++;
				}
				break;
			}
		}																//znalazlem event-> zamienilem atomy miejscami -> policzylem monitory-> uaktualnilem liste eventow
		if( event==target.end() ){
			control_output<<"ERROR in mc::resident(). Problem with event list. Probably inf or nan in calculation of probability. No event was choosen: "<<endl;
			control_output<<"Find event: "<<(*event).first<<" "<<R<<" "<<(*next_event).first<<endl;
			(*event).second.show();
			exit(1);
		}
		double dt =time_increment(SUM);									//increment time
		time += dt;
		opt_equi.add_MCtime(dt);
	}
	}																	// koniec petli for sub step

	return time;
}

//koniec petli resident_time


double RTA_random_alloy(lattice *sample,long number_of_steps,double T, double w1, double w2)
{
//int piec;
		
		double latt_constx = sample->get_latice_const(1,0);
		double latt_consty = sample->get_latice_const(2,0);
		double latt_constz = sample->get_latice_const(3,0);
		
		double boundary_con_x = sample->get_latice_transition(1);
		double boundary_con_y = sample->get_latice_transition(2);
		double boundary_con_z = sample->get_latice_transition(3);

//double beta=1.0/(kB*T);
double *target = new double[8*Vatoms.size()+10];
double time=0.0;
long sasiady[9];
long generator[10];

for (int i=0;i<9;i++)
	sasiady[i]=0;

for (int i=0;i<10;i++)
	generator[i]=0;
	

for(long n=0;n<number_of_steps;n++)
{
	vector <long> vac_to_jump;
	vector <site*> atom_to_jump;
	long itarget=0;
	for(unsigned long i=0;i<=(8*Vatoms.size()+9);i++)
	{target[i]=0.0;}
	
	//set_atoms_list(sample,0);
	//control_output<<"V siz "<<Vatoms.size()<<endl;
	
	vector <long> jumps_sequence;
	
	for(unsigned int i =0;i<Vatoms.size();i++)
	{		
	//	control_output<<"New vacancy nr: "<<endl;
		//control_output<<i<<" adres: "<<Vatoms[i]<<" ";
		//Vatoms[i]->show_site();					// wczytaj wybrana wakancje i jej sasiadow at i en
		//int atom1=Vatoms[i]->get_atom();
		vector <site*> vac_neighbour;
		vector <site*> vac_neighbour_en;
		Vatoms[i]->read_site_neighbours(vac_neighbour,1,0);  //- atomy sasiedzi
		//Vatoms[i]->read_site_neighbours(vac_neighbour_en,0); // - energia sasiedzie
												
		int nr_sasiadow=vac_neighbour.size();
	//	control_output<<"vac_neig at siz "<<vac_neighbour.size()<<endl;
		
		sasiady[nr_sasiadow]++;
		
	//	for(int j =0;j<vac_neighbour.size();j++)			// pokaz sasiadow atomowych
	//	{
	//		control_output<<j<<" adres: "<<vac_neighbour[j]<<" ";
	//		vac_neighbour[j]->show_site();			
	//	}
		
	//	control_output<<"vac_neig en siz "<<vac_neighbour_en.size()<<endl;
	//	for(int j =0;j<vac_neighbour_en.size();j++){
	//			control_output<<j<<" adres: "<<vac_neighbour_en[j]<<" ";
	//			vac_neighbour_en[j]->show_site();
	//		}										// pokaz sasiadow energii
		
		
	/*	
		for(int k =0;k<5;k++)
		{
		//	int o;
			
		//	cout<<"Kbeg "<<k<<endl;	
			long stop = jumps_sequence.size();
			int random_neigh=int(ran01()*vac_neighbour.size());
		//	cout<<random_neigh<<" from: "<<vac_neighbour.size()<<endl;
			
			if(jumps_sequence.size()==0)
			{jumps_sequence.push_back(random_neigh);
		//	cout<<"added0 "<<jumps_sequence.size()<<endl;
			
			}
			else
			{
				int dodaj=0;
				for(long rn=0;rn<stop;rn++)
				{
					
					if(jumps_sequence[rn]==random_neigh)
					{k--;
			//		cout<<random_neigh<<" canceled: "<<k<<endl;
					dodaj=0;
					break;
					}
					else
					{	
						dodaj=1;
					}
				}
				
			if(dodaj)
				{jumps_sequence.push_back(random_neigh);
		//		cout<<"k "<<k<<" added ";
		//		cout<<random_neigh<<" siz "<<jumps_sequence.size()<<endl;
						
			}	
			}
		//	cout<<"Kend "<<k<<endl;
		//	cin>>o;	
		}
			
		for(int kk =0;kk<vac_neighbour.size();kk++)
		{
		//	int o;
			long stop = jumps_sequence.size();
			int dodaj=0;
				for(long rn=0;rn<stop;rn++)
				{
					
					if(jumps_sequence[rn]==kk)
					{
		//			cout<<" canceled: "<<kk<<endl;
					dodaj=0;
					break;
					}
					else
					{	
						dodaj=1;
					}
				}
				
			if(dodaj)
				{jumps_sequence.push_back(kk);
			//	cout<<"k "<<kk<<" added ";
			//	cout<<kk<<" siz "<<jumps_sequence.size()<<endl;
						
			}	
	//	cout<<"Kend "<<kk<<endl;
//			cin>>o;		
		}
		
	*/	
			
		//	control_output<<"Jumps sequence: "<<endl;
		
	//	cout<<" ";   jumps_sequence
		for(unsigned int k =0;k<vac_neighbour.size();k++)			// dla kazdego sasiada atomowego
		{
		
	//		cout<<jumps_sequence[k];		jumps_sequence[k]
		
			vector <site*> jumper_neighbour_en;
			int atom = vac_neighbour[k]->get_atom();		// wczytaj typ atomu sasiada atomwego wakancji
		//	control_output<<"Atom to jump: "<<vac_neighbour[k]<<" ";
			//vac_neighbour[k]->show_site();
			
		//	cout<<"jumper_neig en siz bef "<<jumper_neighbour_en.size()<<endl;
						// wczytaj sasiadow enregii sasiada atomowego
			//vac_neighbour[k]->read_site_neighbours(jumper_neighbour_en,0); // - 0 enegia
			
			
				
	//			control_output<<"jumper_neig en siz: "<<jumper_neighbour_en.size()<<endl;
			
		//	for(int l =0;l<jumper_neighbour_en.size();l++){
		//		control_output<<l<<" adres: "<<jumper_neighbour_en[l]<<" ";
	//			jumper_neighbour_en[l]->show_site();
	//		}									// wypisz sasiadow energii sasiada atomwego wakancji
			
			
		//	cout<<endl;
			//cout<<"Checking energy calculus..."<<endl;
		//	double E1=sample->get_energy(atom,jumper_neighbour_en)+sample->get_energy(0,vac_neighbour_en)-sample->getV(0,atom);
			//cout<<"Eat: "<<sample->get_energy(atom,jumper_neighbour_en);
		//	cout<<" Evac: "<<sample->get_energy(0,vac_neighbour_en)<<" Eat-v: "<<sample->getV(0,atom)<<" E1: "<<E1<<endl;
			//double E2=sample->get_energy(0,jumper_neighbour_en)+sample->get_energy(atom,vac_neighbour_en)+sample->getV(0,atom)-sample->getV(atom,atom)-sample->getV(0,0);
			//cout<<"Eat: "<<sample->get_energy(atom,vac_neighbour_en);
		//	cout<<" Evac: "<<sample->get_energy(0,jumper_neighbour_en)<<" Eat-v: "<<sample->getV(0,atom)<<" 00 "<<sample->getV(0,0)<<" atat "<<sample->getV(atom,atom)<<" E2: "<<E2<<endl;
			
	
			//SPRAWDZ CZY DOBRZE LICZY ENERGIE!!! Sprawdzone
			
			double     barrier=0.0;//sample->get_barrier(*I,*K);
			//int o;
			//cin>>o;
		
			switch(atom)
			{
				case 1 : barrier=w1;break;
				case 2 : barrier=w2;break;
				case 0 : barrier=0.0;break;
			
				
			}  
			    
			 
			
			itarget++;
			target[itarget]=set_prec(target[itarget-1]+barrier)	;	//exp(-beta*((E1+E2)/2.0+barrier-E1));
		//	control_output<<itarget<<" "<<target[itarget]<<" "<<barrier<<endl;
			vac_to_jump.push_back(i);
			atom_to_jump.push_back(vac_neighbour[k]);
		}
	}
	// sprawdzilem wszystkie przeskoki i zapisalem typ przeskoku i jego prawdopodobienstwo
	//control_output<<"Zawartosc talic przeskoku:"<<endl;
	//control_output<<"vac_to_jump siz "<<vac_to_jump.size()<<endl;
		//for(int j =0;j<vac_to_jump.size();j++)
	//	{
	//			control_output<<j<<" "<<vac_to_jump[j]<<" adres: "<<Vatoms[vac_to_jump[j]]<<" ";
	//			Vatoms[vac_to_jump[j]]->show_site();
	//	}
			
	//control_output<<"atom_to_jump siz "<<atom_to_jump.size()<<endl;
	//	for(int j =0;j<atom_to_jump.size();j++)
	//	{
	//			control_output<<j<<" adres: "<<atom_to_jump[j]<<" ";
	//			atom_to_jump[j]->show_site();
	//	}
	
	
	time=time-log(rnd())/target[itarget];
	
	double random_number=rnd(); //	rnd();rownomierny(0.0,1.0)
	double R=random_number*target[itarget];
	//control_output<<"Random number: "<<random_number<<" gives jumps at: "<<R<<"from range 0 to "<<target[itarget]<<endl;
	
	//generator[int(random_number*10)]++;
	
	for(long m=0;m<=itarget-1;m++)			// szukam ktory przeskok bedzie zrealizowany skanujac wszystkie taget
	{
		if(R>target[m] && R<=target[m+1])
		{
			//control_output<<"take jump nr: "<<m<<" target "<<target[m]<<endl;			
			generator[m]++;
			long vac_nr=vac_to_jump[m];
	//		site *adressbufor=0;	
			
			
			// przeliczyc energie czy dobrze liczy
			//teraz trzeba pozmieniac wartosci zmiennych tych sitow na ktore pokazuje wskazniki
			//bufor pokazuje na nowe pozycje atomu (są pod nim teraz zapisane dane wakancji)
			//Vatoms[m] pokazuje na nowa pozycje wakancji (są pod nim teraz zapisane dane atomu)
			
			//tutaj dodane zostaną nowe wartosci przeskoków do dyfuzji
			
			
			//zapamietuje adres oraz parametry wakancji poniewaz pozniej nadpisze ten 
			//wskaznik innym adresem oraz parametrami atomu.
	//		control_output<<endl;
	//		control_output<<endl;
	//		control_output<<"Skok - adres bufor:";
//			adressbufor=Vatoms[vac_nr];
			//long x0=Vatoms[vac_nr]->get_x0();
			//long y0=Vatoms[vac_nr]->get_y0();
			//long z0=Vatoms[vac_nr]->get_z0();
			
			//cout<<adressbufor<<" "<<x0<<" "<<y0<<" "<<z0<<" "<<nr_jumps<<endl;
	
			// wykonuje operacje przeskoku

		//	control_output<<"Atom przed operacja: "<<endl;
		//		control_output<<"adres: "<<atom_to_jump[m]<<" ";
		//		atom_to_jump[m]->show_site();
		//		atom_to_jump[m]->show_site_jumps();
		//		atom_to_jump[m]->show_neigh(1);
		//		atom_to_jump[m]->show_neigh(0);
			
		//	control_output<<"Vacancy przed operacja: "<<endl;
		//		control_output<<"adres: "<<Vatoms[vac_nr]<<" ";
		//		Vatoms[vac_nr]->show_site();
		//		Vatoms[vac_nr]->show_site_jumps();
		//		Vatoms[vac_nr]->show_neigh(1);
		//		Vatoms[vac_nr]->show_neigh(0);
				
		//	control_output<<"Dokonuje przeskoku"<<endl;
		
			long nr_jumps=Vatoms[vac_nr]->get_jumps();
	//		double r2atom=0.0;	//przechowuje kieunki skoku atomu
	//		double r2vac=0.0;
			
			double drxat=atom_to_jump[m]->get_drx();			//DODAC WARUNKI BRZEGOWE
			double xjumper=Vatoms[vac_nr]->get_x();
			double dryat=atom_to_jump[m]->get_dry();
			double yjumper=Vatoms[vac_nr]->get_y();
			double drzat=atom_to_jump[m]->get_drz();
			double zjumper=Vatoms[vac_nr]->get_z();
			//double r2x=(x-x0at)*(x-x0at);
			//double r2y=(y-y0at)*(y-y0at);
			//double r2z=(z-z0at)*(z-z0at);
			//r2atom=atom_to_jump[m]->get_r2()+r2x+r2y+r2z;
		//	cout<<"r0: "<<x0at<<" "<<y0at<<" "<<z0at<<endl;
			double drxva=Vatoms[vac_nr]->get_drx();
			double xvac=atom_to_jump[m]->get_x();
			double dryva=Vatoms[vac_nr]->get_dry();
			double yvac=atom_to_jump[m]->get_y();
			double drzva=Vatoms[vac_nr]->get_drz();
			double zvac=atom_to_jump[m]->get_z();
			//r2x=(x-x0va)*(x-x0va);
			//r2y=(y-y0va)*(y-y0va);
			//r2z=(z-z0va)*(z-z0va);
			//r2vac=Vatoms[vac_nr]->get_r2()+r2x+r2y+r2z;
			double xjump = move(xjumper,xvac,latt_constx,boundary_con_x);
			double yjump = move(yjumper,yvac,latt_consty,boundary_con_y);
			double zjump = move(zjumper,zvac,latt_constz,boundary_con_z);
			// funkcja site::move: oblicza kierunek skoku atomu do vakancji
			//atom_to_jump[m]->r2(_r2,atom_to_jump[m],Vatoms[vac_nr]);		//ZROBIC ZEBY OBLICZAL DYSTANS R^2
			//cout<<"Atom->vacancy: x "<<moves[0]<<" y "<<moves[1]<<" z "<<moves[2]<<endl;
			if(xjump>1 or yjump>1 or zjump >1)
			{
				cout<<"move error "<<xjump<<" "<<yjump<<" "<<zjump<<endl;
				exit(1);
				}
			//w komorce pamieci (Vatoms[vac_nr] - wskaznik do tej komorki) przechowujacej wybrana 
			// zmieniam zmienne wakancje na zmienne atomu. Plus dodatki z dyfuzji
			Vatoms[vac_nr]->set_atom(atom_to_jump[m]->get_atom());
			//Vatoms[vac_nr]->set_r2(r2atom);
			Vatoms[vac_nr]->set_drx(drxat+xjump);
			Vatoms[vac_nr]->set_dry(dryat+yjump);
			Vatoms[vac_nr]->set_drz(drzat+zjump);
			Vatoms[vac_nr]->set_jumps(atom_to_jump[m]->get_jumps()+1);
			//Vatoms[vac_nr]->set_jz(atom_to_jump[m]->get_jz()+moves[2]);
			
			//teraz zmieniam wakancje
			atom_to_jump[m]->set_atom(0);
			//atom_to_jump[m]->set_r2(r2vac);
			atom_to_jump[m]->set_drx(drxva-xjump);
			atom_to_jump[m]->set_dry(dryva-yjump);
			atom_to_jump[m]->set_drz(drzva-zjump);
			atom_to_jump[m]->set_jumps(nr_jumps+1);
				
			//atom_to_jump[m]->set_jz(jz-moves[2]);
			
		//		control_output<<"Atom po operacja: "<<endl;
		//		control_output<<"adres: "<<Vatoms[vac_nr]<<" ";
		//		Vatoms[vac_nr]->show_site();
		//		Vatoms[vac_nr]->show_site_jumps();
		//		Vatoms[vac_nr]->show_neigh(1);
		//		Vatoms[vac_nr]->show_neigh(0);
			
			// nadpisanie wskaznika w Vatoms na nowy miejsce wakancji w pamieci	
			Vatoms[vac_nr]=atom_to_jump[m];
			
		//	control_output<<"Vacancy po operacji: "<<endl;
		//		control_output<<"adres: "<<atom_to_jump[m]<<" ";
		//		atom_to_jump[m]->show_site();
		//		atom_to_jump[m]->show_site_jumps();
		//		atom_to_jump[m]->show_neigh(1);
		//		atom_to_jump[m]->show_neigh(0);
			
			

		//	control_output<<"NEXT JUMP"<<endl;
			// SPRAWDZIC ATOM I WAKANCJE ADRESY,PARAMETRY,SASIADOW DLA 3 SKOKOW - ok
			// dodac obliczanie parametrow
			// czym rozni sie *wskaznik - obiekt? od wskaznik - adres? ??
			
		//		UWAGA
		/*
		 * 
		 * ADRSY ZOSTAWIAM BEZ ZMIAN
		 * ZMIENIAM TYLKO DANE NA KTORE POKAZUJA ADRESY
		 * CZYLI KOPIUJE ZAWARTSCI KOMOREK , WYKONUJE NA NICH OPERACJE I ZAPISUJE
		 * NIE ZMIENIAM ADRESOW, TYLKO ZAWARTOSC NA KTORA POKAZUJA
		 * 
		 * pod adresem nadal bedzie wakancja, tylko bedzie miala nowe wspolrzedne, sasiadow, itd.
		 * ktore beda przepisane z atomu
		 * I analogicznie dla atomu*/	
		}
	}
}
//string name_of_file1="neigh.dat";
string name_of_file2="generator.dat";
//ofstream out_data1(name_of_file1.c_str(),ios :: app);
ofstream out_data2(name_of_file2.c_str(),ios :: app);

//out_data1<<(T+(time/100.0))<<" "<<sasiady[0]<<" "<<sasiady[1]<<" "<<sasiady[2]<<" "<<sasiady[3]<<" "<<sasiady[4]<<" "<<sasiady[5]<<" "<<sasiady[6]<<" "<<sasiady[7]<<" "<<sasiady[8]<<endl;
out_data2<<(T+(time/100.0))<<" "<<generator[0]<<" "<<generator[1]<<" "<<generator[2]<<" "<<generator[3]<<" "<<generator[4]<<" "<<generator[5]<<" "<<generator[6]<<" "<<generator[7]<<" "<<generator[8]<<" "<<generator[9]<<" "<<endl;
delete(target);
return time/100.0;
}



	
double move(double x2, double x1, double latt_const,double boundary_con)
 {
	 //x2 - polozenie wakancji
	 //x1 - polozenie atomu
	 double r2x=0.0;

	if((x2-x1)*(x2-x1)<(latt_const+0.5*latt_const)*(latt_const+0.5*latt_const))	//jesli nie jest na brzegu
	{
		r2x=(x2-x1);	//+ atom ruszyl w prawo , a wakncja w lewo
	}
	else   //jesli jest na brzegu
	{
		if((x2-x1)>0)	//i atom ruszyl w prawo
		{r2x=(x2 - boundary_con -x1);
	//	control_output<<" r2x: "<<x2<<" "<<x1<<" "<<" "<<boundary_con<<" "<<r2x<<endl;
		//int o;
		//cin>>o;
		}
		else
		{
		r2x=(x2 + boundary_con - x1);
	//	control_output<<" r2x: "<<x2<<" "<<x1<<" "<<" "<<boundary_con<<" "<<r2x<<endl;
		//int o;
		//cin>>o;
	}
	}
	return r2x;//dodatnie jesli atom ruszyl w prawo
	
	}

/*-----------------------------------------------------------*/

int save_results(lattice *sample, vector <task> &savings, string output, double a, double b)
{
	//control_output<<"Save: "<<name<<endl;
	
	for(unsigned int i=0;i<savings.size();i++)
	{
		string name="";
		name=savings[i].get_name();
		//if(name==output)
		//{
				
			if(name=="make_pic")
			{	
				vector <double> parameters;
				savings[i].get_parameters(parameters);
				if(parameters.size()!=7){control_output<<"Wrong parameter list in conf.in -> make_pic: 7"<<endl;exit(1);}
				long step_make = parameters[0];
				if(a == 0.0 and b == 0){step_make = 0;}
	    		wektor make_pic_vec_st(parameters[1],parameters[2],parameters[3]);
	    		wektor make_pic_vec(parameters[4],parameters[5],parameters[6]);	
				sample->makepic(b,step_make,make_pic_vec_st,make_pic_vec, output);
			}
			else if(name=="pic_stech")
			{
				vector <double> parameters;
				savings[i].get_parameters(parameters);
				if(parameters.size()!=7){control_output<<"Wrong parameter list in conf.in -> pic_stech: 7"<<endl;exit(1);}
				long step_make = parameters[0];
				if(a == 0.0 and b == 0){step_make = 0;}

		        wektor make_pic_vec_st(parameters[1],parameters[2],parameters[3]);
		        wektor make_pic_vec(parameters[4],parameters[5],parameters[6]);	
				sample->pic_stech(b,step_make,make_pic_vec_st,make_pic_vec, output);
			}
			else if(name=="pic_diff")
			{
				vector <double> parameters;
				savings[i].get_parameters(parameters);
				if(parameters.size()!=7){control_output<<"Wrong parameter list in conf.in -> pic_dif: 7"<<endl;exit(1);}
				long step_make = parameters[0];
				if(a == 0.0 and b == 0){step_make = 0;}

		        wektor make_pic_vec_st(parameters[1],parameters[2],parameters[3]);
		        wektor make_pic_vec(parameters[4],parameters[5],parameters[6]);	
				sample->pic_diff(b,step_make,make_pic_vec_st,make_pic_vec, output);
			}
			else if(name=="Energy")
			{	
				vector <double> parameters;
				savings[i].get_parameters(parameters);
				if(parameters.size()!=1){control_output<<"Wrong parameter list in conf.in -> Energy: 1"<<endl;exit(1);}
				int global_on = parameters[0];
				sample->save_energy(a,b,output,global_on);
			}
			else if(name=="Natoms")
			{
				vector <double> parameters;
				savings[i].get_parameters(parameters);
				if(parameters.size()!=1){control_output<<"Wrong parameter list in conf.in -> Natoms: 1"<<endl;exit(1);}
				int global_on = parameters[0];
				sample->save_Natoms(a,b,output,global_on);
			}
			else if(name=="NandE")
			{
				vector <double> parameters;
				savings[i].get_parameters(parameters);
				if(parameters.size()!=1){control_output<<"Wrong parameter list in conf.in -> EandE: 1"<<endl;exit(1);}
				int global_on = parameters[0];
				sample->save_NandE(a,b,output,global_on);
			}
			else if(name=="SRO")
			{
				sample->save_SRO(a,b,output);
			}
			else if(name=="Rd2")
			{
				sample->save_dR(a,b,output);
			}
			else if(name=="HIST")
			{
				vector <double> parameters;
				savings[i].get_parameters(parameters);
				if(parameters.size()!=4){control_output<<"Wrong parameter list in conf.in -> HIST: 4"<<endl;exit(1);}
				double iod	= parameters[0];
				double ido	= parameters[1];
				double istep= parameters[2];
				double idir = parameters[3];
				
				string name_of_file="";	
				stringstream total(output);
				int word_count=0 ;
				string word;
				while( total >> word ) ++word_count;
				if(word_count == 1)
				{
					name_of_file=output;
				}
				else if(word_count == 2)
				{
					stringstream ss(output);
					int log=0;
					while(ss>>word){
					if(log==0){
						name_of_file=word;log++;}
					else if(log>=1){log++;}
					else {control_output<<"ERROR in lattice::save_Natoms: "<<log<<endl;exit(1);}
					}
				}else if(word_count>2){control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
				else{control_output<<"ERROR in lattice::make_pic->file_name: "<<word_count<<endl;exit(1);}
			
				if( !( opt_equi.save_stauts() ) ){	//init HIST to save
					opt_equi.init_save(iod,ido,istep, idir);
				}
					
				opt_equi.save_hist(a,b,name_of_file);	//save flux fro opcja::HIST
				
			}
			else
			{
				control_output<<"No posibility of doing "<<name<<" to save."<<endl;
			}
	
		};
	
	return 0;
	}

/*-----------------------------------------------------------*/

int execute_task(task &comenda, vector <task> &savings, lattice *sample)
{
	control_output<<" Execute: ";
	comenda.show_task();
	
	string name=comenda.get_name();
	vector <double> parameters;
	comenda.get_parameters(parameters);  
		
	if(name=="INSERT")
	{
		initialize();
		initialize_seed();
		insert_atoms(parameters[0],parameters[1],parameters[2],sample);
		
	}	
			
	if(name=="EXCHANGE")
	{
		initialize();
		initialize_seed();
		exchange_atoms(parameters[0],parameters[1],parameters[2],sample);	
	}	

	if(name=="DIRECT")
	{
		initialize();
		initialize_seed();
		long step=parameters[0];
		long sub_step=parameters[1];
		double T = parameters[2];
		for (long i=1; i<=step;i++){
			direct_exchange(sample, sub_step, T);
			save_results(sample,savings,"dir",sub_step,i);			
		}
	}	
	
	if(name=="WIDOM_RND")
	{
		initialize();
		initialize_seed();
		long step=parameters[0];
		long sub_step=parameters[1];
		double T = parameters[2];	
		for (long i=1; i<=step;i++){
			widom(sample,i,sub_step, T);
			save_results(sample,savings,"wid",sub_step,i);			
		}
	}	

	
	if(name=="WIDOM")
	{
		initialize();
		initialize_seed();
		long step=parameters[0];
		long sub_step=parameters[1];
		double T = parameters[2];
		for (long i=1; i<=step;i++){
			widom_rnd(sample,i,sub_step, T);
			save_results(sample,savings,"wid",sub_step,i);			
		}
	}	

	
	if(name=="SGCMC")
	{
		initialize();
		initialize_seed();
		control_output<<"SGCMC run..."<<endl;
		long main_step=parameters[0];
		long sub_step=parameters[1];
		double Temp=parameters[2];
//		int plik1=parameters[3];	//ile pot. chem. 
//		int plik2=parameters[4];
		//int =parameters[5];
		//double odch=parameters[6];
		//double zmiana=parameters[7];
		//double sr_numer=parameters[8];
		
		//WCZYTAJ chem.in
		ifstream chem("chem.in",ios :: in);	//lista potencjalow
		vector <vector <double> > chem_list;
		
		if( chem.good() )
		{
        std::string napis;
     //   std::cout << "\nZawartosc pliku:" << std::endl;
        while( !chem.eof() )
        {
            getline( chem, napis );
            //control_output<<"String Line: "<< napis << std::endl;
            double data;
       		vector <double> line;
			istringstream string_line(napis);
//			int line_count=0;
			while(string_line>>data)
			{
		//		cout<<typeid(data).name()<<" "<<data<<endl;
				line.push_back(data);
		//		cout<<line.size()<<endl;
			}
			chem_list.push_back(line);
        }
      //  chem_list.erase(chem_list.end()); //usun linie konca pliku dodana do vectora
        chem.close();
		} else {control_output<< "Error! Nie udalo otworzyc sie pliku chem.in!" << std::endl;
		exit(0);
		}
		for(unsigned int i=0; i<chem_list.size();i++)
		{	//cout<<chem_list[i].size()<<endl;
			if((chem_list[i].size() +1) != sample->get_atom_typ_numbers())	// +1 poniewaz wakancja
			{
				control_output<<"ERROR: execute_task: SGCMC in reading chem.in file. ";
				control_output<<"Number of types different than declared in structure.dat file"<<endl;
				control_output<<"chem: "<<(chem_list[i].size()+1)<<" structure: "<<sample->get_atom_typ_numbers()<<endl;
				exit(0);
			}
		}
		//chem_list.size() - ile wierszy jest w pliku (punktow do symulacji)
		//chem_list[i].size() - ile typow atomow (pot. chem) jest zadeklarowanych
		//WCZYTALEM chem.in
		for(unsigned int j=0;j<chem_list.size();j++)
		{
			vector <double> pot_chem;	//chem_list[j]
			pot_chem.clear();
			pot_chem.reserve(5);
			pot_chem.push_back(0.0);
			
			stringstream s;
			string file_name="";

			//tworze nazwe pliku do zapisu
			for(unsigned int element=0;element<chem_list[j].size();element++)	
			{
				double tmp_value = chem_list[j][element];
				//tworzy nazwe pliku do zapisu
				s<<tmp_value;
				file_name=file_name+"u"+s.str();
				s.str("");
				//uzupelniam pot_chem
				pot_chem.push_back(tmp_value);
			}
	
			for(long i=1;i<=main_step;i++)
			{	
				sgcmc(sample,sub_step,Temp,pot_chem);
				save_results(sample,savings,file_name,sub_step,i);
			}
		}
	control_output<<"SGCMC end"<<endl;
	}

	if(name=="RESIDENT_ENERGY")
	{	
		initialize();
		initialize_seed();
		sample->set_atoms_list(Vatoms,0);
		double Time= 0.0;
		long step=0;
		long main_step=parameters[0]+1;
		long sub_step=parameters[1];
		int rob_plik=parameters[2];
		double Temperatura=parameters[3];
	//	double *w1=&parameters[4];
	//	double **w2=&w1;

		int licz_print=0;
		int licz_plik=1;
		int	proc_zad = main_step/10;  
		
		static int files_nrRTAE=0;
		string name_of_file="";


		ifstream file_in("barriers.in",ios :: in);	//lista potencjalow
		//vector <vector <double> > bar_list;
		if( file_in.good() )
		{
        string napis;
     //   std::cout << "\nZawartosc pliku:" << std::endl;
        while( !file_in.eof() )
        {
            getline( file_in, napis );
        //    control_output<<"String Line: "<< napis << std::endl;
            double data;
       		vector <double> line;
			istringstream string_line(napis);
//			int line_count=0;
			while(string_line>>data)
			{
		//		cout<<typeid(data).name()<<" "<<data<<endl;
				line.push_back(data);
		//		cout<<line.size()<<endl;
			}
			simple_bars.push_back(line);
        }
      //  chem_list.erase(chem_list.end()); //usun linie konca pliku dodana do vectora
        file_in.close();
		} else {control_output<< "Error! Nie udalo otworzyc sie pliku barriers.in!" << std::endl;
		exit(0);
		}
		if(simple_bars.size() != sample->get_atom_typ_numbers())
		{
			control_output<<"ERROR: execute_task: RESIDENT_ENERGY in reading barriers.in file. ";
			control_output<<"Number of types different than declared in structure.in file"<<endl;
			control_output<<"bar: "<<(simple_bars.size())<<" structure: "<<sample->get_atom_typ_numbers()<<endl;
			exit(0);
		}
		unsigned int zones = pot.get_coordination_number();
		control_output<<"Barriers initialized ..."<<endl;
		for(unsigned int i=0; i<simple_bars.size();i++)
		{	
			for(unsigned int j=0; j<simple_bars[i].size();j++)
			{
				control_output<<"i "<<i<<"j "<<j<<": "<<simple_bars[i][j]<<endl;
			}
			
			if((simple_bars[i].size()) != zones)
			{
			control_output<<"ERROR: execute_task: RESIDENT_ENERGY in reading barriers.in file. ";
			control_output<<"Number of zones different than declared in energy.in file"<<endl;
			control_output<<"bar: "<<(simple_bars[i].size())<<" energy: "<<zones<<endl;
			exit(0);
			}
		}

		

		for(long j=1;j<main_step;j++)
		{
			
			if(step==0 and Time==0.0)
			{
				files_nrRTAE++;
				stringstream s;
				s<<files_nrRTAE;
				name_of_file=s.str();
			}				
			if(j==(proc_zad*licz_print)+1)
			{
				control_output<<unitbuf<<"#";
				cout<<unitbuf<<"#";
				licz_print++;	
					
				
			}
				
			//if(licz_print>9)
			//{
		//		cout<<endl;
		//		control_output<<endl;
		//	}
			
			Time=Time+residence_time_energy(sample,sub_step,Temperatura,files_nrRTAE);
			step++;
			//control_output<<Time<<" "<<files_nrRTAE<<" "<<files_nrRTAE*step*sub_step<<endl;
			save_results(sample,savings,name_of_file,Time,step);
				
			if(((j*sub_step))==(rob_plik*licz_plik))
			{
				Time=0.0;
				step=0;
				licz_plik++;
				
			}
			
	//		if(licz_plik%rob_plik==0)
		//	{
	//			sample->refresh_structure("1pic0.xyz");
	//		}
		}
	}

	if(name=="RESIDENT_BAR")
	{
		initialize();
		initialize_seed();
		sample->set_atoms_list(Vatoms,0);
		double Time= 0.0;
		long step=0;
		long main_step=parameters[0]+1;
		long sub_step=parameters[1];
		int rob_plik=parameters[2];
		double Temperatura=parameters[3];
	//	double *w1=&parameters[4];
	//	double **w2=&w1;

		int licz_print=0;
		int licz_plik=1;
		int	proc_zad = main_step/10;  
		
		static int files_nrRTAE=0;
		string name_of_file="";

		for(int b1=0;b1<3;b1++){
			for(int b2=0;b2<2;b2++){
				for(int b3=0;b3<2;b3++){
					for(int b4=0;b4<3;b4++){
						for(int b5=0;b5<3;b5++){
							for(int b6=0;b6<3;b6++){
								for(int b7=0;b7<3;b7++){
									for(int b8=0;b8<3;b8++){
										for(int b9=0;b9<3;b9++){
										bar_list[b1][b2][b3][b4][b5][b6][b7][b8][b9]=1000.0;
		}}}}}}}}}

		ifstream file_in("barriers.in",ios :: in);	//lista potencjalow
		//vector <vector <double> > bar_list;
		if( file_in.good() )
		{
        string napis;
//        control_output << "\nZawartosc pliku:" << std::endl;
        while( !file_in.eof() )
        {
            getline( file_in, napis );
//            control_output<<"String Line: "<< napis << std::endl;
            double data;
       		vector <double> line;
			istringstream string_line(napis);
//			int line_count=0;
			while(string_line>>data)
			{
			//	cout<<typeid(data).name()<<" "<<data<<endl;
				line.push_back(data);
			//	cout<<line.size()<<endl;
			}
			if(line.size()==10){
			    bar_list[int(line[0])][int(line[1])][int(line[2])][int(line[3])][int(line[4])][int(line[5])][int(line[6])][int(line[7])][int(line[8])]=double(line[9]);
//			bar_list.push_back(line);
			}
			else if (line.size()!= 0){
			control_output<< "Error! Wrong data format in barriers.in! (size different than 10)" << std::endl;
			control_output<<line.size()<<endl;
			exit(0);
			}			
        }
      //  chem_list.erase(chem_list.end()); //usun linie konca pliku dodana do vectora
        file_in.close();
		} else {control_output<< "Error! Nie udalo otworzyc sie pliku barriers.in!" << std::endl;
		exit(0);
		}


//		for(int b1=0;b1<3;b1++){
//			for(int b2=0;b2<2;b2++){
//				for(int b3=0;b3<2;b3++){
//					for(int b4=0;b4<3;b4++){
//						for(int b5=0;b5<3;b5++){
//							for(int b6=0;b6<3;b6++){
//								for(int b7=0;b7<3;b7++){
//									for(int b8=0;b8<3;b8++){
//										for(int b9=0;b9<3;b9++){
//		control_output<<b1<<" "<<b2<<" "<<b3<<" "<<b4<<" "<<b5<<" "<<b6<<" "<<b7<<" "<<b8<<" "<<b9<<" "<<bar_list[b1][b2][b3][b4][b5][b6][b7][b8][b9]<<endl;
//		}}}}}}}}}
		

		for(long j=1;j<main_step;j++)
		{
			
			if(step==0 and Time==0.0)
			{
				files_nrRTAE++;
				stringstream s;
				s<<files_nrRTAE;
				name_of_file=s.str();
			}				
			if(j==(proc_zad*licz_print)+1)
			{
				control_output<<unitbuf<<"#";
				cout<<unitbuf<<"#";
				licz_print++;	
					
				
			}
				
			//if(licz_print>9)
			//{
		//		cout<<endl;
		//		control_output<<endl;
		//	}
			
			Time=Time+residence_time_energy(sample,sub_step,Temperatura,files_nrRTAE);
			step++;
		//	control_output<<Time<<endl;
			save_results(sample,savings,name_of_file,Time,step);
				
			if(((j*sub_step))==(rob_plik*licz_plik))
			{
				Time=0.0;
				step=0;
				licz_plik++;
				
			}
			
	//		if(licz_plik%rob_plik==0)
		//	{
	//			sample->refresh_structure("1pic0.xyz");
	//		}
		}
		
	
	control_output<<"RTA finished... "<<endl;		
	cout<<"RTA finished... "<<endl;		
	
		
	}


	if(name=="RESIDENT")
	{
		
		save_results(sample,savings,"0",0.0,0);
		initialize();
		initialize_seed();
		
		double Time= 0.0;
		long step=0;
		long main_step=parameters[0]+1;
		long sub_step=parameters[1];
		int rob_plik=parameters[2];
		double Temperatura=parameters[3];

		int licz_print=0;
		int licz_plik=1;
		int	proc_zad = main_step/10;  
		
		static int files_nrRTAE=0;
		string name_of_file="";

		ifstream file_in("barriers.in",ios :: in);	//lista potencjalow

		if( file_in.good() )
		{
        string napis;
     //   std::cout << "\nZawartosc pliku:" << std::endl;
        while( !file_in.eof() )
        {
            getline( file_in, napis );
        //    control_output<<"String Line: "<< napis << std::endl;
            double data;
       		vector <double> line;
			istringstream string_line(napis);
//			int line_count=0;
			while(string_line>>data)
			{
		//		cout<<typeid(data).name()<<" "<<data<<endl;
				line.push_back(data);
		//		cout<<line.size()<<endl;
			}
			simple_bars.push_back(line);
        }
      //  chem_list.erase(chem_list.end()); //usun linie konca pliku dodana do vectora
        file_in.close();
		} else {control_output<< "Error! Nie udalo otworzyc sie pliku barriers.in!" << std::endl;
		exit(0);
		}
		if(simple_bars.size() != sample->get_atom_typ_numbers())
		{
			control_output<<"ERROR: execute_task: RESIDENT_ENERGY in reading barriers.in file. ";
			control_output<<"Number of types different than declared in structure.in file"<<endl;
			control_output<<"bar: "<<(simple_bars.size())<<" structure: "<<sample->get_atom_typ_numbers()<<endl;
			exit(0);
		}
		unsigned int zones = pot.get_coordination_number();
		control_output<<"Barriers initialized ..."<<endl;
		for(unsigned int i=0; i<simple_bars.size();i++)
		{	
			for(unsigned int j=0; j<simple_bars[i].size();j++)
			{
				control_output<<"i "<<i<<"j "<<j<<": "<<simple_bars[i][j]<<endl;
			}
			
			if((simple_bars[i].size()) != zones)
			{
			control_output<<"ERROR: execute_task: RESIDENT_ENERGY in reading barriers.in file. ";
			control_output<<"Number of zones different than declared in energy.in file"<<endl;
			control_output<<"bar: "<<(simple_bars[i].size())<<" energy: "<<zones<<endl;
			exit(0);
			}
		}
		
		sample->set_atoms_list(V_LIST,0);

		for(long j=1;j<main_step;j++)
		{
			
			initialize();
			initialize_seed();
			
			if(step==0 and Time==0.0)
			{
				files_nrRTAE++;
				stringstream s;
				s<<files_nrRTAE<<" "<<rob_plik;
				name_of_file=s.str();
			
			}				
			if(j==(proc_zad*licz_print)+1)
			{
//				control_output<<unitbuf<<"#";
				cout<<unitbuf<<"#";
				licz_print++;	
			}
			
			Time=Time+residence_time(sample,sub_step,Temperatura,files_nrRTAE);
			step++;
			save_results(sample,savings,name_of_file,Time,step);
				
			if(((j))==(rob_plik*licz_plik))
			{
				Time=0.0;
				step=0;
				licz_plik++;
				
			}
			
	//		if(licz_plik%rob_plik==0)
		//	{
	//			sample->refresh_structure("1pic0.xyz");
	//		}
		}
	control_output<<"RTA finished... "<<endl;				
	control_output<<"VAC: "<<V_LIST.size()<<" EVENTS: "<<EVENTS.size()<<endl;
	simple_bars.clear();
	EVENTS.clear();
	V_LIST.clear();
	}  
				
	if(name=="RTA_RAND_ALLOY")
	{	
		initialize();
		initialize_seed();
		sample->set_atoms_list(Vatoms,0);
		double Time= 0.0;
		long step=0;
		long main_step=parameters[0]+1;
		long sub_step=parameters[1];
		int rob_plik=parameters[2];
//		double T=parameters[3];
		double w1=parameters[4];
		double w2=parameters[5];
		
		int licz_print=0;
		int licz_plik=1;
		int	proc_zad = main_step/10;    

		static int files_nrRTAR=0;
		string name_of_file="";
		

		for(long j=1;j<main_step;j++)
		{
			
			if(step==0 and Time==0.0)
			{
				files_nrRTAR++;
				stringstream s;
				s<<files_nrRTAR;
				name_of_file=s.str();
				sample->clear_dR();
			}				
		
					
			if(j==(proc_zad*licz_print)+1)
			{
				control_output<<unitbuf<<"#";
				cout<<unitbuf<<"#";
				licz_print++;	
					
				
			}
				
			//if(licz_print>9)
		//	{
		//		cout<<endl;
		//		control_output<<endl;
		//	}
				
			Time=Time+ RTA_random_alloy(sample,sub_step,Time,w1,w2);
			step++;
		//	control_output<<Time<<endl;
			save_results(sample,savings,"",Time,step);
			
		
			if(((j*sub_step))==(rob_plik*licz_plik))
			{
				Time=0.0;
				step=0;
				licz_plik++;
				
			}		
			
		}
		
	
	control_output<<"RTA_random finished... "<<endl;		
	}
	


	if(name=="VAC_MECH")
	{
		initialize();
		initialize_seed();
		sample->set_atoms_list(Vatoms,0);
		double Time= 0.0;
		long step=0;
		long main_step=parameters[0]+1;
		long sub_step=parameters[1];
		int rob_plik=parameters[2];
		double Temperatura=parameters[3];
		long direct_step=parameters[4];

		int licz_print=0;
		int licz_plik=1;
		int	proc_zad = main_step/10;  
		
		static int files_nrRTAE=0;
		string name_of_file="";

		ifstream file_in("barriers.in",ios :: in);	//lista potencjalow

		if( file_in.good() )
		{
        string napis;
     //   std::cout << "\nZawartosc pliku:" << std::endl;
        while( !file_in.eof() )
        {
            getline( file_in, napis );
        //    control_output<<"String Line: "<< napis << std::endl;
            double data;
       		vector <double> line;
			istringstream string_line(napis);
//			int line_count=0;
			while(string_line>>data)
			{
		//		cout<<typeid(data).name()<<" "<<data<<endl;
				line.push_back(data);
		//		cout<<line.size()<<endl;
			}
			simple_bars.push_back(line);
        }
      //  chem_list.erase(chem_list.end()); //usun linie konca pliku dodana do vectora
        file_in.close();
		} else {control_output<< "Error! Nie udalo otworzyc sie pliku barriers.in!" << std::endl;
		exit(0);
		}
		if(simple_bars.size() != sample->get_atom_typ_numbers())
		{
			control_output<<"ERROR: execute_task: RESIDENT_ENERGY in reading barriers.in file. ";
			control_output<<"Number of types different than declared in structure.in file"<<endl;
			control_output<<"bar: "<<(simple_bars.size())<<" structure: "<<sample->get_atom_typ_numbers()<<endl;
			exit(0);
		}
		unsigned int zones = pot.get_coordination_number();
		control_output<<"Barriers initialized ..."<<endl;
		for(unsigned int i=0; i<simple_bars.size();i++)
		{	
			for(unsigned int j=0; j<simple_bars[i].size();j++)
			{
				control_output<<"i "<<i<<"j "<<j<<": "<<simple_bars[i][j]<<endl;
			}
			
			if((simple_bars[i].size()) != zones)
			{
			control_output<<"ERROR: execute_task: RESIDENT_ENERGY in reading barriers.in file. ";
			control_output<<"Number of zones different than declared in energy.in file"<<endl;
			control_output<<"bar: "<<(simple_bars[i].size())<<" energy: "<<zones<<endl;
			exit(0);
			}
		}


		for(long j=1;j<main_step;j++)
		{
			
			if(step==0 and Time==0.0)
			{
				files_nrRTAE++;
				stringstream s;
				s<<files_nrRTAE<<" "<<rob_plik;
				name_of_file=s.str();
			
			}				
			if(j==(proc_zad*licz_print)+1)
			{
				control_output<<unitbuf<<"#";
				cout<<unitbuf<<"#";
				licz_print++;	
			}
			
			Time=Time+vac_mechanism(sample,sub_step,direct_step,Temperatura,files_nrRTAE);
			step++;
			save_results(sample,savings,name_of_file,Time,step);
				
			if(((j))==(rob_plik*licz_plik))
			{
				Time=0.0;
				step=0;
				licz_plik++;
				
			}
			
	//		if(licz_plik%rob_plik==0)
		//	{
	//			sample->refresh_structure("1pic0.xyz");
	//		}
		}
		simple_bars.clear();
	control_output<<"VAC_MECH finished... "<<endl;				
	}  
	return 0;
}
/*-----------------------------------------------------------*/

int main(int arg, char *argc[])
{
	//inicjalizacja generatora
	cout<<"Adres obiektu potencjal: "<<&pot<<endl;	
	initialize();
	
	// zmienne int
//	double bar1=0;
//	double bar2=0;
	int sizex=10;		//definiuja pamiec
	int sizey=10;
	int sizez=10;
	int ile_struktur=0;
//	int licz_print=1;
	int load_strucrure_flag=0;
	int methods_num=0;
	int saveings_num=0;
//	long main_step=1;
//	long sub_step=1;
//	long proc_zad=0;
//	long par1=0;
//	long par2=0;
//	long par3=0;
//	double Time=0.0;
	double rmin;		// okreslaja obszar miedzy sferami z ktorego brane sa atomy sasiadow
	double rmax;
	
	//zmienne znakowe
	char   name_of_input_file[20]="conf.in";
	char structure_file[20];
	string inf="conf.in";
    string inf1=string(argc[1]);
    string text; // smietnik  
    string pic_flag;
    string method1;
    string method2;
	control_output<<argc[1]<<endl;
	
	//zmienne wektor
	wektor make_pic_vec_st;
	wektor make_pic_vec_ed;
	wektor boundary_con_at;	
	//opcje warunkow brzegowych dla atomow: <0 - wez sasiadow nn spoza obszaru symulacji
	//=0 - wez sasiadow nn ze strefy symulacji bez translacji(warstwa);  
	//>0 - wez sasiadow nn ze strefy symulacji z translacja(bulk);
	wektor boundary_con_en;
	//opcje warunkow brzegowych dla enegii
	wektor st_region;
	wektor end_region;
	wektor st_sim_area;
	wektor end_sim_area;
	//obszar symulacji UWAGA: musi byc <= od fill_vector w pliku structure.dat || x/y/z/left/right/_border w lattice()
	string tab_strct_file_name[15];
	wektor tab_wektorow_st[15];	// - poczatek pudelka
	wektor tab_wektorow_end[15]; // - koniec pudelka
	wektor tab_wektorow_set[15]; // - punkt poczatkowy wczytywania
	int tab_load_lattice[15];
	// obszar konfiguracji kotry bedzie wczytany z pliku. Maksymalnie z 15 plikow mozna wczytac
	// w tej werssji wylaczone
	wektor max_zone;
	// sluzy do sprawdzania czy atom jest na granicy
	if(inf1==inf)
	{
		
// Wczytuje parametry
		
		ifstream input_file(name_of_input_file);
		control_output<<"Reading METHODS..."<<endl;
		vector <task> schedule;
		vector <task> saveings;
		Vatoms.reserve(100000);
		input_file>>text>>methods_num;
		
		for(int i=0;i<methods_num;i++)
		{	
			vector <double> parameters;
			parameters.clear();
			int par_num=0;   
			string method_name=" ";
			input_file>>method_name>>par_num; 
			
			if(method_name == "Save:"){control_output<<"ERROR: bad methods list in conf.in"<<endl; exit(0);}	
			
			for(int j=0;j<par_num;j++)
				{	double par=0;
					input_file>>par;
					parameters.push_back(par);
				}
				
			task thing_to_do(method_name,parameters);	
			schedule.push_back(thing_to_do);	
		}
		
		for(unsigned int i=0;i<schedule.size();i++)
		{
			schedule[i].show_task();
		}

		control_output<<"Reading SAVINGS..."<<endl;
		
		input_file>>text>>saveings_num;
		
		if(text != "Save:"){control_output<<"ERROR: bad methods list in conf.in"<<endl; exit(0);}	

		
		for(int i=0;i<saveings_num;i++)
		{	
			vector <double> parameters;
			parameters.clear();
			int par_num=0;    
			string method_name=" ";
			input_file>>method_name>>par_num; 

			if(method_name == "Options:"){control_output<<"ERROR: bad save list in conf.in"<<endl; exit(0);}	

			
			for(int j=0;j<par_num;j++)
				{	double par=0;
					input_file>>par;
					parameters.push_back(par);
				}
				
			task thing_to_save(method_name,parameters);	
			saveings.push_back(thing_to_save);	
		}

		for(unsigned int i=0;i<saveings.size();i++)
		{
			saveings[i].show_task();
			}
			
		control_output<<"Reading OPTIONS..."<<endl;

		vector <task> options;
		options.reserve(5);
		options.clear();	
		int ile_option=0;
		input_file>>text>>ile_option;

		if(text != "Options:"){control_output<<"ERROR: bad save list in conf.in"<<endl; exit(0);}else{

		for(int i=0;i<ile_option;i++){
			int ile=0;
			input_file>>text>>ile;
			vector <double> parameters;
			for(int j=0;j<ile;j++){	
				double par=0;
				input_file>>par;
				parameters.push_back(par);
			}
			task tmp_option(text,parameters);	
			options.push_back(tmp_option);
		}
		
		for(unsigned int i=0;i<options.size();i++){options[i].show_task();}

		}
		
	
	control_output<<"Reading PARAMETERS..."<<endl;
		input_file>>text>>sizex>>sizey>>sizez;

		if(text != "Memory_Size:"){control_output<<"ERROR: bad reservuar list in conf.in"<<endl; exit(0);}	

		control_output<<"Memory reserved for: "<<sizex<<" "<<sizey<<" "<<sizez<<endl;
		control_output<<"Interacione wector: "<<endl;
		input_file>>text>>rmin>>rmax>>max_zone.x>>max_zone.y>>max_zone.z;
		max_zone.show();
		
	//simulation area, boundary conditions	
		input_file>>text>>st_region.x>>st_region.y>>st_region.z>>end_region.x>>end_region.y>>end_region.z;	
		input_file>>text>>st_sim_area.x>>st_sim_area.y>>st_sim_area.z>>end_sim_area.x>>end_sim_area.y>>end_sim_area.z;
		input_file>>text>>boundary_con_en.x>>boundary_con_en.y>>boundary_con_en.z>>boundary_con_at.x>>boundary_con_at.y>>boundary_con_at.z;
		control_output<<"Sample region:"<<endl;
		st_region.show();
	    end_region.show();
		control_output<<"simulation area:"<<endl;
		st_sim_area.show();
		end_sim_area.show();
		control_output<<"boundary conditions for energy and atom: 0 - non | 1 - nn | normal"<<endl;
		boundary_con_en.show();
		boundary_con_at.show();
		
	//structure creation
		input_file>>text>>ile_struktur>>load_strucrure_flag;
		control_output<<"Bazowe struktury: "<<ile_struktur<<endl;
		for(int i=0;i<ile_struktur;i++)
		{	
			double stpointx,stpointy,stpointz,endpointx,endpointy,endpointz,setpointx,setpointy,setpointz;
			int lattice_nr;
			input_file>>structure_file>>stpointx>>stpointy>>stpointz>>endpointx>>endpointy>>endpointz>>setpointx>>setpointy>>setpointz>>lattice_nr;
			tab_strct_file_name[i]=structure_file;
			tab_wektorow_st[i].x=stpointx;
			tab_wektorow_st[i].y=stpointy;
			tab_wektorow_st[i].z=stpointz;
			tab_wektorow_end[i].x=endpointx;
			tab_wektorow_end[i].y=endpointy;
			tab_wektorow_end[i].z=endpointz;
			tab_wektorow_set[i].x=setpointx;
			tab_wektorow_set[i].y=setpointy;
			tab_wektorow_set[i].z=setpointz;
			tab_load_lattice[i]=lattice_nr;
			control_output<<tab_strct_file_name[i]<<endl;
			control_output<<"st_vec: ";
			tab_wektorow_st[i].show();
			control_output<<"end_vec: ";
			tab_wektorow_end[i].show();
			control_output<<"set_vec: ";
			tab_wektorow_set[i].show();
			control_output<<"load_lattice_nr: "<<tab_load_lattice[i]<<endl;
		}				
	input_file.close();
	// Koniec wczytywania //

	control_output<<"Buduje matrix sitow"<<endl;
 
	//zbudowanie struktury, wczytanie energii
	lattice sample(sizex,sizey,sizez);	
	lattice *_sample;
	_sample=&sample; 
	//	_sample->check_atoms();
	sample.set_alg_objects(EVENTS,simple_bars,pot);
	opt_equi.set_opcja_lattice(_sample);
		
	//wczytuje strukture poczatkowa
	if(load_strucrure_flag){
		control_output<<"Budowanie sample..."<<endl;
		for(int i=0;i<ile_struktur;i++){
			_sample->read_structure(tab_strct_file_name[i],tab_wektorow_st[i],tab_wektorow_end[i],tab_wektorow_set[i], tab_load_lattice[i]);	
		}
		control_output<<"Zbudowano!"<<endl;
	}
		
	//Wczytuje potencjaly
	pot.init(_sample->get_atom_typ_numbers());
	

	//wczytuje do zmiennych w lattice obszar symulacji oraz ustawieienie warunkow brzegowych dla atomow i energii
	_sample->simulation_initialize(rmin,rmax,st_sim_area,end_sim_area,boundary_con_at,boundary_con_en,max_zone, st_region, end_region);

	//wczytuje sasiadow w I strefie zgodnie z podanymi ustawieniami
	_sample->jumps_shell_init();
	_sample->interaction_shell_init();


	//opcja rownowarzenia init	
	control_output<<"Actual options: "<<options.size()<<endl;
	for(unsigned int i=0;i<options.size();i++){
		vector <double> parameters;	parameters.reserve(10);
		options[i].show_task();
		string name=options[i].get_name();
		options[i].get_parameters(parameters);
		opt_equi.execute(name,parameters);
	}
	
	//wykonaj metody	
	control_output<<"Save init..."<<endl;	
	save_results(_sample,saveings,"0",0.0,0);

	control_output<<"Starting schedule"<<endl;
	for(unsigned int i=0;i<schedule.size();i++){
		string name=schedule[i].get_name();
		vector <double> parameters;
		schedule[i].get_parameters(parameters);
		if(name=="LOOP"){
			vector <double> parameters;
			schedule[i].get_parameters(parameters);
			int loop_step = parameters[0];
			int next_methods = parameters[1];
			for(int loop_iter=0;loop_iter<loop_step;loop_iter++){
				for(int method_iter=(i+1);method_iter<=next_methods;method_iter++){
					execute_task(schedule[method_iter],saveings,_sample);
				}
			}
		}
		else{
			execute_task(schedule[i],saveings,_sample);
		}
	}
	
	}
	else
	{
		control_output<<"Error reading conf.in...exiting"<<endl;
		exit(0);
	}

return 0;
}
