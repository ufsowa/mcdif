//opcja_swap//

#include "opcja.h"

void opcja :: execute(string name, vector<double>&parameters){

	if(name=="EQUILIBRATE"){
		this->init_EQ(parameters);
	}
	else if(name=="RESERVUAR"){
		this->init_reservuar(parameters);
	}
	else{
		control_output<<"No posibility of doing "<<name<<" as an option."<<endl;
	}

}


int opcja :: choose_typ(plaster& bin, bool sig){
	
	int A =  bin.size(1);
	int B =  bin.size(2);
	int SIZE = A + B;
	int TYP = -1;
			
//	control_output<<" rozmiar typow 1/2 "<<A<<"/"<<B<<" "<<SIZE<<" ";
	long N=(long)(rnd()*(SIZE))+1;
	if((A > 0) and (0 < N) and (N <= (A+1))){
		TYP=1;
//		control_output<<" nr wylosowanego situ A: "<<N<<" "<<TYP;
	}
	else if((B > 0) and ((A+1) <= N) and (N <= (SIZE+1))){
//		N-=A;
		TYP=2;
//		control_output<<" nr wylosowanego situ B: "<<N<<" "<<TYP;
	}
	else if(sig){control_output<<"Error in opcja::choose_typ: "<<N<<" "<<A<<" "<<B<<" "<<SIZE<<endl;
		control_output<<"Try to remove atom witch does not egzist "<<endl;exit(1);}
	else if(!sig){control_output<<"Warrning in opcja::choose_typ: "<<N<<" "<<A<<" "<<B<<" "<<SIZE<<endl;
		control_output<<"Try to remove atom witch does not egzist "<<endl;}

//	control_output<<N<<" "<<TYP<<endl;
	return TYP;			//TYP < 0 means that there is no atoms in the bin
}


int opcja :: check_stech(double _stech, double _vac, double _size){
	
//	control_output<<"check_stech..."<<endl;
	double status=0;
	//iterowac po equi_curve[stech][vac]
	//szukam wartosci vac dla danego stech
			//przelicz koncentracje na ilosc: 
	double _eq_vac = Ceq_vac(_stech);
	double _err_vac = errCeq_vac(_stech);
	//long size = get_sample_size();
//	int err = (err_vac)*_size;
	double vac=(_vac)*(_size);
//	int eq_vac=(_eq_vac)*(_size);
	double delta_left = ((_eq_vac - _err_vac)*_size);
	double delta_right = ((_eq_vac + _err_vac)*_size);
	if(vac < delta_left){
		status=delta_left-vac;
	}
	else if(vac > delta_right){
		status=delta_right-vac;
	}
	else{
		status=0;
	}	
//	control_output<<" stech "<<_stech<<" vac "<<_vac<<" status "<<status<<endl;
	return round(status);	//if < 0 then must remove vac in plaster
}					//if > 0 then must create vac in plaster

double opcja :: errCeq_vac(double stech){
//	control_output<<"call Ceq_vac("<<stech<<")"<<endl;
	double a=0,b=0,CV=0.0;
	double eq_stech1 = 0.0;
	double eq_stech2 = 0.0;
	double eq_vac1 = 0.0;
	double eq_vac2 = 0.0;
	int size = (equi_curve.size()-1);
//	control_output<<"Ceq_vac ranges"<<endl;
//	for (int i=0; i <= size;i++){
//		eq_stech1 = equi_curve[i][0];
//		eq_vac1 = equi_curve[i][1];
//		control_output<<i<<" "<<eq_stech1<<" "<<eq_vac1<<endl;
//	}
	if(equi_curve[0].size() < 2){control_output<<"ERROR in opcja::errCeg_vac"<<endl;exit(1);}

	if(stech < equi_curve[0][0]){CV=equi_curve[0][2];}
	else if(stech > equi_curve[size][0]){
		//control_output<<"test: "<<equi_curve[size][0]<<endl;
		CV=equi_curve[size][2];}
	else{
	for (int i=1; i <= size;i++){
		eq_stech1 = equi_curve[i-1][0];
		eq_stech2 = equi_curve[i][0];
		if(eq_stech1 <= stech and eq_stech2 >= stech ){	//zakladam ze rosnaco jest stech w pliku
//			control_output<<" for "<<i;
			//fitowanie prostej do 2 punktow i szukanie dokladnej wartosci vac
				eq_vac1 = equi_curve[i-1][2];
				eq_vac2 = equi_curve[i][2];		
				a = (eq_vac2 - eq_vac1)/(eq_stech2 - eq_stech1);
				b = eq_vac1 - eq_stech1*a;
//				control_output<<" a "<<a<<" b "<<b<<" eq_vac "<<(a*stech+b)<<endl;
				break;
		}
	}
		CV=(a*stech+b);
	}
//	control_output<<" a "<<a<<" b "<<b<<" eq_vac "<<CV<<endl;
	return CV;
}
double opcja :: Ceq_vac(double stech){
	
//	control_output<<"call Ceq_vac("<<stech<<")"<<endl;
	double a=0,b=0,CV=0.0;
	double eq_stech1 = 0.0;
	double eq_stech2 = 0.0;
	double eq_vac1 = 0.0;
	double eq_vac2 = 0.0;
	int size = (equi_curve.size()-1);
//	control_output<<"Ceq_vac ranges"<<endl;
//	for (int i=0; i <= size;i++){
//		eq_stech1 = equi_curve[i][0];
//		eq_vac1 = equi_curve[i][1];
//		control_output<<i<<" "<<eq_stech1<<" "<<eq_vac1<<endl;
//	}
	if(equi_curve[0].size() < 1){control_output<<"ERROR in opcja::Ceg_vac"<<endl;exit(1);}
	if(stech < equi_curve[0][0]){CV=equi_curve[0][1];}
	else if(stech > equi_curve[size][0]){
		//control_output<<"test: "<<equi_curve[size][0]<<endl;
		CV=equi_curve[size][1];}
	else{
	for (int i=1; i <= size;i++){
		eq_stech1 = equi_curve[i-1][0];
		eq_stech2 = equi_curve[i][0];
		if(eq_stech1 <= stech and eq_stech2 >= stech ){	//zakladam ze rosnaco jest stech w pliku
//			control_output<<" for "<<i;
			//fitowanie prostej do 2 punktow i szukanie dokladnej wartosci vac
				eq_vac1 = equi_curve[i-1][1];
				eq_vac2 = equi_curve[i][1];		
				a = (eq_vac2 - eq_vac1)/(eq_stech2 - eq_stech1);
				b = eq_vac1 - eq_stech1*a;
//				control_output<<" a "<<a<<" b "<<b<<" eq_vac "<<(a*stech+b)<<endl;
				break;
		}
	}
		CV=(a*stech+b);
	}
//	control_output<<" a "<<a<<" b "<<b<<" eq_vac "<<CV<<endl;
	return CV;
}

bool opcja :: check_rezervuars(int i, int TYP){
	
	int typ_od=TYP,typ_do=TYP;	//,status=0;	//, ile_to_move=-1;
	bool do_move=false;
	long int V=0, Z=0;
		
//	if(TYP==0){
//		for(int t=1;t<reservuars[i].get_size_types();t++){
//		A+=reservuars[i].eq_flux_get(t);}
//		Z=A;
//	}else{
		V = reservuars[i].eq_flux_get(0) + reservuars[i].flux_net_get(0);
		Z=labs(V)+1;
//		}	
	
	
	
	
	if(Z > ROZMIAR[i]){//calkowita zmiana atomow wynosi tyle co 100% jednej plaszczyzny

			if(V>0)	//dV>0 to znaczy ze w rezerwuarze powstala nowa plaszczyzna wakancji -> probka w strone atomow
			{TYP_TO_MOVE=1;}
			else if(V<0)	//dV<0 to znaczy ze w rezerwuarze powstala nowa plaszczyzna atomow -> probka w strone wakancji
			{TYP_TO_MOVE=-1;}
			else{control_output<<"error in opcja::check_reservuars "<<V<<endl;exit(1);}
			do_move=true;
			MOVE_FRAME=true;	//set global variable to true
			REZ_TO_MOVE=i;	//set global var. which rezervuar to move
			//status=1;
			}
	else{	//jesli nie ma nowej plaszczyzny ale zabraklo atomow typu TYP w rezerwuarze		<< nie mozliwy scenariusz bo typ losowany z rezerwuaru
		if(reservuars[i].check(typ_od,typ_do,1)){	//ale jesli brakuje typow
		
		control_output<<"warrning in opcja::check_reservuars -> try remove typ which dose not exist in reservuar "<<TYP<<endl;
			if(TYP>0)	//TYP>0 to znaczy ze w rezerwuarze zabraklo atomow -> probka w strone atomow
			{TYP_TO_MOVE=1;}
			else if(TYP==0)		//TYP==0 to znaczy ze w rezerwuarze zabraklo wakancji -> probka w strone wakancji
			{TYP_TO_MOVE=-1;}
			else{control_output<<"error in opcja::check_reservuars (problem with types < 0) "<<TYP<<endl;exit(1);}
			do_move=true;
			MOVE_FRAME=true;	//set global variable to true
			REZ_TO_MOVE=i;	//set global var. which rezervuar to move
			//status=2;
		}
	}
//	if(MOVE_FRAME or SINGLE){control_output<<"|<"<<TYP<<"|"<<i<<"|"<<V<<"|"<<Z<<"|"<<ROZMIAR[i]<<"|"<<status;}
	return do_move;
}

void opcja :: create_vac(int nr, int ile, bool &FLAG){

	if(MOVE_FRAME or SINGLE){control_output<<" c: "<< ile;}
	bool MOVE = false;
	int rez = -1, j = -1;
	
	for( int i=0; i<(ile);i++){
//		for( int j=1;j<3;j++){		//MOZNA ZROBIC W ZALEZNOSCI OD TYPOW
			long N1,N2;
			site* rnd_vac=0;
			site* rnd_at=0;
			j=choose_typ(BLOKS[nr]);	//losuje typ atomu do wymiany z wakancja
	//		control_output<<"rozmiar typow "<<0<<" w bloku "<<bloks[nr].size(0);
	//		control_output<<" rozmiar typow "<<j<<" w bloku "<<bloks[nr].size(j);
			if(BLOKS[nr].size(j) <= 0){
				control_output<<endl;
				control_output<<"ERROR: in opcja::create_vac -> you want to remove \
element "<<j<<" that does not exist in blok\n	\
Probably error in opcja::init or opcja::reinit_bloks"<<endl; 
				cout<<endl;
				cout<<"ERROR: in opcja::create_vac -> you want to remove \
element "<<j<<" that does not exist in blok\n	\
Probably error in opcja::init or opcja::reinit_bloks"<<endl;exit(0);
			}
			
			while(1){
			N1=(long)(rnd()*(BLOKS[nr].size(j)));
			rnd_at=BLOKS[nr].get_site(j,N1);
			control_output<<" ctyp: "<<rnd_at->get_atom()<<endl;
			
			double d = 0.0, distanceL = 0.0, distanceR = 0.0;
			

			if(BIN_DIRECTION==1)
			{
				d=rnd_at->get_x();
			}
			else if(BIN_DIRECTION==2)
			{
				d=rnd_at->get_y();
			}
			else if(BIN_DIRECTION==3)
			{
				d=rnd_at->get_z();
			}
			else
			{
				cout<<"Wrong direction number in opcja::convert() "<<d<<endl;
				exit(1);
			}
			
			distanceL = abs(d - BIN_ST);
			distanceR = abs(d - BIN_END);
		//	cout<<"Absloute value of distance L/R: "<<distanceL<<"/"<<distanceR<<endl;
			
			if(distanceL > distanceR){
				rez = 1;
			}
			else if(distanceL < distanceR){
				rez = 0;
			}
			else if(distanceL == distanceR){
				double Q = rnd();
				if(Q<=0.5){rez = 0;}else{rez = 1;}
			}
			else{
				cout<<"Wrong direction number in opcja::convert() "<<distanceL<<"/"<<distanceR<<endl;
				exit(1);
			}
			
			MOVE = check_rezervuars(rez,0);	//sprawdz rezerwuwar pod katem dostepnych wakancji
			
			if (!MOVE){
			
			N2=(long)(rnd()*(reservuars[rez].size(0)));
			rnd_vac = reservuars[rez].get_site(0,N2);

//			rnd_vac - pointer to site in reservuar(vac)
//			rnd_at - pointer to site in blok (atom)


//			double old_E=POT.get_energy(rnd_at) + POT.get_energy(rnd_vac);
		//	double old_E=POT.get_energy(rnd_vac);

			rnd_vac->set_atom(j);
				
//

//			double new_E=POT.get_energy(rnd_at) + POT.get_energy(rnd_vac);
		//	double new_E=POT.get_energy(rnd_vac);

		//	double beta=1.0/(kB*T);
		//	double mi=0.0;
		 //   double P1=exp(beta*(mi-(new_E-old_E)));
		//	control_output<<" n: "<<new_E<<" o: "<<old_E<<" "<<(new_E-old_E)<<" "<<beta<<" P: "<<P1<<endl;
		//zmien z powrotem na old_typ jesli zdarzenie to nie zostalo trafione rnd()
		//	if(P1<rnd())
		//	{
		//	rnd_vac->set_atom(0);	
//			rnd_at->set_atom(j);
		//	}
		//	else{
			rnd_at->set_atom(0);
			rnd_vac->reset_site();
			rnd_at->reset_site();
															//kasujemy z listy typow atom z bloku
			BLOKS[nr].delete_site(j,N1);
			BLOKS[nr].add_site(0,rnd_at);
	
												//lista wakancji nie zawiera wakancji -> wymaga odswiezenia 
			reservuars[rez].delete_site(0,N2);
			reservuars[rez].add_site(j,rnd_vac);
			
			Vtoadd.push_back(rnd_at);		// w obszarze symulacji pojawila sie nowa wakancja
			
			break;
		//	}
			
		}else{		//if(MOVE)
								// Jak tak to do_move.
			break;	//wyjdz z petli while
		}
			
		}	//koniec while
//			a co gdy obszar symulacji obejmuje rezerwuary?
	//policz prawdopodobienstwo
//	E1=pot.get_energy(Vatoms[i])+pot.get_energy(vac_neighbour[k])-pot.get_energy(Vatoms[i],vac_neighbour[k]);
//	E2=pot.get_energy(Vatoms[i],atom)+pot.get_energy(vac_neighbour[k],0)-pot.get_energy(Vatoms[i],0,vac_neighbour[k],0)-pot.get_energy(Vatoms[i],atom,vac_neighbour[k],atom) + pot.get_energy(Vatoms[i],atom,vac_neighbour[k],0);
	if(MOVE){control_output<<" c: "<< ile;}
//	if(MOVE_FRAME or SINGLE){control_output<<"||"<<BLOKS[nr].size(j)<<"|"<<reservuars[rez].size(0)<<"|"<<reservuars[rez].size(j)<<"|"<<Vtoadd.size()<<">|";}
	if(MOVE){break;}	//wyjdz z petli for
//	}//for(i=1,2)
//	if(MOVE){break;}	//wyjdz z petli for
	}//for(j=delta)
	
	if(MOVE_FRAME or SINGLE){control_output<<endl;}
	if(MOVE){
		FLAG = true;		//set local FLAG in do_equi_vac
		do_equi_vac();		//rekurencja. Na poczatku sprawdza czy MOVE_FRAME set to TRUE.
	}
}


void opcja :: do_equi_rez(){
	
	//direct exchange dla rezerwuarow	//DODAC PARALLEL
	double beta=1.0/(kB*TEMPERATURE);
	for(unsigned int i=0;i<reservuars.size();i++){	
		long Nsize=reservuars[i].size();
	
	//	cout<<"rez: "<<i<<" size: "<<Nsize<<endl;
		for (long iter=0; iter<DIRECT_STEPS; iter++){
			long N1=(long)(ran01()*Nsize);
			long N2=(long)(ran01()*Nsize);		
			site *site1=0;	
			site *site2=0;
			site1=reservuars[i].get_site(N1);
			site2=reservuars[i].get_site(N2);
			int typ1 = site1->get_atom();
			int typ2 = site2->get_atom();
			double E1=POT.get_energy(site1)+POT.get_energy(site2)-POT.get_energy(site1,site2);
			site1->set_atom(typ2);
			site2->set_atom(typ1);
			double E2=POT.get_energy(site1)+POT.get_energy(site2)-POT.get_energy(site1,site2);		
			double dE=E2-E1,p=0.0,R=0.0;
			if(dE > 0 ){
				p = exp(-beta*dE);
				R = rnd();
				if(R >= p){
				site1->set_atom(typ1);
				site2->set_atom(typ2);
				}
			}
	//	ofstream dir_file("dir_tmp.dat", ios::app); 
	//	dir_file<<i<<" "<<iter<<" "<<E1<<"|-"<<E2<<"|="<<dE<<"|"<<p<<"<?"<<R<<endl;
		}
	}		
}

void opcja :: create_vac_new(int nr, int ile_at, bool &FLAG){

	if(ile_at<0){ile_at=ile_at*-1;}
	unsigned int ile=ile_at;
	if(MOVE_FRAME or SINGLE){control_output<<" c: "<< ile;}
	bool MOVE = false;
	int rez = -1, j = -1;
//	control_output<<"TRYB: "<<TRYB<<endl;
//	vector <int> wybrane_typy; wybrane_typy.reserve(20);
//	for( int i=0; i<(ile);i++){
//		j=choose_typ(BLOKS[nr]);
//		wybrane_typy.push_back(j);
//	}
//	control_output<<"c: "<<nr<<" "<<ile<<": ";
//	for(int i=0;i<wybrane_typy.size();i++){
//	control_output<<wybrane_typy[i];
//	}
//	control_output<<endl;

	for(unsigned int i=0; i<ile;i++){
//		for( int j=1;j<3;j++){		//MOZNA ZROBIC W ZALEZNOSCI OD TYPOW
			long N1,N2;
			site* rnd_vac=0;
			site* rnd_at=0;
			j=choose_typ(BLOKS[nr]);
	//		j=wybrane_typy[i];	//losuje typ atomu do wymiany z wakancja
	//		control_output<<"rozmiar typow "<<0<<" w bloku "<<bloks[nr].size(0);
	//		control_output<<" rozmiar typow "<<j<<" w bloku "<<bloks[nr].size(j);
			if(BLOKS[nr].size(j) <= 0){
				control_output<<endl;
				control_output<<"ERROR: in opcja::create_vac -> you want to remove \
				element "<<j<<" that does not exist in blok\n	\
				Probably error in opcja::init or opcja::reinit_bloks"<<endl; 
				cout<<endl;
				cout<<"ERROR: in opcja::create_vac -> you want to remove \
				element "<<j<<" that does not exist in blok\n	\
				Probably error in opcja::init or opcja::reinit_bloks"<<endl;exit(1);
			}
			
			N1=(long)(rnd()*(BLOKS[nr].size(j)));
			rnd_at=BLOKS[nr].get_site(j,N1);
	//		control_output<<" ctyp: "<<rnd_at->get_atom()<<" "<<nr<<endl;
			
			if(TRYB==2){
				int DIR = decide_direction(rnd_at);	//	-1 left;	+1 right
				vector <site*> migration_path; migration_path.reserve(2000);
				find_migration_path(rnd_at,DIR,migration_path);	
				dislocation_walk(migration_path);
			}
			if(TRYB==1){ //swap
				rez=choose_reservuar(rnd_at);
				MOVE = check_rezervuars(rez,0);	//sprawdz rezerwuwar pod katem dostepnych wakancji
		//		control_output<<"MOVE: "<<MOVE<<" rez: "<<rez<<endl;
				if (!MOVE){
					N2=(long)(rnd()*(reservuars[rez].size(0)));
					rnd_vac = reservuars[rez].get_site(0,N2);
					rnd_vac->set_atom(j);			
					reset_site(rnd_vac);
					reservuars[rez].delete_site(0,N2);
					reservuars[rez].add_site(j,rnd_vac);
					
					rnd_at->set_atom(0);
					reset_site(rnd_at);
					BLOKS[nr].delete_site(j,N1);
					BLOKS[nr].add_site(0,rnd_at);
					BLOKS[nr].prob_update(j,1);
					Vtoadd.push_back(rnd_at);	
				}else{
					//	if(MOVE_FRAME or SINGLE){control_output<<"||"<<BLOKS[nr].size(j)<<"|"<<reservuars[rez].size(0)<<"|"<<reservuars[rez].size(j)<<"|"<<Vtoadd.size()<<">|";}
					break;
				}
					
			}else if(TRYB==0){//convert
				rnd_at->set_atom(0);
				reset_site(rnd_at);
				BLOKS[nr].delete_site(j,N1);
				BLOKS[nr].add_site(0,rnd_at);
				BLOKS[nr].prob_update(j,1);
				Vtoadd.push_back(rnd_at);	
			}
	}//for(j=delta)
	
	if(MOVE_FRAME or SINGLE){control_output<<endl;}
	if(MOVE){
		FLAG = true;		//set local FLAG in do_equi_vac
		do_equi_vac();		//rekurencja. Na poczatku sprawdza czy MOVE_FRAME set to TRUE.
	}
}


int opcja :: choose_reservuar(site* atom){
		
	double d = 0.0, distanceL = 0.0, distanceR = 0.0;
	int rez = -1;
	d=atom->get_position(BIN_DIRECTION);			
	distanceL = abs(d - BIN_ST);
	distanceR = abs(d - BIN_END);
	//	cout<<"Absloute value of distance L/R: "<<distanceL<<"/"<<distanceR<<endl;
		
	if(distanceL > distanceR){
		rez = 1;
	}else if(distanceL < distanceR){
		rez = 0;
	}else if(distanceL == distanceR){
		double Q = rnd();
		if(Q<=0.5){rez = 0;}else{rez = 1;}
	}else{
		cout<<"Wrong direction number in opcja::convert() "<<distanceL<<"/"<<distanceR<<endl;
		exit(1);
	}
	return rez;
}

void opcja :: do_equi_vac(){

	bool LOCAL_MOVE = false;
	if(MOVE_FRAME or SINGLE){	
	control_output<<"Rownowarze: "<<LOCAL_MOVE<<"|"<<MOVE_FRAME<<endl;
	refresh(1);}
	
	if(MOVE_FRAME){
		move_frame();	//do_move resetuje global MOVE_FRAME
	}
	
//	find_matano_plane();

//	control_output<<"Rownowaga: "<<endl;
	for (unsigned int i=0;i<BLOKS.size();i++){
	//	bloks[i].calc_stech();
		if(MOVE_FRAME or SINGLE)
		{control_output<<i;}
		
		double stech=BLOKS[i].get_stech();
		double vac=BLOKS[i].get_vac();
		double size=BLOKS[i].size();
		int delta_vac = check_stech(stech,vac,size);	//zwracac ile wakancji remove/create
	//w zaleznosci od wyniku wykonac remove/create
		//delta_vac if < 0 then must remove vac in plaster
		//delta_vac if > 0 then must create vac in plaster
		
//		control_output<<i<<" "<<stech<<" "<<vac<<" "<<size<<" "<<delta_vac<<" "<<endl;
//lokalna zmienna MOVE zeby ominuac pozostale bloki
		
		if(delta_vac < -1){
			remove_vac_new(i,-delta_vac, LOCAL_MOVE);
		}
		else if (-1 <= delta_vac and delta_vac <= 1)
		{
			if(MOVE_FRAME or SINGLE){	
			control_output<<" do nothing"<<endl;}
		}
		else if ( delta_vac > 1 )
		{
			create_vac_new(i,delta_vac, LOCAL_MOVE);
		}
		else{
			cout<<"Error with equilibrate vacancy opcja::do_equi_vac()"<<endl;
			exit(1);
		}
		
		
		
		if(LOCAL_MOVE){
			control_output<<" po: "<<i<<"|"<<LOCAL_MOVE<<"|"<<MOVE_FRAME<<endl;
			break;}	//jesli bylo do_move przestan rownowazyc pozostale bloki
	}	//for bloks
	
//	refresh(1);
//	control_output<<"Koniec Rownowaga!"<<endl;
			if(!LOCAL_MOVE){
			if(MOVE_FRAME or SINGLE){
			control_output<<" po: all|"<<LOCAL_MOVE<<"|"<<MOVE_FRAME<<endl;}}
	//uaktulanic liste wakancji w  resident !!!!!!!!!!!! UWAGA !!!!!!!!!	
}

void opcja :: equilibrate(vector <site*> &kontener){
	
	refresh(0);
	do_equi_vac();
	do_equi_rez();
	refresh_sim_area(kontener);

	
}

void opcja :: flux_net_add(double pos_V, double pos_A, int typ, vector<plaster>& layer){
	int id_A=-2,id_V=-2;							
	//znajdz index bin dla x_A oraz x_V
	for(unsigned int i=0; i< layer.size(); i++){
		if(layer[i].get_st() <= pos_A  and pos_A < layer[i].get_end()){
			id_A=i;
		}
		if(layer[i].get_st() <= pos_V  and pos_V < layer[i].get_end()){
			id_V=i;
		}
	}
	if(id_V >=0 and (id_V != id_A)){
		//vakancja pojawia sie w plastrze id_V oraz atom znika z plastra id_V
		layer[id_V].flux_net_delta(typ, 0);
		layer[id_V].flux_net_delta(0, 1);
	}
	if(id_A >=0 and (id_V != id_A)){
		//vakancja znika z plastra id_A oraz atom pojawia sie w plastrze id_A
		layer[id_A].flux_net_delta(typ, 1);
		layer[id_A].flux_net_delta(0, 0);
	}	
}

/*
void opcja :: flux_add(double pos_V, double pos_A, int typ, vector<plaster>& layer){
	int id_A=-2,id_V=-2;
	double dir = pos_V - pos_A;	//if > 0 - wakancja skoczyla z lewej strony na prawa. Strumien dla wakancji jest +
								//													  Strumien dla atomu jest -
	//znajdz index bin dla x_A oraz x_V
	for(unsigned int i=0; i< layer.size(); i++){
		if(layer[i].get_st() <= pos_A  and pos_A < layer[i].get_end()){
			id_A=i;
		}
		if(layer[i].get_st() <= pos_V  and pos_V < layer[i].get_end()){
			id_V=i;
		}
	}
	if(id_V >=0 and (id_V != id_A)){
		//vakancja pojawia sie w plastrze id_V oraz atom znika z plastra id_V
	//	layer[id_V].flux_net_delta(typ, 0);
		layer[id_V].jump_occured();
		if(dir > 0){	//wakancja wplynela z lewej strony; atom wyplyna w lewo
			layer[id_V].flux_net_delta(0, 1);
			layer[id_V].flux_net_delta(typ, 0);
		}
		if(dir < 0){	//wakancja wplynela z prawej strony; atom wyplynal w prawo
			layer[id_V].eq_flux_delta(0, 0);
			layer[id_V].eq_flux_delta(typ, 1);
		}

	}
	if(id_A >=0 and (id_V != id_A)){
		//vakancja znika z plastra id_A oraz atom pojawia sie w plastrze id_A
	//	layer[id_A].flux_net_delta(typ, 1);
		layer[id_A].jump_occured();
		if(dir < 0){	//atom wplynal z lewej strony; wakancja wyplynela w lewo
			layer[id_A].flux_net_delta(typ, 1);
			layer[id_A].flux_net_delta(0, 0);
		}
		if(dir > 0){	//atom wplynal z prawej strony; wakancja wyplynela w prawo
			layer[id_A].eq_flux_delta(typ, 0);
			layer[id_A].eq_flux_delta(0, 1);
		}

	}	
}
*/

void opcja :: flux_add(site* VAC, site* ATO, vector<plaster>& layer){
	int id_A=-2,id_V=-2;				
	unsigned int DIR = (layer[0]).get_direction();
	double pos_V = VAC->get_position(DIR);
	double pos_A = ATO->get_position(DIR);
	double dir = pos_V - pos_A;											//if > 0 - wakancja skoczyla z lewej strony na prawa. Strumien dla wakancji jest +

	//znajdz index bin dla x_A oraz x_V
	id_V=VAC->get_hist_index();
	id_A=ATO->get_hist_index();

	int typ = ATO->get_atom();
//	for(unsigned int i=0; i< layer.size(); i++){
//		if(layer[i].get_st() <= pos_A  and pos_A < layer[i].get_end()){
//			id_A=i;
//		}
//		if(layer[i].get_st() <= pos_V  and pos_V < layer[i].get_end()){
//			id_V=i;
//		}
//	}

	if(id_V >=0 and (id_V != id_A)){
		//vakancja pojawia sie w plastrze id_V oraz atom znika z plastra id_V		
		if(dir > 0){	//wakancja wplynela z lewej strony; atom wyplyna w lewo
			for(int ID=id_V; ID>id_A; ){
				layer[ID].jump_occured_dislocation();
				layer[ID].flux_net_delta(0, 1);
				layer[ID].flux_net_delta(typ, 0);
				ID--;
			}
		}
		if(dir < 0){	//wakancja wplynela z prawej strony; atom wyplynal w prawo
			for(int ID=id_V; ID<id_A; ){
				layer[ID].jump_occured_dislocation();
				layer[ID].eq_flux_delta(0, 0);
				layer[ID].eq_flux_delta(typ, 1);
				ID++;
			}
		}	
	}
	if(id_A >=0 and (id_V != id_A)){
		//vakancja pojawia sie w plastrze id_V oraz atom znika z plastra id_V
		if(dir < 0){	//wakancja wplynela z lewej strony; atom wyplyna w lewo
			for(int ID=id_A ;ID>id_V ;){
				layer[ID].jump_occured_dislocation();
				layer[ID].flux_net_delta(typ, 1);
				layer[ID].flux_net_delta(0, 0);
				ID--;
			}
		}
		if(dir > 0){	//wakancja wplynela z prawej strony; atom wyplynal w prawo
			for(int ID=id_A ;ID<id_V ;){
				layer[ID].jump_occured_dislocation();
				layer[ID].eq_flux_delta(typ, 0);
				layer[ID].eq_flux_delta(0, 1);
				ID++;
			}
		}	
	}


}

void opcja :: flux_add_dislocation(site* VAC, site* ATO, vector<plaster>& layer){
	int id_A=-2,id_V=-2;				
	unsigned int DIR = (layer[0]).get_direction();
	double pos_V = VAC->get_position(DIR);
	double pos_A = ATO->get_position(DIR);
	double dir = pos_V - pos_A;											//if > 0 - wakancja skoczyla z lewej strony na prawa. Strumien dla wakancji jest +

	//znajdz index bin dla x_A oraz x_V
	id_V=VAC->get_hist_index();
	id_A=ATO->get_hist_index();

	int typ = ATO->get_atom();
//	for(unsigned int i=0; i< layer.size(); i++){
//		if(layer[i].get_st() <= pos_A  and pos_A < layer[i].get_end()){
//			id_A=i;
//		}
//		if(layer[i].get_st() <= pos_V  and pos_V < layer[i].get_end()){
//			id_V=i;
//		}
//	}

	if(id_V >=0 and (id_V != id_A)){
		//vakancja pojawia sie w plastrze id_V oraz atom znika z plastra id_V		
		if(dir > 0){	//wakancja wplynela z lewej strony; atom wyplyna w lewo
			for(int ID=id_V; ID>id_A; ){
				layer[ID].jump_occured_dislocation();
				layer[ID].prob_hist_l(0, 1);
				layer[ID].prob_hist_l(typ, 0);
				ID--;
			}
		}
		if(dir < 0){	//wakancja wplynela z prawej strony; atom wyplynal w prawo
			for(int ID=id_V; ID<id_A; ){
				layer[ID].jump_occured_dislocation();
				layer[ID].prob_hist_r(0, 0);
				layer[ID].prob_hist_r(typ, 1);
				ID++;
			}
		}	
	}
	if(id_A >=0 and (id_V != id_A)){
		//vakancja pojawia sie w plastrze id_V oraz atom znika z plastra id_V
		if(dir < 0){	//wakancja wplynela z lewej strony; atom wyplyna w lewo
			for(int ID=id_A ;ID>id_V ;){
				layer[ID].jump_occured_dislocation();
				layer[ID].prob_hist_l(typ, 1);
				layer[ID].prob_hist_l(0, 0);
				ID--;
			}
		}
		if(dir > 0){	//wakancja wplynela z prawej strony; atom wyplynal w prawo
			for(int ID=id_A ;ID<id_V ;){
				layer[ID].jump_occured_dislocation();
				layer[ID].prob_hist_r(typ, 0);
				layer[ID].prob_hist_r(0, 1);
				ID++;
			}
		}	
	}


}

void opcja :: find_matano_plane(){
	double x0 = BIN_ST + fabs(  (BIN_ST-BIN_END)/2.0 ) ;
//	cout<<BIN_DIRECTION<<" "<<BIN_ST<<" "<<BIN_END<<" "<< x0<<endl;


	vector <double> X;
	vector <double> Y;
	
	double SUM=0;
	typedef vector <plaster>::iterator it2pl;
	for( it2pl IT=HIST.begin();IT!=HIST.end();++IT){
		double y = IT->flux_net_get(0);
		double x1=IT->get_st();
		double x2=IT->get_end();
		double x= (x2+x1)/2.0;
//		cout<<"OPCJA: "<<x1<<" "<<x2<<" "<<x<<" "<<y<<endl;
		X.push_back(x);
		Y.push_back(y);
	}

		
	integral_data(X, Y, x0);
	
	
	
}

void opcja :: call_flux(site* vac_after_jump,site* atom_after_jump){

	double x_A = atom_after_jump->get_position(BIN_DIRECTION);
	double x_V = vac_after_jump->get_position(BIN_DIRECTION);
	int typ=atom_after_jump->get_atom();
	flux_net_add(x_V,x_A,typ,BLOKS);
	flux_net_add(x_V,x_A,typ,reservuars);
	flux_add(vac_after_jump,atom_after_jump,HIST);
}

void opcja :: call_flux_dislocation(site* vac_after_jump,site* atom_after_jump){
	flux_add_dislocation(vac_after_jump,atom_after_jump,HIST);	
}

void opcja :: set_opcja_lattice(lattice *sample){
	
	SAMPLE=sample;
	BIN_ATOMS_TYP=sample->get_atom_typ_numbers();
}

void opcja :: init_EQ(vector <double> &parameters ){
	
	if(parameters.size()!=7){cout<<"ERROR in opcja::init. Wrong parameters list in conf.in"<<endl;exit(1);}

	//dla -1: tylko monitoruj	+biny -rezerwuary +flux	+flux_eq -do_equi
	//dla 0 nie rob nic:		-biny -rezerwuary -flux -flux_eq -do_equi
	//dla 1 do_eq				+biny +rezerwuary +flux +flux_eq +do_equi
	int tryb = parameters[0];
	if(tryb > 0){		
	control_output<<"Init equilibrations..."<<endl;


	int tryb = parameters[0];

	int co_ile = parameters[1];
	double od_kod = parameters[2];
	double do_kod = parameters[3]; 
	double co_ile_bin = parameters[4]; 
	int bin_direction = parameters[5]; 
	int if_move = parameters[6];

	TRYB=tryb;	
	EQ_STEP=co_ile;			
	BIN_ST=od_kod;
	BIN_END=do_kod;
	BIN_SIZE=co_ile_bin;
	BIN_DIRECTION=bin_direction;
	ST_VOL=BIN_ST;
	END_VOL=BIN_END;
	
	MOVE_SIM_REGION = (if_move != 0);
	wektor a(0.0,0.0,0.0);
	del_L_sim=a;
	del_R_sim=a;
	
	cout<<"Adres potencjal w opcja: "<<&POT<<endl;
	cout<<"Adres bariery w opcja: "<<&BARRIERS<<endl;
	cout<<"Adres sample w opcja: "<<SAMPLE<<endl;
						
	build_bins(BLOKS,"block");				//zbuduj puste plasterki i zainicjuj parametry
	if(SAMPLE==0){control_output<<"ERROR in opcja::init_EQ. SAMPLE adress not set -> opcja::set_opcja_lattice"<<endl;exit(1);}
	SAMPLE->get_sites(BLOKS);		//wywolac funkcje z lattice ktora wczyta sity do plasterkow
//	show();
	
	//przeliczam stech i vac dla plasterkow
	for (unsigned int i=0; i < BLOKS.size(); i++){
		control_output<<"do init_calc in blok: "<<i<<" "; 
		BLOKS[i].init_calc(1);
	}
	//wczytuje stech z pliku
	if(EQ_STEP>0){
		string file_name="stech_curve.in";
		read_file(file_name);
		}
	}	
}


void opcja :: build_bins(vector<plaster>& layer, string name){	
	layer.clear();
	int ile = (int((BIN_END - BIN_ST)/BIN_SIZE));
	double miarka[ile];
	miarka[0]=BIN_ST;
	for(int ii=1; ii<=ile; ii++){miarka[ii] = miarka[ii-1] + BIN_SIZE;}
	//for(int ij=0; ij<ile; ij++){control_output<<ij<<" "<<miarka[ij]<<endl;}

	for (int i=0;i<ile;i++)
	{
		plaster tmp(2000,BIN_ATOMS_TYP,BIN_DIRECTION,i,miarka[i],miarka[i+1],name);
		layer.push_back(tmp);
	}
//	control_output<<ile<<" "<<layer.size()<<endl;
}	


void opcja :: init_reservuar(vector <double> &parameters){
	

	if(SAMPLE == NULL){
		control_output<<"Error: pointer do sample not initialized -> opcja :: init_reservuar"<<endl;
		exit(0);
	}
	
	if(parameters.size()!=5){control_output<<"ERROR in opcja::init_reservuars. Wrong parameters list in conf.in"<<endl;exit(1);}
	
	//liczba sitow w plastrze komorek elementarnych
	double LATT = SAMPLE->get_latice_const(BIN_DIRECTION, 0);
	plaster tmp_rozmiar(2000,BIN_ATOMS_TYP,BIN_DIRECTION,0,0,LATT, "tmp");
	SAMPLE->get_sites(tmp_rozmiar);
	long rozmiar = (tmp_rozmiar.size()/2)*(parameters[2]);
	ROZMIAR.push_back(rozmiar);				//liczba sitow o ktore moze zmienic sie rezerwuar (domyslnie jedna plaszczyzna)

	double od_kod = parameters[0]; 
	double width = parameters[1];
	double o_ile = parameters[2];		//-> init_reservuars [2] to move frame lenght
	int direction =parameters[3];
	long direct_step = parameters[4];
	double do_kod = od_kod + width;
	
	DIRECT_STEPS=direct_step;
	if( !(DIRECT_STEPS > 0) ){control_output<<"WARRNING in opcja::init_reservuar. DIRECT_STEP equal to: "<<DIRECT_STEPS<<endl;}

	vector <double> par;
	par=parameters;
	reservuars_par.push_back(par);

	plaster tmp(2000,BIN_ATOMS_TYP,BIN_DIRECTION,(ROZMIAR.size()-1),od_kod,do_kod, "rez");
	SAMPLE->get_sites(tmp);
	control_output<<"do init_calc in reservuar: "; 
	tmp.init_calc(1);
	reservuars.push_back(tmp);
	
	ST_VOL=min(od_kod,ST_VOL);
	END_VOL=max(do_kod,END_VOL);

}


void opcja :: move_frame(){
	
	control_output<<"Przesowam ramke: "<<REZ_TO_MOVE<<"|"<<TYP_TO_MOVE<<endl;
//	cout<<"Przesowam ramke, press any kay..."<<endl;
//	int o;
//	cin>>o;
	if(REZ_TO_MOVE<0){
		control_output<<"ERROR: opcja:move_frame() REZ: "<<REZ_TO_MOVE<<endl;
		exit(1);
		}
	if(TYP_TO_MOVE==0){
		control_output<<"ERROR: opcja:move_frame() TYP: "<<TYP_TO_MOVE<<endl;
		exit(1);
		}
	save_write();
	reinit_reservuars(REZ_TO_MOVE,TYP_TO_MOVE);
	reinit_bloks();

	//reset global FLAGS	REZ_TO_MOVE is reset in refresh_sim_area
	
	REZ_TO_MOVE = -1; 	
	TYP_TO_MOVE = 0;
	
};

void opcja :: refresh(int on){

	//przeliczam stech i vac dla plasterkow
	for (unsigned int i=0; i < BLOKS.size(); i++){
				BLOKS[i].init_calc(on);
	}
	
	for (unsigned int i=0; i < reservuars.size(); i++){
			reservuars[i].init_calc(on);

	}
	
	for (unsigned int i=0; i < HIST.size(); i++){
			HIST[i].init_calc(on);

	}
	
}



void opcja :: reinit_reservuars(int nr, int typ){
		
	double odkod = reservuars_par[nr][0]; 
	double stwidth = reservuars_par[nr][1]; 
	double dokod = odkod + stwidth;
	
	double ile = reservuars_par[nr][2];		//-> init_reservuars [2] to move frame lenght
	int direction =int(reservuars_par[nr][3]);
	if(direction < 1 or direction > 4){control_output<<"ERROR: wrong direction: "<< direction<<" must be 1,2,3"<<endl; exit(0);}
	double s[3]={0.0,0.0,0.0};
	double e[3]={0.0,0.0,0.0};
	
	//Identyfikuje po ktorej stronie lezy rezerwuwar
	int POZ = 0;	//-1 lewa ;+1 prawa
	if(odkod < BIN_ST){ POZ=1;}	//to rez jest po lewej. Atomy sa po prawej, a wakancje po lewej.
	if(odkod >= BIN_END){ POZ=-1;}	//to rez jest po prawej. Atomy sa po lewej, a wakancje po prawej.
		
	//W ktora strone chcesz przesunac? do atomow czy do wakancji?	
	int DIR = 0;	
	if(typ < 0){ DIR=-1;}	//Chce w strone wakancji
	if(typ > 0){ DIR=1;}	//Chce w strone atomow
	///////////////////
	//		POZ	 POZ
	//	DIR  11  1-1	
	//	DIR -11 -1-1
	///////////////////
	
	//Przyjeto staly width rezerwuaru
	if(POZ==1 and DIR == 1){	//jestem po lewej i chce do atomow (prawa scianke ruszam)
	reservuars_par[nr][0] = reservuars_par[nr][0] + ile*POZ*DIR;
	//Blocks przesuwam
	BIN_ST += ile*POZ*DIR;	//przesowam biny poczatek w prawo
	s[(direction-1)] += ile*POZ*DIR;
	}
	if(POZ==1 and DIR == -1){	//jestem po lewej i chce do wakancji (lewa scianke ruszam)
	reservuars_par[nr][0] = reservuars_par[nr][0] + ile*POZ*DIR;	//WARUNEK BRZEGOWY!!!!!!!!!!
	//Blocks przesuwam
	BIN_ST += ile*POZ*DIR;	//przesowam biny poczatek w prawo
	s[(direction-1)] += ile*POZ*DIR;
	}

	if(POZ==-1 and DIR == -1){	//jestem po prawej i chce do wakancji (prawa scianke ruszam)
	reservuars_par[nr][0] = reservuars_par[nr][0] + ile*POZ*DIR;	//WARUNEK BRZEGOWY!!!!!!!!!!
	BIN_END += ile*POZ*DIR;
	e[(direction-1)] += ile*POZ*DIR;
	}
	if(POZ==-1 and DIR == 1){	//jestem po prawej i chce do atomow (lewa scianke ruszam)
	reservuars_par[nr][0] = reservuars_par[nr][0] + ile*POZ*DIR;
	//Blocks przesuwam
	BIN_END += ile*POZ*DIR;	//przesowam biny koniec w lewo
	e[(direction-1)] += ile*POZ*DIR;
	}

//nadpisane	
	double od_kod = reservuars_par[nr][0]; 
	double width = reservuars_par[nr][1];
	double do_kod = od_kod + width;
//	int typy_atomow=SAMPLE->get_atom_typ_numbers();
	ST_VOL=min(od_kod,ST_VOL);
	END_VOL=max(do_kod,END_VOL);
	
	plaster tmp(2000,BIN_ATOMS_TYP,BIN_DIRECTION,nr,od_kod,do_kod, "rez");
	SAMPLE->get_sites(tmp);
	control_output<<"reinit reservuar: "<<nr<<" "<<odkod<<"|"<<stwidth<<"|"<<dokod<<"||"<<od_kod<<"|"<<width<<"|"<<do_kod<<" "<<direction<<endl;
	control_output<<"do init_calc in reservuar: "; 
	//przerzuc zostawione atomy do nowego rezerwuaru (swap atom<->wakancja)
	if(DIR==1){
	tmp.init_calc(1);
	tmp.swap(reservuars[nr],1);
	}
	tmp.init_calc(1);	
	reservuars[nr]=tmp;

	wektor a(s[0],s[1],s[2]);
	wektor b(e[0],e[1],e[2]);
	del_L_sim += a;
	del_R_sim += b;

			
}

void opcja :: reinit_bloks(){
	build_bins(BLOKS);				
	SAMPLE->get_sites(BLOKS);

	control_output<<"do init_calc in blok: "<<endl;
	for (unsigned int i=0; i < BLOKS.size(); i++){
		control_output<<i<<" "; 
		BLOKS[i].init_calc(1);
	}
	//show();
}

void opcja :: remove_vac(int b, int vac, bool &FLAG){
	
	if(MOVE_FRAME or SINGLE)
	{control_output<<" r: "<< vac;}
	
	bool MOVE = false;
	int rez = -1, j=-1;
	
	
	for( int i=0; i<(vac);i++){
//		for( int j=1;j<3;j++){		//MOZNA ZROBIC W ZALEZNOSCI OD TYPOW
			long N1,N2;
			site* rnd_vac=0;
			site* rnd_at=0;
			j = choose_typ(BLOKS[b]);	//losuje typ atomu do wymiany z wakancja z blokow

			if(BLOKS[b].size(0) <= 0){
				control_output<<endl;
				control_output<<"ERROR: in opcja::remove_vac -> you want to remove \
element "<<0<<" that does not exist in blok\n	\
Probably error in opcja::init or opcja::reinit_bloks"<<endl; 
				cout<<endl;
				cout<<"ERROR: in opcja::remove_vac -> you want to remove \
element "<<0<<" that does not exist in blok\n	\
Probably error in opcja::init or opcja::reinit_bloks"<<endl;exit(0);
			}			
			
			while(1){
			N1=(long)(rnd()*(BLOKS[b].size(0)));
			
			rnd_vac=BLOKS[b].get_site(0,N1);

	//		control_output<<" jej adres: "<<rnd_vac;//<<" "<<N2<<" "<<rnd_vac2;


			double d = 0.0, distanceL = 0.0, distanceR = 0.0;
			rez = 0;

			if(BIN_DIRECTION==1)
			{
				d=rnd_vac->get_x();
			}
			else if(BIN_DIRECTION==2)
			{
				d=rnd_vac->get_y();
			}
			else if(BIN_DIRECTION==3)
			{
				d=rnd_vac->get_z();
			}
			else
			{
				cout<<"Wrong direction number in opcja::convert() "<<d<<endl;
				exit(1);
			}
			
			distanceL = abs(d - BIN_ST);
			distanceR = abs(d - BIN_END);
	//		cout<<"Absloute value of distance L/R: "<<distanceL<<"/"<<distanceR<<endl;

//PENTLA FOR PO REZERWUARACH - funkcja czyta pozycje plastra			
			if(distanceL > distanceR){
				rez = 1;
			}
			else if(distanceL < distanceR){
				rez = 0;
			}
			else if(distanceL == distanceR){
				double Q = rnd();
				if(Q<=0.5){rez = 0;}else{rez = 1;}
			}
			else{
				cout<<"Wrong direction number in opcja::convert() "<<distanceL<<"/"<<distanceR<<endl;
				exit(1);
			}
			
			MOVE = check_rezervuars(rez,j);	//sprawdz rezerwuwar
			
			if (!MOVE){

			
			N2=(long)(ran01()*(reservuars[rez].size(j)));
//			control_output<<" nr wylosowanego situ z listy "<<N2;
			rnd_at = reservuars[rez].get_site(j,N2);


	//		double old_E=POT.get_energy(rnd_at);

//			rnd_at->set_atom(0);
				
//			double new_E=POT.get_energy(rnd_at) + POT.get_energy(rnd_vac);
	//		double new_E=POT.get_energy(rnd_vac);

	//		double beta=1.0/(kB*T);
	//		double mi=0.0;
	//	    double P1=exp(beta*(mi-(new_E-old_E)));
	//		control_output<<" n: "<<new_E<<" o: "<<old_E<<" "<<(new_E-old_E)<<" "<<beta<<" P: "<<P1<<endl;
	//	//zmien z powrotem na old_typ jesli zdarzenie to nie zostalo trafione rnd()
	//		if(P1<rnd())
	//		{
	//		rnd_at->set_atom(j);	
//			rnd_at->set_atom(j);
	//		}
	//		else{
	
			rnd_at->set_atom(0);
			rnd_vac->set_atom(j);
	
			rnd_vac->reset_site();	
			rnd_at->reset_site();
												//kasujemy z listy typow vakancje w bloku
			BLOKS[b].delete_site(0,N1);
			BLOKS[b].add_site(j,rnd_vac);
	
												//lista wakancji zawiera atomy -> wymaga odswiezenia 
			reservuars[rez].delete_site(j,N2);
			reservuars[rez].add_site(0,rnd_at);
			
	//		Vtoadd.push_back(rnd_at);		// po za obszarem symulacji pojawila sie nowa wakancja
//			a co gdy obszar symulacji obejmuje rezerwuary?	!!!!!	NIE MA TAKIEJ OPCJI !!!!!!!
//			na wszelki wypadek dodam wszystkie nowe wakancje
// 			w refresh Vatoms sprawdze ktora wakancji z Vtoadd nalezy do obszaru symulacji
//			i te dodam do listy

//			a co z wakancjami ktore zostaly zastapion atomami?
//			nic, w f. refresh Vatoms sprawdze ktore sity to atomy i je usune z listy
			
//			control_output<<" adres situ "<<rnd_at<<endl;
			break;
		//	}
		}else{
			break;	//wyjdz z petli while
			}
		}	//koniec while
	//policz prawdopodobienstwo
//	E1=pot.get_energy(Vatoms[i])+pot.get_energy(vac_neighbour[k])-pot.get_energy(Vatoms[i],vac_neighbour[k]);
//	E2=pot.get_energy(Vatoms[i],atom)+pot.get_energy(vac_neighbour[k],0)-pot.get_energy(Vatoms[i],0,vac_neighbour[k],0)-pot.get_energy(Vatoms[i],atom,vac_neighbour[k],atom) + pot.get_energy(Vatoms[i],atom,vac_neighbour[k],0);

//		if(MOVE_FRAME or SINGLE){control_output<<"||"<<BLOKS[b].size(0)<<"|"<<reservuars[rez].size(0)<<"|"<<reservuars[rez].size(j)<<"|"<<Vtoadd.size()<<">|";}
		
		if(MOVE){break;}	
	//	}
	//	if(MOVE){break;}	
	}//wyszedlem z wszystkich pentli
	if(MOVE_FRAME or SINGLE){control_output<<endl;}
	if(MOVE){
		
		FLAG = true;		//set local FLAG in do_equi_vac
		do_equi_vac();		//rekurencja. Na poczatku sprawdza czy MOVE_FRAME set to TRUE.
	}
}

bool opcja :: check_x_belonging_volume(double x){
				   
	if( (set_prec(x) >= set_prec(get_start_volume()) ) and ( set_prec(x) < set_prec(get_end_volume()) ) ){
				return true;
	}		
	return false;
}

int opcja :: decide_direction(site *node){	//do zmiany na kryterium granicy faz			UWAGA UWAGA UWAGA		<----------------------------------------!!!!!!!!!!!!!!!!!!!

	//pobrac strumien lewy i prawy
	unsigned int id=node->get_hist_index();
	double left_area = -HIST[id].net_flux_get(0);
	double right_area = -HIST[id].eq_flux_get(0);

	//policz prace dla skoku w lewo/prawo
	double W_left = -left_area/Actual_MCtime;
	double W_right = right_area/Actual_MCtime;
	if(Actual_MCtime==0){
		W_left = -left_area/1;
		W_right = right_area/1;
	}

	W_left *= 0.5;
	W_right *= 0.5;

	
	//policz prawdopodobienstwo
	double p_left = exp(W_left/(kB*TEMPERATURE));
	double p_right = exp(W_right/(kB*TEMPERATURE));
	double SUM = p_left + p_right;
	
	double shot = rnd()*SUM;

	vector <double> target;
	target.clear(); target.reserve(10);
	target.push_back(0.0);
	target.push_back(p_left);
	target.push_back(SUM);


	int move = 0;
	int event = -1;
	control_output<<endl;
	for (int i = 1; i < (target.size()); i++){
		control_output<<" "<<W_left<<" "<<W_right<<" "<<i<<" "<<target[i-1]<<" "<<target[i]<<" "<<shot<<" "<<kB*TEMPERATURE<<endl;
		if( (target[i-1] <= shot) and (shot < target[i]) ){
			control_output<<"take event: "<<i<<" "<<target[i-1]<<" "<<target[i]<<" "<<shot<<endl;
			//break;
			event=i;
		}
	}
	
	if(event == 1){move = -1;}
	else if (event == 2){move = 1;}
	else{
		control_output<<"ERROR in opcja::decide_direction(): "<<event<<" "<<p_left<<" "<<p_right<<" "<<shot<<endl; exit(1);
	} 
	
	return move;

}

void opcja :: cal_angles(site *node, wektor &main, vector <site*> &wynik_at, vector <site*> &wynik_vac){
	
	vector <site*> neighs;													//dla noda pobierz sasiadow
	vector <site*>::iterator atom;
	int dir = get_direction();
	
//	node->show_site();
	wektor r0 = node->get_position();
	double x0 = r0[dir];
	bool IN_VOLUME = check_x_belonging_volume(r0[dir]);

	if(IN_VOLUME){

	bool IN_REZ = false, ON_WALL = false;	
	vector <plaster>::iterator it2REZ;
	for( it2REZ = reservuars.begin(); it2REZ != reservuars.end(); ++it2REZ){
		double start = it2REZ->get_st();
		double koniec = it2REZ->get_end();
		if( x0 > start and x0 < koniec ){								//node is in reservuar
			IN_REZ=true;
			break;
		}
		if( x0 == get_start_volume() or x0 == get_end_volume() ){								//node is in reservuar
			IN_REZ=true;
			break;
		}
	}
	
	node->read_site_neighbours(neighs, 1, 0);
	double bufor=0;
	for(atom=neighs.begin();atom!=neighs.end();++atom){

		wektor r1 = (*atom)->get_position();							//	r1.show();	
		wektor r = r1 - r0;
		wektor PB = ( r > ( SAMPLE->get_PB()*0.5 ) );					//warunki brzegowe
		
		r = r - ( PB * ( r.wersor() ) * ( SAMPLE->get_PB() ) );			//	r.show();	PB.show();
		r1 = r0 + r;
	//	wektor chPB = ( r.wersor() ) * ( SAMPLE->get_PB() );	chPB.show();
	//	wektor addPB = PB * chPB;	addPB.show();
	//	r = r - addPB;	r.show();
	//	r0.show();r1.show();r.show();
		double x0 = r1[dir];
		bool IN_VOLUME = check_x_belonging_volume(r1[dir]);

		if(IN_VOLUME){


		double cosL = main.cos_AB(r);
																		//	control_output<<"kat: "<<cosL<<" "<<rad2degree(cosL)<<endl;
		if( !IN_REZ ){
			if(cosL > 0){														// biore katy tylko zgodne z kierunkiem ruchu dyslokacji		//wybierz skok z katem najblizszym 0. cos -> 1
				int typ = (*atom)->get_atom();
				if( typ>0 ){														//biore tylko atomy
					if(cosL > bufor){											//biore tylko atomy z najmniejszym katem
						bufor = cosL;
						wynik_at.clear();
						wynik_at.push_back( (*atom) );
					}else if(cosL == bufor){									//jesli sa dwa o takim samym kacie?? //dodaj do listy 
						wynik_at.push_back( (*atom) );
					}
				}else{
					wynik_vac.push_back( (*atom) );
				}
			}
		}else{
			if(cosL > 0){
				wynik_at.push_back( (*atom) );
			}
			wynik_vac.clear();
		}
	}
		
	}		//koniec for neighs

	}else{																	//jesli site jest poza VOLUME (sim_area + reservuars), do not count it.
		control_output<<"ERROR: opcja:: cal_angle(). site outside Volume"<<endl;
		node->show_site(); exit(1);
	}
	
//	control_output<<"After cal_angle: "<<wynik_at.size()<<" "<<wynik_vac.size() <<endl;

}

void opcja :: cal_angles_strong(site *node, wektor &main, vector <site*> &wynik_at, vector <site*> &wynik_vac){
	
	wynik_at.clear();
	wynik_vac.clear();
	int dir = get_direction();
//	node->show_site();
	wektor r0 = node->get_position();
	double x0 = r0[dir];
	bool IN_VOLUME = check_x_belonging_volume(r0[dir]);

	if(IN_VOLUME){

	bool IN_REZ = false, ON_WALL = false;	
	vector <plaster>::iterator it2REZ;
	for( it2REZ = reservuars.begin(); it2REZ != reservuars.end(); ++it2REZ){
		double start = it2REZ->get_st();
		double koniec = it2REZ->get_end();
		if( x0 >= start and x0 <= koniec ){								//node is in reservuar
			IN_REZ=true;
			break;
		}	}
	
	vector <site*> neighs;													//dla noda pobierz sasiadow
	vector <site*>::iterator atom;	
	node->read_site_neighbours(neighs, 1, 0);
	double bufor=0;
	for(atom=neighs.begin();atom!=neighs.end();++atom){

		wektor r1 = (*atom)->get_position();							//	r1.show();	
		wektor r = r1 - r0;
		wektor PB = ( r > ( SAMPLE->get_PB()*0.5 ) );					//warunki brzegowe
		
		r = r - ( PB * ( r.wersor() ) * ( SAMPLE->get_PB() ) );			//	r.show();	PB.show();
		r1 = r0 + r;
	//	wektor chPB = ( r.wersor() ) * ( SAMPLE->get_PB() );	chPB.show();
	//	wektor addPB = PB * chPB;	addPB.show();
	//	r = r - addPB;	r.show();
	//	r0.show();r1.show();r.show();
		double x0 = r1[dir];
		bool IN_VOLUME = check_x_belonging_volume(r1[dir]);

		if(IN_VOLUME){

		double cosL = main.cos_AB(r);
																		//	control_output<<"kat: "<<cosL<<" "<<rad2degree(cosL)<<endl;
		if( !IN_REZ ){
			if(cosL >= 0){														// biore katy tylko zgodne z kierunkiem ruchu dyslokacji		//wybierz skok z katem najblizszym 0. cos -> 1
				int typ = (*atom)->get_atom();
				if( typ>0 ){														//biore tylko atomy
					if(cosL > bufor){											//biore tylko atomy z najmniejszym katem
						bufor = cosL;
						wynik_at.clear();
						wynik_at.push_back( (*atom) );
					}else if(cosL == bufor){									//jesli sa dwa o takim samym kacie?? //dodaj do listy 
						wynik_at.push_back( (*atom) );
					}
				}else{
					wynik_vac.push_back( (*atom) );
				}
			}
		}else{
			if(cosL >= 0){
				wynik_at.push_back( (*atom) );
			}
			wynik_vac.clear();
		}
	}
		
	}		//koniec for neighs

	}else{																	//jesli site jest poza VOLUME (sim_area + reservuars), do not count it.
		control_output<<"ERROR: opcja:: cal_angle(). site outside Volume"<<endl;
		node->show_site(); exit(1);
	}
	
	//control_output<<"After cal_angle_strong: "<<wynik_at.size()<<" "<<wynik_vac.size() <<endl;
}


void opcja :: find_migration_path(site *first_node,int DIR, vector <site*> &migration_path){
	
	wektor kierunek;
	if(DIR==1){ kierunek(1.0,0.0,0.0);}else if(DIR== -1){kierunek(-1.0,0.0,0.0);}
	else{control_output<<"ERROR in opcja""find_migration_path(). Wrong direction: "<<DIR<<endl;exit(1);}
	kierunek.show();
	
	site* node=first_node;
	migration_path.push_back( first_node );
	unsigned int TRY_WALL = 0;
	bool WALL = false;
	do{		
		vector <site*> site_bufor;
		vector <site*> vac_bufor;
		if( !WALL ){
			cal_angles(node, kierunek, site_bufor, vac_bufor);
		}else{
			wektor oposite_kier = kierunek*(-1.0);	
			cal_angles_strong(node, oposite_kier, site_bufor, vac_bufor);
		}
		
		if(site_bufor.size() > 0 ){											// jesli jestes w sim_area to masz tylko atomy w site_buf oraz tylko vac w vac_buf
			int rndIndex = rand() % site_bufor.size();						// jesli jestes w reservuar to masz atomy i vacancy w site i pusty vac_buf
			int typ = (site_bufor[rndIndex])->get_atom();					//losuj sita
			if(typ==0 ){													//jesli jestes w rezerwuar i jesli vacanc to zakoncz szukanie siezki
				if(first_node->get_atom()>0){								
					migration_path.push_back( (site_bufor[rndIndex]) );
		//			(site_bufor[rndIndex])->show_site();					
				}											
				node=site_bufor[rndIndex];									//w node jest vacans stop while
				break;
			}
			else if (typ > 0){												//jesli atom to powtorz procedure dla nowo wybranego situ				
				migration_path.push_back( (site_bufor[rndIndex]) );
				node=site_bufor[rndIndex];										
			}
			else{	control_output<<"ERROR in opcja::find_migration_path(). \
				Atom type < 0: "<<typ<<endl;exit(1);}
			TRY_WALL = 0; WALL = false;
		}
		else if(site_bufor.size() == 0 and vac_bufor.size()>0){				//Jesli w sim_area natrafisz na faze wakancyjna, por, klaster wakancji, to zakoncz sciezke.
			int rndIndex = rand() % vac_bufor.size();						
			int typ = (vac_bufor[rndIndex])->get_atom();
			if(typ==0){
				if(first_node->get_atom()>0){								
					migration_path.push_back( (vac_bufor[rndIndex]) );		//zapisz adres vacancy do sciezki jesli poczatek byl atom. Create vac.
		//			(vac_bufor[rndIndex])->show_site();
				}
				node=vac_bufor[rndIndex];									//w node jest vacans stop while	
				break;
			}
			else{
				control_output<<"ERROR in opcja::find_migration_path(). \
				Wrong atom type: "<<typ<<endl;exit(1);
			}
			TRY_WALL = 0; WALL = false;
		}								
		else if(site_bufor.size() == 0 and vac_bufor.size() == 0){			//Jesli path doszla do sciany?? Nie ma ani atomu ani wakancji do skoku. TO odbij wektor kierunek.
		//	control_output<<"WARRNING in opcja::find_migration_path(). \
			Path reached a wall. No atoms or vacancy available to contiune path."<<endl;//exit(1);		
			WALL = true;TRY_WALL++;
			if(TRY_WALL>3){
				control_output<<"ERROR in opcja::find_migration_path(). \
				Node is oscilating: "<<TRY_WALL<<" "<<site_bufor.size()<<" "<<vac_bufor.size() <<endl;
				exit(1);
				}
		}else{
			control_output<<"ERROR in opcja::find_migration_path(). \
			Undefined condition: "<<site_bufor.size()<<" "<<vac_bufor.size() <<endl;
			exit(1);
		}	
		
	//migration_path.back()->show_site();	
	}while(node->get_atom()>0);		
	
	//migration_path.back()->show_site();	
	//node->show_site();
																		//UWAGA: na perkolacje!					-> odbija wektor kierunek gdy dojdzie do sciany reserwuaru.
																		//UWAGA: A co gdy mam juz wszedzie rownowage i w jednym miejscu musze dodac/odjac wakancje
																		// +/- jeden rodzaj atomu w kazdym plasterku nic nie zmienia (wciaz w granicy sredniej stechiometrii)
																		// +/- jedna wakancja w danym plasterku też nic nie zmienia (wciaz w granicy sredniej koncentracji)
																		//Aktualnie algorytm omija wakancje w sim_area. Chyba ze zabraknie mu mozliwosci, to wybiera wakancje.
}

void opcja :: dislocation_walk(vector <site*> &path){
	
	control_output<<"Path size: "<<path.size()<<endl;
	site* first = path.front();
	site* last = path.back();
	
	vector <site*>::iterator point;
	//first->show_site();
	//last->show_site();

	//sciezka ma na poczatku wakancje (na koncu atom) -> remove wakancja
	if( first->get_atom() == 0 ){	//forrward vacancy out
		vector<site*>::iterator i = path.begin(); ++i;
		for ( ; i != path.end(); ++i ){
			vector<site*>::iterator prev = i; --prev;
			virtual_jump_vac_atom( (*prev), *i);	
			(*prev)->show_site();
			(*i)->show_site();
		}	
		reset_site( *(--i) );
		Vtoadd.push_back( *i );	
		(*i)->show_site();

	}else if( last->get_atom()== 0 ){	//backward vacancy in
		vector<site*>::reverse_iterator i = path.rbegin(); ++i;
		for ( ; i != path.rend(); ++i ){
			vector<site*>::reverse_iterator prev = i; --prev;
			virtual_jump_vac_atom( (*prev), *i);
			(*prev)->show_site();
			(*i)->show_site();
		}	
		reset_site( *(--i) );
		Vtoadd.push_back( *i );	
		(*i)->show_site();
	}else{
		control_output<<"ERROR in opcja::dislocation_walk()"<<endl; exit(1);
	}
}

void opcja :: virtual_jump_vac_atom( site* VAC, site* ATOM){

	//one direction exchange of sites. Virtual jump of vac to atom.
//	control_output<<"Przed:"<<endl;
//	VAC->show_site();
//	ATOM->show_site();
	if(VAC->get_atom() != 0){control_output<<"ERROR in lattice::exchange_sites(). Type of vacancy not 0: "<<VAC->get_atom()<<endl;exit(1);}
	if(ATOM->get_atom() <= 0){control_output<<"ERROR in lattice::exchange_sites(). Type of atom not >0: "<<ATOM->get_atom()<<endl;exit(1);}

	update_opcja(VAC,0);
	update_opcja(ATOM,0);

	site bufor(VAC);
	VAC->change_to( *ATOM );
	SAMPLE->update_events( VAC );
	ATOM->change_to(bufor);

	update_opcja(VAC,1);
	update_opcja(ATOM,1);

	call_flux_dislocation(ATOM,VAC);
//	control_output<<"Po:"<<endl;
//	VAC->show_site();
//	ATOM->show_site();	
}

void opcja :: update_opcja( site* node, bool status){
//	control_output<<"\nupdate_opcja: ";								//refresh opcja::fields: BLOCKS, HIST, REZ
	int ID_B = node->get_block_index();
	int ID_H = node->get_hist_index();
	int ID_R = node->get_rez_index();
//	control_output<<ID_B<<" "<<ID_H<<" "<<ID_R<<endl;
	if(ID_B >=0){
//		control_output<<"in block ";
		BLOKS[ID_B].update_plaster(node,status);
	}
	if(ID_H >=0){
//		control_output<<"in hist "<<endl;;
//		HIST[ID_H].update_hist(node,status);							//will change HIST, which keep data to print. Not to equlibrate.
	}																	//if ON, then it means that flux comming from dislocation movement 
	if(ID_R >=0){														// is added to flux of particles. 	
//		control_output<<"in rez ";
		reservuars[ID_R].update_plaster(node,status);
	}
}

void opcja :: remove_vac_new(int b, int ile_vac, bool &FLAG){
	
	if(ile_vac<0){ile_vac=ile_vac*-1;}
	unsigned int vac = ile_vac;
	if(MOVE_FRAME or SINGLE){control_output<<" r: "<< vac;}
	bool MOVE = false;
	int rez = -1;
//	control_output<<"TRYB: "<<TRYB<<endl;
	vector <int> wybrane_typy; wybrane_typy.reserve(20);


//	control_output<<"r: "<<b<<" "<<vac<<": ";
//	for(int i=0;i<wybrane_typy.size();i++){
//		control_output<<wybrane_typy[i];
//	}
//	control_output<<endl;

	for(unsigned int i=0; i<vac;i++){
	//		for( int j=1;j<3;j++){		//MOZNA ZROBIC W ZALEZNOSCI OD TYPOW
			long N1,N2;
			site* rnd_vac=0;
			site* rnd_at=0;
//			j = wybrane_typy[i];	//losuje typ atomu do wymiany z wakancja z blokow

			if(BLOKS[b].size(0) <= 0){
				control_output<<endl;
				control_output<<"ERROR: in opcja::remove_vac -> you want to remove \
				element "<<0<<" that does not exist in blok\n	\
				Probably error in opcja::init or opcja::reinit_bloks"<<endl; 
				cout<<endl;
				cout<<"ERROR: in opcja::remove_vac -> you want to remove \
				element "<<0<<" that does not exist in blok\n	\
				Probably error in opcja::init or opcja::reinit_bloks"<<endl;exit(0);
			}			
			
			N1=(long)(rnd()*(BLOKS[b].size(0)));	
			rnd_vac=BLOKS[b].get_site(0,N1);
						
			if(TRYB==2){
				int DIR = -decide_direction(rnd_vac);	//	-1 left;	+1 right
				vector <site*> migration_path; migration_path.reserve(2000);
				find_migration_path(rnd_vac,DIR,migration_path);	
				dislocation_walk(migration_path);							
			}
			else if(TRYB==1){ //swap
				rez=choose_reservuar(rnd_vac);
				int j = choose_typ(reservuars[rez],false);
				if(j<0){j=j*-1;}
				MOVE = check_rezervuars(rez,j);	//sprawdz rezerwuwar pod katem dostepnych atomow
		//			control_output<<"MOVE: "<<MOVE<<" rez: "<<rez<<endl;
				if (!MOVE){
					N2=(long)(rnd()*(reservuars[rez].size(j)));
					rnd_at = reservuars[rez].get_site(j,N2);
		//			control_output<<" rtyp: "<<rnd_at->get_atom()<<" "<<b<<endl;
					rnd_at->set_atom(0);
					reset_site(rnd_at);
					reservuars[rez].delete_site(j,N2);
					reservuars[rez].add_site(0,rnd_at);
					reservuars[rez].prob_update(j,0);

					rnd_vac->set_atom(j);
					reset_site(rnd_vac);	
					BLOKS[b].delete_site(0,N1);
					BLOKS[b].add_site(j,rnd_vac);
					BLOKS[b].prob_update(j,0);
				}else{
					//control_output<<" r: "<< vac;
					//		if(MOVE_FRAME or SINGLE){control_output<<"||"<<BLOKS[b].size(0)<<"|"<<reservuars[rez].size(0)<<"|"<<reservuars[rez].size(j)<<"|"<<Vtoadd.size()<<">|";}
					break;
				}
			}
			else if(TRYB==0){
				int j = choose_typ(BLOKS[b]);
				rnd_vac->set_atom(j);
				reset_site(rnd_vac);	
				BLOKS[b].delete_site(0,N1);
				BLOKS[b].add_site(j,rnd_vac);
				BLOKS[b].prob_update(j,0);
			}
			else{
				control_output<<"ERROR in opcja::remove_vac_new(). Wrong TRYB: "<<TRYB<<endl;
				exit(1);
			}
	}//end of for j
	if(MOVE_FRAME or SINGLE){control_output<<endl;}
	if(MOVE){		
		FLAG = true;		//set local FLAG in do_equi_vac
		do_equi_vac();		//rekurencja. Na poczatku sprawdza czy MOVE_FRAME set to TRUE.
	}
}


void opcja :: refresh_sim_area(vector <site*> &kontener){
	
	if(!MOVE_FRAME){		//nie bylo przesuwania rezerwuarow i blokow
		if(SINGLE){control_output<<"refresh_vac_vector "<<kontener.size()<<endl;}
		int typ=-1;
		bool log=false;
		vector <site* > tmp;
		tmp.reserve(10000);
		tmp.clear();
		for (unsigned int i=0; i < kontener.size(); i++){
			typ=kontener[i]->get_atom();
			log=SAMPLE->check_site_belonging_to_sim_area(kontener[i]);
			if((typ==0) and log ){tmp.push_back(kontener[i]);}else{
				kontener[i]->show_site();
			}
		}
		if(SINGLE){control_output<<"|>"<<tmp.size()<<endl;}
		//teraz trzeba dodac nowe wakancje z Vtoadd
		int count_vac_ok=0;
		for (unsigned int i=0; i < Vtoadd.size(); i++){
			typ=Vtoadd[i]->get_atom();
			log=SAMPLE->check_site_belonging_to_sim_area(Vtoadd[i]);
			if( log and (typ == 0) ){tmp.push_back(Vtoadd[i]);count_vac_ok++;}
//		{control_output<<"ERROR: opcja::refresh_vac_vector, atoms in Vtoadd: "<<log<<"/"<<typ<<endl; 
			Vtoadd[i]->show_site();
//			exit(1);}
		}
		if(SINGLE){control_output<<"|+"<<count_vac_ok;}
		Vtoadd.clear();
		kontener.clear();
		//przepisz i nadaj Vindexy sitom
		for (int i=0; i < tmp.size(); i++){
			tmp[i]->set_vindex(i);
			kontener.push_back(tmp[i]);
		} 
		if(SINGLE){control_output<<"|= "<<kontener.size()<<endl;}
	}else{		//bylo przesowanie blokow i rezerwuwarow
		control_output<<"opcja::refresh_vac_vector "<<kontener.size()<<endl;
		if(MOVE_SIM_REGION){	//przesunac obszar symulacji
			SAMPLE -> reinit_sim_area(del_L_sim, del_R_sim);			
		}
		SAMPLE -> set_atoms_list(kontener,0);
		//resetuj parametry kontrolne przesowania
		wektor a(0.0,0.0,0.0);
		del_L_sim=a;
		del_R_sim=a;
		Vtoadd.clear();
		//clear eq_flux in plasters
	}
	if(SINGLE or MOVE_FRAME){
		for (unsigned int i=0; i < BLOKS.size(); i++){
	//	control_output<<"do init_calc in blok: "<<i<<" "; 
		BLOKS[i].init_calc(1);
		}
		
		for (unsigned int i=0; i < reservuars.size(); i++){
	//	control_output<<"do init_calc in res: "<<i<<" "; 
		reservuars[i].init_calc(1);
		}
	}
	
	if(MOVE_FRAME){MOVE_FRAME = false;}
	if(SINGLE){SINGLE = true;}
}

void opcja :: reset_site(site *sajt){

	SAMPLE->update_events(sajt);
	sajt->reset_site();
	
}

void opcja :: read_file(string filename){
		//WCZYTAJ chem.in
		ifstream chem(filename.c_str(),ios :: in);	//lista potencjalow
		//vector <vector <double> > chem_list;
		
		if( chem.good() )
		{
        string napis;
     //   std::cout << "\nZawartosc pliku:" << std::endl;
        while( !chem.eof() )
        {
            getline( chem, napis );
            //control_output<<"String Line: "<< napis << std::endl;
            double data;
       		vector <double> line;
			istringstream string_line(napis);
	//		int line_count=0;
			while(string_line>>data)
			{
		//		cout<<typeid(data).name()<<" "<<data<<endl;
				line.push_back(data);
		//		cout<<line.size()<<endl;
			}

			if( napis != ""){
				equi_curve.push_back(line);}
        
        }
      //  chem_list.erase(chem_list.end()); //usun linie konca pliku dodana do vectora
        chem.close();
		} else {control_output<< "Error! Nie udalo otworzyc sie pliku "<<filename<<endl;
		exit(0);
		}
	//	for(int i=0; i<equi_curve.size();i++)
	//	{	//cout<<chem_list[i].size()<<endl;
			//if((chem_list[i].size() +1) != sample->get_atom_typ_numbers())	// +1 poniewaz wakancja
	//		{
	//			control_output<<"ERROR: execute_task: SGCMC in reading chem.in file. ";
	//			control_output<<"Number of types different than declared in structure.dat file"<<endl;
			//	control_output<<"chem: "<<(chem_list[i].size()+1)<<" structure: "<<sample->get_atom_typ_numbers()<<endl;
	//			exit(0);
	//		}
	//	}
		//chem_list.size() - ile wierszy jest w pliku (punktow do symulacji)
		//chem_list[i].size() - ile typow atomow (pot. chem) jest zadeklarowanych
		//WCZYTALEM chem.in
}

void opcja :: save_hist(double MCtime, double MCstep, string fileout){
	
	SAVE_file = fileout;
	SAVE_MCstep = MCstep;
	SAVE_MCtime = MCtime;
	save_call();
	save_write();
	
}

void opcja :: init_save(double iod,double ido,double istep, double direction){
	if(!SAVE_BUILDED){
		double tmp_st=0;
		double tmp_end=0;
		double tmp_size=0;
		double tmp_dir=0;
		
		if(TRYB>0){
			tmp_st=BIN_ST;
			tmp_end=BIN_END;
			tmp_size=BIN_SIZE;
			tmp_dir=BIN_DIRECTION;
		}
		
		BIN_ST=iod;
		BIN_END=ido;
		BIN_SIZE=istep;
		BIN_DIRECTION=direction;
		build_bins(HIST,"hist");		
		
		if(TRYB>0){
			BIN_ST=tmp_st;
			BIN_END=tmp_end;
			BIN_SIZE=tmp_size;
			BIN_DIRECTION=tmp_dir;
		}
		SAMPLE->get_sites(HIST);
		for (unsigned int i=0; i < HIST.size(); i++){
			control_output<<"do init_calc in hist: "<<i<<" "; 
			HIST[i].init_calc(1);
		}
		SAVE_BUILDED=true;
	}
}

void opcja :: save_call(){
	
	refresh();
	vector<plaster>& wsk2hist=HIST;
	vector<plaster>& wsk2blok=BLOKS;
	vector<plaster>& wsk2rez=reservuars;
	for(unsigned int i =0;i<wsk2hist.size();i++){
		wsk2hist[i].cumulate();	
		}	
	for(unsigned int i =0;i<wsk2blok.size();i++){
		wsk2blok[i].cumulate();	
		}	
	for(unsigned int i =0;i<wsk2rez.size();i++){
		wsk2rez[i].cumulate();	//dodaje aktualna wartosc parametrow plastra do avg_par's
		}
}


void opcja :: save_write(){
	
	ofstream data_bin;
	ofstream data_hist;	

	vector<plaster>& wsk2hist=HIST;
	vector<plaster>& wsk2blok=BLOKS;
	vector<plaster>& wsk2rez=reservuars;
	vector<double> results;
	results.reserve(30);
	
	if(wsk2blok.size()>0){
		string name=SAVE_file+"bin.dat";
		data_bin.open(name.c_str(),ios::app);
	}
	if(wsk2hist.size()>0){
		string name=SAVE_file+"hist.dat";
		data_hist.open(name.c_str(),ios::app);	
	}
	
	for(unsigned int i =0;i<wsk2blok.size();i++){
		wsk2blok[i].call_avg(results);
		data_bin<<SAVE_MCstep<<" "<<SAVE_MCtime;
		for(unsigned int j=0;j<results.size();j++){
			data_bin<<" "<<results[j];
			}
			data_bin<<endl;
		
	}
	for(unsigned int i =0;i<wsk2rez.size();i++){
		wsk2rez[i].call_avg(results);
		data_bin<<SAVE_MCstep<<" "<<SAVE_MCtime;
		for(unsigned int j=0;j<results.size();j++){
			data_bin<<" "<<results[j];
			}
			data_bin<<endl;	
	}	
	data_bin<<endl;
	for(unsigned int i =0;i<wsk2hist.size();i++){
		wsk2hist[i].call_avg(results);
		data_hist<<SAVE_MCstep<<" "<<SAVE_MCtime;
		for(unsigned int j=0;j<results.size();j++){
			data_hist<<" "<<results[j];
			}
			data_hist<<endl;	
	}	
	data_hist<<endl;
}


void opcja :: show(){

	for (unsigned int i=0;i<BLOKS.size();i++){
		for (unsigned int j=0;j<BLOKS[i].size();j++){
			
			control_output<<i<<" "<<j<<" "<<BLOKS[i][j]<<endl;
			BLOKS[i][j]->show_site();
		}
	}
}



