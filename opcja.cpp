//opcja_swap//

#include "opcja.h"

void opcja :: execute(string name, vector<double>&parameters){

	if(name=="EQUILIBRATE"){
		this->init_EQ(parameters);
	}else if(name=="RESERVUAR"){
		this->init_reservuar(parameters);
	}else if(name=="L_SIM_B"){
		SAMPLE->init_sim_boundary(parameters);
	}else{
		control_output<<"No posibility of doing "<<name<<" as an option."<<endl;
	}

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

bool opcja :: check_rezervuars(site* first, site* &last){
	if(DEBUG){control_output<<"check reservuars st";}
	site* actual = last;
	int TYP = first->get_atom();	
	bool do_move=false;
	
	unsigned int nr_rez = actual->get_rez_index();

	if(nr_rez<0){
		control_output<<"ERROR: opcja::check_rez(). Site not in the reservuar"<<endl;
		actual->show_site();
		control_output<<"Wrong definition of res."<<endl;
		exit(1);
	}

	if(TYP==0){
		int at_typ = (reservuars[nr_rez]).choose_typ();
		if(at_typ <= 0){
			control_output<<"warrning in opcja::check_reservuars -> try remove typ which dose not exist in reservuar "<<at_typ<<endl;
			TYP_TO_MOVE=1;
			do_move=true;
			MOVE_FRAME=true;											//set global variable to true
			REZ_TO_MOVE=nr_rez;											//set global var. which rezervuar to move			
		}else{
			last = (reservuars[nr_rez]).choose_atom(at_typ);
		}
	}else if(TYP>0){
		if(reservuars[nr_rez].check(0,0,1)){
			control_output<<"warrning in opcja::check_reservuars -> try remove typ which dose not exist in reservuar "<<0<<endl;
			TYP_TO_MOVE=-1;
			do_move=true;
			MOVE_FRAME=true;											//set global variable to true
			REZ_TO_MOVE=nr_rez;											//set global var. which rezervuar to move
		}else{
			//choose vacancy from rez
			last = (reservuars[nr_rez]).choose_atom(0);
			}	
	}else{
		control_output<<"ERROR: opcja::check_rez(). Wrong type: "<<endl;
		first->show_site();
		exit(1);
	}		
	if(DEBUG){control_output<<"check reservuars end"<<endl;}

	return do_move;
}

bool opcja :: check_rez_dN(){
	
	bool do_move=false;
	long int V=0, Z=0;
	if(DEBUG or DEBUG_SMALL){control_output<<"check_rez_dN:->";}

//	if( (i<0) and (i>=reservuars.size()) ){
//		control_output<<"ERROR: opcja::check_rez(): 209. Wrong rez: "<<i<<endl;
//		exit(1);}
		
	for(unsigned int i =0;i<reservuars.size();i++){
		V = reservuars[i].eq_flux_get(0) + reservuars[i].flux_net_get(0);
		Z=labs(V);

		if(Z > ROZMIAR[i]){												//calkowita zmiana atomow wynosi tyle co 100% jednej plaszczyzny
			if(V>0){													//dV>0 to znaczy ze w rezerwuarze powstala nowa plaszczyzna wakancji -> probka w strone atomow
				TYP_TO_MOVE=1;
			}else if(V<0){													//dV<0 to znaczy ze w rezerwuarze powstala nowa plaszczyzna atomow -> probka w strone wakancji
				TYP_TO_MOVE=-1;
			}else{control_output<<"error in opcja::check_reservuars "<<V<<endl;
				exit(1);
			}

			do_move=true;
			MOVE_FRAME=true;	//set global variable to true
			NEW_PLANE=true;
			REZ_TO_MOVE=i;	//set global var. which rezervuar to move
			break;
		}
	}

	if(DEBUG or DEBUG_SMALL){control_output<<"|->check_rez_dN";}

	return do_move;
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
	}else{	//jesli nie ma nowej plaszczyzny ale zabraklo atomow typu TYP w rezerwuarze		<< nie mozliwy scenariusz bo typ losowany z rezerwuaru
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

void opcja :: do_equi_rez(){
	if(DEBUG or DEBUG_SMALL){	control_output<<"do_equi_rez:->";}

	//direct exchange dla rezerwuarow	//DODAC PARALLEL
	double beta=1.0/(kB*TEMPERATURE);
	for(unsigned int i=0;i<reservuars.size();i++){	
		long Nsize=reservuars[i].size();
	
	//	cout<<"rez: "<<i<<" size: "<<Nsize<<endl;
		for (long iter=0; iter<DIRECT_STEPS; iter++){
			long N1=(long)(rnd()*Nsize);
			long N2=(long)(rnd()*Nsize);		
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
				if(R > p){
				site1->set_atom(typ1);
				site2->set_atom(typ2);
				}
			}
	 		
			SAMPLE->update_events( site1 );
			SAMPLE->update_events( site2 );
			if(site1->events_size() > 0){
				Vtoadd.insert(site1);
			}
			if(site2->events_size() > 0){
				Vtoadd.insert(site2);
			}
		
	//	ofstream dir_file("dir_tmp.dat", ios::app); 
	//	dir_file<<i<<" "<<iter<<" "<<E1<<"|-"<<E2<<"|="<<dE<<"|"<<p<<"<?"<<R<<endl;
		}
	}		
	if(DEBUG or DEBUG_SMALL){	control_output<<"|->do_equi_rez"<<endl;}

}
  
site* opcja :: get_node(int in_bin, bool create, int for_rez, long int &nr_site){
	
	site* node=0; int j=-1;
	vector <int> exclude;exclude.reserve(10);
	if(create){
		j=(reservuars[for_rez]).choose_typ();
	}else{
		j=0;
	}
	if(DEBUG_SMALL){control_output<<"get_node:->";}

	if(DEBUG){
		control_output<<"get node "<< in_bin<<" "<<for_rez<<" "<<nr_site<<" "<<create<<" ";
		control_output<<node<<" "<<j<<endl;
	}
	
	bool avail =  ( (BLOKS[in_bin]).size(j) > 0);
	while(!avail){
		if(j>0){																	//take another schoot exlusivly
			exclude.push_back(j);
			j=(reservuars[for_rez]).choose_typ(exclude);	
			if(j<=0){
				control_output<<"WARRNING: in opcja::get_node: types in rez "<<for_rez<<" do not match types in bin "<<in_bin<<endl;
				control_output<<"Type choosen proportional to bin."<<endl;
				j=(BLOKS[in_bin]).choose_typ();	
			}
			avail = (BLOKS[in_bin].size(j) > 0);
		}else{
			control_output<<"ERROR: in opcja::get_node -> you want to remove type that does not exist in blok\n"<<endl;
			control_output<<create<<" "<<j<<" "<<avail<<" "<< (BLOKS[in_bin]).size(j)<<endl;
			(reservuars[for_rez]).show();
			(BLOKS[in_bin]).show();
			exit(1);
		}
	}

	node = (BLOKS[in_bin]).choose_atom(j);
	
	if(DEBUG){
		control_output<<"get node end "<< in_bin<<" "<<for_rez<<" "<<nr_site<<" "<<create<<" ";
		control_output<<node<<" "<<j<<" "; node->show_site();
	}
	int for_typ=node->get_atom();
	if( for_typ==0 and create){control_output<<"ERROR:opcja::get_node:641: create  vac, but trying remove vac from system: "<<for_typ<<endl; 
		control_output<<create<<" "<<j<<" "<<nr_site<<" ";node->show_site();
		(reservuars[for_rez]).show();
		(BLOKS[in_bin]).show();
		exit(1);}
	if( for_typ>0 and !create){control_output<<"ERROR:opcja::get_node:346: remove  vac, but trying remove atom from system: "<<for_typ<<endl; 
		control_output<<create<<" "<<j<<" "<<nr_site<<endl;
		(reservuars[for_rez]).show();
		(BLOKS[in_bin]).show();
		exit(1);}
	
	if(DEBUG_SMALL){control_output<<"|->get_node";}

			
return node;	
}

double opcja :: call_total_flux(double range, unsigned int nr_bin){

	if(DEBUG_SMALL){control_output<<"call_tot_flux:->";}

	long SUM=0;long iter=0;
//	unsigned int st = nr_bin - int(range/2.0);
	
	//using bin to calculate over hist??
	//use position of the bin, terate over hist, check position of hist if it is within range, then cumulate.
	//if position of the hist is outside bin then cumulate twice time in oposite direction
	// http://www.cprogramming.com/c++11/c++11-ranged-for-loop.html
	//actually total flux is calculated 
	
	vector <plaster>::iterator P=HIST.begin();
	
	for(; P!=HIST.end(); ++P){
		SUM += P->net_flux_get(0); iter++;
	}
	if(DEBUG_SMALL){control_output<<"|->call_tot_flux";}	
	if(iter==0){iter=1;}
	return SUM/iter;
}

site* opcja :: source_sink_localize(int in_bin, bool create, int &from_rez, long int &nr_site, int &in_dir){

	double displace = 0;
	double X0=( (BLOKS[in_bin]).get_st() + (BLOKS[in_bin]).get_end() )/2.0;
	double C = (BLOKS[in_bin]).get_stech();
//	control_output<<"Decide for: "<<X0<<" "<<bin<<" "<<C<<" ";

	double maxY=0, norma =1.0, min=1.0, max=0.0; 
	int REZ=-1, sign = 1;
	long N=-1;
	site* node =0;
	vector <unsigned int> mykey; mykey.reserve(10);
	
	if(DEBUG_SMALL){control_output<<"sink_loc:->";}
	if(DEBUG){	control_output<<"sink loc "<< in_bin<<" "<<from_rez<<" "<<in_dir<<" "<<nr_site<<" "<<N<<" "<<create<<" ";
	control_output<<maxY<<" "<<REZ<<" "<<N<<" "<<node<<endl;}

	for( unsigned rez = 0; rez < reservuars.size(); rez++){				//find maximum and min stechiometry range
		double CR = (reservuars[rez]).get_stech();
		if(CR >= max){max=CR;}
		if(CR <= min){min=CR;}
	}
	norma=fabs(max - min);
	if(norma > 1.0){
			control_output<<"ERROR:opcja::source_localize():400 "<<min<<" "<<max<<" "<<norma<<endl;
			exit(1);
	}
																		//separate diffusion zone
	for( unsigned rez = 0; rez < reservuars.size(); rez++){
		
		if( ((reservuars[rez]).get_st() <= X0) and (X0 <= reservuars[rez].get_end()) ){
			control_output<<"ERROR:opcja::source_localize():-> bad definition of rezervuars: "<<X0<<endl;
			(reservuars[rez]).show();
			(BLOKS[in_bin]).show();
			exit(1);
		}
		double CR = (reservuars[rez]).get_stech();
		double Y = 1.0 - fabs( CR - C);	
		if( (Y >= maxY - 0.5*norma) and (Y <= maxY + 0.5*norma) ){				
			mykey.push_back(rez);					
		}else if(Y > maxY + 0.5*norma){
			maxY=Y;
			mykey.clear();
			mykey.push_back(rez);		
		}else if (Y < maxY - 0.5*norma){
			continue;
		}else{
			control_output<<"ERROR:opcja::source_localize():-> bad definition of rezervuars:419"<<X0<<endl;
			(reservuars[rez]).show();
			(BLOKS[in_bin]).show();
			exit(1);
		}
	}																	//		control_output<<"	"<<REZ<<" "<<Y<<" "<<minY<<" "<<XL<<" "<<XP<<" "<<X<<" "<<minX<<endl;
	if(DEBUG){	control_output<<"sink loc "<<maxY<<" "<<mykey.size()<<":";
	for(unsigned int i=0; i<mykey.size();i++){
		control_output<<" "<<mykey[i];
	}control_output<<endl;}
																		
	if(mykey.size() == 1){												//when out off diffusion zone, get nearest rez
		REZ=mykey[0];
		if(DEBUG){		control_output<<"sink loc 1 "<<maxY<<" "<<REZ<<" "<<displace<<endl;}
	
	}else if(mykey.size() > 1){											//when in diffusion zone, test local net flux
		double Ft = call_total_flux(10, in_bin);									//change to temporary flux (time) and local (position). Add rounding of flux to include fluctuations
		
		if(!create){sign=-1;}													//if vacancy is to be removed -> change direction of movement
	
		//build list of rez
		typedef vector <pair <double,double> > lista; double sum=0;
		lista tmp_rta;
		for(unsigned int i=0; i<mykey.size();i++){
			double left = (reservuars[REZ]).get_st() - X0;
			double right = (reservuars[REZ]).get_end() - X0;
			double dx = fabs(left) > fabs(right) ? right : left;
			double P = exp(-(sign*dx*Ft)/(TEMPERATURE*kB));				
			tmp_rta.push_back(make_pair(sum,mykey[i]));sum += P;
		}
		tmp_rta.push_back(make_pair(sum,-1));

		//select rez
		double R=rnd()*sum; 
		lista::iterator event=tmp_rta.begin();
		lista::iterator next_event=tmp_rta.begin();
		for( ++next_event ; next_event != tmp_rta.end(); ++event, ++next_event){	
			double Lvalue = (*event).first;
			double Rvalue = (*next_event).first;	
			if( R>=Lvalue and R < Rvalue){
				REZ=(*event).second;
			}
		}		
	}else{
		control_output<<"ERROR: opcja::source_sink_localize:469"<<endl;exit(1);
	}

	//select initial node from bin
	if(DEBUG){	control_output<<"sink loc 6 "<<N<<" "<<REZ<<" "<<displace<<" "<<endl;}
	if(REZ>=0 and ( N < 0 and node == 0 )){
		node=get_node(in_bin,create,REZ,N);
	}else{
		control_output<<"ERROR: opcja::source_sink_localize:476"<<endl;exit(1);
	}
	
	//calculate direction
	double wal_l = (reservuars[REZ]).get_st();
	double wal_r = (reservuars[REZ]).get_end();
	if(wal_l < X0 and X0 < wal_r){
		control_output<<"ERROR: opcja::source_sink_localize:483"<<endl;exit(1);
	}
	double left = wal_l - X0;
	double right = wal_r - X0;
	displace = fabs(left) > fabs(right) ? right : left;				//take shorter distance to the wall
																	
	from_rez=REZ;
	nr_site=N;
	in_dir=displace;
	if(DEBUG){	control_output<<"sink loc end "<< in_bin<<" "<<from_rez<<" "<<in_dir<<" "<<nr_site<<" "<<N<<" "<<create<<" ";
	control_output<<maxY<<" "<<REZ<<" "<<N<<" "<<node<<endl;}
	if(DEBUG_SMALL){control_output<<"|->sink_loc";}

	return node;
}


void opcja :: source_sink_act(int in_bin, int ile_at, bool &FLAG){

	if(MOVE_FRAME or SINGLE){control_output<<" c: "<< ile_at;}
	bool MOVE = false, CREATE=true;
	if(DEBUG_SMALL){control_output<<"sink_act:->";}
	if(DEBUG){control_output<<"sink act "<<in_bin<<" "<<ile_at<<" "<<FLAG<<" "<<MOVE<<" "<<CREATE<<endl;}

	if(ile_at < 0){
		CREATE=false;
		ile_at=abs(ile_at);
	}else if (ile_at > 0){
		CREATE=true;
		ile_at=abs(ile_at);
	}else{
		return ;
	}
	if(DEBUG){control_output<<"sink act "<<ile_at<<" "<<CREATE<<endl;}

	for(int i=0; i<ile_at;i++){
		int from_rez = -1, in_dir = 0, for_typ=-1, to_typ=-1;
		site* AT1=0; site* AT2=0;
		long N1=-1,N2=-1;
	if(DEBUG){	control_output<<"sink act "<<" "<<i<<" "<< in_bin<<" "<<from_rez<<" "<<in_dir<<" "<<N1<<" "<<N2<<" "<<for_typ<<" "<<to_typ<<" ";
				control_output<<AT1<<" "<<AT2<<" "<<endl;}


		if(TRYB==2){													//dislocation move
			AT1=source_sink_localize(in_bin, CREATE, from_rez, N1, in_dir);		
			for_typ=AT1->get_atom();	
			if(for_typ==0 and CREATE){control_output<<"ERROR:opcja::source_sink_act:630: create  vac, but remove vac from system type: "<<for_typ<<endl; exit(1);}
			if(for_typ>0 and !CREATE){control_output<<"ERROR:opcja::source_sink_act:631: remove  vac, but move atom from sysstem type: "<<for_typ<<endl; exit(1);}
			if(DEBUG){	control_output<<"sink act "<<i<<" "<<from_rez<<" "<<in_dir<<" "<<for_typ<<" "<<to_typ<<" "<<AT1<<" "<<AT2<<" ";
			control_output<<N1<<" "<<N2<<" "<<CREATE<<endl;}

			vector <site*> migration_path; migration_path.reserve(2000);
			MOVE = find_migration_path(AT1,in_dir,migration_path);	
			
			if(!MOVE){
				if(!migration_path.empty()){
					int tmp_typ=(migration_path.front())->get_atom();	
					if(tmp_typ==0 and CREATE){control_output<<"ERROR:opcja::source_sink_act:639: create  vac, but remove vac from system type: "<<tmp_typ<<endl; exit(1);}
					if(tmp_typ>0 and !CREATE){control_output<<"ERROR:opcja::source_sink_act:640: remove  vac, but move atom from sysstem type: "<<tmp_typ<<endl; exit(1);}
				
					dislocation_walk(migration_path);
					MOVE=check_rez_dN();
					if(MOVE){
						break;
					}
				}
			}else{
				break;
			}
		}else if(TRYB==1){ 												//swap
			AT1=source_sink_localize(in_bin, CREATE, from_rez, N1, in_dir);
			for_typ=AT1->get_atom();		

			if(CREATE){
				to_typ=0;
			}else{
				to_typ = (reservuars[from_rez]).choose_typ();
			}
			MOVE = check_rezervuars(from_rez,to_typ);
			if (!MOVE){
				N2=(long)(rnd()*(reservuars[from_rez].size(to_typ)));
				AT2 = (reservuars[from_rez]).get_site(to_typ,N2);

				AT2->set_atom(for_typ);			
				reset_site(AT2);
				reservuars[from_rez].delete_site(to_typ,N2);
				reservuars[from_rez].add_site(for_typ,AT1);
				
				AT1->set_atom(to_typ);
				reset_site(AT1);
				BLOKS[in_bin].delete_site(for_typ,N1);
				BLOKS[in_bin].add_site(to_typ,AT1);
				BLOKS[in_bin].prob_update(for_typ,1);
				Vtoadd.insert(AT1);	
			}else{
				//	if(MOVE_FRAME or SINGLE){control_output<<"||"<<BLOKS[nr].size(j)<<"|"<<reservuars[rez].size(0)<<"|"<<reservuars[rez].size(j)<<"|"<<Vtoadd.size()<<">|";}
				break;
			}
		}else if(TRYB==0){												//convert
			if(CREATE){
				to_typ=0;
				for_typ = (BLOKS[in_bin]).choose_typ();
				N1=(long)(rnd()*(BLOKS[in_bin].size(for_typ)));
				AT1 = (BLOKS[in_bin]).get_site(for_typ,N1);
			}else{
				for_typ=0;
				to_typ = (BLOKS[in_bin]).choose_typ();
				N1=(long)(rnd()*(BLOKS[in_bin].size(for_typ)));
				AT1 = (BLOKS[in_bin]).get_site(for_typ,N1);
			}
			AT1->set_atom(to_typ);
			reset_site(AT1);
			BLOKS[in_bin].delete_site(for_typ,N1);
			BLOKS[in_bin].add_site(to_typ,AT1);
			BLOKS[in_bin].prob_update(for_typ,1);
			Vtoadd.insert(AT1);	
		}
	}																	// end for(j=delta)
	
	if(MOVE_FRAME or SINGLE){control_output<<endl;}
	if(MOVE){
		FLAG = true;													//set local FLAG in do_equi_vac
		do_equi_vac();													//recurency. At the beginnig check if MOVE_FRAME is set to TRUE.
	}
	if(DEBUG or DEBUG_SMALL){	control_output<<"|->sink_act"<<endl;}


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

	if(DEBUG or DEBUG_SMALL){control_output<<"do_equi:->"<<endl;}

	bool LOCAL_MOVE = false;
	if(MOVE_FRAME or SINGLE){	
	control_output<<"Rownowarze: "<<LOCAL_MOVE<<"|"<<MOVE_FRAME<<endl;
	refresh(1);}
	
	if(MOVE_FRAME){
		move_frame();	//do_move resetuje global MOVE_FRAME
	}
	
	for (unsigned int i=0;i<BLOKS.size();i++){
		if( 1 ){
	//	bloks[i].calc_stech();
		if(MOVE_FRAME or SINGLE)
		{control_output<<i;}
		
		double stech=BLOKS[i].get_stech();
		double vac=BLOKS[i].get_vac();
		double size=BLOKS[i].size();
		int delta_vac = check_stech(stech,vac,size);	//zwracac ile wakancji remove/create
	if(DEBUG){	control_output<<"do_equi "<<i<<" "<<stech<<" "<<vac<<" "<<size<<" "<<delta_vac<<" "<<endl;}
		
		if(delta_vac < 0){
			source_sink_act(i, delta_vac, LOCAL_MOVE);
		}
		else if (delta_vac == 0)
		{
			if(MOVE_FRAME or SINGLE){	
	if(DEBUG){			control_output<<" do nothing"<<endl;}}
		}
		else if ( delta_vac > 0 )
		{
			source_sink_act(i,delta_vac, LOCAL_MOVE);
		}
		else{
			cout<<"Error with equilibrate vacancy opcja::do_equi_vac()"<<endl;
			exit(1);
		}
		
		
		
		if(LOCAL_MOVE){
	if(DEBUG){			control_output<<" po: "<<i<<"|"<<LOCAL_MOVE<<"|"<<MOVE_FRAME<<endl;}
			break;}	//jesli bylo do_move przestan rownowazyc pozostale bloki
	}	//if exclude block
	}	//for bloks
	
//	refresh(1);
//	control_output<<"Koniec Rownowaga!"<<endl;
			if(!LOCAL_MOVE){
			if(MOVE_FRAME or SINGLE){
			control_output<<" po: all|"<<LOCAL_MOVE<<"|"<<MOVE_FRAME<<endl;}}
	//uaktulanic liste wakancji w  resident !!!!!!!!!!!! UWAGA !!!!!!!!!	
	if(DEBUG or DEBUG_SMALL){	control_output<<"|->do_equi"<<endl;}


}

void opcja :: equilibrate(){
	
	refresh(0);
	do_equi_vac();
	do_equi_rez();
	refresh_vac_list();
	
//	int eventy = 0;
//	for (std::set<site*>::iterator it=VAC_LIST.begin(); it!=VAC_LIST.end(); ++it){
//		(*it)->show_site();
//		eventy += (*it)->events_size();
//	}
//	control_output<<"Print VAC_LIST: "<<VAC_LIST.size()<<"/"<<eventy<<endl;
	
}

void opcja :: flux_net_add(double pos_V, double pos_A, int typV, int typA, vector<plaster>& layer){
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
		layer[id_V].flux_net_delta(typA, 0);
		layer[id_V].flux_net_delta(typV, 1);
	}
	if(id_A >=0 and (id_V != id_A)){
		//vakancja znika z plastra id_A oraz atom pojawia sie w plastrze id_A
		layer[id_A].flux_net_delta(typA, 1);
		layer[id_A].flux_net_delta(typV, 0);
	}	
}

void opcja :: flux_add(site* VAC, site* ATO, vector<plaster>& layer){

	if(SAVE_BUILDED){
	int id_A=-2,id_V=-2;				
	unsigned int DIR = (layer[0]).get_direction();
	double pos_V = VAC->get_position(DIR);								//vacancy after jump
	double pos_A = ATO->get_position(DIR);								//atom after jump
	double dir = SAMPLE->move(pos_V,pos_A,DIR);							// + if vac moved to right
	//move takes position of vac before jump and position of atom before jump
	//if dir > 0 - wakancja skoczyla z lewej strony na prawa. Strumien dla wakancji jest +

	id_V=VAC->get_hist_index();
	id_A=ATO->get_hist_index();

	int typV = VAC->get_atom();
	int typA = ATO->get_atom();
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
			if(id_V >= id_A){
				for(int ID=id_V; ID>id_A; ){
					layer[ID].jump_occured();
					layer[ID].flux_net_delta(typV, 1);
					layer[ID].flux_net_delta(typA, 0);
					ID--;
				}
			}else{
				layer[id_V].jump_occured();
				layer[id_V].flux_net_delta(typV, 1);
				layer[id_V].flux_net_delta(typA, 0);
			}
		}
		if(dir < 0){	//wakancja wplynela z prawej strony; atom wyplynal w prawo
			if(id_V <= id_A){
				for(int ID=id_V; ID<id_A; ){
					layer[ID].jump_occured();
					layer[ID].eq_flux_delta(typV, 0);
					layer[ID].eq_flux_delta(typA, 1);
					ID++;
				}
			}else{
					layer[id_V].jump_occured();
					layer[id_V].eq_flux_delta(typV, 0);
					layer[id_V].eq_flux_delta(typA, 1);		
			}
		}	
	}
	if(id_A >=0 and (id_V != id_A)){
		//vakancja pojawia sie w plastrze id_V oraz atom znika z plastra id_V
		if(dir < 0){	//wakancja w lewo strony; atom wyplyna w prawo
			if(id_A >= id_V){
				for(int ID=id_A ;ID>id_V ;){
					layer[ID].jump_occured();
					layer[ID].flux_net_delta(typA, 1);
					layer[ID].flux_net_delta(typV, 0);
					ID--;
				}
			}else{
				layer[id_A].jump_occured();
				layer[id_A].flux_net_delta(typA, 1);
				layer[id_A].flux_net_delta(typV, 0);				
			}
		}
		if(dir > 0){	//wakancja wplynela z lewej strony; atom wyplynal w lewo
			if(id_A <= id_V){
				for(int ID=id_A ;ID<id_V ;){
					layer[ID].jump_occured();
					layer[ID].eq_flux_delta(typA, 0);
					layer[ID].eq_flux_delta(typV, 1);
					ID++;
				}
			}else{
				layer[id_A].jump_occured();
				layer[id_A].eq_flux_delta(typA, 0);
				layer[id_A].eq_flux_delta(typV, 1);
			}
		}	
	}
}

}

void opcja :: flux_add_dislocation(site* VAC, site* ATO, vector<plaster>& layer){
	int id_A=-2,id_V=-2;				
	unsigned int DIR = (layer[0]).get_direction();
	double pos_V = VAC->get_position(DIR);
	double pos_A = ATO->get_position(DIR);
	double dir = SAMPLE->move(pos_V,pos_A,DIR);											//if > 0 - wakancja skoczyla z lewej strony na prawa. Strumien dla wakancji jest +

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

void opcja :: find_interface(){

	if(0){
	double x0 = BIN_ST + fabs(  (BIN_ST-BIN_END)/2.0 ) ;
//	cout<<BIN_DIRECTION<<" "<<BIN_ST<<" "<<BIN_END<<" "<< x0<<endl;

	

	vector <double> X;
	vector <double> Y;
	
	typedef vector <plaster>::iterator it2pl;
	for( it2pl IT=reservuars.begin();IT!=reservuars.end();++IT){
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
	
}

void opcja :: call_flux(site* vac_after_jump,site* atom_after_jump){

	if(SAVE_BUILDED){
		flux_add(vac_after_jump,atom_after_jump,HIST);
	}

	if(EQ_BUILDED){	
		double x_A = atom_after_jump->get_position(BIN_DIRECTION);
		double x_V = vac_after_jump->get_position(BIN_DIRECTION);
		int typA=atom_after_jump->get_atom();
		int typV=vac_after_jump->get_atom();
		flux_net_add(x_V,x_A,typV,typA,BLOKS);
		flux_net_add(x_V,x_A,typV,typA,reservuars);

		MOVE_FRAME=check_rez_dN();
		if(MOVE_FRAME){
			move_frame();	
		}	
	}
}

void opcja :: call_flux_dislocation(site* vac_after_jump,site* atom_after_jump){
	flux_add_dislocation(vac_after_jump,atom_after_jump,HIST);	
}

void opcja :: set_opcja_lattice(lattice *sample){
	
	SAMPLE=sample;
	BIN_ATOMS_TYP=sample->get_atom_typ_numbers();
}

void opcja :: init_EQ(vector <double> &parameters ){
	
	if(parameters.size()!=6){cout<<"ERROR in opcja::init. Wrong parameters list in conf.in"<<endl;exit(1);}

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

	TRYB=tryb;	
	EQ_STEP=co_ile;			
	BIN_ST=od_kod;
	BIN_END=do_kod;
	BIN_SIZE=co_ile_bin;
	BIN_DIRECTION=bin_direction;
	ST_VOL=BIN_ST;
	END_VOL=BIN_END;

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
	
	EQ_BUILDED=true;	
}


void opcja :: build_bins(vector<plaster>& layer, string name){	
	
	for ( vector<plaster>::iterator it2pl = layer.begin(); it2pl != layer.end(); ++it2pl){
		it2pl->reset_indexes();
	}	
	layer.clear();
	
	int ile = (int((BIN_END - BIN_ST)/BIN_SIZE));
	double miarka[ile];
	miarka[0]=BIN_ST;
	for(int ii=1; ii<=ile; ii++){miarka[ii] = miarka[ii-1] + BIN_SIZE;}
	
//	for(int ij=0; ij<ile; ij++){control_output<<ij<<" "<<miarka[ij]<<endl;}

	
	for (int i=0;i<ile;i++){
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
//	double o_ile = parameters[2];		//-> init_reservuars [2] to move frame lenght
//	int direction =parameters[3];
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
	reinit_reservuars(REZ_TO_MOVE,TYP_TO_MOVE);
	reinit_bloks();

	//reset global FLAGS	REZ_TO_MOVE is reset in refresh_sim_area
	
	REZ_TO_MOVE = -1; 	
	TYP_TO_MOVE = 0;
	NEW_PLANE = false;
	refresh_simarea();
	
};

void opcja :: refresh(int on){

	if(DEBUG or DEBUG_SMALL){control_output<<"refresh st ";}
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
	if(DEBUG or DEBUG_SMALL){control_output<<"->refresh end"<<endl;}
	
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
	
	(reservuars[nr]).reset_indexes();
	
	plaster tmp(2000,BIN_ATOMS_TYP,BIN_DIRECTION,nr,od_kod,do_kod, "rez");
	SAMPLE->get_sites(tmp);
	control_output<<"reinit reservuar and volume: "<<ST_VOL<<" "<<END_VOL<<" ";
	control_output<<nr<<" "<<odkod<<"|"<<stwidth<<"|"<<dokod<<"||"<<od_kod<<"|"<<width<<"|"<<do_kod<<" "<<direction<<endl;

	control_output<<"do init_calc in reservuar: "; 
	//przerzuc zostawione atomy do nowego rezerwuaru (swap atom<->wakancja)
	if(DIR==1){
		tmp.init_calc(1);
		swap(reservuars[nr],tmp,1);
	}
	tmp.init_calc(1);
	if(!NEW_PLANE){
		tmp.copy_fluxes(reservuars[nr]);
	}
	reservuars[nr]=tmp;

	wektor a(s[0],s[1],s[2]);
	wektor b(e[0],e[1],e[2]);
	del_L_sim += a;
	del_R_sim += b;			

	ST_VOL=od_kod;END_VOL=do_kod;
	for(unsigned int my_rez=0; my_rez < reservuars.size();my_rez++){
		double tmp_s=reservuars[my_rez].get_st();
		double tmp_e=reservuars[my_rez].get_end();
		ST_VOL=min(tmp_s,ST_VOL);
		END_VOL=max(tmp_e,END_VOL);
	}
	control_output<<"New volume: "<<ST_VOL<<" "<<END_VOL<<" "; 

}

void opcja :: reinit_bloks(){
	build_bins(BLOKS,"block");				
	SAMPLE->get_sites(BLOKS);
	control_output<<"do init_calc in blok: "<<endl;
	for (unsigned int i=0; i < BLOKS.size(); i++){
		control_output<<i<<" "; 
		BLOKS[i].init_calc(1);
	}
	//show();
}

bool opcja :: check_x_belonging_volume(double x){
				   
	if( (set_prec(x) >= set_prec(get_start_volume()) ) and ( set_prec(x) < set_prec(get_end_volume()) ) ){
				return true;
	}
	//else{
	//	control_output<<"WARRNING:opcja::check_volume():-> site out of volume: "<<x<<" "<<ST_VOL<<" "<<END_VOL<<endl;	
	//}
	return false;
}

void opcja :: cal_angles(site *node, wektor &main, vector <site*> &wynik_at, vector <site*> &wynik_vac){
	
	vector <site*> neighs;													//dla noda pobierz sasiadow
	vector <site*>::iterator atom;
	int dir = get_direction();
	
//	control_output<<"weak ";node->show_site();
	wektor r0 = node->get_position();
	bool IN_VOLUME = check_x_belonging_volume(r0[dir]);

	if(IN_VOLUME){

	bool IN_REZ = false;	

//	vector <plaster>::iterator it2REZ;
//	for( it2REZ = reservuars.begin(); it2REZ != reservuars.end(); ++it2REZ){
//		double start = it2REZ->get_st();
//		double koniec = it2REZ->get_end();
//		if( x0 > start and x0 < koniec ){								//node is in reservuar inside walls
//			IN_REZ=true;
//			break;
//		}
//		if( x0 == get_start_volume() or x0 == get_end_volume() ){								//node is in reservuar on external wall
//			IN_REZ=true;
//			break;
//		}
//	}

	if( node->get_rez_index() >=0 ){
		IN_REZ=true;
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
//	control_output<<"str "; node->show_site();
	wektor r0 = node->get_position();

	bool IN_VOLUME = check_x_belonging_volume(r0[dir]);

	if(IN_VOLUME){

	bool IN_REZ = false;

//	vector <plaster>::iterator it2REZ;
//	for( it2REZ = reservuars.begin(); it2REZ != reservuars.end(); ++it2REZ){
//		double start = it2REZ->get_st();
//		double koniec = it2REZ->get_end();
//		if( x0 >= start and x0 <= koniec ){								//node is in reservuar
//			IN_REZ=true;
//			break;
//		}	}
	
	if( node->get_rez_index() >=0 ){
		IN_REZ=true;
	}
	
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


bool opcja :: find_migration_path(site *first_node,int DIR, vector <site*> &migration_path){

	if(DEBUG or DEBUG_SMALL){control_output<<"find_mig:->";}

	bool MOVE_MIG = false;
	wektor kierunek;
	if(DIR>0){ 
		kierunek(1.0,0.0,0.0);
	}else if(DIR < 0){
		kierunek(-1.0,0.0,0.0);
	}else if(DIR== 0){
		kierunek(0.0,0.0,0.0);
	}else{
		control_output<<"ERROR in opcja""find_migration_path(). Wrong direction: "<<DIR<<endl;exit(1);
	}
//	kierunek.show();

	site* node=first_node;
	migration_path.push_back( first_node );
	bool WALL = false;
	
	if(DIR != 0){
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
		}else if(site_bufor.size() == 0 and vac_bufor.size()>0){				//Jesli w sim_area natrafisz na faze wakancyjna, por, klaster wakancji, to zakoncz sciezke.
			int rndIndex = rand() % vac_bufor.size();						
			int typ = (vac_bufor[rndIndex])->get_atom();
			if(typ==0){
				if(first_node->get_atom()>0){								
					migration_path.push_back( (vac_bufor[rndIndex]) );		//zapisz adres vacancy do sciezki jesli poczatek byl atom. Create vac.					
				}
				if(DEBUG){								//w node jest vacans stop while	
				control_output<<"Cluster Vac in sample for: "; node->show_site();
				for(unsigned int i = 0; i<migration_path.size();i++){
					migration_path[i]->show_site();
				}
				control_output<<"Vac from cluster "; node->show_site();
				node->show_neigh(1);
				}
				node=vac_bufor[rndIndex];	
				
				break;
			}
			else{
				control_output<<"ERROR in opcja::find_migration_path(). \
				Wrong atom type: "<<typ<<endl;exit(1);
			}
		}else if(site_bufor.size() == 0 and vac_bufor.size() == 0){			//Jesli path doszla do sciany?? Nie ma ani atomu ani wakancji do skoku. TO wylosuj site z rezerwuara. Jesli nie ma to przesun sample.
			if(DEBUG){control_output<<"WARRNING in opcja::find_migration_path(). ";
			control_output<<"Path reached a wall. No atoms or vacancy available to contiune path: "<<endl;}		
	
				MOVE_MIG=check_rezervuars(first_node, node);
				if(MOVE_MIG){
					migration_path.clear();
					break;
				}else{
					if(node != migration_path.back()){							
						migration_path.push_back( node );					//node is atom if first_node is vac. 
					}
					break;												//node is vac if first_node is atom. 
				}
		}else{
			control_output<<"ERROR in opcja::find_migration_path(). \
			Undefined condition: "<<site_bufor.size()<<" "<<vac_bufor.size() <<endl;
			exit(1);
		}	
		
	//migration_path.back()->show_site();	
	}while(node->get_atom()>0);
	}
		
	if(DEBUG){control_output<<"find mig. pth: "<<migration_path.size()<<endl;}

	if( (migration_path.back() == first_node) and (migration_path.size()==1) ){		//first_node (vacancy) is next to cluster of vacancies				
		int nr_rez = node->get_rez_index();
		if(nr_rez >=0 ){															//cluster of vacancies is located in rezervuar
			MOVE_MIG=check_rezervuars(first_node, node);
			if(MOVE_MIG){
				migration_path.clear();
			}else{
				if(node != migration_path.back()){							
					migration_path.push_back( node );					//node is atom if first_node is vac. 
				}else{
					control_output<<"ERROR:opcja::find_migration_path(): 1599"<<endl;exit(1); 
				}
			}
		}else{
			migration_path.clear();										//do nothing. Allow cluster of vacancies in sim area exists.
		}
	}
	
	if(DEBUG or DEBUG_SMALL){control_output<<"|->find_mig";}
	return MOVE_MIG;
}

void opcja :: dislocation_walk(vector <site*> &path){
	
	if(DEBUG or DEBUG_SMALL){control_output<<"dis_walk:->";}

	if(path.size()>0){
	if(DEBUG){	control_output<<"Path size... "<<path.size()<<endl;}
	site* first = path.front();
	site* last = path.back();
	int typ1 = first->get_atom();
	int typ2 = last->get_atom();
	
	if(path.size()==1){
		control_output<<"ERROR in opcja::dislocation_walk(): 1619"<<endl; exit(1);
	}
	vector <site*>::iterator point;
	//first->show_site();
	//last->show_site();

	//sciezka ma na poczatku wakancje (na koncu atom) -> remove wakancja
	if( typ1 == 0 and typ2 > 0){	//forrward vacancy out
		vector<site*>::iterator i = path.begin(); ++i;
		for ( ; i != path.end(); ++i ){
			vector<site*>::iterator prev = i; --prev;
			virtual_jump_vac_atom( (*prev), *i);	
//			(*prev)->show_site();
//			(*i)->show_site();
		}	
		reset_site( *(--i) );
		Vtoadd.insert( *i );	
//		(*i)->show_site();
	}else if( typ2 == 0 and typ1 > 0){	//backward vacancy in
		vector<site*>::reverse_iterator i = path.rbegin(); ++i;
		for ( ; i != path.rend(); ++i ){
			vector<site*>::reverse_iterator prev = i; --prev;
			virtual_jump_vac_atom( (*prev), *i);
//			(*prev)->show_site();
//			(*i)->show_site();
		}	
		reset_site( *(--i) );
		Vtoadd.insert( *i );	
//		(*i)->show_site();
	}else{
		control_output<<"ERROR in opcja::dislocation_walk(): "<<typ1<<" "<<typ2<<endl; 
		first->show_site();
		last->show_site();
		exit(1);
	}
	}
	if(DEBUG or DEBUG_SMALL){control_output<<"|->dis_walk";}

}

void opcja :: virtual_jump_vac_atom( site* VAC, site* ATOM){

	//one direction exchange of sites. Virtual jump of vac to atom.
//	control_output<<"Przed: "; VAC->show_site(); 
	//control_output<<"Przed: "; ATOM->show_site();
	if(VAC->get_atom() != 0){control_output<<"ERROR in opcja::virtual_jump(). Type of vacancy not 0: "<<VAC->get_atom()<<endl;exit(1);}
	if(ATOM->get_atom() <= 0){control_output<<"ERROR in opcja::virtual_jump(). Type of atom not >0: "<<ATOM->get_atom()<<endl;exit(1);}

	update_opcja(VAC,0);
	update_opcja(ATOM,0);

	site bufor(VAC);
	VAC->change_to( *ATOM );
	SAMPLE->update_events( VAC );
	ATOM->change_to(bufor);

	update_opcja(VAC,1);
	update_opcja(ATOM,1);

	call_flux_dislocation(ATOM,VAC);
//	control_output<<"Po: ";	VAC->show_site();	
//	control_output<<"Po: ";	ATOM->show_site();	
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

void opcja :: refresh_vac_list(){

	if(DEBUG  or DEBUG_SMALL){control_output<<"refresh_vac:->";}
	if(!MOVE_FRAME){
		if(SINGLE){
			if(DEBUG){	control_output<<"refresh_vac_vector "<<VAC_LIST.size()<<endl;}}
		SAMPLE->update_vac_list(Vtoadd, VAC_LIST);
		if(SINGLE){	if(DEBUG){control_output<<"|= "<<VAC_LIST.size()<<endl;}}
	}
	if(DEBUG  or DEBUG_SMALL){control_output<<"|->refresh_vac"<<endl;}

}

void opcja :: refresh_simarea(){
	
	if(MOVE_FRAME){														//bylo przesowanie blokow i rezerwuwarow
		control_output<<"opcja::refresh_sim_area "<<VAC_LIST.size()<<endl;
	
		bool NEW_SIM_AREA = SAMPLE->reinit_sim_area(del_L_sim, del_R_sim, VAC_LIST);			
		if(! NEW_SIM_AREA){
			refresh_vac_list();
		}
		//resetuj parametry kontrolne przesowania
		wektor a(0.0,0.0,0.0);
		del_L_sim=a;
		del_R_sim=a;
		Vtoadd.clear();
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

void opcja :: swap(plaster &source, plaster &destination, int FLAG){	//move all atoms (no vacancies) from source to destination
	
	control_output<<"swap:->";
	int count=0;
	unsigned int dir = source.get_direction();
	int P0 = destination.get_st();
	int P1 = destination.get_end();
	
	if(FLAG){control_output<<source.size()<<" "<<destination.size()<<endl;
		source.show();
		destination.show();

	}
	
	for(unsigned int i=0;i<source.size();i++){
		
		site* rnd_vac=0;
		site* rnd_at=0;
		rnd_at = source[i];
		double x = rnd_at->get_position(dir);
				
																		//	if(FLAG){control_output<<i<<" "<<P0<<" "<<x<<" "<<P1<<endl;}
		if( !( (P0 <= x) and (x < P1) ) ){		
			int typ = rnd_at->get_atom();
																		//		if(FLAG){control_output<<source.ref_to_sites[i]<<" "<<typ<<endl;}
			if(typ>0){
				long N =(long)(rnd()*(destination.size(0)) );			//		if(FLAG){control_output<<size(0)<<" "<<N<<endl;}
				rnd_vac=destination.get_site(0,N);
				double xv = rnd_vac->get_position(dir); bool common = false;
				if( (source.get_st() <= xv) and (xv<source.get_end()) ){
					common = true;
				}		
				count++;												//		if(FLAG){control_output<<"delete"<<endl;}
				destination.delete_site(0,N);
				if(common){
					source.plaster_delete_site( rnd_vac );				
				}
				rnd_vac->set_atom(typ);
				destination.add_site(typ,rnd_vac);	
				if(common){
					source.add_site(typ,rnd_vac);		
				}
				
				source.plaster_delete_site( rnd_at );							//no need because later destroyed
				rnd_at->set_atom(0);
				source.add_site(0,rnd_at);		
				
				SAMPLE->update_events(rnd_vac);
				SAMPLE->update_events(rnd_at);	
			}
		}
	}
	if(FLAG){
		source.show();
		destination.show();
	}
//	if sth wrong then you will see some atoms which stay out of sample
	for(unsigned int i=0;i<source.size();i++){
		site* rnd_at=0;
		rnd_at = source[i];
		double x = rnd_at->get_position(dir);				
		int typ = rnd_at->get_atom();
		if( (P0 <= x) and (x < P1) ){
			continue;
		}else{
			if(typ > 0){
				control_output<<"ERROR:opcja::swap():no all atoms moved to new rez"<<endl;
				source.show();
				destination.show();
				rnd_at->show_site();
				exit(1);
			}	
		}
	}

	control_output<<"|->swap";
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
	results.reserve(50);
	
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



