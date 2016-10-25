#include "site.h"



pairjump :: pairjump(){
		bariera=1000.0;
		vac_to_jump=NULL;
		atom_to_jump=NULL;
		E1=0.0;
		E2=0.0;
		bar=1000.0;
		}
pairjump :: pairjump(site* from, site* to,double e1,double e2,double bar_con,double bar_act){
		bariera=bar_act;
		E1=e1;
		E2=e2;
		bar=bar_con;
		vac_to_jump=from;
		atom_to_jump=to;
		}
		
pairjump :: ~pairjump(){}
	
double pairjump :: get_barier() const {
	return bariera;
}

double pairjump :: get_e1() const {
	return E1;
}

double pairjump :: get_e2() const {
	return E2;
}

site* pairjump :: get_vac_to_jump() const {
	return vac_to_jump;
}

site* pairjump :: get_atom_to_jump() const {
	return atom_to_jump;
}
	
void pairjump :: show() const {
	control_output<<"Skok: "<<E1<<" "<<E2<<" "<<bar<<" "<<bariera<<endl;
	vac_to_jump->show_site();
	atom_to_jump->show_site();
}
