#include "site.h"



pairjump :: pairjump(){
		bariera=1000.0;
		vac_to_jump=NULL;
		atom_to_jump=NULL;
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
	
double pairjump :: get_barier(){
	return bariera;
}

site* pairjump :: get_vac_to_jump(){
	return vac_to_jump;
}

site* pairjump :: get_atom_to_jump(){
	return atom_to_jump;
}
	
void pairjump :: show(){
	control_output<<"Skok: "<<E1<<" "<<E2<<" "<<bar<<" "<<bariera<<endl;
	vac_to_jump->show_site();
	atom_to_jump->show_site();
}
