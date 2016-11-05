
class site;

class pairjump
{
private:
double bariera,E1,E2,bar;
site* vac_to_jump;
site* atom_to_jump;

public:

pairjump();
pairjump(site* from, site* to,double e1,double e2,double bar_con,double bar_act);
~pairjump();
double get_barier() const;
double get_e1() const;
double get_e2() const;
site* get_vac_to_jump() const;
site* get_atom_to_jump() const ;
void show() const;

};
