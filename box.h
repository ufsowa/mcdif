#ifndef MYMATH_H
#define MYMATH_H
#include "mymath.h"
#endif


using namespace std;

class box
{
	private:
	vector <site> sity;
	long box_id;
	
	public:
	
	box(){
		sity.reserve(4);
//		cout<<"Pojemonosc tablicy rezerwuje: "<<sity.capacity()<<" "<<&sity<<" ";
		};
	box(long _id){
		sity.reserve(4);
		box_id=_id;
//		cout<<"Pojemonosc tablicy rezerwuje: "<<sity.capacity()<<" "<<&sity<<" ";
		};	
	~box(){};
	
	box(const box &i)
	{	
		sity.clear();
		sity.reserve(4);
		sity=i.sity;
		box_id=i.box_id;
	};
	
	box& operator= (const box &i)
	{
		if(this != &i)
		{
			sity.clear();
			sity.reserve(4);
			sity=i.sity;
			return *this;
		}
		return *this;
	};
	
	void set_boxid(long iter)
	{
		box_id=iter;
		};
	
	long get_box_id()
	{
		return box_id;
	};
	
	void clear_sity()
	{
	//	control_output<<"sity_size "<<sity.size()<<endl;
		sity.clear();
		
	//	control_output<<"clear sity_size "<<sity.size()<<endl;
		}
	
	void put_site(site &Site)
	{
	//	control_output<<"Dodalem sita bdfb do boxa "<<sity.size()<<endl;
	//	cout<<"Adres boxu: "<<this<<", tablicy sitow w boxie: "<<&sity<<" pojemonosc tablicy przed: "<<sity.capacity()<<endl;
	//	cout<<" Adres dodawanego sita w put_site: "<<&Site<<endl;
		int add =1;  
	//	int addd=0;
		vector <site> ::iterator I; 
		
		if(sity.size()==0)
		{sity.push_back(Site);}
		else
		{
		
		//for(int j=0;j<sity.size();j++)
		
		for(I=sity.begin();I!=sity.end();I++)
		{
		//	control_output<<"W for"<<endl;
			//sity[j]
			if ( site(*I) == site(Site) )
			{
				//control_output<<"replace:"<<endl;
				
				//for(int i=0;i<sity.size();i++)
				//{
				//	sity[i].show_site();
				//}
				
				sity.erase(I);
				//sity.erase(sity.begin()+j);
				//control_output<<"errase:"<<endl;
				//for(int i=0;i<sity.size();i++)
				//{
				//	sity[i].show_site();
				//}
				sity.push_back(Site);
				//control_output<<"adding:"<<endl;
				//for(int i=0;i<sity.size();i++)
				//{
				//	sity[i].show_site();
				//}
				add=0;
				
			}
			
		
		}
		
		if(add)
		{
	//		control_output<<"addd"<<endl;
			sity.push_back(Site);}
		
		}
	//	cout<<"Adres tablicy sitow w boxie po: "<<&sity<<" pojemonosc tablicy po: "<<sity.capacity()<<endl;
	//	cout<<"Adresy zapisanych sitow w boxie: "<<endl; 
	//	for(int i=0;i<sity.size();i++)
	//	{
	//		cout<<i<<" "<<&sity[i]<<" ";
	//		sity[i].show_site();
	//	}
	};
	
		void put_atom(site &Site)
	{
	//	control_output<<"Dodalem sita bdfb do boxa "<<sity.size()<<endl;
		
		int add =1;  
//		int addd=0;
//		int atom=Site.get_atom();
		vector <site> ::iterator I; 
		
		if(sity.size()==0)
		{
			sity.push_back(Site);
		}
		else
		{
		
		//for(int j=0;j<sity.size();j++)
		
		for(I=sity.begin();I!=sity.end();I++)
		{
//			control_output<<"W for"<<endl;
			//sity[j]
			if ( site(*I) == site(Site) )
			{
		//		cout<<"replace site at: ";
		//		I->show_site();
				
				//for(int i=0;i<sity.size();i++)
				//{
				//	sity[i].show_site();
				//}
				
//				I->set_atom(atom);
				*I=Site;
				
		//		cout<<"to site: ";
		//		I->show_site();
				
				//sity.erase(sity.begin()+j);
				//control_output<<"errase:"<<endl;
				//for(int i=0;i<sity.size();i++)
				//{
				//	sity[i].show_site();
				//}
				//sity.push_back(Site);
				//control_output<<"adding:"<<endl;
				//for(int i=0;i<sity.size();i++)
				//{
				//	sity[i].show_site();
				//}
				add=0;		
			}
		}
		
		if(add)
		{
	//		control_output<<"addd"<<endl;
			sity.push_back(Site);
		}
	}
		
		
	};
	
	
	void show_sity()
	{
		for(unsigned int i=0;i<sity.size();i++)
		{
		control_output<<i<<" "<<&sity[i]<<" ";
		sity[i].show_site();
		}
	};
	
	void get_sity_in_box(vector <site*> &tmp_atom_list)
	{	//int o;
		//cin>>o;	
	//	control_output<<"sity size: "<<sity.size()<<endl;
	//	cin>>o;
		for(unsigned int i=0;i<sity.size();i++)
		{
		site* wsk_sita=0;
		wsk_sita=&sity[i];
	//	control_output<<"Getting sites at: "<<wsk_sita<<endl;
		tmp_atom_list.push_back(wsk_sita);
		
		}
		
		};
	
};
