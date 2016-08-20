#ifndef MYMATH_H
#define MYMATH_H
#include "mymath.h"
#endif


class task
{
	
	private:
	string task_name;
	vector <double> task_parameters; 
	//long task_main_step;
	
	
	public:
	task(string name, vector <double> &parameters)
	{
		task_name=name;
		task_parameters=parameters;
	};
		
	~task()
	{};
	
	void show_task()
	{
		control_output<<task_name<<" par: ";
		for(unsigned int i=0;i<task_parameters.size();i++)
		{
			control_output<<task_parameters[i]<<" ";
			}
		control_output<<endl;	
	};
	
	string get_name()
	{
		return task_name;
		};
	
	void get_parameters(vector <double> &parameters)
	{
		parameters=task_parameters;
		};
		
	long get_main_step()
	{
		return task_parameters[0];
		//return task_main_step;
		};	
		
	};
