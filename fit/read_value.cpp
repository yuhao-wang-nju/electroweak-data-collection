#include<fstream>
#include<string>
#include<iostream>
#include<vector>
#include"TMatrix.h"
#include"TArrayD.h"
using namespace std;

void read_value(const char* filepwd,std::vector<string>& v_observable,std::vector<Double_t>& v_value,std::vector<Double_t>& v_error)
{

	ifstream inputfile(filepwd);
	string line;

	if (!inputfile){cout<<"no such file"<<endl;}
	if (inputfile)
	{
		int n_line=0;
		int position;		
		Int_t num_observable=0;

		while (getline(inputfile,line))
		{
			n_line++;
			stringstream ss(line);

			position = line.find("//");
            if(position==0 || line.size()==0){continue;}

            string str;
            int col_index = 0;
            while(getline(ss,str,',')){
                col_index++;
                stringstream iss(str);
                switch (col_index)
                {
                case 1:
                    v_observable.push_back(str);
                    break;
                case 2:
                    
                    Double_t value;
                    iss>>value;
                    v_value.push_back(value);
                    break;
                case 3:
                    
                    Double_t error;
                    iss>>error;
                    v_error.push_back(error);
                default:
                    break;
                }
            }
		
			
		}

	
		//cout<<v_observable.at(19)<<" "<<v_value.at(19)<<" "<<v_error.at(19)<<endl;

	}
}