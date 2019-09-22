#include<fstream>
#include<string>
#include<iostream>
#include<vector>
#include"TMatrix.h"
#include"TMatrixD.h"
#include"TArrayD.h"
using namespace std;

void read_matrix(const char* filepwd,vector<string>& v_horizontal,vector<string>& v_vertical,TMatrixD& a)
{

	ifstream inputfile(filepwd);
	string line;

	if (!inputfile){cout<<"no such file"<<endl;}
	if (inputfile)
	{
		int n_line=0;
		bool horizontal_op = false;
		bool vertical_op = false;
		bool matrix_found = false;
		int position;
		

		TArrayD data;
		
		Int_t num_horizontal=0;
		Int_t num_vertical=0;
		Int_t matrix_row_counter = 0;
		Int_t matrixelementcounter = 0;

		while (getline(inputfile,line))
		{
			n_line++;
			stringstream ss(line);
			//cout<<n_line<<endl;
			position = line.find("//");
			//cout<<position<<endl;
			if(position==0){
				if(line.find("horizontal")!=line.npos){horizontal_op=true;}
				if(line.find("vertical")!=line.npos){vertical_op=true;}
				if(line.find("Matrix")==2){matrix_found=true;}
			}
			
			if(horizontal_op & (position!=0)){
				string str;
				
				while(getline(ss,str,',')){
					num_horizontal++;
					v_horizontal.push_back(str);					
				}

				horizontal_op =false;
				//cout<<num_horizontal<<endl;
			}

			if(vertical_op & (position!=0)){
				string str;
				
				while(getline(ss,str,',')){
					num_vertical++;
					v_vertical.push_back(str);
				}

				vertical_op =false;
				//cout<<num_vertical<<endl;
			}

			data.Set(num_horizontal*num_vertical);
			
			if(matrix_found & (position!=0)){
				string str;
				Int_t matrix_col_counter = 0;
				while(getline(ss,str,',')){
					matrix_col_counter++;
					stringstream iss(str);
					Double_t matrixelement;
					iss>>matrixelement;
					
					data[matrixelementcounter]=matrixelement;
					matrixelementcounter++;
				}
				if(matrix_col_counter==0){matrix_found=false;continue;}
				matrix_row_counter++;
			}
			
		
			
		}
		a.Use(num_vertical,num_horizontal,data.GetArray());
	
		//cout<<a.operator()(0,0)<<endl;
		for(int i =0;i<a.GetNrows();i++){
        for(int j =0;j<a.GetNcols();j++){cout<<a.operator()(i,j)<<" ";}
        cout<<endl;
    }
	cout<<&a<<endl;
	}
}