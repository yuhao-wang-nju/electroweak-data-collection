#include "TMatrixDBase.h"
#include "TMatrixT.h"
#include "TMatrixDSparse.h"
#include "TMatrixDEigen.h"
#include"TArrayI.h"
#include "TDecompSVD.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include <string>
#include <vector>
#include <sstream>
#include "TMatrixD.h"
#include "TROOT.h"
using namespace std;

void read_value(const char* filepwd,std::vector<string>& v_observable,std::vector<Double_t>& v_value,std::vector<Double_t>& v_error)
{
	v_observable.resize(0);
	v_value.resize(0);
	v_error.resize(0);
	
	ifstream inputfile(filepwd);
	string line;

	if (!inputfile){cout<<"no such file named: "<<filepwd<<endl;}
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

void read_matrix(const char* filepwd,vector<string>& v_horizontal,vector<string>& v_vertical,TMatrixD& mtx)
{
	v_horizontal.resize(0);
	v_vertical.resize(0);
	ifstream inputfile(filepwd);
	string line;
	//cout<<"Hello?"<<endl;
	if (!inputfile){cout<<"no such file named: "<<filepwd<<endl;}
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
		mtx.ResizeTo(num_vertical,num_horizontal);
		mtx.SetMatrixArray(data.GetArray());
		//cout<<mtx.IsValid()<<endl;
		//cout<<a.operator()(0,0)<<endl;
		// for(int i =0;i<a.GetNrows();i++){
        // for(int j =0;j<a.GetNcols();j++){cout<<a.operator()(i,j)<<" ";}
        // cout<<endl;
    //}
	//cout<<&a<<endl;
	}
}

template <class _Tp, class _Alloc>
int vector_compare(vector<_Tp, _Alloc>& v_1,vector<_Tp, _Alloc>& v_2 ){
//This function will compare if two vector has the same content
//if they are all the same, it will return 0;
//if they have different size, it will return -1
//if they have different element, it will return the position of first difference apperars

    bool isitsame=true;
    int i=0;
    if(v_1.size()!=v_2.size()){
        cout<<"two vector do not have same size!"<<"("<<v_1.size()<<"and"<<v_2.size()<<")"<<endl;
        i=-1;
        isitsame=false;
    }
    else{
        while(i<v_1.size() && isitsame){
            isitsame=(v_1.at(i)==v_2.at(i));
            i++;
        }
        if(!isitsame){
        cout<<"The two vector do not have same element, begin at position"<<i<<endl;
        //cout<<v_1[i]<<" "<<v_2[i]<<endl;
    }
    }   
    if(isitsame){i=0;}
    return i;
    
}

template <class Element>
void combine_matrix(TMatrixTBase<Element>& matrix_1,  TMatrixTBase<Element>& matrix_2 ){
	Int_t nrow_1=matrix_1.GetNrows();
	Int_t ncol_1=matrix_1.GetNcols();
	Int_t nrow_2=matrix_2.GetNrows();
	Int_t ncol_2=matrix_2.GetNcols();

	matrix_1.ResizeTo(nrow_1+nrow_2,ncol_1+ncol_2);
	matrix_1.SetSub(nrow_1,ncol_1,matrix_2);
}

void chi2_fit(TMatrixD& mtx_op_obs,TMatrixD& mtx_cor,TMatrixD& y_forfit,TMatrixD& theta,TMatrixD& cov_mtx){
    TMatrixD mtx_cor_inv(TMatrixD::kInverted,mtx_cor);
    TMatrixD m2(mtx_op_obs*mtx_cor_inv,TMatrixD::kMultTranspose,mtx_op_obs);
    //m2.Print(); 
    TDecompSVD msvd(m2);
    msvd.Decompose();
    cov_mtx.ResizeTo(mtx_op_obs.GetNrows(),mtx_op_obs.GetNrows());
    cov_mtx=msvd.Invert();
    
    theta.ResizeTo(mtx_op_obs.GetNrows(),1);
    theta=cov_mtx*mtx_op_obs*mtx_cor_inv*y_forfit;
    
}


void selectrow(TMatrixD& aftermtx,TMatrixD& beforemtx,std::vector<Int_t>& v_selectrow){
    aftermtx.ResizeTo(v_selectrow.size(),beforemtx.GetNcols());
    for(int i =0;i<v_selectrow.size();i++){
        TMatrixDRow(aftermtx,i)=TMatrixDRow(beforemtx,v_selectrow[i]);
    }
}
void selectcol(TMatrixD& aftermtx,TMatrixD& beforemtx,std::vector<Int_t>& v_selectcol){
    aftermtx.ResizeTo(beforemtx.GetNrows(),v_selectcol.size());
    for(int i =0;i<v_selectcol.size();i++){
        TMatrixDColumn(aftermtx,i)=TMatrixDColumn(beforemtx,v_selectcol[i]);
    }
}

template <class _Tp, class _Alloc>
void vectorselect(vector<_Tp, _Alloc>&v_after,vector<_Tp, _Alloc>&v_before,vector<int>& v_select){
	v_after.resize(0);
	for(int i =0;i<v_select.size();i++){
		v_after.push_back(v_before.at(v_select[i]));
	}
	
}

