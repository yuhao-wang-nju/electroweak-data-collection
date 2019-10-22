#ifndef DEFINE
#define DEFINE

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
#include <iomanip>
#include "TMatrixD.h"
#include "TROOT.h"
//#include "read_matrix.cpp"
#include "function_for_fit.cpp"
using namespace std;


class EWFit
{
    public:
    
    vector<string> v_operator;
    vector<string> v_observable;
    TMatrixD correlation_mtx;
    TMatrixD contribution_mtx;
    vector<double_t> v_smvalue;
    vector<double_t> v_smerror;
    vector<double_t> v_expvalue;
    vector<double_t> v_experror;
    string energyregion;
    Int_t statuscode;

    TMatrixD exp_sm_diff;
    TMatrixD LStheta;
    TMatrixD covariance_matrix;
    TMatrixD LSerror;

    string folderpwd;
    string correlation_pwd="Correlation_matrix.txt";
    string contribution_pwd="Operator-Observable-alphascheme-Matrix.txt";
    string smpwd="SM_prediction.txt";
    string exppwd="Experimental_value.txt";

    int readall(){
        //This function will read all the file under the folder, save it to the class EWFit and will return a status code
        //status code=11111, no folder pwd is assigned
        //status code=10000, correlation matrix is not a square matrix
        //status code=01000, num of observable in contribution matrix is smaller than in correlation matrix
        //status code=00100, num of observable in contribution matrix is bigger than in correlation matrix, it will replace the v_observable in struct
        //status code=00010, SM prediction observable is not as same as in contribution matrix
        //status code=00001, experimental observable is not as same as in contribution matrix
        //The final status code is the bitwise AND result of all the above status code
        
        int kk=0b00000;
        cout<<endl<<"#########################"<<endl<<"Begin Reading Folder!"<<endl;
    if(folderpwd==""){cout<<"please assign the folder pwd first"<<endl;kk=0b11111;}
    if(folderpwd!=""){
        cout<<"folder: "<<folderpwd<<" found! Reading..."<<endl;
        stringstream ssnew;
        string str;
        vector<string> v_tmp;
        ssnew<<folderpwd<<correlation_pwd;
        str=ssnew.str();
        cout<<"correlation matrix Reading..."<<endl;
        const char* file_pwd=str.c_str();
        read_matrix(file_pwd,v_observable,v_tmp,correlation_mtx);
        int k1=vector_compare(v_observable,v_tmp);
        if(k1!=0){k1=0b10000;}
        
        cout<<"contribution matrix Reading..."<<endl;
        ssnew.str("");
        ssnew<<folderpwd<<contribution_pwd;
        str=ssnew.str();
        file_pwd=str.c_str();
        read_matrix(file_pwd,v_tmp,v_operator,contribution_mtx);
        int k2=vector_compare(v_observable,v_tmp);
        if(k2>0){k2=0b01000;}
        int sizediff=k2*(v_observable.size()-v_tmp.size());
        if(sizediff>0){v_observable.resize(0);v_observable.swap(v_tmp);k2=0b00100;}
        //cout<<v_operator.size()<<endl;

        cout<<"Standard Model Prediction Reading..."<<endl;  
        ssnew.str("");
        ssnew<<folderpwd<<smpwd;
        str=ssnew.str();
        file_pwd=str.c_str();
        read_value(file_pwd,v_tmp,v_smvalue,v_smerror);
        int k3=vector_compare(v_observable,v_tmp);
        //cout<<v_tmp.at(1)<<endl;
        if(k3!=0){k3=0b00010;}
        
        cout<<"Experimental Value Reading..."<<endl;
        ssnew.str("");
        ssnew<<folderpwd<<exppwd;
        str=ssnew.str();
        file_pwd=str.c_str();
        read_value(file_pwd,v_tmp,v_expvalue,v_experror);
        int k4=vector_compare(v_observable,v_tmp);      
        if(k4!=0){k4=0b00001;}
        else{exp_sm_diff.ResizeTo(v_observable.size(),1);
        for(int i =0;i<v_observable.size();i++){
            exp_sm_diff(i,0)=v_expvalue.at(i)-v_smvalue.at(i);
        }
        }
        //cout<<v_expvalue.size()<<endl;
        kk=k1|k2|k3|k4;
    }
    cout<<"Reading Ends!"<<endl;
    cout<<"Matrix Reading Status Code:"<<bitset<5>(kk)<<endl;
    cout<<"############################"<<endl<<endl;
    statuscode=kk;
    return kk;
}

void complete_correlation_withI(){
    //This function will complete the correlation which is smaller than observable size with Identity matrix
    if(statuscode!=0b00100&&statuscode!=0b00000){cout<<"There is other problem in the EWfit data, Please check first!";}
    if(statuscode==0b00100){
        Int_t size_I=contribution_mtx.GetNcols()-correlation_mtx.GetNcols();
        TMatrixDSparse unit1(size_I,size_I);
        TArrayI row(size_I),col(size_I);
        for(int i =0;i<size_I;i++){row[i]=col[i]=i;}
        TArrayD data(size_I);
        data.Reset(1.);
        unit1.SetMatrixArray(size_I,row.GetArray(),col.GetArray(),data.GetArray());

        combine_matrix(correlation_mtx,unit1);
        statuscode=0b000000;
    }
}

int check_size(){
    
    int i=v_observable.size();
    int j=v_operator.size();

    int k=(i==correlation_mtx.GetNrows())&&(i==correlation_mtx.GetNcols())&&(i==contribution_mtx.GetNcols());
    int q=(i==v_smvalue.size())&&(i==v_expvalue.size());
    int p=(j==contribution_mtx.GetNrows());

    return k&&p&&q;

}

void print_size(){
    cout<<"Size Information:"<<endl;
    cout<<"Correlation Matrix:"<<correlation_mtx.GetNrows()<<"*"<<correlation_mtx.GetNcols()<<endl;
    cout<<"Contribution Matrix:"<<contribution_mtx.GetNrows()<<"*"<<contribution_mtx.GetNcols()<<endl;
    cout<<"exp_sm_diff vector:"<<exp_sm_diff.GetNrows()<<"*"<<exp_sm_diff.GetNcols()<<endl;

}

void LSFit(){
    //This function will do the Least Square Fit and save the result to class attributes
    int p=check_size();
    if(p==1){chi2_fit(contribution_mtx,correlation_mtx,exp_sm_diff,LStheta,covariance_matrix);
    LSerror.ResizeTo(LStheta.GetNrows(),1);
    for(int i =0;i<LStheta.GetNrows();i++){
        LSerror(i,0)=sqrt(covariance_matrix(i,i));
    }}
    else{cout<<"It can not fit due to the incompatible size! Please check"<<endl;}
}

EWFit SetSub(vector<int>& v_op,vector<int>& v_ob){
    EWFit Sub;
    TMatrixD subcorr1,subcon1;
    
    selectrow(subcorr1,correlation_mtx,v_ob);
    selectcol(Sub.correlation_mtx,subcorr1,v_ob);
    
    selectrow(subcon1,contribution_mtx,v_op);
    selectcol(Sub.contribution_mtx,subcon1,v_ob);

    selectrow(Sub.exp_sm_diff,exp_sm_diff,v_ob);
    
    vectorselect(Sub.v_operator,v_operator,v_op);
    vectorselect(Sub.v_observable,v_observable,v_ob);
    vectorselect(Sub.v_expvalue,v_expvalue,v_ob);
    vectorselect(Sub.v_experror,v_experror,v_ob);
    vectorselect(Sub.v_smvalue,v_smvalue,v_ob);
    vectorselect(Sub.v_smerror,v_smerror,v_ob);

    

    int k=Sub.check_size();
    if(k==0){cout<<"There may be some mistake in the Sub!";}
    else{Sub.statuscode=0b00000;}
    return Sub;

}
EWFit SetSub(vector<string>&v_op, vector<string>& v_ob){
    vector<int> v_1;
    vector<int> v_2;
    for(int i=0;i<v_ob.size();i++){
        for(int j =0;j<v_observable.size();j++){
            //cout<<v_ob[i]<<v_observable[j]<<endl;
            if(v_ob.at(i)==v_observable.at(j)){v_1.push_back(j);}
        }
    }
    for(int i=0;i<v_op.size();i++){
        for(int j =0;j<v_operator.size();j++){
            if(v_op.at(i)==v_operator.at(j)){v_2.push_back(j);}
        }
    }
 //cout<<v_2.size()<<v_1.size()<<endl;
    EWFit Sub=SetSub(v_2,v_1);
    return Sub;
}

template <typename _Tp,unsigned N,unsigned N2>
EWFit SetSub(_Tp(& arrop)[N],_Tp(& arrob)[N2]){
    vector<_Tp> v_op,v_ob;
    v_op.resize(sizeof(arrop)/sizeof(_Tp));
    for(int i=0;i<v_op.size();i++){v_op[i]=arrop[i];}
    v_ob.resize(sizeof(arrob)/sizeof(_Tp));
    for(int i=0;i<v_ob.size();i++){v_ob[i]=arrob[i];}
    EWFit Sub=SetSub(v_op,v_ob);
    return Sub;
}

TEllipse* GetEllipse(){
    TEllipse* elp;
    if(covariance_matrix.GetNcols()!=2){cout<<"It is not a two-dimension fit, please check the size!"<<endl;}
    else{
        TMatrixDEigen cov_mtxe(covariance_matrix);
        TMatrixD cov_egvalue=cov_mtxe.GetEigenValues();
        TMatrixD cov_egvector=cov_mtxe.GetEigenVectors();
        elp=new TEllipse(LStheta(0,0),LStheta(1,0),sqrt(cov_egvalue(0,0)),sqrt(cov_egvalue(1,1)),0,360,atan2(cov_egvector(1,0),cov_egvector(0,0))/3.1415926*180);
    }
    return elp;
}

EWFit RuleOut0Contribution(){
    vector<Int_t> v_nonzerorow;
    vector<Int_t> v_zerorow;
    vector<Int_t> v_allcol;
    for(int i=0;i<contribution_mtx.GetNrows();i++){
        TVectorD vrow=TMatrixDRow(contribution_mtx,i);
        if(vrow.NonZeros()!=0){v_nonzerorow.push_back(i);}
        else {v_zerorow.push_back(i);}
    }
    for(int i=0;i<contribution_mtx.GetNcols();i++){v_allcol.push_back(i);}
    //cout<<v_allcol.size()<<endl;
    EWFit nonzero=SetSub(v_nonzerorow,v_allcol);

    return nonzero;
}

EWFit CombineEWFit(EWFit info2){
    vector<int> rearrange;
    vector<int> info2ob;
    EWFit Combine;
    bool compatible=true;
    
    for(int i=0;i<v_operator.size();i++){
        bool Isexist=false;
        for(int j=0;j<info2.v_operator.size();j++){
            if(v_operator.at(i)==info2.v_operator.at(j)){rearrange.push_back(i);Isexist=true;}
        }
        if(!Isexist){cout<<"The two energy do not have same operators, Please check!"<<endl;}
        compatible*=Isexist;
    }
    for(int i=0;i<info2.v_observable.size();i++){info2ob.push_back(i);}

    EWFit reinfo2=info2.SetSub(rearrange,info2ob);
    Combine.contribution_mtx.ResizeTo(v_operator.size(),v_observable.size()+info2.v_observable.size());
    Combine.contribution_mtx.SetSub(0,0,contribution_mtx);
    Combine.contribution_mtx.SetSub(0,v_observable.size(),info2.contribution_mtx);
    
    Combine.correlation_mtx.ResizeTo(v_observable.size()+info2.v_observable.size(),v_observable.size()+info2.v_observable.size());
    Combine.correlation_mtx.SetSub(0,0,correlation_mtx);
    Combine.correlation_mtx.SetSub(v_observable.size(),v_observable.size(),info2.correlation_mtx);
    
    Combine.v_observable.insert(Combine.v_observable.end(),v_observable.begin(),v_observable.end());
    Combine.v_observable.insert(Combine.v_observable.end(),info2.v_observable.begin(),info2.v_observable.end());
    
    Combine.v_operator=v_operator;
    
    Combine.v_smvalue.insert(Combine.v_smvalue.end(),v_smvalue.begin(),v_smvalue.end());
    Combine.v_smvalue.insert(Combine.v_smvalue.end(),info2.v_smvalue.begin(),info2.v_smvalue.end());

    Combine.v_smerror.insert(Combine.v_smerror.end(),v_smerror.begin(),v_smerror.end());
    Combine.v_smerror.insert(Combine.v_smerror.end(),info2.v_smerror.begin(),info2.v_smerror.end());

    Combine.v_expvalue.insert(Combine.v_expvalue.end(),v_expvalue.begin(),v_expvalue.end());
    Combine.v_expvalue.insert(Combine.v_expvalue.end(),info2.v_expvalue.begin(),info2.v_expvalue.end());
    
    Combine.v_experror.insert(Combine.v_experror.end(),v_experror.begin(),v_experror.end());
    Combine.v_experror.insert(Combine.v_experror.end(),info2.v_experror.begin(),info2.v_experror.end());
    Combine.exp_sm_diff.ResizeTo(v_observable.size()+info2.v_observable.size(),1);
    Combine.exp_sm_diff.SetSub(0,0,exp_sm_diff);
    Combine.exp_sm_diff.SetSub(v_observable.size(),0,info2.exp_sm_diff);
    
    
    return Combine;
}
void print_result(){
    cout<<endl<<"#####################"<<endl<<"Fitting result:"<<endl;
    if(LStheta.GetNrows()==0){cout<<"Need to fit first!"<<endl;}
    else{
        for(int i=0;i<LStheta.GetNrows();i++){
            cout<<left<<setw(4)<<v_operator.at(i)<<":";
            cout<<left<<setw(10)<<LStheta(i,0)<<"+-";
            cout<<LSerror(i,0)<<endl;
        }
    }
}





};


#endif