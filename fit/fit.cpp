
#include "read_matrix.cpp"
#include "read_value.cpp"
#include "TMatrixDBase.h"
void fit(){
    const char* filepwd = "/home/allenwang/Documents/Files/EWPM/analysis_tools/electro-weak-github/Z_pole/Correlation_matrix.txt";
    vector<string> v_cor_horizontal;
    vector<string> v_cor_vertical;
    TMatrixD cor_matrix;
    read_matrix(filepwd,v_cor_horizontal,v_cor_vertical,cor_matrix);

    
    //cout<<&cor_matrix<<endl;
    //cout<<cor_matrix.GetNrows()<<" "<<cor_matrix.GetNcols()<<endl;
    //cout<<cor_matrix.operator()(0,0)<<endl;
    const char* smpwd = "../Z_pole/Experimental_value.txt";
    vector<string> v_smobservable;
    vector<Double_t> v_smvalue;
    vector<Double_t> v_smerror;
    read_value(smpwd,v_smobservable,v_smvalue,v_smerror);

    const char* op_ob_pwd= "/home/allenwang/Documents/Files/EWPM/analysis_tools/electro-weak-github/Z_pole/Operator-Observable-alphascheme-Matrix.txt";
    vector<string> v_op_ob_horizontal;
    vector<string> v_op_ob_vertical;
    TMatrixD op_ob_matrix;
    read_matrix(op_ob_pwd,v_op_ob_horizontal,v_op_ob_vertical,op_ob_matrix);
    cout<<op_ob_matrix.GetNrows()<<" "<<op_ob_matrix.GetNcols()<<endl;
    
    //cout<<op_ob_matrix.operator()(0,3)<<endl;
    // for(int i =0;i<op_ob_matrix.GetNrows();i++){
    //     for(int j =0;j<op_ob_matrix.GetNcols();j++){cout<<op_ob_matrix.operator()(i,j)<<" ";}
    //     cout<<endl;
    // }
    

    Int_t num_without_cor = op_ob_matrix.GetNcols()-cor_matrix.GetNcols();
    TMatrixD a(op_ob_matrix.GetNcols(),op_ob_matrix.GetNcols());
    a.Zero();
    cout<<a.GetNrows()<<a.GetNcols()<<endl;
    //Create a unit matrix for the observables without correlation
    TMatrixDSparse unit1(num_without_cor,num_without_cor);
    TArrayI row(num_without_cor),col(num_without_cor);
    for (Int_t i = 0; i < num_without_cor; i++) row[i] = col[i] = i;
    TArrayD data(num_without_cor); data.Reset(1.);
    unit1.SetMatrixArray(num_without_cor,row.GetArray(),col.GetArray(),data.GetArray());
    // for(int i =0;i<num_without_cor;i++){
    //     for(int j =0;j<num_without_cor;j++){cout<<unit1.operator()(i,j);}
    //     cout<<endl;
    // }

    a.SetSub(0,0,cor_matrix);
    a.SetSub(cor_matrix.GetNcols(),cor_matrix.GetNcols(),unit1);
    for(int i =0;i<a.GetNrows();i++){
        for(int j =0;j<a.GetNcols();j++){cout<<a.operator()(i,j)<<" ";}
        cout<<endl;
    }
}