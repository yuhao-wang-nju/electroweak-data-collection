This is the document of the fit function
A class named EWFit is created to help do the fit
A ROOT version of 6.18 is needed, older version may not be supported
Least Square Fit Method is used

Here is some information of the class EWFit:

class EWFit
{
public Attributes
    vector<string> v_operator;  //name list of the operator
    vector<string> v_observable;    //name list of the observable
    TMatrixD correlation_mtx;   //corrlation matrix between observable
    TMatrixD contribution_mtx;  //contribution matrix(operator-observable)
    vector<double_t> v_smvalue; //Standard Model prediction of observable
    vector<double_t> v_smerror; //SM prediction error
    vector<double_t> v_expvalue;    //Experimental value of observable
    vector<double_t> v_experror;    //Experimental error of observable
    string energyregion;    //energy region information
    Int_t statuscode;   //status code returned in the readall function

    TMatrixD exp_sm_diff;   //Experimental Value minus SM value(used in fit)
    TMatrixD LStheta;   //the predictor in Least Square Fit Method(fit result)
    TMatrixD covariance_matrix; //Covariance Matrix of the fit result
    TMatrixD LSerror;   //error of the fit result

    string folderpwd;   //The folder path of the data
    string correlation_pwd="Correlation_matrix.txt";    //correlation matrix data file name
    string contribution_pwd="Operator-Observable-alphascheme-Matrix.txt";   //contribution matrix data file name
    string smpwd="SM_prediction.txt";   //SM prediction data file name
    string exppwd="Experimental_value.txt"; //Experimental value data file name

public Member Function
    int readall()
        //This function will read all the file under the folder, save it to the class EWFit and will return a status code
        //status code=11111, no folder pwd is assigned
        //status code=10000, correlation matrix is not a square matrix
        //status code=01000, num of observable in contribution matrix is smaller than in correlation matrix
        //status code=00100, num of observable in contribution matrix is bigger than in correlation matrix, it will replace the v_observable in struct
        //status code=00010, SM prediction observable is not as same as in contribution matrix
        //status code=00001, experimental observable is not as same as in contribution matrix
        //The final status code is the bitwise AND result of all the above status code
    int check_size()
        //This function will check if the size of the matrixs and vectors are compatible for fit, 1 for compatible, 0 for not
    void print_size()
        //This function will print the size of CorrelationMatrix,ContributionMatrix and exp_sm_diff (num_rows*num_cols)
    void complete_correlation_withI()
        //This function will complete the correlation matrix which is smaller than observable size with Identity matrix
    EWFit RuleOut0Contrinbution()
        //This function will rule out the operator which do not have any contribution to any observable and return a new EWFit
    EWFit CombinEWFit(EWFit info2)
        //This function will combine two EWFit class together, with the info2 attached after the applied EWFit, the operator number of two EWFit should have same size, the order of operator will follow the applied EWFit and the correlation matrix should have been filled before this function is used
    void LSFit()
        //This function will do the Least Square Fit and save the result to class attributes
    EWFit SetSub(vector<_Tp>& v_op,vector<_Tp>& v_ob)
    EWFit SetSub(_Tp(& arrop)[N],_Tp(& arrob)[N2])
        //This funtion will return a Sub EWFit class of the current one
        //v_op is the vector of operator is selected, v_ob is the observable
        //two kind of type is supported: int and string(you can use the position/name of the operator/observable)
    TEllipse* GetEllipse()
        //If the fit uses only two operator, this function could return a Ellipse, with the first operator as x-axis and second as y-axis
    void print_result()
        //print the fitting result

}

Other function in the function_for_fit.cpp

    void read_value(const char* filepwd,std::vector<string>& v_observable,std::vector<Double_t>& v_value,std::vector<Double_t>& v_error)
        //This function will read the vector-like value(such as SM prediction) file and save the result to the input vectors
    void read_matrix(const char* filepwd,vector<string>& v_horizontal,vector<string>& v_vertical,TMatrixD& mtx)
        //This function will read the matrix-like value(such as correlation matrix) file and save the result to the input TMatrixD, with two vector as the name list of matrix 
    int vector_compare(vector<_Tp, _Alloc>& v_1,vector<_Tp, _Alloc>& v_2 )
        //This function will compare two vector if they are the same and return value i, i=-1 for they do not have same size, i>0 is that they do not have same element, the first difference appear at position i, i=0 for all the same
    void combine_matrix(TMatrixTBase<Element>& matrix_1,  TMatrixTBase<Element>& matrix_2 )
        //This function will combine two matrix diagonally and saved to matrix_1
    void selectrow(TMatrixD& aftermtx,TMatrixD& beforemtx,std::vector<Int_t>& v_selectrow)
        //This function will select the row of beforemtx according to the row_no provided in the vector v_selectrow and saved to aftermtx
    void selectcol(TMatrixD& aftermtx,TMatrixD& beforemtx,std::vector<Int_t>& v_selectcol)
        //This function will select the col of beforemtx
    void vectorselect(vector<_Tp, _Alloc>&v_after,vector<_Tp, _Alloc>&v_before,vector<int>& v_select)
        //This function will select the element in v_before and saved to the v_after according to the position provided in the vector v_select

