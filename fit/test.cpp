#include "classdefine.h"
//#include "function_for_fit.cpp"

void test(){
    EWFit zpole;
    zpole.folderpwd="/Users/yuhaowang/Documents/electroweak-data-collection-fit/Z_pole/";
    zpole.readall();
    zpole.complete_correlation_withI();
    zpole.print_size();

    EWFit wpole;
    wpole.folderpwd="/Users/yuhaowang/Documents/electroweak-data-collection-fit/W_pole/";
    wpole.readall();
    wpole.complete_correlation_withI();
    wpole.print_size();

    EWFit z_w_combine=zpole.CombineEWFit(wpole);
    z_w_combine.print_size();
    EWFit z_w_ruleout0=z_w_combine.RuleOut0Contribution();
    z_w_ruleout0.print_size();
    z_w_ruleout0.LSFit();
    z_w_ruleout0.print_result();

    string op_1[]={"cpWB","cHQP"};
    string ob_1[]={"sigma_had","AFB_b"};
    string ob_2[]={"Gamma_W","Mw"};
    string ob_3[]={"sigma_had","AFB_b","Gamma_W","Mw"};

    EWFit SubFit1=z_w_ruleout0.SetSub(op_1,ob_1);
    EWFit SubFit2=z_w_ruleout0.SetSub(op_1,ob_2);
    EWFit SubFit3=z_w_ruleout0.SetSub(op_1,ob_3);

    SubFit1.LSFit();
    TEllipse* elp1=SubFit1.GetEllipse();
    SubFit2.LSFit();
    TEllipse* elp2=SubFit2.GetEllipse();
    SubFit3.LSFit();
    TEllipse* elp3=SubFit3.GetEllipse();

    
    elp1->SetFillColorAlpha(2,.3);
    elp2->SetFillColorAlpha(3,.3);
    elp3->SetFillColorAlpha(4,.7);

    double xmin=-5;
    double xmax=5;
    double ymin=-5;
    double ymax=5;

    TCanvas *mycanvas=new TCanvas();
    TGaxis *axisx=new TGaxis(xmin,0,xmax,0,xmin,xmax,510,"");
    TGaxis *axisy=new TGaxis(0,ymin,0,ymax,ymin,ymax,510,"");
    axisx->SetTitle(op_1[0].c_str());
    axisy->SetTitle(op_1[1].c_str());
    mycanvas->Range(xmin,ymin,xmax,ymax);
    elp1->Draw();
    elp2->Draw();
    elp3->Draw();
    axisx->Draw();
    axisy->Draw();



}
