#include"iostream" 
#include"fstream"
#include"string"
#include"TFile.h"
#include"TGraphErrors.h"
#include"TF1.h"
#include"TCanvas.h"
#include"geometric.h"

int main(int argc,char **argv)
{  
  std::string outfile = argv[1];
  outfile = outfile + ".csv";
  std::ofstream ofs(outfile);
  std::string str_header ="Serial,Nominal_HV,TT_Nominal,TT_1150V"; 
  ofs << str_header<<std::endl;  

  for(int lp_file = 0;lp_file < argc -1;lp_file++ )
    {
      TFile* infile  = TFile::Open(argv[lp_file+1]);
      TString str_fn = argv[lp_file+1];
      int fn_init = str_fn.First("AA");
      TString Serial = str_fn(fn_init,6);
      TGraphErrors* gr_tt  = (TGraphErrors*)infile->Get("Transittime_Graph");
      TGraphErrors* gr_gain= (TGraphErrors*)infile->Get("grph_Gain");
      TF1*   func_tt   = gr_tt->GetFunction("func_fit");
      TF1*   func_gain = gr_gain->GetFunction("func_fit");
      float nom_HV = func_gain->GetX(40000,800,1500);
      float nom_TT = func_tt->Eval(nom_HV);
      float fix_TT = func_tt->Eval(1150);
      //double* x_tt = gr_tt->GetY();
      double* y = gr_gain->GetY();      
      double* x = gr_tt->GetX();
      for(int i=0;i<5;i++)
	{
	  if(y[i]>40000)
	    {
	      nom_HV = xlinear(x[i],y[i],x[i-1],y[i-1],40000);
	    }
	}
      nom_TT = func_tt->Eval(nom_HV);
      std::cout<<"Nominal HV is           : "<<nom_HV<<std::endl;
      std::cout<<"Nominal Transittime  is : "<<nom_TT<<std::endl;
      std::cout<<"Transittime @ 1150V  is : "<<fix_TT<<std::endl;
      ofs<<Serial.Data()<<","<<nom_HV<<","<<nom_TT<<","<<fix_TT<<std::endl;
    }
  return 0;
}
