/////////////////////////////////////////////
// File include                            //
/////////////////////////////////////////////
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TString.h"
#include "TMath.h"
#include "string"
#include "iostream"
#include <sys/stat.h>
#include <sys/types.h>
#include <complex>
#include <unistd.h>

/////////////////////////////////////////////
//User's include file                      //
/////////////////////////////////////////////
#include <global.h>
#include <pre_analysis.h> 
#include <chrg.h>
#include <name_edit.h>

/////////////////////////////////////////////
// Main Program                            //
/////////////////////////////////////////////
int main( int argc, char **argv)
{
  TGraphErrors* grph_volt = new TGraphErrors();
  float v_gain[argc];//={};
  float g_err[argc];//={};
  float v_widl[argc];//={};
  float v_widu[argc];//={};
  float gw_err[argc];//={};
  for(int lp_file = 1;lp_file < argc;lp_file++)
    {
      std::string filename = argv[lp_file];      
      TFile *openfile = TFile::Open(argv[lp_file]);
      ///////////In <name_edit.h>,Get PMT Serial and Measurement Date///////
      TString Serial = GetSerial(filename);
      TGraphErrors* grph_gain  = (TGraphErrors*)(openfile->Get("grph_Gain"));
      TGraphErrors* grph_width = (TGraphErrors*)(openfile->Get("grph_Width"));
      TF1* func_gain = grph_gain->GetFunction("func_fit");
      std::cout<<__LINE__<<std::endl;
      //func_gain->SetName("func_gain");
        std::cout<<__LINE__<<std::endl;
      TF1* func_width= grph_width->GetFunction("func_fit");
      std::cout<<__LINE__<<std::endl;
      //////////////////////////////////////////////////////////////////////
      float a      = func_gain->GetParameter(0)            ;
      float b      = (func_gain->GetParameter(1))*0.75     ;
      float a_err  = func_gain->GetParError(0)             ;
      float b_err  = (func_gain->GetParError(1))*0.75             ;
      std::cout<<__LINE__<<std::endl;
      //v_gain[lp_file] = func_gain->GetX(40000.)            ;std::cout<<__LINE__<<std::endl;      
      v_gain[lp_file] = 350 + std::pow( (40000./a) , (1./b) )                          ;
      a_err           = (1./b)*std::pow( (40000./a) , (1./b - 1.) )*a_err              ;
      b_err           = std::pow( (40000./a) , (1./b) )*log(40000./a)*(1./(b*b))*b_err ;
      g_err[lp_file ] = sqrt(a_err*a_err + b_err*b_err)                                ;
      float chi_g = func_gain->GetChisquare();
      float chi_w = func_width->GetChisquare();
      float wa = func_width->GetParameter(0)            ;
      float wb = func_width->GetParameter(1)            ;
      float wc = func_width->GetParameter(2)            ;
      float wa_err = func_width->GetParError(0)         ;
      float wb_err = func_width->GetParError(1)         ;
      float wc_err = func_width->GetParError(2)         ;
      //v_widl[lp_file] = func_width->GetX(3.0);
      v_widl[lp_file] =wa/(3.0-wc)+wb                  ; 
      wa_err          =wa_err/(3.0-wc)         ; 
      wc_err          =wa*wc_err/((3.0-wc)*(3.0-wc))  ; 
      gw_err[lp_file] =sqrt(wa_err*wa_err + wb_err*wb_err + wc_err*wc_err);
      std::cout<<Serial<<std::endl;
      std::cout<<"40,000 Gain Voltage is "<< v_gain[lp_file]<<std::endl
	       <<"  Error is             "<<g_err[lp_file ]<<std::endl
	       <<"  Chi square is        "<<chi_g<<std::endl
	       <<"                       "<<func_gain->GetNDF()<<std::endl;   
      std::cout<<"    3.0 ns  Voltage is "<< v_widl[lp_file]<<std::endl
	       <<"  Error is             "<< gw_err[lp_file]<<std::endl
	       <<"  Chi square is        "<<chi_w<<std::endl   
      	       <<"                       "<<func_width->GetNDF()<<std::endl;   
      grph_volt->SetPoint(lp_file,v_gain[lp_file],v_widl[lp_file]);      
      openfile->Close();
    }
  TFile *savefile = TFile::Open("Sample.root","RECREATE");
  grph_volt->Write();
  savefile->Close();
  return 0;
}

