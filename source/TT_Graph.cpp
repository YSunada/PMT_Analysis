#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TGraphErrors.h"
#include <sys/stat.h>
#include <unistd.h>
#include "iostream"
#include "stdlib.h"
#include "name_edit.h"
#include "geometric.h"
#include "TF1.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TSystem.h"
//#include "pre_analysis.h"

float smpl_frq = 5.0; //Sampling Speed of DRS4 [GHz]
int pls_wdth = 5;//About Pulse Width [ns] FWFM
int nSlice = 1024;//Number of DRS4 Channel
int k_num_slice =nSlice;
int k_edg_cut =20;

#include "getopt.h"
#include "chrg.h"
int main (int argc,char **argv)
{
  TCanvas* cnv = new TCanvas("ccc","ccc",600,600);
  std::string str_fname ="NULL";
  //Option Control//////////////////////////////////////////////////////////////
  int opt;
  opterr = 0;
  bool fflag= false;
  bool hflag= false;
  while ((opt = getopt(argc,argv,"f:h")) != -1)
    {
      switch (opt)
	{
	case 'f':
	  str_fname = optarg;
	  fflag =true;
	  break;
	case 'h':
	  hflag =true;
	  break;
	}
    }
  
  std::cout<<"line : "<<__LINE__<<std::endl;
  ////////////////////////////////////////////////////////////////////////////////
  
  //Option Control//////////////////////////////////////////////////////////////
  if(hflag || argc <2)
    {
      std::cout<<"/////////////////////////////////////////////////////////////////"<<std::endl
	       <<"  Option List                "<<std::endl
	       <<"  -f  : input filename (needed)  "<<std::endl
	       <<"/////////////////////////////////////////////////////////////////"<<std::endl;
      return -1;
    }

 
  ////////////////////////////////////////////////////////////////////////////////
  TString fname = str_fname;
  TFile* infile = TFile::Open(fname);

  //std::string str_fname =infilename.Data();
  int  fn = str_fname.rfind("TT");
    std::cout<<__LINE__<<std::endl
	     <<fn<<std::endl
    	     <<str_fname<<std::endl;
  TString str_pmt  = str_fname.substr(fn+2,6);
  TString str_spm  = str_fname.substr(fn+11,6);
  TString str_date = str_fname.substr(fn+18,12);
  TFile* ofile=new TFile("PMT_Result/"+str_pmt+"_" +str_date +".root","UPDATE");
  TGraphErrors* grp_TT  = new TGraphErrors();
  TGraphErrors* grp_TTS = new TGraphErrors();
  TString str_vol[5] = {"1000","1100","1200","1300","1400"};
  std::cout<<__LINE__<<std::endl
    	   <<str_pmt<<std::endl
    	   <<str_spm<<std::endl
    	   <<str_date<<std::endl;
  TTree* works=new TTree();
  for(int i=0;i<5;i++)
    {
      TString treename =str_pmt+"_"+str_spm+"_TreeTT_Fil5_"+str_vol[i]+"V_0";
      std::cout<<treename<<std::endl;
      works  = (TTree*)infile->Get(treename);
      TH1F* hst = new TH1F("hst","",400,10,30);
      works->Draw("ptime_01-rtime_02>>hst","","",works->GetEntries(),0);
      float mean = hst->GetMean();
      float sigma = hst->GetRMS();
      hst->Fit("gaus","","",mean-sigma,mean+sigma);
      mean = hst->GetFunction("gaus")->GetParameter(1);
      sigma = hst->GetFunction("gaus")->GetParameter(2);
      hst->Fit("gaus","","",mean - 3.0*sigma,mean + 3.0*sigma);
      mean = hst->GetFunction("gaus")->GetParameter(1);
      sigma = hst->GetFunction("gaus")->GetParameter(2);
      float mean_err  = hst->GetFunction("gaus")->GetParError(1);
      float sigma_err = hst->GetFunction("gaus")->GetParError(2);
      grp_TT->SetPoint(i,str_vol[i].Atoi(),mean);
      grp_TT->SetPointError(i,0,mean_err);
      grp_TTS->SetPoint(i,str_vol[i].Atoi(),sigma);
      grp_TTS->SetPointError(i,0,sigma_err);
      hst->SetName("Transitime_"+str_vol[i]+"V");
      hst->SetTitle(str_pmt+"  Transitime @ "+str_vol[i]+"V");
      ofile->cd();
      hst->Write();
      hst->Delete();
    }
  TF1* func_TT = new TF1("func_fit","[0]*TMath::Power(x,-0.5)+[3]",1000,1400);
  grp_TT->Fit(func_TT);
  grp_TT->SetName("Transittime_Graph");
  grp_TT->SetTitle(str_pmt + "  Transittime : HV");
  grp_TT->Write();
  /*
  cnv->cd();
  grp_TT->Draw("AP");
  cnv->Modified();
  cnv->Update();
  gSystem->ProcessEvents();
  std::getchar();
  */
  //return 0;  
  // Save Data /////////////////////////////////////////////////////////////////////////////////
  ofile->Close();
  std::cout<<"line : "<<__LINE__<<std::endl;  
  //mkdir("PMT_Result",0777);
 
  
  return 0;
}

