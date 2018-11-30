
//Name      : F_Factor.cpp
//Created by: G.Nishiyama
//Contents  : This program write 1p.e. charge distrebution and evaluate F_Facor.


////////////////////////////////////////////////////////////////////////////////
// ABBREVIATION RULE                                                          //
// after pulse -> ap							      //
// canvas      -> cnvs							      //
// histogram   -> hst							      //
// loop        -> lp                                                          //
// measurement -> msr							      //
// number      -> num                                                         //
// threshold   -> trsld							      //
// valuable    -> val							      //
// waveform    -> wform                                                       //
////////////////////////////////////////////////////////////////////////////////

#include "TApplication.h"
#include "TString.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include <cmath>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// DEFINITION OF VALUABLES                                                    //
// k_max_seq        : Time range what you want to analyse [us]                //
// k_num_ch         : Number of the DRS4 ch                                   //
// k_lp_start_slice : Number of start fill                                    //
// k_lp_extry_max   : Number of analysis data                                 //
// k_cell_cut_ch1   : Cell number for noise reduction in ch 1                 //
// k_cell_cut_ch2   : Cell number for noise reduction in ch 2                 //
// k_cell_cut_ch3   : Cell number for noise reduction in ch 3                 //
// k_cell_cut_ch4   : Cell number for noise reduction in ch 4                 //
// k_cut_swich      : Switch for noise reduction [true/cut,false/without]     //
// k_trsld_cut      : Threshold of the noise reduction range [mV]             //
// k_trsld_upper    : Upper limit of the noise reduction range [mV]           //
// msr_vals         : Measurement data group                                  //
// tree_ap          : Tree which contain measurement data                     //
// hst_wform        : Histogram for waveform                                  //
// time_msr         : date of measurement [MM/dd/hh/mm]                       //
// dhm_msr          : month.day.hour.minutes                                  // 
////////////////////////////////////////////////////////////////////////////////


//---- config gloabal valiable ----
Int_t nbin = 100;
Int_t XMinRange = -50;
Int_t XMaxRange = 200;
Float_t data_time[1024];
Float_t data_wform0[1024];   // negative
Float_t data_wform1[1024];   // positive
Float_t data_wform2[1024];   // negative
Float_t data_wform3[1024];   // positive
Int_t   data_stopcell = 0;   //stopcell
TString TreeName = "TreeOnePhoto_0";


//---- define function ----
Float_t GetCharge0();
Float_t GetCharge1();
void ShowUsage();


//---- Main Program ----
TH1F* OnePhoto(TFile* InFile,int ch,float *ope,float *operr)
{
  //  if (argc != 3) {
  // ShowUsage();
  // return -1;
  //}
  
  //---- read input file ----
  //  TFile *InFile = TFile::Open(filename);
  TTree *InTree = (TTree*)InFile->Get(TreeName);
  //  InTree->Print();
  cout << "" << endl;
  cout << "InFile has been read." << endl;

  if (ch == 0) {
    cout << "0 selected (ch1,ch2)" << endl;
    InTree->SetBranchAddress("time",data_time);
    InTree->SetBranchAddress("wform0",data_wform0);
    InTree->SetBranchAddress("wform1",data_wform1);
  }
  
  else if (ch == 1) {
    cout << "1 selected (ch3,ch4)" << endl;
    InTree->SetBranchAddress("time",data_time);
    InTree->SetBranchAddress("wform2",data_wform2);
    InTree->SetBranchAddress("wform3",data_wform3);
  }
  
  InTree->SetBranchAddress("stopcell",&data_stopcell);
  //int a = 0;
  //char** b;
  //TApplication theApp("App",&a,b);
  
  //---- init hist ----
  TH1F *hst = new TH1F("hst","1p.e. Charge distribution",nbin,XMinRange,XMaxRange);
  TCanvas *cnvs = new TCanvas;
  
  int EventNum = InTree->GetEntries();
  
  for (int nEvent=0; nEvent<EventNum; nEvent++) {
    InTree->GetEntry(nEvent);

    Float_t charge;
    
    if (ch == 0) {
      charge = GetCharge0();
    }
    else if (ch == 1) {
      charge = GetCharge1();
    }
    
    hst->Fill(charge);
    
    if (nEvent%5000 == 0) {
      cout << nEvent << "th Event has been finished" << endl;
    }
  }

  //hst->Draw();
  hst->SetXTitle("Charge[mV ns]");
  hst->SetYTitle("[events]");
  cnvs->SetLogy();

  
  // Float_t oneMax = 0;
  // Int_t oneN;
    
  // if (ch == 0) {
  //   for (int sliceNum=20;sliceNum<1024;sliceNum++) {
  //     if (oneMax < data_wform1[sliceNum]-data_wform0[sliceNum]) {
  // 	oneMax = data_wform1[sliceNum]-data_wform0[sliceNum];
  // 	oneN = sliceNum;
  //     }
  //   }
  // }

  // if (ch == 1) {
  //   for (int sliceNum=20;sliceNum<1024;sliceNum++) {
  //     if (oneMax < data_wform3[sliceNum]-data_wform2[sliceNum]) {
  // 	oneMax = data_wform3[sliceNum]-data_wform2[sliceNum];
  // 	oneN = sliceNum;
  //     }
  //   }
  // }
  
  
  TF1 *func = new TF1("func","gaus(0)+gaus(3)",XMinRange,XMaxRange);
  //func->SetParameter(0,10000);
  //  func->SetParLimits(0,9000,20000);
  //func->SetParameter(1,0);
  // func->SetParLimits(1,-10,10);
  // func->SetParameter(2,20);
  //func->SetParameter(3,oneMax);
  // func->SetParLimits(3,50,500);
  //func->SetParameter(4,oneN);
  // func->SetParLimits(4,20,70);
  //func->SetParameter(5,40);
  // func->SetParLimits(5,10,100);
  // hst->Fit("func","","",XMinRange,XMaxRange);

  float pdstl = hst->GetBinCenter(hst->GetMaximumBin());
  std::cout<<"Pedestal is "<<pdstl<<std::endl;
  hst->Fit("gaus","","",pdstl - 5.,pdstl + 5.);
  TF1 *func_pdstl = hst->GetFunction("gaus");
  double pdstl_par[3]={0};//constant ,mean ,sigma
  func_pdstl->GetParameters(pdstl_par);
  TGraph *grph_onepht = new TGraph();
  for(int lp_bin = 0 ; lp_bin < nbin ; lp_bin++)
    {
      float x = hst->GetBinCenter(lp_bin+1);
      float y = hst->GetBinContent(lp_bin+1);
      std::cout<<"x is "<<x<<"   : y is "<<y<<"    : f "<<func_pdstl->Eval(x)<<std::endl;
      grph_onepht->SetPoint(lp_bin, x, y - func_pdstl->Eval(x));
    }
  grph_onepht->Fit("gaus","","",3.*pdstl_par[2],XMaxRange);
  grph_onepht->Draw();
  //theApp.Run();
  double onepht_par[3]={0};//constant ,mean ,sigma
  grph_onepht->GetFunction("gaus")->GetParameters(onepht_par);
  //(TF1*)(grph_onepht->Get("gaus"))->GetParameters(onepht_par);
  TF1 *func_fit = new TF1("func_fit","gaus(0)+gaus(3)",XMinRange,XMaxRange);
  func_fit->SetParameter(0,pdstl_par[0]);
  func_fit->SetParameter(1,pdstl_par[1]); 
  func_fit->SetParameter(2,pdstl_par[2]); 
  func_fit->SetParameter(3,onepht_par[0]);
  func_fit->SetParameter(4,onepht_par[1]);
  func_fit->SetParameter(5,onepht_par[2]);
  cout<<"const : "<<onepht_par[0]<<endl;
  cout<<"mean  : "<<onepht_par[1]<<endl;
  cout<<"sigma : "<<onepht_par[2]<<endl;
  hst->Fit("func_fit","","",XMinRange,XMaxRange);
  
  gStyle->SetOptFit(1111);
  
  Double_t FitPar[6];
  Double_t F_Factor;
  func_fit->GetParameters(FitPar);
  cout << "mean = " << FitPar[4] << endl;
  cout << "sigma = " << FitPar[5] << endl;
  F_Factor = sqrt(1+pow(FitPar[5]/FitPar[4],2.0));
  cout << "F-Factor = " << F_Factor << endl;

  *ope = func_fit->GetParameter(4);
  *operr = func_fit->GetParError(4);
  return hst;
}


Float_t GetCharge0()
{
  Float_t max = 0;
  Int_t maxN = 0;
  Float_t a = 0;
  Float_t b = 0;
  Float_t gt = 0;
  Float_t v = 0;
  Int_t minN = 0;
  Float_t min = 2000;
  
  for (int sliceNum=1; sliceNum<1023; sliceNum++)
    {
      if (max < (data_wform1[sliceNum-1]-data_wform0[sliceNum-1])+(data_wform1[sliceNum]-data_wform0[sliceNum])+(data_wform1[sliceNum+1]-data_wform0[sliceNum+1]))
	{
	  max = (data_wform1[sliceNum-1]-data_wform0[sliceNum-1])+(data_wform1[sliceNum]-data_wform0[sliceNum])+(data_wform1[sliceNum+1]-data_wform0[sliceNum+1]);
	  maxN = sliceNum;
	}
    }
  
  if (100<maxN && maxN<1000)
    {
      for (int sliceNum=maxN-10; sliceNum < maxN+11; sliceNum++)
	{
	  a += data_wform1[sliceNum]-data_wform0[sliceNum];
	  b += data_time[sliceNum]*(data_wform1[sliceNum]-data_wform0[sliceNum]);
	}
    }
  else
    {
      a = 1;
      b = 10000;
    }
  gt = b/a;
  
  if (50 < gt && gt <500)
    {      
      for (int sliceNum=0; sliceNum<1024; sliceNum++)
	{
	  if (abs(data_time[sliceNum]-gt) < min)
	    {
	      min = abs(data_time[sliceNum]-gt);
	      minN = sliceNum;
	      
	    }
	}
      
      for (int sliceNum=minN-3; sliceNum<minN+4; sliceNum++)
	{
	  if( (sliceNum + data_stopcell)%1024 == 175 )
	    {
	      v = 150;
	      break;
	    }
	  v += data_wform1[sliceNum]-data_wform0[sliceNum];
	}
    }
  else
    {
      for (int sliceNum=250; sliceNum<257; sliceNum++)
	{
	  if( (sliceNum + data_stopcell)%1024 == 175 )
	    {
	      v = 150;
	      break;
	    }
	  v += data_wform1[sliceNum]-data_wform0[sliceNum];
	}
    }
  return v;
}

Float_t GetCharge1()
{
  Float_t max = 0;
  Int_t maxN = 0;
  Float_t a = 0;
  Float_t b = 0;
  Float_t gt = 0;
  Float_t v = 0;
  Int_t minN = 0;
  Float_t min = 2000;
  
  for (int sliceNum=1; sliceNum<1023; sliceNum++)
    {
      if (max < (data_wform3[sliceNum-1]-data_wform2[sliceNum-1])+(data_wform3[sliceNum]-data_wform2[sliceNum])+(data_wform3[sliceNum+1]-data_wform2[sliceNum+1]))
	{
	  max = (data_wform3[sliceNum-1]-data_wform2[sliceNum-1])+(data_wform3[sliceNum]-data_wform2[sliceNum])+(data_wform3[sliceNum+1]-data_wform2[sliceNum+1]);
	  maxN = sliceNum;
	}
    }
  
  if (100<maxN && maxN<1000)
    {
      for (int sliceNum=maxN-10; sliceNum<maxN+11; sliceNum++)
	{
	  a += data_wform3[sliceNum]-data_wform2[sliceNum];
	  b += data_time[sliceNum]*(data_wform3[sliceNum]-data_wform2[sliceNum]);
	}
    }
  else
    {
      a = 1;
      b = 10000;
    }
  gt = b/a;
  
  if (50 < gt && gt <500)
    {
      for (int sliceNum=0; sliceNum<1024; sliceNum++)
	{
	  if (abs(data_time[sliceNum]-gt) < min)
	    {
	      min = abs(data_time[sliceNum]-gt);
	      minN = sliceNum;
	    }
	}
      for (int sliceNum=minN-3; sliceNum<minN+4; sliceNum++)
	{
	  v += data_wform3[sliceNum]-data_wform2[sliceNum];
	}
    }
  else
    {
      for (int sliceNum=250; sliceNum<257; sliceNum++)
	{
	  v += data_wform3[sliceNum]-data_wform2[sliceNum];
	}
    }
  return v;
}

void ShowUsage()
{
  cout << "Usage : F_Factor FileName 0(1,2ch)or1(3,4ch)" << endl;
}
