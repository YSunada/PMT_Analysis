
#include "TString.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include <cmath>

#include <pre_analysis.h>
#include <chrg_rng.h>
#include <chrg.h>

/////////////////////////////////////////////
// Main Program                            //
/////////////////////////////////////////////


TH1F* OnePhoto(TTree* works, int ch)
{
  std::cout<<"One Photo Analysis Start !!!"<<std::endl;
  int k_event_1pe = works->GetEntries();
  std::cout<<__LINE__<<std::endl;
  
  TH1F* hst_onepht_avrg = new TH1F();
  float wform_cg[k_event_1pe] ;
  float offset[k_event_1pe]   ;

  hst_onepht_avrg = pre_analysis(works,ch,k_event_1pe,offset,wform_cg);
  hst_onepht_avrg->Draw();
  std::cout<<__LINE__<<std::endl;
  float charge[k_event_1pe] ;
  float chrg_max  = 0;
  float chrg_min  = 0;
  int data_stpcll = 0;
  int   event_max = 0;
  int   event_min = 0;
   
  float data_nwform[k_num_slice] = {0};
  float data_pwform[k_num_slice] = {0};
  float data_time[k_num_slice]   = {0};
  std::cout<<__LINE__<<std::endl;
  works->SetBranchAddress("time",data_time);
  works->SetBranchAddress("stopcell",&data_stpcll);      
  if(ch==0)
    {
      works->SetBranchAddress("wform0",data_nwform);
      works->SetBranchAddress("wform1",data_pwform);
    }
  else if(ch==1)
    {
      works->SetBranchAddress("wform2",data_nwform);
      works->SetBranchAddress("wform3",data_pwform);
    }
  std::cout<<__LINE__<<std::endl;
  for(int lp_event = 0; lp_event < k_event_1pe; lp_event ++)
    {
      works->GetEntry(lp_event);
      
      float wform[k_num_slice] = {0};

      for(int lp_slice = 0; lp_slice < k_num_slice; lp_slice ++)
	{
	  wform[lp_slice] = data_pwform[lp_slice] - data_nwform[lp_slice] - offset[lp_event];
	}
      
      charge[lp_event] = chrg(ch,wform,data_time,wform_cg[lp_event] -3. ,
			      wform_cg[lp_event] +3., data_stpcll);
      if(charge[lp_event] != -7777.)
	{
	  if(lp_event == 0)
	    {
	      chrg_min = charge[lp_event];
	      chrg_max = charge[lp_event];
	    }
	  
	  if(charge[lp_event] > chrg_max)
	    {	
	      chrg_max = charge[lp_event];
	    }
	  
	  if(charge[lp_event] < chrg_min)
	    {
	      chrg_min = charge[lp_event];
	    }
	}

      if(TMath::Abs(charge[lp_event]) > 1000)
	{
	  std::cout<<"charge : "<<charge[lp_event]<<std::endl
		   <<"offset : "<<offset[lp_event]<<std::endl
		   <<"gc     : "<<wform_cg[lp_event]<<std::endl;
	}

    }

  TH1F *ChargeDist = new TH1F("","",
			      int(TMath::Sqrt(float(k_event_1pe))),chrg_min - 10.,chrg_max + 10.);
  
  for(int lp_event=0;lp_event<k_event_1pe;lp_event++)
    {
      if(charge[lp_event] != -7777.)
	{
	  ChargeDist->Fill(charge[lp_event]);
	}
    }
  
  float pdstl = ChargeDist->GetBinCenter(ChargeDist->GetMaximumBin());
  std::cout<<"in OnePhoto.h : "<<__LINE__<<std::endl;
  double pdstl_par[3] =  {0};//constant ,mean ,sigma
  ChargeDist->Fit("gaus","","",pdstl - 5.,pdstl + 5.);
 std::cout<<"in OnePhoto.h : "<<__LINE__<<std::endl; 
 TF1 *func_pdstl     =  ChargeDist->GetFunction("gaus");
 std::cout<<"in OnePhoto.h : "<<__LINE__<<std::endl;
  func_pdstl         ->  GetParameters(pdstl_par); 


  TGraph *grph_onepht =  new TGraph(); 
  for(int lp_bin = 0 ; lp_bin < int(TMath::Sqrt(float(k_event_1pe))) ; lp_bin++)
    {
      float x = ChargeDist->GetBinCenter(lp_bin+1);
      float y = ChargeDist->GetBinContent(lp_bin+1);
       grph_onepht->SetPoint(lp_bin, x, y - func_pdstl->Eval(x));
    }
  
  double onepht_par[3]={0};//constant ,mean ,sigma
  grph_onepht->Fit("gaus","","",3.*pdstl_par[2],chrg_max);
  grph_onepht->GetFunction("gaus")->GetParameters(onepht_par);

  TF1 *func_fit = new TF1("func_fit","gaus(0)+gaus(3)",chrg_min,chrg_max);
  func_fit->SetParameter(0,pdstl_par[0]);
  func_fit->SetParameter(1,pdstl_par[1]); 
  func_fit->SetParameter(2,pdstl_par[2]); 
  func_fit->SetParameter(3,onepht_par[0]);
  func_fit->SetParameter(4,onepht_par[1]);
  func_fit->SetParameter(5,onepht_par[2]);
  ChargeDist->Fit("func_fit","","",chrg_min,chrg_max);

  return ChargeDist;
}
