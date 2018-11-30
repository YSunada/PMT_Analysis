

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
#include "TGraph.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include <cmath>

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


/////////////////////////////////////////////
// Define constant                         //
/////////////////////////////////////////////

const int k_rng_cg            =10;    //Search range for center of gracity
const int k_num_slice         =1024;  //Total number of drs4 channnel
const int k_edg_cut           =10;    //Deviation of the center of gravity
const float k_chrg_per_up     =0.9;   //
const float k_chrg_per_dwn    =0.1;   //  
const float k_chrg_per =k_chrg_per_up - k_chrg_per_dwn; //
const float k_echrg = 1.602176e-19;   //Elementary Charge
const float k_impd  = 1.2e3;          //Impedance of PACTA
const float pls_wid = 5.0;             //FWFM pulse width [ns] 
//float smpl_frq = 5.0; //Sampling Speed of DRS4 [GHz]
//gSystem->Load("../include/chrg_rng.h")
//gSystem->Load("../include/chrg.h")
#include <pre_analysis.h>
#include <chrg_rng.h>
#include <chrg.h>
/////////////////////////////////////////////
// Main Program                            //
/////////////////////////////////////////////

int main (int argc,char **argv)
{
  TString filename = argv[1];
  int ch = atoi(argv[2]);
  TFile* openfile = TFile::Open(filename);
  TApplication theApp("theApp",&argc,argv);
  //  void OnePhoto(TFile* openfile,int ch){
 
  TTree *works    = (TTree*)openfile->Get("TreeOnePhoto_0");
  TTree *dark     = (TTree*)openfile->Get("TreeDark_0");
  int k_event = works->GetEntries();
  float wform_cg[k_event] ;
  float offset[k_event]   ;
  float dwform_cg[k_event] ;
  float doffset[k_event]   ;
  float charge[k_event] ;
  float dcharge[k_event] ;
  float chrg_max =0;
  float chrg_min =0;
  float dchrg_max =0;
  float dchrg_min =0;
  float gc_av =0;
  float dgc_av =0;
  int   event_max=0;
  int   event_min=0;
  int   data_stopcell=0;
  int   dark_stopcell=0;
  float peak_time=0;
  int peak_bin=0;
  float data_nwform[k_num_slice] = {0};
  float data_pwform[k_num_slice] = {0};
  float data_time[k_num_slice]   = {0};
  float dark_nwform[k_num_slice] = {0};
  float dark_pwform[k_num_slice] = {0};
  float dark_time[k_num_slice]   = {0};
  
  works->SetBranchAddress("time",data_time);
  dark ->SetBranchAddress("time",dark_time);
  works->SetBranchAddress("stopcell",&data_stopcell);
  dark ->SetBranchAddress("stopcell",&dark_stopcell);
  
  
  if(ch==0)
    {
      works->SetBranchAddress("wform0",data_nwform);
      works->SetBranchAddress("wform1",data_pwform);
      dark ->SetBranchAddress("wform0",dark_nwform);
      dark ->SetBranchAddress("wform1",dark_pwform);
    }
  else if(ch==1)
    {
      works->SetBranchAddress("wform2",data_nwform);
      works->SetBranchAddress("wform3",data_pwform);
      dark ->SetBranchAddress("wform2",dark_nwform);
      dark ->SetBranchAddress("wform3",dark_pwform);
    }

  TH1F* hst = new TH1F("","",k_num_slice,0,(int)k_num_slice/pls_wid); 
  for(int lp_event=0;lp_event<k_event;lp_event++)
    {
      works->GetEntry(lp_event);
      for(int lp_slice=k_edg_cut;lp_slice<k_num_slice-k_edg_cut;lp_slice++)
	{
	  if(data_pwform[lp_slice] - data_nwform[lp_slice] > 0)
	  hst->Fill(data_time[lp_slice],data_pwform[lp_slice] - data_nwform[lp_slice]);
	}
    }
  
  peak_bin = hst->GetMaximumBin();
  peak_time= hst->GetBinCenter(peak_bin);
 
  for(int lp_event=0;lp_event<k_event;lp_event++)
    {
      int cnt=0;
      works->GetEntry(lp_event);
      dark ->GetEntry(lp_event);
      offset[lp_event] = 0;
      doffset[lp_event] = 0;
      for(int lp_slice=k_edg_cut;data_time[lp_slice]<peak_time-pls_wid;lp_slice++)
	{
	  cnt++;
	  offset[lp_event]  +=  data_pwform[lp_slice] - data_nwform[lp_slice];
	  doffset[lp_event] +=  dark_pwform[lp_slice] - dark_nwform[lp_slice];
	}
      offset[lp_event] = offset[lp_event]/(float)cnt;
      doffset[lp_event] = doffset[lp_event]/(float)cnt;      
      cnt =0;

    }
  //Search Gravity Center
  int cnt =0;
  int dcnt=0;
  for(int lp_event=0;lp_event<k_event;lp_event++)
    {
      works->GetEntry(lp_event);
      dark ->GetEntry(lp_event);
      float timevolt = 0;
      float voltsum = 0;
      float dtimevolt = 0;
      float dvoltsum = 0;
      for(int lp_slice = peak_bin - smpl_frq*pls_wid;lp_slice<peak_bin + smpl_frq*pls_wid;
	  lp_slice++)
	{
	  
	  timevolt += (data_pwform[lp_slice] - data_nwform[lp_slice] - offset[lp_event])*data_time[lp_slice];
	  voltsum  += data_pwform[lp_slice] - data_nwform[lp_slice]- offset[lp_event];
	  	  dtimevolt += (dark_pwform[lp_slice] - dark_nwform[lp_slice] - doffset[lp_event])*dark_time[lp_slice];
	  dvoltsum  += dark_pwform[lp_slice] - dark_nwform[lp_slice]- doffset[lp_event];
	}      
      wform_cg[lp_event] =timevolt/voltsum;
      dwform_cg[lp_event] =dtimevolt/dvoltsum;
      
      if(wform_cg[lp_event] > peak_time+k_rng_cg || wform_cg[lp_event] < peak_time-k_rng_cg)
	{
	  wform_cg[lp_event]=k_num_slice*2;
	}
      else
	{
	  gc_av += wform_cg[lp_event];
	  cnt ++ ;
	}

      if(dwform_cg[lp_event] > peak_time+k_rng_cg || dwform_cg[lp_event] < peak_time-k_rng_cg)
	{
	  dwform_cg[lp_event]=k_num_slice*2;
	}
      else
	{
	  dgc_av += dwform_cg[lp_event];
	  dcnt ++ ;
	}
    }
  gc_av = gc_av/(float)(cnt);
  dgc_av = dgc_av/(float)(dcnt);
  std::cout<<"gcav : "<<gc_av<<std::endl
	   <<"dgcav : "<<dgc_av<<std::endl;
    
  for(int lp_event = 0; lp_event < k_event; lp_event ++)
    {
      works->GetEntry(lp_event);
      dark ->GetEntry(lp_event);
      float wform[k_num_slice] = {0};
      float dwform[k_num_slice] = {0};
      if(wform_cg[lp_event] == k_num_slice*2)
	{
	  wform_cg[lp_event] = gc_av;
	}
     
      if(dwform_cg[lp_event] == k_num_slice*2)
	{
	  dwform_cg[lp_event] = dgc_av;
	}
      
      for(int lp_slice = 0; lp_slice < k_num_slice; lp_slice ++)
	{
	  wform[lp_slice] = data_pwform[lp_slice] - data_nwform[lp_slice] - offset[lp_event];
	  dwform[lp_slice] = dark_pwform[lp_slice] - dark_nwform[lp_slice] - doffset[lp_event];
	}
      charge[lp_event] = chrg(ch,wform,data_time, wform_cg[lp_event] -3. , wform_cg[lp_event] +3.,data_stopcell);
      dcharge[lp_event] = chrg(ch,dwform,dark_time, dwform_cg[lp_event] -3. , dwform_cg[lp_event] +3.,dark_stopcell);
	    //dcharge[lp_event] = chrg(ch,dwform,dark_time, peak_time-3., peak_time+3.,dark_stopcell);


      
      if(wform_cg[lp_event] != 2048)
	{
	  if(lp_event == 0)
	    {
	      chrg_min = charge[lp_event];
	      chrg_max = charge[lp_event];
	    }
	  
	  if(charge[lp_event] > chrg_max)
	    {
	      chrg_max = charge[lp_event];
	      event_max= lp_event;
	    }
	  
	  if(charge[lp_event] < chrg_min)
	    {
	      chrg_min = charge[lp_event];
	      event_min= lp_event;
	    }
	}
      if(dwform_cg[lp_event] != 2048)
	{
	  if(lp_event == 0)
	    {
	      dchrg_min = dcharge[lp_event];
	      dchrg_max = dcharge[lp_event];
	    }
	  
	  if(charge[lp_event] > dchrg_max)
	    {
	      dchrg_max = dcharge[lp_event];
	    }
	  
	  if(charge[lp_event] < dchrg_min)
	    {
	      dchrg_min = dcharge[lp_event];
	    }
	}      
    }
  //int bin_num = int(TMath::Sqrt(float(k_event)));
  //chrg_min =0.;
  float bin_num = (chrg_max-chrg_min)/0.1;
  TH1F *ChargeDist = new TH1F("","",(int)bin_num,chrg_min,chrg_max);
  TH1F *ChargeDist2 = new TH1F("","",(int)(bin_num),chrg_min,chrg_max);
  TH1F *DarkDist = new TH1F("","",(int)bin_num,chrg_min,chrg_max); 
  TGraph *grph1   = new TGraph() ;
  TGraph *grph2   = new TGraph() ;
  std::cout<<"min is : "<<event_min<<" : "<<chrg_min<<"  GC :"<<wform_cg[event_min]<<std::endl
	   <<"max is : "<<event_max<<" : "<<chrg_max<<"  GC :"<<wform_cg[event_max]<<std::endl;
  
  for(int lp_event=0;lp_event<k_event;lp_event++)
    {
      ChargeDist->Fill(charge[lp_event]) ;
      DarkDist  ->Fill(dcharge[lp_event]);      
    }
  int max_bin = ChargeDist->GetMaximumBin();
  int dmax_bin = DarkDist->GetMaximumBin();
  
  float single_max=ChargeDist->GetBinContent(max_bin-1)
    +ChargeDist->GetBinContent(max_bin)
    +ChargeDist->GetBinContent(max_bin+1);

  float dark_max =DarkDist->GetBinContent(dmax_bin-1)
    +DarkDist->GetBinContent(dmax_bin)
    +DarkDist->GetBinContent(dmax_bin+1);
  
  DarkDist->Scale(single_max/dark_max);

  
  TCanvas* cnvs1 =new TCanvas("ALL","ALL",600,600);
  cnvs1->cd();
  ChargeDist->Draw();
  DarkDist->SetLineColor(3);
  DarkDist->Draw("SAME");
  TCanvas* cnvs2 =new TCanvas("Dark","Dark",600,600);
  cnvs2->cd();
  //  DarkDist->Fit("gaus","","",0.,10.);
  DarkDist->Draw();

  ChargeDist2->Add(ChargeDist,DarkDist,1,-1.);
  ChargeDist2->Rebin(30);
  
  for(int lp_bin = 0; lp_bin<bin_num/10;lp_bin++)
    {
      if(ChargeDist2->GetBinContent(lp_bin)<0 || ChargeDist2->GetBinCenter(lp_bin)<0 )
	{
	  ChargeDist2->SetBinContent(lp_bin,0.);
	}
    }
  
  int bin_num2 = ChargeDist2->GetNbinsX();
  float two_phe[bin_num2];
  float sum = DarkDist->Integral();
  TH1F *ChargeDist3 = new TH1F("","",sum,chrg_min,chrg_max);
  ChargeDist3->SetLineColor(3);
  for(int lp_x=0;lp_x<bin_num2;lp_x++)
    {
      two_phe[lp_x+1] =0;
      for(int lp_y=0;lp_y<lp_x;lp_y++)
	{
	  two_phe[lp_x+1] += (ChargeDist2->GetBinContent(lp_x-lp_y))
	    *(ChargeDist2->GetBinContent(lp_y+1))/(2*sum); 
	}
      float x = ChargeDist2->GetBinCenter(lp_x+1);
      ChargeDist3->Fill(x,two_phe[lp_x+1]);
      std::cout<<two_phe[lp_x+1]<<std::endl;
    }
 
  TCanvas* cnvs3 =new TCanvas("SinglePhoto1","SinglePhoto1",600,600);
  cnvs3->cd();
  ChargeDist2->Draw();
  ChargeDist3->Draw("SAME");
  theApp.Run();
  return 0;
}

