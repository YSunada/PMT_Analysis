#ifndef INCLUDED_PREANALISYS_SUNA
#define INCLUDED_PREANALISYS_SUNA

#include "TTree.h"
#include "TH1F.h"
#include <global.h>
float smpl_frq = 5.0; //Sampling Speed of DRS4 [GHz]

TH1F* pre_analysis(TTree* works,int ch ,float* offset ,float* wform_cg)
{
  //  const int g_chrg_rng = 10;//Search Range for Gravity Center
  TH1F *hst1 = new TH1F("","",g_num_slice,0,g_num_slice/smpl_frq);
  float time[g_num_slice]={0};
  float ngtv[g_num_slice]={0};
  float pstv[g_num_slice]={0};
  float gc_av = 0;
  float peak = 0;
  int  peak_bin = 0; 
  int  cnt  = 0;

  //Data Set
  int event_num =  works->GetEntries();
  works->SetBranchAddress("time",time);  
  if(ch==0)
    {
      works->SetBranchAddress("wform0",ngtv);
      works->SetBranchAddress("wform1",pstv);
    }
  else if(ch==1)
    {
      works->SetBranchAddress("wform2",ngtv);
      works->SetBranchAddress("wform3",pstv);
    }

  //Peak Search
  for(int lp_event=0;lp_event<event_num;lp_event++)
    {
      works->GetEntry(lp_event);
      for(int lp_slice=g_edg_cut;lp_slice<g_num_slice - g_edg_cut;lp_slice++)
	{
	  if(pstv[lp_slice]-ngtv[lp_slice]>0)
	    {
	      hst1->Fill(time[lp_slice],pstv[lp_slice]-ngtv[lp_slice]);	  
	    }
	}
    }
  peak_bin =hst1->GetMaximumBin();
  peak = hst1->GetBinCenter(peak_bin); 
  std::cout<<"peak Bin Number: "<<peak_bin<<std::endl;
  std::cout<<"In pre_analysis.h  :  "<<__LINE__<<std::endl;
 
  //Offset Search
  for(int lp_event=0;lp_event<event_num;lp_event++)
    {
      offset[lp_event] = 0;
      works->GetEntry(lp_event);
      for(int lp_slice=g_edg_cut; lp_slice < peak_bin - g_chrg_rng*smpl_frq ; lp_slice++)
	{
	  offset[lp_event] += (pstv[lp_slice]-ngtv[lp_slice]);
	  cnt++;
	}
      offset[lp_event] = offset[lp_event]/float(cnt);
      cnt = 0;
    }

  //Search Gravity Center
  for(int lp_event=0;lp_event<event_num;lp_event++)
    {
      works->GetEntry(lp_event);
      float timevolt = 0;
      float voltsum = 0;
      for(int lp_slice = peak_bin - g_chrg_rng*smpl_frq;lp_slice<peak_bin + g_chrg_rng*smpl_frq;lp_slice++)
	{
	  timevolt += (pstv[lp_slice] - ngtv[lp_slice] - offset[lp_event])*time[lp_slice];
	  voltsum  += pstv[lp_slice] - ngtv[lp_slice]- offset[lp_event];	  
	}
      wform_cg[lp_event] =timevolt/voltsum;
      
      if(wform_cg[lp_event] > peak+g_chrg_rng || wform_cg[lp_event] < peak-g_chrg_rng)
	{
	  wform_cg[lp_event]=g_num_slice*2;
	}
      else if(wform_cg[lp_event] < peak+g_chrg_rng && wform_cg[lp_event] > peak-g_chrg_rng)
	{
	  gc_av += wform_cg[lp_event] ;
	  cnt ++ ;
	}
    }
    for(int lp_event = 0 ;lp_event<event_num ;lp_event++ )
    {
      works->GetEntry(lp_event);
      if(wform_cg[lp_event] > g_num_slice || wform_cg[lp_event] < 0)
	{
	  wform_cg[lp_event] = gc_av;
	}
    }
    return hst1;
}

#endif
