#ifndef INCLUDED_CHRG_SUNA
#define INCLUDED_CHRG_SUNA

#include <global.h>
#include <geometric.h>

float chrg(int ch,float* wform, float* time, float lft, float rght, int stpcll/*,int drs4*/)
{
  float chrg = 0;
  float  wform_lft  = 0;
  float  wform_rght = 0;
  int num1 = 0;
  int num2 = 0;

  
for(int lp_slice = g_edg_cut; lp_slice < g_num_slice-g_edg_cut; lp_slice ++)
    {
      if((num1 == 0) && (time[lp_slice] > lft ))
	{
	  num1 = lp_slice;
	}

      if((num2 == 0) && (time[lp_slice] > rght ))
	{
	  num2 = lp_slice;
	}
    }
  for(int lp_slice = num1 ; lp_slice < num2-1 ; lp_slice ++)
    {
      if(ch ==0)
	{
	  if((lp_slice + stpcll)%1024 == 175 || (lp_slice + stpcll)%1024 == 158 ||
	     (lp_slice + stpcll)%1024 == 460 || (lp_slice + stpcll)%1024 == 430 || 
	     (lp_slice + stpcll)%1024 == 637 )
	    {
	      wform[lp_slice] = (wform[lp_slice-1]+wform[lp_slice+1])/2.;
	    }
	}
      else if(ch ==1)
	{
	  if((lp_slice + stpcll)%1024 == 430 || (lp_slice + stpcll)%1024 == 677  ||
	     (lp_slice + stpcll)%1024 == 677)
	    {
	      wform[lp_slice] = (wform[lp_slice-1]+wform[lp_slice+1])/2.;
	    }

	}
       chrg += Trapezoid(wform[lp_slice],wform[lp_slice+1],
      		time[lp_slice+1] - time[lp_slice ]);
    }
  wform_lft  = ylinear(time[num1-1],wform[num1-1],time[num1],wform[num1],lft);
  wform_rght = ylinear(time[num2-1],wform[num2-1],time[num2],wform[num2],rght);
  chrg += Trapezoid(wform[num1],wform_lft,time[num1]-lft);
  chrg += Trapezoid(wform[num2],wform_rght,rght-time[num2]);
  return chrg;
}



TH1F* ChrgDist(TTree* tree,int ch,float peak,float*offset,int n=1)
{
  int event_num = tree->GetEntries();
  int stopcell =0;
  float time[g_num_slice];
  float ngtv[g_num_slice];
  float pstv[g_num_slice];
  float charge[event_num];
  float chrg_max = -99999.;
  float chrg_min =  99999.;
  tree->SetBranchAddress("time",time);
  tree->SetBranchAddress("stopcell",&stopcell);  
  if(ch==0)
    {
      tree->SetBranchAddress("wform0",ngtv);
      tree->SetBranchAddress("wform1",pstv);
    }
  else if(ch==1)
    {
      tree->SetBranchAddress("wform2",ngtv);
      tree->SetBranchAddress("wform3",pstv);
    }
  
  for(int lp_event=0;lp_event<event_num;lp_event++)
    {
      tree->GetEntry(lp_event);
      float wform[g_num_slice] = {0};
      for(int lp_slice = 0; lp_slice < g_num_slice; lp_slice ++)
	{
	  wform[lp_slice] = pstv[lp_slice] - ngtv[lp_slice] - offset[lp_event];
	}
        charge[lp_event] = chrg(ch,wform,time,peak - g_chrg_rng,
					  peak + g_chrg_rng,stopcell);
	if(charge[lp_event] > chrg_max)
	  {
	    chrg_max = charge[lp_event];
	  }
	if(charge[lp_event] < chrg_min)
	  {
	    chrg_min = charge[lp_event];
	  }
    }
  TH1F* hst_chrg = new TH1F("Charge Distribution","chrg_dist",
			    (chrg_max - chrg_min)/0.1,chrg_min,chrg_max);
  for(int lp_event=0;lp_event<event_num;lp_event++)
    {
      hst_chrg->Fill(charge[lp_event]);
    }
  hst_chrg->Rebin(n);
  hst_chrg->SetXTitle("Charge (mV*ns)");
  hst_chrg->SetYTitle("count");
  hst_chrg->GetXaxis()->CenterTitle();
  hst_chrg->GetYaxis()->CenterTitle();
  hst_chrg->GetYaxis()->SetTitleOffset(1.2); 
  return hst_chrg;
}

#endif
