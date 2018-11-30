#ifndef INCLUDED_WIDTH_SUNA
#define INCLUDED_WIDTH_SUNA

#include <global.h>
#include <geometric.h>

float width(int ch,float* wform, float* time,int stpcll)
{
  int num1 = 0;
  int num2 = 0;
  int peak_bin=0; 
  float peak=0;  

  for(int lp_slice = g_edg_cut; lp_slice < g_num_slice-g_edg_cut; lp_slice ++)
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
    }

  for(int lp_slice = g_edg_cut; lp_slice < g_num_slice-g_edg_cut; lp_slice ++)
    {
      if(peak<wform[lp_slice])
	{
	  peak = wform[lp_slice];
	  peak_bin = lp_slice;
	}
    }

  for(int lp_slice = peak_bin; lp_slice < g_num_slice-g_edg_cut; lp_slice--)
    {
      if(wform[lp_slice]<peak/2.)
	{
	  num1=lp_slice;
	  break;
	}

    }

  for(int lp_slice = peak_bin; lp_slice < g_num_slice-g_edg_cut; lp_slice++)
    {
      if(wform[lp_slice]<peak/2.)
	{
	  num2=lp_slice;
	  break;
	}      
    }  
  float  wform_lft  = xlinear(time[num1],wform[num1],time[num1+1],wform[num1+1],peak/2.);
  float  wform_rght = xlinear(time[num2-1],wform[num2-1],time[num2],wform[num2],peak/2.);
  if(peak<6.0)
    {
      return -1;
    }
  else
    {
      return (wform_rght - wform_lft);
    }
}



TH1F* WidthDist(TTree* tree,int ch,float*offset)
{
  int event_num = tree->GetEntries();
  int stopcell =0;
  float time[g_num_slice];
  float ngtv[g_num_slice];
  float pstv[g_num_slice];
  float wid[event_num];
  float wid_max = -99999.;
  float wid_min =  99999.;
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
        wid[lp_event] = width(ch,wform,time,stopcell);
	if(wid[lp_event] > wid_max)
	  {
	    wid_max = wid[lp_event];
	  }
	if(wid[lp_event] < wid_min && wid[lp_event] >0)
	  {
	    wid_min = wid[lp_event];
	  }
    }
  TH1F* hst_wid = new TH1F("Width Distribution","width_dist",
			   TMath::Sqrt(event_num),wid_min,wid_max);
  for(int lp_event=0;lp_event<event_num;lp_event++)
    {
      if(wid[lp_event]>0)
	{
	  hst_wid->Fill(wid[lp_event]);
	}
    }
  hst_wid->SetXTitle("Width FWHM (ns)");
  hst_wid->SetYTitle("count");
  hst_wid->GetXaxis()->CenterTitle();
  hst_wid->GetYaxis()->CenterTitle();
  hst_wid->GetYaxis()->SetTitleOffset(1.2); 
  return hst_wid;
}

#endif
