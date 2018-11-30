#ifndef INCLUDED_CHRG_RNG_SUNA
#define INCLUDED_CHRG_RNG_SUNA

#include<geometric.h>

void chrg_rng(TH1F* AverageWform,float* lft,float* rght)
{

  //Convert histgram to array
  float wform[k_num_slice] = {0};
  float time[k_num_slice]  = {0};      
  for(int lp_slice=k_edg_cut;lp_slice<k_num_slice - k_edg_cut ;lp_slice++)
    {
      wform[lp_slice] = AverageWform->GetBinContent(lp_slice+1);
      time[lp_slice]  = AverageWform->GetXaxis()->GetBinCenter(lp_slice+1);
    }

  int peak = AverageWform->GetMaximumBin();

  float timevolt = 0;
  float voltsum  = 0;
  float chrgsum  = 0;  
  float chrg_max = 0;
  float grvty_cnt= 0;
  float chrg[k_num_slice] = {0};

  for(int lp_slice=k_edg_cut+1;lp_slice<k_num_slice - k_edg_cut;lp_slice++)
    {
      chrgsum += Trapezoid(wform[lp_slice],wform[lp_slice-1],time[lp_slice]-time[lp_slice-1]);
      
      if(chrg_max<chrgsum)
	{
	  chrg_max=chrgsum;
	}
      chrg[lp_slice]=chrgsum;
      
      if((lp_slice >= peak - k_rng_cg) && (lp_slice<=peak + k_rng_cg))
	{    
	  voltsum +=wform[lp_slice];
	  timevolt +=time[lp_slice] * wform[lp_slice];
	}
    }
  grvty_cnt = timevolt/voltsum;


  int num1 = 0;

  for(int lp_slice=k_edg_cut;lp_slice<k_num_slice-k_edg_cut;lp_slice++)
    {
      if(chrg[lp_slice]>chrg_max*k_chrg_per_dwn)
	{
	  num1=lp_slice;
	  //	  std::cout<<"left range is "<<num1<<std::endl;
	  break;
	}
    }

  *lft = xlinear(time[num1],chrg[num1],time[num1-1],chrg[num1-1],chrg_max*k_chrg_per_dwn)-grvty_cnt;

  
  for(int lp_slice=k_edg_cut;lp_slice<k_num_slice-k_edg_cut;lp_slice++)
    {
      if(chrg[lp_slice]>chrg_max*k_chrg_per_up)
	{
	  num1=lp_slice;
	  //	  std::cout<<"right range is "<<num1<<std::endl;
	  break;
	}
    }

  *rght = xlinear(time[num1],chrg[num1],time[num1-1],chrg[num1-1],chrg_max*k_chrg_per_up)-grvty_cnt;

  
}


#endif
