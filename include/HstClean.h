#include "TH1F.h"
#include "TString.h"

TH1F* HstClean(TH1F* hst)
{
  TString name = hst->GetName();
  name =name+ "Clean";
  TH1F* hst_rtn = (TH1F*)hst->Clone(name);
  int bin_num = hst_rtn->GetNbinsX();
  float under=0;
  for(int i=0;i<bin_num;i++)
    {
      float a = hst_rtn->GetBinContent(i+1);
      if(a<0)
	{
	  under+=a;
	}
    }
  int bin_num_edg=0;
  for(int i=0;i<bin_num;i++)
    {
      float a = hst_rtn->GetBinContent(i+1);
      if(a>0)
	{
	  under+=a;
	}
      if(under>0)
	{
	  bin_num_edg=i+1;
	  break;
	}
    }
  
  for(int i=1;i<bin_num_edg;i++)
    {
      hst_rtn->SetBinContent(i,0.);
      hst_rtn->SetBinError(i,0.);
    } 
  hst_rtn->SetBinContent(bin_num_edg,under);

  return hst_rtn;
}
