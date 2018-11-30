#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include <sys/stat.h>
#include <unistd.h>
#include "iostream"
#include "fstream"
#include "stdlib.h"
#include "name_edit.h"
#include "geometric.h"

float smpl_frq = 5.0; //Sampling Speed of DRS4 [GHz]
int pls_wdth = 5;//About Pulse Width [ns] FWFM
int nSlice = 1024;//Number of DRS4 Channel
int k_num_slice =nSlice;
int k_edg_cut =20;

#include "chrg.h"
int main (int argc,char **argv)
{ 
  std::ofstream ofs(argv[2]);
  int lp_ch = atoi(argv[1]);

  for(int lp_file=3;lp_file<argc;lp_file++)
    {
      TString infilename =argv[lp_file];
      std::string date = argv[lp_file];
      date = date.substr(date.rfind("/")+7,6);
      float ttime = (atof(date.substr(0,2).c_str()) -17.)*3600.;
      ttime += atof(date.substr(2,2).c_str())*60. + atof(date.substr(4,2).c_str());
      TString treename ="TreeStabilityTest_0";
      TString ofilename ="NULL";
      float peak = 0;
      int  peak_bin = 0; 
      int  cnt  = 0;
      
      TFile* infile = TFile::Open(infilename);
      TTree* works  = (TTree*)infile->Get(treename);
      
      
      TH1F *hst2 = new TH1F("","",nSlice,0,nSlice/smpl_frq);
      TH1F *hst3 = new TH1F("","",nSlice,0,nSlice/smpl_frq);
      TString Pmt  = "NULL";
      TString Date = GetDate(infilename);
      std::cout<<"line : "<<__LINE__<<std::endl;
      //  TFile* infile = TFile::Open(infilename);
      std::cout<<"line : "<<__LINE__<<std::endl;
      std::cout<<treename<<std::endl;
      std::cout<<"line : "<<__LINE__<<std::endl;
      int    nEvent = works->GetEntries();
      std::cout<<"line : "<<__LINE__<<std::endl;
      float data_time[nSlice];
      float nwform[nSlice]     ;
      float pwform[nSlice]     ;
      float wform_cg[nEvent] ;
      float gc_av=0          ;
      std::cout<<"line : "<<__LINE__<<std::endl;
      float offset[nEvent] ;
      float pvalue[nEvent] ;
      float pbin[nEvent] ;
      bool  pulse[nEvent];
      std::cout<<"line : "<<__LINE__<<std::endl;
      
      //Data Set///////////////////////////////////////////////////////////////////////
      if(lp_ch == 0)
	{
	  Pmt = GetSerial1(infilename);
	}
      else if(lp_ch == 1)
	{
	  Pmt = GetSerial2(infilename);
	}      
      std::cout<<"line : "<<__LINE__<<std::endl;
      int data_stopcell = 0;
      works->SetBranchAddress("time",data_time);
      works->SetBranchAddress("stopcell",&data_stopcell);  
      if(lp_ch==0)
	{
	  works->SetBranchAddress("wform0",nwform);
	  works->SetBranchAddress("wform1",pwform);
	  std::cout<<"Selected_PMT_is,"<<Pmt<<std::endl;
	}
      else if(lp_ch==1)
	{
	  works->SetBranchAddress("wform2",nwform);
	  works->SetBranchAddress("wform3",pwform);
	  std::cout<<"Selected_PMT_is,"<<Pmt<<std::endl;
	} 
      std::cout<<"line : "<<__LINE__<<std::endl;
      //////////////////////////////////////////////////////////////////////////////////
      
      
      //Peak Search/////////////////////////////////////////////////////////////////////
      TH1F *hst1 = new TH1F("","",nSlice,0,nSlice/smpl_frq);
      for(int lp_event=0;lp_event<nEvent;lp_event++)
	{
	  works->GetEntry(lp_event);
	  pvalue[lp_event] = 0; 
	  for(int lp_slice=k_edg_cut;lp_slice<nSlice - k_edg_cut;lp_slice++)
	    {
	      if(pwform[lp_slice]-nwform[lp_slice]>3)
		{
		  hst1->Fill(data_time[lp_slice],pwform[lp_slice]-nwform[lp_slice]);	  
		}
	      if(pwform[lp_slice]-nwform[lp_slice]>pvalue[lp_event])
		{
		  pvalue[lp_event] = pwform[lp_slice]-nwform[lp_slice];
		  pbin[lp_event]   = lp_slice;
		}
	    }
	  
	  peak_bin =hst1->GetMaximumBin();
	  peak = hst1->GetBinCenter(peak_bin);
	} 
      std::cout<<"line : "<<__LINE__<<std::endl;
      /////////////////////////////////////////////////////////////////////////////////////
      
      
      //Offset Search//////////////////////////////////////////////////////////////////////
      for(int lp_event=0;lp_event<nEvent;lp_event++)
	{
	  works->GetEntry(lp_event);
	  offset[lp_event]= 0;
	  for(int lp_slice=k_edg_cut; lp_slice < peak_bin - pls_wdth*smpl_frq/2.; lp_slice++)
	    {
	      offset[lp_event] += (pwform[lp_slice]-nwform[lp_slice]);
	      cnt++;
	    }
	  offset[lp_event] = offset[lp_event]/float(cnt);
	}
      cnt = 0;
      
      
      ///////////////////////////////////////////////////////////////////////////////////////
      
      
      //Search Gravity Center/////////////////////////////////////////////////////////////////////////
      for(int lp_event=0;lp_event<nEvent;lp_event++)
	{
	  works->GetEntry(lp_event);
	  pulse[lp_event] = false;
	  float timevolt = 0;
	  float voltsum = 0;
	  for(int lp_slice = peak_bin - pls_wdth*smpl_frq/2;
	      lp_slice<peak_bin + pls_wdth*smpl_frq/2.;
	      lp_slice++)
	    {
	      timevolt += (pwform[lp_slice] - nwform[lp_slice] - offset[lp_event])*data_time[lp_slice];
	      voltsum  += pwform[lp_slice] - nwform[lp_slice]- offset[lp_event];	  
	    }
	  wform_cg[lp_event] =timevolt/voltsum;
	  if(wform_cg[lp_event] <= peak+pls_wdth/2. && wform_cg[lp_event] >= peak-pls_wdth/2.)
	    {
	      gc_av += wform_cg[lp_event] =timevolt/voltsum;
	      cnt ++ ;
	    }
	}
      std::cout<<cnt<<std::endl;
      gc_av = gc_av/float(cnt);
      for(int lp_event=0;lp_event<nEvent;lp_event++)
	{
	  if(wform_cg[lp_event] <= peak+pls_wdth/2. || wform_cg[lp_event] >= peak-pls_wdth/2.)
	    {
	      wform_cg[lp_event] =gc_av;
	    }
	}
      
      
      std::cout<<"line : "<<__LINE__<<std::endl;
      ///////////////////////////////////////////////////////////////////////////////////////////////
      
      
      
      //Create Average Waveform /////////////////////////////////////////////////////////////////////
      for(int lp_event = 0 ;lp_event<nEvent ;lp_event++ )
	{
	  works->GetEntry(lp_event);
	  if(offset[lp_event]>1000)
	    {	  
	      continue;
	    }
	  
	  for(int lp_slice = k_edg_cut;lp_slice < nSlice - k_edg_cut ;lp_slice ++)
	    {	  
	      hst2 -> Fill(data_time[lp_slice] - wform_cg[lp_event] + gc_av,
			   pwform[lp_slice] - nwform[lp_slice] - offset[lp_event]);
	      hst3 -> Fill(data_time[lp_slice] - wform_cg[lp_event] + gc_av);
	    }
	}
      hst2->Divide(hst2,hst3);
      ///////////////////////////////////////////////////////////////////////////////////////////////
      std::cout<<"#################################################"<<std::endl;
      ofs<<ttime<<" "<<hst2->GetMaximum()<<std::endl;
      std::cout<<"#################################################"<<std::endl;
      
      

      
      
      
      
      
      // Save Data /////////////////////////////////////////////////////////////////////////////////
      
      ///////////////////////////////////////////////////////////////////////////////////////////////
    }      
  return 0;
}

