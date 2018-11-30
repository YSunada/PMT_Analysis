#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include <sys/stat.h>
#include <unistd.h>
#include "iostream"
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
  TString infilename ="NULL";
  TString treename ="NULL";
  TString ofilename ="NULL";
  float peak = 0;
  int  peak_bin = 0; 
  int  cnt  = 0;
  int  ch = -1;
  
  
  //Option Control//////////////////////////////////////////////////////////////
  int opt;
  opterr = 0;
  bool cflag = false;
  bool pflag = false;
  bool hflag = false;
  bool eflag = false;
  bool oflag = false;
  bool Iflag = false;
  int xxx=0;
  while ((opt = getopt(argc,argv,"f:t:c:p:w:s:heoI")) != -1)
    {
      switch (opt)
	{
	case 'f':
	  infilename = optarg;
	  break;
	case 't':
	  treename = optarg;
	  break;
	case 'c':
	  ch = atoi(optarg);
	  if(ch>1 || ch <0)
	    {
	      ch =0;
	    }
	  cflag = true;
	  std::cout<<ch<<std::endl;
	  break;
	case 'p':
	  peak = atof(optarg);
	  pflag = true;
	  break;
	case 'w':
	  pls_wdth = atof(optarg);
	  break;
	case 's':
	  smpl_frq = atof(optarg);
	  break;
	case 'h':
	  hflag = true;
	  break;
	case 'e':
	  eflag = true;
	  break;
	case 'o':
	  std::cout<<"option o "<<std::endl;
	  oflag = true;
	  break;
	case 'I':
	  Iflag = true;
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
	       <<"  -t  : input treename (needed)  "<<std::endl
	       <<"  -c  : input DRS4 channnel 0 or 1 (default both chnnel)"<<std::endl
	       <<"  -p  : input peak time (this option skip peak seach)"<<std::endl
	       <<"  -w  : input about pulse width [ns] (default 5.0 [ns])"<<std::endl
	       <<"  -s  : input DRS4 sampling speed [GHz] (default 5.0 [GHz])"<<std::endl
	       <<"  -e  : this option cut error capasita "<<std::endl
      	       <<"  -o  : this option skip offset search  "<<std::endl
	       <<"  -I  : this option define integration　range arround average peak  "<<std::endl;
      std::cout<<"/////////////////////////////////////////////////////////////////"<<std::endl;
      return -1;
    }

  ////////////////////////////////////////////////////////////////////////////////

  TFile* infile = TFile::Open(infilename);
  TTree* works  = (TTree*)infile->Get(treename);

  for(int lp_ch = 0;lp_ch<2;lp_ch++)
    {
      if(cflag)
	{
	  lp_ch = ch;
	}

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
  
      //Noise Check/////////////////////////////////////////////////////////////////////  
      TH2F* Noise = new TH2F("","",nSlice,0,200,200,-20,20);
      for(int lp_event=0;lp_event<nEvent;lp_event++)
	{
	  works->GetEntry(lp_event);
	  for(int lp_slice=0;lp_slice<nSlice;lp_slice++)
	    {
	      if(eflag)
		{
		  if(lp_ch ==0)
		    {
		      if(
			 (lp_slice + data_stopcell)%1024 == 175  ||
			 (lp_slice + data_stopcell)%1024 == 32   ||
			 (lp_slice + data_stopcell)%1024 == 356  	        
			 //(lp_slice + data_stopcell)%1024 == 1016 ||
			 //(lp_slice + data_stopcell)%1024 == 158  ||
			 //(lp_slice + data_stopcell)%1024 == 460  ||
			 //(lp_slice + data_stopcell)%1024 == 677  || 
			 //(lp_slice + data_stopcell)%1024 ==  637
			 )
			{
			  continue;
			}
		    }
		  else if(lp_ch ==1)
		    {
		      if(
			 (lp_slice + data_stopcell)%1024 == 32   ||
			 (lp_slice + data_stopcell)%1024 == 677  ||
			 (lp_slice + data_stopcell)%1024 == 428  
			 //(lp_slice + data_stopcell)%1024 == 430  ||
			 //(lp_slice + data_stopcell)%1024 == 689  ||
			 //(lp_slice + data_stopcell)%1024 == 191  ||
			 //(lp_slice + data_stopcell)%1024 == 269  ||
			 //(lp_slice + data_stopcell)%1024 == 815  ||
			 //(lp_slice + data_stopcell)%1024 == 806  ||    
			 //(lp_slice + data_stopcell)%1024 == 1005 ||
			 //(lp_slice + data_stopcell)%1024 == 1008 || 
			 //(lp_slice + data_stopcell)%1024 == 1010
			 )
			{
			  continue;
			}
		    
		    }        
		}
	      Noise->Fill(data_time[lp_slice],pwform[lp_slice]-nwform[lp_slice]);
	    }
	}
  
      //Peak Search/////////////////////////////////////////////////////////////////////
      TH1F *hst1 = new TH1F("","",nSlice,0,nSlice/smpl_frq);
      if(!pflag)
	{
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
	    }
	  peak_bin =hst1->GetMaximumBin();
	  peak = hst1->GetBinCenter(peak_bin);
	} 
      else
	{
	  peak_bin = hst1->FindBin(peak);
	}
      
      //  std::cout<<"peak_Bin_Number:,"<<peak_bin<<std::endl;
      std::cout<<"line : "<<__LINE__<<std::endl;
      /////////////////////////////////////////////////////////////////////////////////////
  
  
      //Offset Search//////////////////////////////////////////////////////////////////////
      for(int lp_event=0;lp_event<nEvent;lp_event++)
	{
	  works->GetEntry(lp_event);
	  offset[lp_event]= 0;
	  if(oflag != true)
	    {
	      for(int lp_slice=k_edg_cut; lp_slice < peak_bin - pls_wdth*smpl_frq/2.; lp_slice++)
		{
		  offset[lp_event] += (pwform[lp_slice]-nwform[lp_slice]);
		  cnt++;
		}
	      offset[lp_event] = offset[lp_event]/float(cnt);
	    }
	      cnt = 0;
	}
  
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
  

      //Search Charge ////////////////////////////////////////////////////////////////////////////////
      float charg[nEvent];
      float wform[nSlice];
      for(int lp_event=0;lp_event<nEvent;lp_event++)
	{
	  float peak_pos = 0;
	  works->GetEntry(lp_event);
	  charg[lp_event] = 0;
	  if(offset[lp_event]>1000)
	    {
	      std::cout<<"line : "<<__LINE__<<std::endl;
	      continue;
	    }
	  if(Iflag)
	    {
	      peak_pos = peak;
	    }
	  else
	    {
	      peak_pos = wform_cg[lp_event];
	    }	  
	  for(int lp_slice=k_edg_cut; lp_slice < nSlice - k_edg_cut; lp_slice++)      
	    {
	      wform[lp_slice]=pwform[lp_slice]-nwform[lp_slice]-offset[lp_event];
	    }
	  charg[lp_event] = chrg(ch,wform,data_time, peak_pos - 5.,
				  peak_pos + 5.,data_stopcell);
	  if(TMath::Abs(charg[lp_event]) >100)
	    {
	      std::cout<<"charge : "<<charg[lp_event]<<std::endl
		       <<"offset : "<<offset[lp_event]<<std::endl
		       <<"gc     : "<<wform_cg[lp_event]<<std::endl;
	    }
	}

      //Search Width ////////////////////////////////////////////////////////////////////////////////
      float width[nEvent]; 
      for(int lp_event=0;lp_event<nEvent;lp_event++)
	{
	  works->GetEntry(lp_event);
	  float left=0;
	  float rght=0;
	  int   lbin=0;
	  int   rbin=0;
	  float phalf = (pvalue[lp_event]-offset[lp_event])/2.; 
	  for(int lp_slice = pbin[lp_event];
	      pwform[lp_slice]-nwform[lp_slice]-offset[lp_event]>phalf;lp_slice++)
	    {
	      rbin = lp_slice;
	    }
	  for(int lp_slice = pbin[lp_event];
	      pwform[lp_slice]-nwform[lp_slice]-offset[lp_event]>phalf;lp_slice--)
	    {
	      lbin = lp_slice;
	    }
	  left = xlinear(data_time[lbin],pwform[lbin]-nwform[lbin]-offset[lp_event],
			 data_time[lbin-1],pwform[lbin-1]-nwform[lbin-1]-offset[lp_event],phalf);
	  rght = xlinear(data_time[rbin],pwform[rbin]-nwform[rbin]-offset[lp_event],
			 data_time[rbin+1],pwform[rbin+1]-nwform[rbin+1]-offset[lp_event],phalf);
	  width[lp_event] = rght - left;
	  if(TMath::Abs(width[lp_event]>100))
	    {
	      std::cout<<lbin<<" : "<<rbin<<" : "<<phalf<<" : "<<pbin[lp_event]<<std::endl; 
	    }
	}
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
  
      // Save Data /////////////////////////////////////////////////////////////////////////////////
      std::cout<<"line : "<<__LINE__<<std::endl;  
      mkdir("PMT_Result",0777);
      TFile* ofile=new TFile("PMT_Result/pre"+Pmt +"_"+Date+".root","UPDATE");
      if(ofile->Get(treename))
	{
	  ofile->Delete(treename+";1");
	  std::cout<<"Tree Delete"<<std::endl;
	}
      if(ofile->Get(treename+"_AverageWform"))
	{
	  ofile->Delete(treename+"_AverageWform"+";1");
	  std::cout<<"Hist Delete"<<std::endl;
	}
      if(eflag)
	{
	  if(ofile->Get(treename+"_NoiseCheck_CapCut"))
	    {
	      ofile->Delete(treename+"_NoiseCheck_CapCut"+";1");
	      std::cout<<"Hist Delete"<<std::endl;
	    }
	  
	}
      else
	{
	  if(ofile->Get(treename+"_NoiseCheck"))
	    {
	      ofile->Delete(treename+"_NoiseCheck"+";1");
	      std::cout<<"Hist Delete"<<std::endl;
	    }
	}

      std::cout<<"line : "<<__LINE__<<std::endl;
      ofile->Close();
      ofile->Delete();
      std::cout<<"line : "<<__LINE__<<std::endl;
      TFile* wfile=new TFile("PMT_Result/pre"+Pmt +"_"+Date+".root","UPDATE");
      TTree* otree=new TTree(treename,treename);
      float ctrd =0;
      float ofst =0;
      float peakbin =0;
      float peakvalue =0;
      float wid=0;
      float charge=0;
      otree->Branch("cntrd",&ctrd,"ctrd/F");
      otree->Branch("offset",&ofst,"offset/F");
      otree->Branch("pbin",&peakbin,"peakbin/F");
      otree->Branch("pvalue",&peakvalue,"peakvalue/F");
      otree->Branch("width",&wid,"wid/F");
      otree->Branch("chrg",&charge,"charge/F");
      std::cout<<nEvent<<std::endl;
      std::cout<<wform_cg[100]<<" : "<<offset[100]<<std::endl;
      for(int lp_event=0;lp_event<nEvent;lp_event++)
	{
	  //if(pulse[lp_event])
	  {
	    ctrd = wform_cg[lp_event];
	    ofst = offset[lp_event];
	    peakbin = pbin[lp_event];
	    peakvalue=pvalue[lp_event];
	    wid = width[lp_event];
	    charge =charg[lp_event]; 
	    otree->Fill();
	  }
	}
      hst2->SetXTitle("time[ns]");
      hst2->SetYTitle("Pulse Height [mV]"); 
      hst2->GetXaxis()->CenterTitle();
      hst2->GetYaxis()->CenterTitle();
      hst2->GetYaxis()->SetTitleOffset(1.2); 
      hst2->SetName(treename+"_AverageWform");
      hst2->SetTitle(Pmt+"_"+Date+"_"+treename);
      hst2->Write();

      Noise->SetXTitle("time[ns]");
      Noise->SetYTitle("Pulse Height [mV]"); 
      Noise->GetXaxis()->CenterTitle();
      Noise->GetYaxis()->CenterTitle();
      Noise->GetYaxis()->SetTitleOffset(1.2);
      if(eflag)
	{
	  Noise->SetName(treename+"_NoiseCheck_CapCut"); 
	}
      else
	{
	  Noise->SetName(treename+"_NoiseCheck");
	}
      Noise->SetTitle(Pmt+"_"+Date+"_"+treename);
      Noise->Write();    
      wfile->Write();
      wfile->Close();
      if(cflag)
	{
	  break;
	}
    }
  infile->Close();
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  return 0;
}

