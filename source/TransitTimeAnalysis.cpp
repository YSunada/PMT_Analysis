#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TGraphErrors.h"
#include <sys/stat.h>
#include <unistd.h>
#include "iostream"
#include "stdlib.h"
#include "name_edit.h"
#include "geometric.h"
#include "TF1.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TSystem.h"
//#include "pre_analysis.h"

float smpl_frq = 5.0; //Sampling Speed of DRS4 [GHz]
int pls_wdth = 5;//About Pulse Width [ns] FWFM
int nSlice = 1024;//Number of DRS4 Channel
int k_num_slice =nSlice;
int k_edg_cut =20;

#include "getopt.h"
#include "chrg.h"
int main (int argc,char **argv)
{ 
  TString infilename ="NULL";
  TString treename ="NULL";
  TString ofilename ="NULL";
  TString nfilename ="NULL";

  //Option Control//////////////////////////////////////////////////////////////
  int opt;
  opterr = 0;
  int Fm0 =0;
  int Fm1 =0;
  int fflag = 0;
  int tflag = 0;

  bool Fflag = false;
  bool hflag= false;
  bool nflag= false;
  while ((opt = getopt(argc,argv,"f:t:s:F:n:h")) != -1)
    {
      switch (opt)
	{
	case 'f':
	  infilename = optarg;
	  fflag =1;
	  break;
	case 't':
	  treename = optarg;
	  tflag = 1;
	  break;
	case 's':
	  smpl_frq = atof(optarg);
	  break;
	case 'F':
	  Fm0 = atoi(optarg)/10;
	  Fm1 = atoi(optarg)%10;
	  Fflag = true;
	  break;
	case 'N':
	  nfilename = optarg;
	  nflag = true;
	  break;
	case 'h':
	  hflag = true;
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
	       <<"  -s  : input DRS4 sampling speed [GHz] (default 5.0 [GHz])"<<std::endl
	       <<"  -F  : This option skipps waveform fitting for decide rise time   "<<std::endl
	       <<"      : Please input 00 ,01, 10 or 11                                        "<<std::endl
	       <<"      : 00 means channnel 0 connected PMT and channnel connected 1 PMT "<<std::endl
	       <<"      : 01 means channnel 0 connected PMT and channnel connected 1 SiPM "<<std::endl
	       <<"/////////////////////////////////////////////////////////////////"<<std::endl;
      return -1;
    }

 
  ////////////////////////////////////////////////////////////////////////////////

  TFile* infile = TFile::Open(infilename);
  TFile* nfile  = new TFile();
  std::cout<<"line : "<<__LINE__<<std::endl;
  if(nflag){nfile = TFile::Open(nfilename);}
  std::cout<<"line : "<<__LINE__<<std::endl;
  TTree* works  = (TTree*)infile->Get(treename);
  TString strPmtName0;
  TString strPmtName1;
  TString Date = GetDate(infilename);
  int    nEvent = works->GetEntries();
  if(nEvent>50000)  
    nEvent = 50000;
  std::cout<<"nEvent is "<<nEvent<<std::endl;
  float data_time[nSlice];
  float wformP0[nSlice]     ;
  float wformN0[nSlice]     ;
  float wformP1[nSlice]     ;
  float wformN1[nSlice]     ;
  float pdst0[nSlice]       ;
  float pdst1[nSlice]       ;
  float cerr0[nSlice]       ;
  float cerr1[nSlice]       ;

  std::cout<<"line : "<<__LINE__<<std::endl;
 
  //Data Set///////////////////////////////////////////////////////////////////////
  strPmtName0 = GetSerial1(infilename);
  strPmtName1 = GetSerial2(infilename);

  std::cout<<"line : "<<__LINE__<<std::endl;
  int data_stopcell = 0;
  works->SetBranchAddress("time",data_time);
  works->SetBranchAddress("stopcell",&data_stopcell);  
  
  works->SetBranchAddress("wform0",wformN0);
  works->SetBranchAddress("wform1",wformP0);

  works->SetBranchAddress("wform2",wformN1);
  works->SetBranchAddress("wform3",wformP1);  

  if(nflag)
    {
      TH2F* nhst0 = (TH2F*)nfile->Get("NoiseDist0");
      TH2F* nhst1 = (TH2F*)nfile->Get("NoiseDist1");
      for(int i=0;i<nSlice;i++)
	{
	  pdst0[i] = (float)(nhst0->ProjectionY("",i+1,i+1,"")->GetMean());
	  pdst1[i] = (float)(nhst1->ProjectionY("",i+1,i+1,"")->GetMean());
	  cerr0[i] = (float)(nhst0->ProjectionY("",i+1,i+1,"")->GetRMS());
	  cerr1[i] = (float)(nhst1->ProjectionY("",i+1,i+1,"")->GetRMS());
	}
    }
  else
    {
 for(int i=0;i<nSlice;i++)
	{
	  pdst0[i] = 0;
	  pdst1[i] = 0;
	  cerr0[i] = 0.6;
	  cerr1[i] = 0.6;
	}
    }
  std::cout<<"line : "<<__LINE__<<std::endl;
  //////////////////////////////////////////////////////////////////////////////////
  
  
  
  //Peak Search/////////////////////////////////////////////////////////////////////
  TH1F *hst1 = new TH1F("","",nSlice,0,nSlice/smpl_frq);
  TH1F *hst0 = new TH1F("","",nSlice,0,nSlice/smpl_frq);
  float peak0[nEvent];
  float peak1[nEvent];
  float pTime0[nEvent] ;
  float pTime1[nEvent] ;
  int pbin0[nEvent] ;
  int pbin1[nEvent] ; 

  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      works->GetEntry(lp_event);
      for(int lp_slice=k_edg_cut;lp_slice<nSlice - k_edg_cut;lp_slice++)
	{
	  if((wformP0[lp_slice]-wformN0[lp_slice]) > 0.0)
	    hst0->Fill(data_time[lp_slice],wformP0[lp_slice]-wformN0[lp_slice]);
	  
	  if((wformP1[lp_slice]-wformN1[lp_slice]) > 0.0)
	    hst1->Fill(data_time[lp_slice],wformP1[lp_slice]-wformN1[lp_slice]);
	}
    }
float pktm0 =  hst0->GetBinCenter(hst0->GetMaximumBin());
float pktm1 =  hst1->GetBinCenter(hst1->GetMaximumBin());

  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      float p0 = 0;
      float p1 = 0;
      works->GetEntry(lp_event);
      for(int lp_slice=k_edg_cut;lp_slice<nSlice - k_edg_cut;lp_slice++)
	{
	  if((data_time[lp_slice]>pktm0-10) && (data_time[lp_slice]<pktm0+10)){
	    if(wformP0[lp_slice]-wformN0[lp_slice]>p0)
	      {
		p0 = wformP0[lp_slice]-wformN0[lp_slice];
		peak0[lp_event] = wformP0[lp_slice]-wformN0[lp_slice];
		pTime0[lp_event] = data_time[lp_slice];
	      pbin0[lp_event]=lp_slice;
	      }
	  }
	  if((data_time[lp_slice]>pktm1-10) && (data_time[lp_slice]<pktm1+10))	    
	    {
	    if(wformP1[lp_slice]-wformN1[lp_slice]>p1)
	      {
		p1 = wformP1[lp_slice-1]-wformN1[lp_slice-1] + wformP1[lp_slice]-wformN1[lp_slice];
		peak1[lp_event] = wformP1[lp_slice]-wformN1[lp_slice];
		pTime1[lp_event] = data_time[lp_slice];
		pbin1[lp_event]=lp_slice;
	      }
	    }
	  //	  std::cout<<"peak bins are "<<pbin0[lp_event]<<" "<<pbin1[lp_event]<<std::endl;
	}
    }
  std::cout<<"line : "<<__LINE__<<std::endl;

  //Offset Search//////////////////////////////////////////////////////////////////////
  float offset0[nEvent];
  float offset1[nEvent];
  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      works->GetEntry(lp_event);
      offset0[lp_event] =0;
      offset1[lp_event] =0;
      int cnt0 = 0;
      int cnt1 = 0;
      for(int lp_slice=k_edg_cut; lp_slice < pbin0[lp_event] - 10*smpl_frq; lp_slice++)
	{
	  offset0[lp_event] += (wformP0[lp_slice]-wformN0[lp_slice]);
	  cnt0++;
	}
      for(int lp_slice=k_edg_cut; lp_slice < pbin1[lp_event] - 10*smpl_frq; lp_slice++)
	{
	  offset1[lp_event] += (wformP1[lp_slice]-wformN1[lp_slice]);
	  cnt1++;
	}
      //      std::cout<<"offset is "<<offset1[lp_event]/float(cnt1)<<std::endl;
      offset0[lp_event] = offset0[lp_event]/float(cnt0);
      offset1[lp_event] = offset1[lp_event]/float(cnt1);
      if(nflag){
	offset0[lp_event] = 0;
	offset1[lp_event] = 0;
      }    
    }

  std::cout<<"line : "<<__LINE__<<std::endl;  
  
  //// Rise Time Search ////////////
  float rTime0[nEvent];
  float dTime0[nEvent];
  int flag0 =0;

  float rTime1[nEvent];
  float dTime1[nEvent];
  int flag1 =0;
  
  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      works->GetEntry(lp_event);
      rTime0[lp_event] = 0;
      rTime1[lp_event] = 0;
      for(int lp_slice = pbin0[lp_event];lp_slice>k_edg_cut;lp_slice--)
	{
	  if(wformP0[lp_slice] - wformN0[lp_slice] - offset0[lp_event] < (peak0[lp_event] - offset0[lp_event])/10.)
	    {
	      if(wformP0[lp_slice - 1] - wformN0[lp_slice - 1]< wformP0[lp_slice] - wformN0[lp_slice])
		{
		  rTime0[lp_event] = xlinear(data_time[lp_slice],wformP0[lp_slice] - wformN0[lp_slice] - offset0[lp_event],
				   data_time[lp_slice+1],wformP0[lp_slice+1] - wformN0[lp_slice+1] - offset0[lp_event],
				   (peak0[lp_event] - offset0[lp_event])/10.);
		  break;
		}
	    }
	}
      for(int lp_slice = pbin0[lp_event];lp_slice<nSlice - k_edg_cut;lp_slice++)
	{
	  if(wformP0[lp_slice] - wformN0[lp_slice] - offset0[lp_event] < (peak0[lp_event] - offset0[lp_event])/2.)
	    {
	      //if(wformP0[lp_slice - 1] - wformN0[lp_slice - 1] - offset0[lp_event] < wformP0[lp_slice] - wformN0[lp_slice] - offset0[lp_event])
		{
		  dTime0[lp_event] = xlinear(data_time[lp_slice],wformP0[lp_slice] - wformN0[lp_slice] - offset0[lp_event],
						 data_time[lp_slice-1],wformP0[lp_slice-1] - wformN0[lp_slice-1] - offset0[lp_event],
						 (peak0[lp_event] - offset0[lp_event])/2.);
		  break;
		}
		dTime0[lp_event] = rTime0[lp_event]+5;
	    }
	}
      for(int lp_slice = pbin1[lp_event];lp_slice>k_edg_cut;lp_slice--)
	{
	  if(wformP1[lp_slice] - wformN1[lp_slice] - offset1[lp_event] < (peak1[lp_event] - offset1[lp_event])/10.)
	    {
	      if(wformP1[lp_slice - 1] - wformN1[lp_slice - 1] - offset1[lp_event] < wformP1[lp_slice] - wformN1[lp_slice] - offset1[lp_event])
		{
		  rTime1[lp_event] = xlinear(data_time[lp_slice],wformP1[lp_slice] - wformN1[lp_slice] - offset1[lp_event],
				   data_time[lp_slice+1],wformP1[lp_slice+1] - wformN1[lp_slice+1] - offset1[lp_event],
				   (peak1[lp_event] - offset1[lp_event])/10.);
		  break;	
		}
	    }
	}
      for(int lp_slice = pbin1[lp_event];lp_slice < nSlice-k_edg_cut;lp_slice++)
	{
	  if(wformP1[lp_slice] - wformN1[lp_slice] - offset1[lp_event] < (peak1[lp_event] - offset1[lp_event])/2.)
	    {
	      //if(wformP1[lp_slice - 1] - wformN1[lp_slice - 1] > wformP1[lp_slice] - wformN1[lp_slice])
		{
		  dTime1[lp_event] = xlinear(data_time[lp_slice],wformP1[lp_slice] - wformN1[lp_slice] - offset1[lp_event],
						 data_time[lp_slice-1],wformP1[lp_slice-1] - wformN1[lp_slice-1] - offset1[lp_event],
						 (peak1[lp_event] - offset1[lp_event])/2.);
		  break;
		}
		dTime1[lp_event] = rTime1[lp_event]+5;
	    }
	}
    }

 
  /////////////////////////////////////////////////////////////////////////////////////
  
  
  //Search Gravity Center/////////////////////////////////////////////////////////////////////////
  float  wform_cg0[nEvent];
  int cgbin0[nEvent];  float gc_av0=0;
  int cnt0=0;
  float  wform_cg1[nEvent];
  int cgbin1[nEvent];
  float gc_av1=0;
  int cnt1=0;
  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      works->GetEntry(lp_event);
      float timevolt0 = 0;
      float voltsum0 = 0;
      float timevolt1 = 0;
      float voltsum1 = 0;
      for(int lp_slice = pbin0[lp_event] - pls_wdth*smpl_frq/2;
	  lp_slice<pbin0[lp_slice] + pls_wdth*smpl_frq/2.;
	  lp_slice++)
	{
	  timevolt0 += (wformP0[lp_slice] - wformN0[lp_slice] - offset0[lp_event])*data_time[lp_slice];
	  voltsum0  += wformP0[lp_slice] - wformN0[lp_slice]- offset0[lp_event];	  
	}
      wform_cg0[lp_event] =timevolt0/voltsum0;
      if((wform_cg0[lp_event] <= pTime0[lp_event]+pls_wdth/2.) && (wform_cg0[lp_event] >= pTime0[lp_event]-pls_wdth/2.))
	{
	  gc_av0 += wform_cg0[lp_event] =timevolt0/voltsum0;
	  cnt0 ++ ;
	}

      for(int lp_slice = pbin1[lp_event] - pls_wdth*smpl_frq/2;
	  lp_slice<pbin1[lp_event] + pls_wdth*smpl_frq/2.;
	  lp_slice++)
	{
	  timevolt1 += (wformP1[lp_slice] - wformN1[lp_slice] - offset1[lp_event])*data_time[lp_slice];
	  voltsum1  += wformP1[lp_slice] - wformN1[lp_slice]- offset1[lp_event];	  
	}
      wform_cg1[lp_event] =timevolt1/voltsum1;
      if((wform_cg1[lp_event] <= pTime1[lp_event]+pls_wdth/2.) && (wform_cg1[lp_event] >= pTime1[lp_event]-pls_wdth/2.))
	{
	  gc_av1 += wform_cg1[lp_event] =timevolt1/voltsum1;
	  cnt1 ++ ;
	}
      
    }

  gc_av0 = gc_av0/float(cnt0);
  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      if((wform_cg0[lp_event] <= pTime0[lp_event]+pls_wdth/2.) || (wform_cg0[lp_event] >=pTime0[lp_event] -pls_wdth/2.))
	{
	  wform_cg0[lp_event] =gc_av0;
	}
    }

  
  gc_av1 = gc_av1/float(cnt1);
  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      if((wform_cg1[lp_event] <= pTime1[lp_event]+pls_wdth/2.) || (wform_cg1[lp_event] >=pTime0[lp_event] -pls_wdth/2.))
	{
	  wform_cg1[lp_event] =gc_av1;
	}
    }
  

  std::cout<<"line : "<<__LINE__<<std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////////////
  


  //search rise time /////////////////////////////////////////////////////////////////////////////
  TApplication app("app",&argc,argv);
  TCanvas* c = new TCanvas("ccc","ccc",800,400);
  c->Divide(2,1);
  float fit0[nEvent][5];
  float fit1[nEvent][5];
  float fitC0[nEvent];
  float fitC1[nEvent];
  float fitN0[nEvent];
  float fitN1[nEvent];
  TGraphErrors* grp0= new TGraphErrors();
  TGraphErrors* grp1= new TGraphErrors();
  //  TGraph* grp0 = new TGraph();
  //  TGraph* grp1 = new TGraph();
  std::cout<<"line : "<<__LINE__<<std::endl;
  
  TF1* f0;
  TF1* f1;
    std::cout<<"line : "<<__LINE__<<std::endl;
  
  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      works->GetEntry(lp_event);
      if(!Fflag)break;
      int n =0;
      int m =0;
      //for(int lp_slice= 0; data_time[lp_slice] < (wform_cg0[lp_event]+10.); lp_slice++)
      for(int lp_slice= k_edg_cut; lp_slice< nSlice - k_edg_cut; lp_slice++)      
	{
	  //if(data_time[lp_slice] > (wform_cg0[lp_event]-10.))
	  //if((wform_cg0[lp_event]-5. < data_time[lp_slice] )&&(data_time[lp_slice] < wform_cg0[lp_event]+1.))
	    {
	      grp0->SetPoint(n,data_time[lp_slice],wformP0[lp_slice] - wformN0[lp_slice] -offset0[lp_event] - pdst0[lp_slice]);	  
	      grp0->SetPointError(n,0,cerr0[lp_slice]);
	      n++;
	    }
	}
      
      //for(int lp_slice= pbin1[lp_event] - 20; lp_slice < (pbin1[lp_event] + 5); lp_slice++)
      for(int lp_slice= k_edg_cut; lp_slice< nSlice - k_edg_cut; lp_slice++)      
	{
	  //if((rTime1[lp_event]-10 <data_time[lp_slice])&&(data_time[lp_slice]<rTime1[lp_event]+10))
	    {
	      grp1->SetPoint(m,data_time[lp_slice],wformP1[lp_slice] - wformN1[lp_slice] -offset1[lp_event] - pdst1[lp_slice]);	  
	      grp1->SetPointError(m,0,cerr1[lp_slice]);
	      m++;
	    }
	}
      
      /////ch0 Fitting /////////////////
      if(Fm0 == 1)
	{	  
	  f0 = new TF1("f0","[0] + [1]/(exp(([2]-x)/[3]) + exp((x-[2])/[4]))",0,200);
	  f0->SetParameters(0,peak0[lp_event]*1.2,pTime0[lp_event] - 0.9 ,0.2,18);
	  grp0->Fit("f0","Q","",wform_cg0[lp_event]-15,wform_cg0[lp_event]+10);
	  for(int i = 0;i<5;i++){
	    fit0[lp_event][i] = f0->GetParameter(i);
	  }
	}
      else
	{
	  if(lp_event==0)std::cout<<"ch0 Fit by gaus function "<<std::endl;
	  f0 = new TF1("f0","gaus",0,200);
	  f0->SetParameters(peak0[lp_event],wform_cg0[lp_event]-0.5,0.9);
	  f0->SetParLimits(1,wform_cg0[lp_event] -10.,wform_cg0[lp_event]+10);
	  f0->SetParLimits(2,0.1,2.5);
	  grp0->Fit("f0","Q","",wform_cg0[lp_event]-15,wform_cg0[lp_event]+5);
	  float mean = f0->GetParameter(1);
	  float sigma = f0->GetParameter(2);
	  grp0->Fit("f0","Q","",mean - 2*sigma,mean + 0.8*sigma);
	  for(int i = 0;i<3;i++){
	    fit0[lp_event][i] = f0->GetParameter(i);
	  }
	  fit0[lp_event][3] = 0;
	  fit0[lp_event][4] = 0;
	}

      /////ch1 Fitting /////////////////
      if(Fm1 == 1)
	{
	  f1 = new TF1("f1","[0] + [1]/(exp(([2]-x)/[3]) + exp((x-[2])/[4]))",0,200);
	  f1->SetParameters(0,wform_cg1[lp_event]*1.2,wform_cg1[lp_event] - 0.9 ,0.2,18);
	  grp1->Fit("f1","Q","",wform_cg1[lp_event]-15,wform_cg1[lp_event] +10);	
	  //grp1->Fit("f1","Q");
	  for(int i = 0;i<5;i++){
	    fit1[lp_event][i] = f1->GetParameter(i);
	  }
	}
      else
	{
	  if(lp_event==0) std::cout<<"ch 1 Fit by gaus function "<<std::endl;
	  f1 = new TF1("f1","gaus",0,200);
	  f1->SetParameters(peak1[lp_event],wform_cg1[lp_event]-0.5,0.9);
	  f1->SetParLimits(1,wform_cg1[lp_event] -10.,wform_cg1[lp_event]+10.);
	  f1->SetParLimits(2,0.1,2.5);
	  grp1->Fit("f1","Q","",wform_cg1[lp_event]-15,wform_cg1[lp_event]+5);
	  float mean = f1->GetParameter(1);
	  float sigma = f1->GetParameter(2);
	  grp1->Fit("f1","Q","",mean - 2*sigma,mean + 0.8*sigma);
	  for(int i = 0;i<3;i++){
	    fit1[lp_event][i] = f1->GetParameter(i);
	  }
	  fit1[lp_event][3] = 0;
	  fit1[lp_event][4] = 0;
	}
      

      fitC0[lp_event] = f0->GetChisquare();
      fitC1[lp_event] = f1->GetChisquare();
      fitN0[lp_event] = f0->GetNDF();
      fitN1[lp_event] = f1->GetNDF();
      

      peak0[lp_event]  = f0->GetMaximum();
      pTime0[lp_event] = f0->GetX(peak0[lp_event]);
      rTime0[lp_event]  = f0->GetX(peak0[lp_event]/10.,0,pTime0[lp_event]);

      peak1[lp_event]  = f1->GetMaximum();
      pTime1[lp_event] = f1->GetX(peak1[lp_event]);
      rTime1[lp_event]  = f1->GetX(peak1[lp_event]/10.,0,pTime1[lp_event]);


      

      if(lp_event%(nEvent/20) == 0 )
	{
	  std::string title = "Title";	  
	  c->cd(1);
	  grp0->SetTitle(title.data());
	  grp0->GetXaxis()->SetRangeUser(pTime0[lp_event]-25,pTime0[lp_event]+15);
	  grp0->Draw("AP");
	  c->cd(2);
	  grp1->SetTitle(title.data());
	  grp1->GetXaxis()->SetRangeUser(pTime1[lp_event]-15,pTime1[lp_event]+15);
	  grp1->Draw("AP");
	  c->Modified();
	  c->Update();
	  gSystem->ProcessEvents();
	  //std::getchar(); 	  
	  //c->WaitPrimitive();
	  //std::cin<<
	  std::cout<<"*************** rise time and fall time **********************"<<std::endl;
	  std::cout<<rTime0[lp_event]<<"     "<<dTime0[lp_event]<<std::endl
		   <<rTime1[lp_event]<<"     "<<dTime1[lp_event]<<std::endl;
	  //app.Run();
	  std::cout<<"**************************************************************"<<std::endl;
	}
      f0->Delete();
      f1->Delete();
      grp0->Clear();
      grp1->Clear();

    }
  
  std::cout<<"line : "<<__LINE__<<std::endl;  
  
  // Save Data /////////////////////////////////////////////////////////////////////////////////
  std::cout<<"line : "<<__LINE__<<std::endl;  
  mkdir("PMT_Result",0777);
  TFile* ofile=new TFile("PMT_Result/TT"+strPmtName0 +"and"+ strPmtName1 +"_"+Date+".root","UPDATE");
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
  /*
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
  */
  
  std::cout<<"line : "<<__LINE__<<std::endl;
  ofile->Close();
  ofile->Delete();
  std::cout<<"line : "<<__LINE__<<std::endl;
  TFile* wfile=new TFile("PMT_Result/TT"+strPmtName0 +"and"+ strPmtName1 +"_"+Date+".root","UPDATE");
  TTree* otreeP=new TTree(strPmtName0+"_"+strPmtName1+"_"+treename,strPmtName0+"_"+strPmtName1+"_"+treename);
  //TTree* otreeP=new TTree(strPmtName0+treename,strPmtName0+treename);
  //TTree* otreeS=new TTree(strPmtName1+treename,strPmtName1+treename);
  float ofstP =0;
  int peakbinP =0;
  float peakvalueP =0;
  float peakTimeP =0;
  float NdofP=0;
  float chi2P =0;
  float parP[5]={0};
  float RiseP=0;
  float rtP=0;
  
  float ofstS =0;
  int peakbinS =0;
  float peakvalueS =0;
  float peakTimeS =0;
  float chi2S =0;
  float parS[5]={0};
  float NdofS=0;
  float RiseS=0;
  float rtS=0;  

  otreeP->Branch("offset_01",&ofstP,"ofstP/F");
  otreeP->Branch("pbin_01",&peakbinP,"pbinP/I");
  otreeP->Branch("pvalue_01",&peakvalueP,"peakvalueP/F");
  otreeP->Branch("ptime_01",&peakTimeP,"peakTimeP/F");
  otreeP->Branch("rtime_01",&rtP,"rtP/F");

  
  otreeP->Branch("offset_02",&ofstS,"ofstS/F");
  otreeP->Branch("pbin_02",&peakbinS,"pbinS/I");
  otreeP->Branch("pvalue_02",&peakvalueS,"peakvalueS/F");
  otreeP->Branch("ptime_02",&peakTimeS,"peakTimeS/F");
  otreeP->Branch("rtime_02",&rtS,"rtS/F");

  if(Fflag){
    otreeP->Branch("fit_chi2_01",&chi2P,"chi2P/F");
    otreeP->Branch("fit_Ndof_01",&NdofP,"NdofP/F");
    otreeP->Branch("fit_par_01",parP,"parP[5]/F");
    otreeP->Branch("fit_chi2_02",&chi2S,"chi2S/F");
    otreeP->Branch("fit_Ndof_02",&NdofS,"NdofS/F");
    otreeP->Branch("fit_par_02",parS,"parS[5]/F");
  }
  
  std::cout<<nEvent<<std::endl;
  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      //if(pulse[lp_event])
      {
	ofstP = offset0[lp_event];
	peakbinP = pbin0[lp_event];
	peakvalueP=peak0[lp_event];
	peakTimeP=pTime0[lp_event];
	rtP=rTime0[lp_event];




	ofstS = offset1[lp_event];
	peakbinS = pbin1[lp_event];
	peakvalueS=peak1[lp_event];
	peakTimeS=pTime1[lp_event];
	rtS=rTime1[lp_event];



	if(Fflag){
	  NdofP = fitN0[lp_event];
	  chi2P = fitC0[lp_event];
	  parP[0] = fit0[lp_event][0];
	  parP[1] = fit0[lp_event][1];
	  parP[2] = fit0[lp_event][2];
	  parP[3] = fit0[lp_event][3];
	  parP[4] = fit0[lp_event][4];
	  NdofS = fitN1[lp_event];
	  chi2S = fitC1[lp_event];
	  parS[0] = fit1[lp_event][0];
	  parS[1] = fit1[lp_event][1];
	  parS[2] = fit1[lp_event][2];
	  parS[3] = fit1[lp_event][3];
	  parS[4] = fit1[lp_event][4];
	}
	
	otreeP->Fill();
	//otreeS->Fill();
      }
    }

  wfile->Write();
  wfile->Close();
  infile->Close();
  
  return 0;
}

