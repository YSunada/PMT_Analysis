#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TGraph.h"
#include <sys/stat.h>
#include <unistd.h>
#include "iostream"
#include "stdlib.h"
#include "name_edit.h"
#include "geometric.h"
#include "TF1.h"
#include "TApplication.h"
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
  //Option Control//////////////////////////////////////////////////////////////
  int opt;
  opterr = 0;
  int Sch = 0;
  int Pch = 0;
  int Sflag = 0;
  int Pflag = 0;
  int fflag = 0;
  int tflag = 0;

  bool Fflag = false;
  bool hflag= false;
  while ((opt = getopt(argc,argv,"f:t:s:S:P:Fh")) != -1)
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
	case 'S':
	  Sch = atof(optarg);
	  Sflag = 1;
	  break;
	case 'P':
	  Pflag = atof(optarg);
	  Pflag = 1;
	  break;
	case 'F':
	  Fflag = true;
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
	//	       <<"  -P  : input DRS4 channnel conected PMT  ,0 or 2   "<<std::endl
	       <<"  -S  : input DRS4 channnel conected SiPM ,0 to 3   "<<std::endl
	       <<"  -F  : This option skipps waveform fitting for decide rise time   "<<std::endl;
	std::cout<<"/////////////////////////////////////////////////////////////////"<<std::endl;
      return -1;
    }

 
  ////////////////////////////////////////////////////////////////////////////////

  TFile* infile = TFile::Open(infilename);
  TTree* works  = (TTree*)infile->Get(treename);
  TString strSi1 = "NULL";
  TString strSi2= "NULL";
  TString Date = GetDate(infilename);
  int    nEvent = works->GetEntries();
  std::cout<<"nEvent is "<<nEvent<<std::endl;
  float data_time[nSlice];
  float Si1P[nSlice]     ;
  float Si1N[nSlice]     ;
  float Si2P[nSlice]     ;
  float Si2N[nSlice]     ;
  std::cout<<"line : "<<__LINE__<<std::endl;
 
  //Data Set//////////////////////////////////////////////////////////////////////

  int data_stopcell = 0;
  works->SetBranchAddress("time",data_time);
  works->SetBranchAddress("stopcell",&data_stopcell);  

  works->SetBranchAddress("wform0",Si1N);
  works->SetBranchAddress("wform1",Si1P);
  works->SetBranchAddress("wform2",Si2N);
  works->SetBranchAddress("wform3",Si2P);

  strSi1 = "SiPM01";
  strSi2 = "SiPM02";
  std::cout<<"line : "<<__LINE__<<std::endl;
 



  //Peak Search/////////////////////////////////////////////////////////////////////
  TH1F *hst1 = new TH1F("","",nSlice,0,nSlice/smpl_frq);
  float peakS1[nEvent];
  float peakS2[nEvent];
  float pTimeS1[nEvent] ;
  float pTimeS2[nEvent] ;
  int pbinS1[nEvent] ;
  int pbinS2[nEvent] ;
  float riseS1[nEvent];
  float riseS2[nEvent];
  float rTimeS1[nEvent];
  float rTimeS2[nEvent];
  std::cout<<"nEvent is "<<nEvent<<std::endl;
  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      works->GetEntry(lp_event);
      peakS1[lp_event] = 0; 
      peakS2[lp_event] = 0;
      float p = 0;
      float s = 0;
      float rs1 = 0;
      float rs2 = 0;
      for(int lp_slice=k_edg_cut;lp_slice<nSlice - k_edg_cut;lp_slice++)
	{
	  //// Peak Time Search ////////////
	 
	  if(Si1P[lp_slice-1]-Si1N[lp_slice-1] + Si1P[lp_slice]-Si1N[lp_slice]
	     + Si1P[lp_slice+1]-Si1N[lp_slice+1] > p)
	    {
	      p = Si1P[lp_slice-1]-Si1N[lp_slice-1] + Si1P[lp_slice]-Si1N[lp_slice]
		+ Si1P[lp_slice+1]-Si1N[lp_slice+1];
	      peakS1[lp_event] = Si1P[lp_slice]-Si1N[lp_slice];
	      pTimeS1[lp_event] = data_time[lp_slice];
	      pbinS1[lp_event]=lp_slice;
	    }
	  if(Si2P[lp_slice-1]-Si2N[lp_slice-1] + Si2P[lp_slice]-Si2N[lp_slice]
	     + Si2P[lp_slice+1]-Si2N[lp_slice+1] > s)
	    {
	      s = Si2P[lp_slice-1]-Si2N[lp_slice-1] + Si2P[lp_slice]-Si2N[lp_slice]
		+ Si2P[lp_slice+1]-Si2N[lp_slice+1];
	      peakS2[lp_event] = Si2P[lp_slice]-Si2N[lp_slice];
	      pTimeS2[lp_event] = data_time[lp_slice];
	      pbinS2[lp_event]=lp_slice;
	    }
	}
	  //// Rise Time Search ////////////
      for(int lp_slice=k_edg_cut;lp_slice<pbinS1[lp_event];lp_slice++)
	{ 
	  float rsb1 = ((Si1P[lp_slice]-Si1N[lp_slice])- (Si1P[lp_slice-1]-Si1N[lp_slice-1]))/
	    (data_time[lp_slice] - data_time[lp_slice-1]);
	  float rsf1 = ((Si1P[lp_slice+1]-Si1N[lp_slice+1])- (Si1P[lp_slice]-Si1N[lp_slice]))/
	    (data_time[lp_slice+1] - data_time[lp_slice]);
	  float rsff1 = ((Si1P[lp_slice+2]-Si1N[lp_slice+2])- (Si1P[lp_slice+1]-Si1N[lp_slice+1]))/
	    (data_time[lp_slice+2] - data_time[lp_slice+1]);

	  if((rsff1+rsf1+rsb1 > rs1) && (rsf1>0 && rsb1 >0))
	    {
	      rs1 = rsb1+rsf1;
	      riseS1[lp_event] = rsf1;
	      rTimeS1[lp_event] = data_time[lp_slice];
	    } 
	}

      for(int lp_slice=k_edg_cut;lp_slice<pbinS2[lp_event];lp_slice++)
	{ 
	  float rsb2 = ((Si2P[lp_slice]-Si2N[lp_slice])- (Si2P[lp_slice-1]-Si2N[lp_slice-1]))/
	    (data_time[lp_slice] - data_time[lp_slice-1]);
	  float rsf2 = ((Si2P[lp_slice+1]-Si2N[lp_slice+1])- (Si2P[lp_slice]-Si2N[lp_slice]))/
	    (data_time[lp_slice+1] - data_time[lp_slice]);
	  float rsff2 = ((Si2P[lp_slice+2]-Si2N[lp_slice+2])- (Si2P[lp_slice+1]-Si2N[lp_slice+1]))/
	    (data_time[lp_slice+2] - data_time[lp_slice+1]);

	  if((rsff2+rsf2+rsb2 > rs2) && (rsf2>0 && rsb2 >0))
	    {
	      rs2 = rsb2+rsf2;
	      riseS2[lp_event] = rsf2;
	      rTimeS2[lp_event] = data_time[lp_slice];
	    } 
	}	  
      
    }

  
 
  /////////////////////////////////////////////////////////////////////////////////////
  
  
  //Offset Search//////////////////////////////////////////////////////////////////////
  float offsetS1[nEvent];
  float offsetS2[nEvent];
  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      works->GetEntry(lp_event);
      offsetS1[lp_event] =0;
      offsetS2[lp_event] =0;
      int cntS1 = 0;
      int cntS2 = 0;
      for(int lp_slice=k_edg_cut; lp_slice < pbinS1[lp_event] - 10*smpl_frq; lp_slice++)
	{
	  offsetS1[lp_event] += (Si1P[lp_slice]-Si1N[lp_slice]);
	  cntS1++;
	}
      for(int lp_slice=k_edg_cut; lp_slice < pbinS2[lp_event] - 10*smpl_frq; lp_slice++)
	{
	  offsetS2[lp_event] += (Si2P[lp_slice]-Si2N[lp_slice]);
	  cntS2++;
	}


      offsetS1[lp_event] = offsetS1[lp_event]/float(cntS1);
      offsetS2[lp_event] = offsetS2[lp_event]/float(cntS2);
    }
  


  //serch rise time /////////////////////////////////////////////////////////////////////////////


  /*
    //TApplication app("app",&argc,argv);
  float fitP[nEvent][5];
  float fitS[nEvent][5];
  float fitCP[nEvent];
  float fitCS[nEvent];
  float fitNP[nEvent];
  float fitNS[nEvent];
  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      if(!Fflag)break;
      int n =0;
      int m =0;
      TGraph* grpPmt  = new TGraph();
      TGraph* grpSiPM = new TGraph();
      TF1* f0 = new TF1("f0","[0] + [1]/(exp(([2]-x)/[3]) + exp((x-[2])/[4]))",0,200);
      TF1* f1 = new TF1("f1","[0] + [1]/(exp(([2]-x)/[3]) + exp((x-[2])/[4]))",0,200);
      for(int lp_slice=(pbinP[lp_event] - 10.*smpl_frq)/1; lp_slice < pbinP[lp_event] + 10.*smpl_frq; lp_slice++)      
	{	  
	  grpPmt->SetPoint(lp_slice,data_time[lp_slice],PmtP[lp_slice] - PmtN[lp_slice] -offsetP[lp_event]);	  
	}
      for(int lp_slice=(pbinP[lp_event] - 10.*smpl_frq)/1; lp_slice < pbinP[lp_event] + 10.*smpl_frq; lp_slice++)
	{	  
	  grpSiPM->SetPoint(lp_slice-k_edg_cut,data_time[lp_slice],SiPM[lp_slice]  -offsetS[lp_event]);
	}
      f0->SetParameters(0,peakP[lp_event],pTimeP[lp_event],1,5);
      f1->SetParameters(0,peakS[lp_event],pTimeS[lp_event],0.3,15);
      grpPmt->Fit("f0","Q");
      grpSiPM->Fit("f1","Q");
      
      for(int i = 0;i<5;i++){
	fitP[lp_event][i] = f0->GetParameter(i);
	fitS[lp_event][i] = f1->GetParameter(i);
      }
      fitCS[lp_event] = f0->GetChisquare();
      fitCP[lp_event] = f1->GetChisquare();
      fitNS[lp_event] = f0->GetNDF();
      fitNP[lp_event] = f1->GetNDF();
      
      
      f1->Delete();
      f0->Delete();
      grpPmt->Delete();
      grpSiPM->Delete();
      //std::cout<<"line : "<<__LINE__<<std::endl;
      if(lp_event%(nEvent/100) == 0 )
	std::cout<<lp_event<<"/"<<nEvent<<"  Fitting completed "<<std::endl;
    }
  
  std::cout<<"line : "<<__LINE__<<std::endl;  
  
  */

  // Save Data /////////////////////////////////////////////////////////////////////////////////
  std::cout<<"line : "<<__LINE__<<std::endl;  
  mkdir("PMT_Result",0777);
  TFile* ofile=new TFile("PMT_Result/SiPMTTCal_"+Date+".root","UPDATE");
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
  TFile* wfile=new TFile("PMT_Result/SiPMTTCal_"+Date+".root","UPDATE");
  TTree* otreeP=new TTree(strSi1+treename,strSi1+treename);
  TTree* otreeS=new TTree(strSi2+treename,strSi2+treename);
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

  otreeP->Branch("offset",&ofstP,"ofstP/F");
  otreeP->Branch("pbin",&peakbinP,"pbinP/I");
  otreeP->Branch("pvalue",&peakvalueP,"peakvalueP/F");
  otreeP->Branch("ptime",&peakTimeP,"peakTimeP/F");

  otreeP->Branch("rise",&RiseP,"RiseP/F");
  otreeP->Branch("risetime",&rtP,"rtP/F");

  
  otreeS->Branch("offset",&ofstS,"ofstS/F");
  otreeS->Branch("pbin",&peakbinS,"pbinS/I");
  otreeS->Branch("pvalue",&peakvalueS,"peakvalueS/F");
  otreeS->Branch("ptime",&peakTimeS,"peakTimeS/F");
  otreeS->Branch("rise",&RiseS,"RiseS/F");
  otreeS->Branch("risetime",&rtS,"rtS/F");

  /*
  if(Fflag){
    otreeP->Branch("fit_chi2",&chi2P,"chi2P/F");
    otreeP->Branch("fit_Ndof",&NdofP,"NdofP/F");
    otreeP->Branch("fit_par",parP,"parP[5]/F");
    otreeS->Branch("fit_chi2",&chi2S,"chi2S/F");
    otreeS->Branch("fit_Ndof",&NdofS,"NdofS/F");
    otreeS->Branch("fit_par",parS,"parS[5]/F");
  }
  */

  std::cout<<nEvent<<std::endl;
  for(int lp_event=0;lp_event<nEvent;lp_event++)
    {
      //if(pulse[lp_event])
      {
	ofstP = offsetS1[lp_event];
	peakbinP = pbinS1[lp_event];
	peakvalueP=peakS1[lp_event];
	peakTimeP=pTimeS1[lp_event];

	RiseP   = riseS1[lp_event];
	rtP  = rTimeS1[lp_event];


	ofstS = offsetS2[lp_event];
	peakbinS = pbinS2[lp_event];
	peakvalueS=peakS2[lp_event];
	peakTimeS=pTimeS2[lp_event];

	RiseS   = riseS2[lp_event];
	rtS  = rTimeS2[lp_event];
	/*
	if(Fflag){
	  NdofP = fitNP[lp_event];
	  chi2P = fitCP[lp_event];
	  parP[0] = fitP[lp_event][0];
	  parP[1] = fitP[lp_event][1];
	  parP[2] = fitP[lp_event][2];
	  parP[3] = fitP[lp_event][3];
	  parP[4] = fitP[lp_event][4];
	  NdofS = fitNS[lp_event];
	  chi2S = fitCS[lp_event];
	  parS[0] = fitS[lp_event][0];
	  parS[1] = fitS[lp_event][1];
	  parS[2] = fitS[lp_event][2];
	  parS[3] = fitS[lp_event][3];
	  parS[4] = fitS[lp_event][4];
	}
	*/
	otreeP->Fill();
	otreeS->Fill();
      }
    }

  wfile->Write();
  wfile->Close();
  infile->Close();
  
  return 0;
}

