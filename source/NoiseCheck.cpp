#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "name_edit.h"
#include "iostream"
#include <sys/stat.h>
#include <unistd.h>
#include "stdlib.h"

float smpl_frq = 5.0; //Sampling Speed of DRS4 [GHz]
int pls_wdth = 5;//About Pulse Width [ns] FWFM
int nSlice = 1024;//Number of DRS4 Channel
int k_num_slice =nSlice;
int k_edg_cut =20;

#include "getopt.h"
#include "chrg.h"


void ShowHelp();

int main (int argc,char **argv)
{
  TString infilename = "NULL" ;
  TString treename   = "NULL" ;

  int opt;
  opterr = 0;
  while ((opt = getopt(argc,argv,"f:t:h")) != -1)
    {
      switch (opt)
	{
	case 'f':
	  infilename = optarg;
	  break;
	case 't':
	  treename = optarg;
	  break;
	case 'h':
	  ShowHelp();
	  return -1;
	  break;
	}
    }
  if(argc <2)
    {
      ShowHelp();
      return -1;
    }

  TFile* infile = TFile::Open(infilename);
  TTree* works  = (TTree*)infile->Get(treename);
  TString strPmtName0 = GetSerial1(infilename);
  TString strPmtName1 = GetSerial2(infilename);
  TString Date = GetDate(infilename);
  int    nEvent = works->GetEntries();
  
  TH2F* hst0 = new TH2F("NoiseDist0","Noise Distribution" + strPmtName0,1024,0,200,400,-20,20);
  TH2F* hst1 = new TH2F("NoiseDist1","Noise Distribution" + strPmtName1,1024,0,200,400,-20,20);
  
  works->Draw("wform1-wform0:time>>NoiseDist0","","",nEvent,0);
  works->Draw("wform3-wform2:time>>NoiseDist1","","",nEvent,0);
  //  hst0->SetName(treename+"_ch0");  
  //  hst1->SetName(treename+"_ch1");  

  ////// Save Data /////////////
  mkdir("PMT_Result",0777);  
  TFile* ofile=new TFile("PMT_Result/Noise"+strPmtName0 +"and"+ strPmtName1 +"_"+Date+".root","UPDATE");
  hst0->Write();
  hst1->Write();
  ofile->Write();
  //////////////////////////////

  infile->Close();
  ofile->Close();

  return 0;

}


void ShowHelp(){
      std::cout<<"/////////////////////////////////////////////////////////////////"<<std::endl
	       <<"  Option List                "<<std::endl
	       <<"  -f  : input filename (needed)  "<<std::endl
	       <<"  -t  : input treename (needed)  "<<std::endl
	       <<"  -h  : Show this message   "<<std::endl
	       <<"/////////////////////////////////////////////////////////////////"<<std::endl;
}






