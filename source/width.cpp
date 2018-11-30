
/////////////////////////////////////////////
// File include                            //
/////////////////////////////////////////////

#include "iostream"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "string"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TString.h"
#include "TMath.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <complex>
#include <unistd.h>
//#include <regex>

/////////////////////////////////////////////
// Define constant                         //
/////////////////////////////////////////////
#include <global.h>
#include <pre_analysis.h> 
#include <width.h>
#include <name_edit.h>

/////////////////////////////////////////////
// Main Program                            //
/////////////////////////////////////////////

int main( int argc, char **argv)
{
  int opt;
  opterr = 0;
  int cflag = -1;
  while ((opt = getopt(argc,argv,"c:")) != -1)
    {
      switch (opt)
	{
	case 'c':
	  if(atoi(optarg) == 0 || atoi(optarg) == 1)
	    {
	      cflag = atoi(optarg);
	    }
	  break;
	}
    }

  for(int lp_file = 1;lp_file < argc;lp_file++)
    {
      std::cout <<std::endl<<"   "<<lp_file<<"/"<<argc-1<<" "
		<<"analysis start"<<std::endl<<std::endl
		<<argv[lp_file]<<std::endl;
      std::string filename = argv[lp_file];
      TString DataType = GetDataType(filename);
 
      if(DataType != "w")
	{
	  //	  std::cout<<"////////////////////////ERROR/////////////////////////"<<std::endl
	  //	   <<std::endl
	  //	   <<filename<<" don't have OnPhoto and Gain Data !!"<<std::endl
	  //	   <<std::endl
	  //	   <<"//////////////////////////////////////////////////////"<<std::endl;
	  continue;
	}
      
      mkdir("PMT_Result",0777);  
      
      TFile *openfile = TFile::Open(argv[lp_file]);

      for(int ch = 0;ch<2;ch++)
	{
	  if(cflag != -1)
	    {
	      ch = cflag;
	    }
		  
	  ///////////In <name_edit.h>,Get PMT Serial and Measurement Date///////
	  TString Serial = "NULL";
	  TString Date = "NULL";
	  GetSerialAndDate(&Serial,&Date,ch,filename);
	  //////////////////////////////////////////////////////////////////////
	  
	  //TApplication theApp("App",&argc,argv);

	  TGraphErrors  *grph = new TGraphErrors;  
	  TH1F* Width[6];
	  TString Voltage[6] = {"900","1000","1100","1200","1300","1400"};
	  float voltage[6]={900,1000,1100,1200,1300,1400};
	  TFile *savefile = new TFile("PMT_Result/Wid"+Serial+"_"+Date+".root","recreate");
	  for(int lp_volt = 5; lp_volt >= 0; lp_volt --)
	    {
	      //      TString treename ="TreeWidth" + Voltage[lp_volt] + "V_0";
	      TString treename ="TreeWidth" + Voltage[lp_volt] + "V_0"; 
	      if(lp_volt==6)	      
		{
		  treename ="TreeWidth1100Vhigh_0"; 
		}	
	      TTree *works  = (TTree *)openfile->Get(treename);
	      /////In <pre_analysis.h>, Make Average Waveform and //////////////////////
	      /////search each event's offset and gravity center   /////////////////////
	      int k_num_event = works->GetEntries();	      
	      float wform_cg[k_num_event] ;
	      float offset[k_num_event] ;     	      
	      TH1F*  AverageWform = pre_analysis(works,ch,offset,wform_cg);
	      //theApp.Run();
	      //////////////////////////////////////////////////////////////////////////
	      
	      //Search Width ///////////////////////////////////////////////////////////
	      Width[lp_volt]= WidthDist(works,ch,offset);
	      Width[lp_volt]->SetTitle("Width "+Voltage[lp_volt]);
	      Width[lp_volt]->SetName("Width_"+Voltage[lp_volt]);
	      Width[lp_volt]->Write();
	      grph->SetPoint(lp_volt,voltage[lp_volt],Width[lp_volt]->GetMean());
	      grph->SetPointError(lp_volt,0,Width[lp_volt]->GetMeanError());
	    }
	  TF1 *func_fit =new TF1("func_fit","[0]/(x-[1]) + [2]",900,1400);
	  grph->SetTitle("Width");
	  grph->SetName("grph_Width");
	  grph->SetMarkerSize(1);
	  grph->GetYaxis()->SetTitle("Width FWHM (ns)");
	  grph->GetXaxis()->SetTitle("Supply Voltage (V)");
	  grph->GetYaxis()->CenterTitle();
	  grph->GetXaxis()->CenterTitle();
	  grph->GetYaxis()->SetTitleOffset(1.3);
	  grph->SetMarkerStyle(4);
	  grph->SetMarkerSize(1);
	  grph->Fit("func_fit");
	  grph->Write();
	  savefile->Close();
	  if(cflag != -1)
	    {
	      break;
	    }
	}
      openfile->Close();
    }
  return 0;
}

