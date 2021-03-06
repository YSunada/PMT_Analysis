/////////////////////////////////////////////
// File include                            //
/////////////////////////////////////////////
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TString.h"
#include "TMath.h"
#include "string"
#include "iostream"
#include <sys/stat.h>
#include <sys/types.h>
#include <complex>
#include <unistd.h>

/////////////////////////////////////////////
//User's include file                      //
/////////////////////////////////////////////
#include <global.h>
#include <pre_analysis.h> 
#include <chrg.h>
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
      std::string filename = argv[lp_file];      
      if(GetDataType(filename) != "OG")
	{
	  //continue;
	}
      TFile *openfile = TFile::Open(argv[lp_file]);      
      mkdir("PMT_Result",0777);        

      //for(int ch = 0;ch<2;ch++)
      for(int ch = 0;ch<1;ch++)
	{
	  if(cflag != -1)
	    {
	      //ch = cflag;
	    }
	  
	  ///////////In <name_edit.h>,Get PMT Serial and Measurement Date///////
	  TString Serial = "NULL";
	  TString Date = "NULL";
	  GetSerialAndDate(&Serial,&Date,ch,filename);
	  //////////////////////////////////////////////////////////////////////
	  std::cout<<__LINE__<<std::endl;
	  TGraphErrors  *grph_volt_chrg = new TGraphErrors;  
	  TGraphErrors  *grph_volt_gain = new TGraphErrors;  
	  TString Voltage[5] = {"1000","1100","1200","1300","1400"};
	  TH1F* hst_chrg[5];
	  TFile *savefile = new TFile("PMT_Result/Chrg"+Serial+"_"+Date+".root","RECREATE");
	  std::cout<<__LINE__<<std::endl;
	  for(int lp_volt = 0; lp_volt <5; lp_volt ++)
	    {
	      /////In <pre_analysis.h>, Make Average Waveform and //////////////////////
	      /////search each event's offset and gravity center   /////////////////////
	      //TString treename ="TreeMultiGain" + Voltage[lp_volt] + "V_0"; 
	      TString treename ="TreeTT_Fil5_" + Voltage[lp_volt] + "V_0"; 
	      std::cout<<__LINE__<<std::endl;
	      TTree *works  = (TTree*)openfile->Get(treename);
	      std::cout<<__LINE__<<std::endl;
	      int event_num = works->GetEntries();
	      std::cout<<__LINE__<<std::endl;
	      float wform_cg[event_num];
	      std::cout<<__LINE__<<std::endl;
	      float offset[event_num];
	      std::cout<<__LINE__<<std::endl;
	      TH1F* AverageWform = pre_analysis(works,ch,offset,wform_cg);
	      std::cout<<__LINE__<<std::endl;
	      //////////////////////////////////////////////////////////////////////////

	      std::cout<<__LINE__<<std::endl;	      
	      float peak_time = AverageWform->GetBinCenter(AverageWform->GetMaximumBin());
	      hst_chrg[lp_volt] = ChrgDist(works,ch,peak_time,offset,50);
	      hst_chrg[lp_volt]->SetName("ChargeDist"+Voltage[lp_volt]);
	      hst_chrg[lp_volt]->SetTitle("Charge Distribution "+Voltage[lp_volt]+"V_"+Serial+"_"+Date);
	      hst_chrg[lp_volt]->SetName("ChargeDist"+Voltage[lp_volt]);
	      hst_chrg[lp_volt]->Write();	      
	      float chrg_mean = hst_chrg[lp_volt]->GetMean();
	      float chrg_err  = hst_chrg[lp_volt]->GetMeanError();
	      std::cout<<Voltage[lp_volt]+"V Charge is "<< chrg_mean<<"mV*ns"<<std::endl;
	      grph_volt_chrg ->SetPoint(lp_volt,1000.+100.*float(lp_volt),chrg_mean);
	      grph_volt_chrg ->SetPointError(lp_volt,0,chrg_err);

	      ///////////////////////////////////////////////////
	      //                Saving Data                    //
	      ///////////////////////////////////////////////////
	      AverageWform->SetXTitle("time(ns)");
	      AverageWform->SetYTitle("pulse height(mV/event)");
	      AverageWform->GetXaxis()->CenterTitle();
	      AverageWform->GetYaxis()->CenterTitle();
	      AverageWform->GetYaxis()->SetTitleOffset(1.2);
	      AverageWform->SetName("AverageWform"+Voltage[lp_volt]+"V");
	      AverageWform->SetTitle("Average Wform "+Voltage[lp_volt]+"V "+Serial+"_"+Date);
	      AverageWform->Write();
	    }     

	  ///////////////////////////////////////////////////
	  //                Saving Data                    //
	  ///////////////////////////////////////////////////
	  grph_volt_chrg->SetTitle("Charge - Suply Voltage "+Serial+"_"+Date);
	  grph_volt_chrg->SetName("grph_Charge");
	  grph_volt_chrg->SetMarkerSize(1);
	  grph_volt_chrg->GetYaxis()->SetTitle("Charge (mV*ns)");
	  grph_volt_chrg->GetXaxis()->SetTitle("Supply Voltage (V)");
	  grph_volt_chrg->GetYaxis()->CenterTitle();
	  grph_volt_chrg->GetXaxis()->CenterTitle();
	  grph_volt_chrg->GetYaxis()->SetTitleOffset(1.3);
	  TF1 *HV_Chrg =new TF1("HV_Chrg","[0]*((x-350.0))**[1]",0,2000);
	  double a =4.80e-10;
	  HV_Chrg->SetParameters(a,4.29);
	  grph_volt_chrg ->Fit("HV_Chrg");
	  grph_volt_chrg->Write();
	  HV_Chrg->Write();     
	  savefile->Close();

	  delete grph_volt_chrg;
	  if(cflag != -1)

	    {
	     break;
	    }
	}
      openfile->Close();
	
    }
  return 0;
}

