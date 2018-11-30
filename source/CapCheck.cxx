#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH2F.h"
#include "iostream"
#include "TMath.h"
int  ybin = 400;
float ylow =-20;
float yup  = 20;
int xbin = 1024;
float xlow = 0.5;
float xup  = 1024.5;
int nslice = 0;

void CapCheck(TString filename,TString treename)
{
  //  TApplication theApp("theApp",)
  TFile* file = TFile::Open(filename);
  TTree* tree = (TTree*)(file->Get(treename));
  int numeve  = tree->GetEntries();
  float wform0[numeve];
  float wform1[numeve];
  float wform2[numeve];
  float wform3[numeve];
  float d_time[numeve];
  int   stopcell = 0;
  tree->SetBranchAddress("time",d_time);
  tree->SetBranchAddress("stopcell",&stopcell);
  tree->SetBranchAddress("wform0",wform0);
  tree->SetBranchAddress("wform1",wform1);
  tree->SetBranchAddress("wform2",wform2);
  tree->SetBranchAddress("wform3",wform3);
  TH2F* hst0 = new TH2F("hst0","ch 0;cpasita number;pulse height",
			xbin,xlow,xup,ybin,ylow,yup);
  TH2F* hst1 = new TH2F("hst1","ch 1;cpasita number;pulse height",
			xbin,xlow,xup,ybin,ylow,yup);
  TH2F* hst2 = new TH2F("hst2","ch 2;cpasita number;pulse height",
			xbin,xlow,xup,ybin,ylow,yup);
  TH2F* hst3 = new TH2F("hst3","ch 3;cpasita number;pulse height",
			xbin,xlow,xup,ybin,ylow,yup);

  TH1F* hsta = new TH1F("hst0","ch 0;cpasita number",xbin,xlow,xup);
  TH1F* hstb = new TH1F("hst1","ch 1;cpasita number",xbin,xlow,xup);
  
  for(int lp_eve = 0;lp_eve<numeve;lp_eve ++)
    {
      tree->GetEntry(lp_eve);
      for(int lp_cap = 20;lp_cap<1000;lp_cap++)
	{
	  /*
	  hst0->Fill((lp_cap+stopcell)%1024,wform0[lp_cap]);
	  hst1->Fill((lp_cap+stopcell)%1024,wform1[lp_cap]);
	  hst2->Fill((lp_cap+stopcell)%1024,wform2[lp_cap]);
	  hst3->Fill((lp_cap+stopcell)%1024,wform3[lp_cap]);
	 */
	  if(TMath::Abs(wform1[lp_cap]-wform0[lp_cap])>5)
	    {
	      hsta->Fill((lp_cap+stopcell)%1024);
	    }
	  if(TMath::Abs(wform3[lp_cap]-wform2[lp_cap])>5)
	    {
	      hstb->Fill((lp_cap+stopcell)%1024);
	    }
	}
    }
  /*
  TCanvas* cvs0 = new TCanvas("cvs0","cvs0",600,600);
  cvs0->cd();
  hst0->Draw("colz");
  TCanvas* cvs1 = new TCanvas("cvs1","cvs1",600,600);
  cvs1->cd();
  hst1->Draw("colz");
  TCanvas* cvs2 = new TCanvas("cvs2","cvs2",600,600);
  cvs2->cd();
  hst2->Draw("colz");
  TCanvas* cvs3 = new TCanvas("cvs3","cvs3",600,600);
  cvs3->cd();
  hst3->Draw("colz");
  */

  TCanvas* cvs0 = new TCanvas("cvs0","cvs0",600,600);
  cvs0->cd();
  hsta->Draw();
  TCanvas* cvs1 = new TCanvas("cvs1","cvs1",600,600);
  cvs1->cd();
  hstb->Draw();

}
