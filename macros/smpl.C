
void smpl(TString name)
{
  TFile* file = TFile::Open(name+".root","recreate");
  TCanvas* cnv0 =new TCanvas(name+"_cnv0",name+"_cnv0",400,400);
  TCanvas* cnv1 =new TCanvas(name+"_cnv1",name+"_cnv1",400,400);
  TCanvas* cnv2 =new TCanvas(name+"_cnv2",name+"_cnv2",400,400);
  TGraphErrors* grp = new TGraphErrors();
  TH1F* Trn[5] ;//= new TH1F("","",80,15,23);
  grp->SetName(name);
  grp->GetXaxis()->SetTitle("Voltage");
  grp->GetYaxis()->SetTitle("transitTime (ns)");
  //cnv->Divide(2,2);
  cnv0->cd();
  TString str[5] = {"1000","1050","1100","1150","1200"};
  int col[5] = {kRed,kBlue,kOrange,kGreen,kBlack};
  bool flag = true; 
  for(int V=0;V<5;V++)
    {
      Trn[V] = new TH1F("","",80,15,23);
      TString trname = "TransitTime";
      TTree* w = (TTree*)_file0->Get("Tree"+trname+str[V]+"V_0");
      TTree* s = (TTree*)_file1->Get("SiPM01Tree"+trname+str[V]+"V_0");
      TTree* p = (TTree*)_file1->Get("AA6811Tree"+trname+str[V]+"V_0");
      float dtime[1024];
      float sp[1024];
      float sn[1024];
      float pp[1024];
      float pn[1024];
      float pt;
      float rt;
      TGraph* grpP0 = new TGraph();
      TGraph* grpP1 = new TGraph();
      TGraph* grpS0 = new TGraph();
      TGraph* grpS1 = new TGraph();
      TGraph* RvsP = new TGraph();
      
      int nEv = 3001;
      int nSl = 1024;
      
      w->SetBranchAddress("time",dtime);
      w->SetBranchAddress("wform0",pn);
      w->SetBranchAddress("wform1",pp);
      w->SetBranchAddress("wform2",sn);
      w->SetBranchAddress("wform3",sp);
      s->SetBranchAddress("risetime",&rt);
      p->SetBranchAddress("ptime",&pt);
      
      for(int i=0;i<nEv;i++)
	{
	  
	  w->GetEntry(i);
	  s->GetEntry(i);
	  p->GetEntry(i);
	  /*
	  if(rt>100)
	    {
	      for(int n=0;n<nSl;n++)
		{
		  grpS0->SetPoint(n+i*3001,dtime[n],sp[n]-sn[n]);
		}
	    }
	  else
	    {
	      for(int n=0;n<nSl;n++)
		{
		  grpS1->SetPoint(n+i*300,dtime[n],sp[n]-sn[n]);
		}
	    }
	  if(pt>121)
	    {
	      for(int n=0;n<nSl;n++)
		{
		  grpP0->SetPoint(n+i*3001,dtime[n],pp[n]-pn[n]);
		}
	    }
	  else
	    {
	      for(int n=0;n<nSl;n++)
		{
		  grpP1->SetPoint(n+i*300,dtime[n],pp[n]-pn[n]);
		}
	    }
      */
	  RvsP->SetPoint(i,rt,pt);
	  Trn[V]->Fill(pt-rt);
	}
      
      Trn[V]->GetXaxis()->SetTitle("Transit Time (ns)");
      Trn[V]->GetYaxis()->SetTitle("Events");
      Trn[V]->SetTitle(name+str[V]+"V");
      Trn[V]->SetLineColor(col[V]);
      Trn[V]->SetFillColor(col[V]);
      Trn[V]->SetFillStyle(3003);
      Trn[V]->SetName(name+str[V]+"V");
      Trn[V]->SetTitle(name+str[V]+"V");
      if(flag)
	{
	  std::cout<<"#################"<<std::endl;
	  Trn[V]->Draw();
	  cnv1->cd();
	  Trn[V]->Draw();
	  flag = false;
	}
      Trn[V]->Draw("SAME");
      Trn[V]->Write();
      grp->SetPoint(V,1000+V*50 +10,Trn[V]->GetMean());
      grp->SetPointError(V,0,Trn[V]->GetStdDev());
    }
  cnv2->cd();
  grp->Draw("AP");
  grp->Write();
  cnv0->Write();
  //file->Write();
  file->Close();
}














