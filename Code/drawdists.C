#include "TApplication.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TRandom2.h"
#include "TMath.h"
#include <utility>  // for std::pair
#include <cstdio>
#include <iostream>
#include "TGraph.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include "stdlib.h"
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include <TMinuit.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TString.h>
#include <TAxis.h>
#include <TLine.h>
#include <TFile.h>
#include <algorithm>

using namespace std;

int drawdists()
{
  TFile *f0 = new TFile("test.root","new");
  TCanvas *c1 = new TCanvas("c1","c1");
  TCanvas *c2 = new TCanvas("c2","c2");
  int powers[3] = {-1,4,8};
  gStyle->SetOptStat(0);
  for(int l=0; l<3; l++)
    {
      TFile *f1 = new TFile(("./compiled/Data/Pyth_10_90_Train_Trees"+to_string(powers[l])+".root").c_str(), "read");
      TFile *f2 = new TFile(("./compiled/Data/Pyth_40_60_Test_Trees"+to_string(powers[l])+".root").c_str(), "read");
      TTree *t1 = f1->Get<TTree>("FastJet_Tree");
      TTree *t2 = f2->Get<TTree>("FastJet_Tree");
      float pt1[20];
      float pt2[20];
      int jetn1;
      int jetn2;
      t1->SetBranchAddress("jet_pt", pt1);
      t2->SetBranchAddress("jet_pt", pt2);
      t1->SetBranchAddress("jet_n", &jetn1);
      t2->SetBranchAddress("jet_n", &jetn2);
      TH1D *h[2];
      h[0] = new TH1D("hist1","hist1",4000,40,60);
      h[1] = new TH1D("hist2","hist2",4000,40,60);
      for(int i=0; i<t1->GetEntries(); ++i)
	{
	  t1->GetEntry(i);
	  for(int j=0; j<jetn1; j++)
	    {
	      if(pt1[j] > 10 && pt1[j] < 90) h[0]->Fill(pt1[j]);
	    }
	}
      for(int i=0; i<jetn2; ++i)
	{
	  t2->GetEntry(i);
	  for(int j=0; j<jetn2; j++)
	    {
	      if(pt2[j] > 40 && pt2[j] < 60 && pt2[j] != 52) h[1]->Fill(pt2[j]);
	    }
	}
      c1->cd();
      string title("10-90 GeV Train Truth Jet p_{T} Distribution; p_{T} Bias = " + to_string(powers[l]));
      h[0]->SetTitle(title.c_str());
      //h[0]->GetYaxis()->SetRangeUser(1,100000);
      c1->SetLogy();
      h[0]->Draw();
      c1->SaveAs(("10_90_bias" + to_string(powers[l]) + ".pdf").c_str());
      c2->cd();
      string title2("40-60 GeV Test Truth Jet p_{T} Distribution; p_{T} Bias = " + to_string(powers[l]));
      h[1]->SetTitle(title2.c_str());
      //h[1]->GetYaxis()->SetRangeUser(1,100000);
      c2->SetLogy();
      h[1]->Draw("");
      c2->SaveAs(("40_60_bias" + to_string(powers[l]) + ".pdf").c_str());
      c1->Clear();
      c2->Clear();
      f0->cd();
      h[0]->Write("test.root");
      h[1]->Write("test.root");
    }
  return 0;
}
