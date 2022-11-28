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

int basicdrawtruth()
{
  string names[4] = {"Pyth_10_90_Train_Trees-1.root","Pyth_10_90_Train_Trees4.root","Pyth_10_90_Train_Trees8.root","Pyth_40_60_Test_Trees-1.root"};
  string titles[4] = {"Truth p_{T} Distribution All Jets","Truth p_{T} Distribution 10-90 GeV Train Bias 4 All Jets","Truth p_{T} Distribution 10-90 GeV Train Bias 8 All Jets","Truth p_{T} Distribution 40-60 GeV Train Bias 0 All Jets"};
  string types[4] = {"10-90 GeV, 0 Bias","10-90 GeV, p_{T}^{4} Bias","10-90 GeV, p_{T}^{8} Bias","40-60 GeV, 0 Bias"};
  gStyle->SetOptStat(0);
  TH1F* pthist[4];
  TCanvas *c1 = new TCanvas("c1","c1");
  auto legend = new TLegend(0.75,0.75,1.0,1.0);
  for(int i=0; i<4; i++)
    {
      c1->SetLogy();
      if(i<3)
	{
	  pthist[i] = new TH1F("temp",titles[i].c_str(), 40, 10, 90);
	}
      else
	{
	  pthist[i] = new TH1F("temp",titles[i].c_str(),10, 40,60);
	}
      string name("./compiled/Data/" + names[i]);
      TFile* file = TFile::Open(name.c_str());
      int jetn;
      float jetpt[100];
      TTree* tree = (TTree*) file->Get("FastJet_Tree");
      tree->SetBranchAddress("jet_n", &jetn);
      tree->SetBranchAddress("jet_pt", jetpt);
      legend->AddEntry(pthist[i],types[i].c_str(),"p");
      for(int j=0; j<tree->GetEntries(); ++j)
	{
	  tree->GetEntry(j);
	  for(int k=0; k<jetn; ++k)
	    {
	      
	      pthist[i]->SetMarkerColor(i+1);
	      pthist[i]->SetMarkerStyle(i+39);
	      if(i==3)
		{
		  pthist[i]->Fill(jetpt[k],0.0045);
		}
	      else
		{
		  pthist[i]->Fill(jetpt[k]);
		}
	    }
	}
      pthist[i]->GetYaxis()->SetRangeUser(1,1000000);
      pthist[i]->SetLineWidth(5);
      pthist[i]->GetYaxis()->SetTitle("Counts");
      pthist[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      pthist[i]->Draw("P0 HIST SAME");
    }
  legend->Draw();
  c1->SaveAs("Truthptdistalljets.pdf");
  return 0;
}
