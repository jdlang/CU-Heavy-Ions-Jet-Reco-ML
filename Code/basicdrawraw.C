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

int basicdrawraw()
{
  string names[4] = {"Pyth_10_90_Train_Trees-1.root","Pyth_10_90_Train_Trees4.root","Pyth_10_90_Train_Trees8.root","Pyth_40_60_Test_Trees-1.root"};
  string titles[4] = {"Uncorrected p_{T} Distribution 10-90 GeV Train Bias 0 All Jets","Uncorrected p_{T} Distribution 10-90 GeV Train Bias 4 All Jets","Uncorrected p_{T} Distribution 10-90 GeV Train Bias 8 All Jets","Uncorrected p_{T} Distribution 40-60 GeV Train Bias 0 All Jets"};
  gStyle->SetOptStat(0);
  TH1F* pthist[4];
  for(int i=0; i<4; i++)
    {
      TCanvas *c1 = new TCanvas("c1","c1");
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
      TTree* tree = (TTree*) file->Get("Otherfj_tree");
      tree->SetBranchAddress("jet_n1", &jetn);
      tree->SetBranchAddress("jet_pt", jetpt);
      for(int j=0; j<tree->GetEntries(); ++j)
	{
	  tree->GetEntry(j);
	  for(int k=0; k<jetn; ++k)
	    {
	      pthist[i]->Fill(jetpt[k]);
	    }
	}
      pthist[i]->GetYaxis()->SetRangeUser(1,1000000);
      pthist[i]->SetLineWidth(5);
      pthist[i]->GetYaxis()->SetTitle("Counts");
      pthist[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      pthist[i]->Draw();
      c1->SaveAs(("rawptdist" + to_string(i) + "alljets.pdf").c_str());
    }
  
  return 0;
}
