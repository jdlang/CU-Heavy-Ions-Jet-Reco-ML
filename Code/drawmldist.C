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

int drawmldist()
{
  string names[5] = {"ML_Prep_10_90_Train-1.root","ML_Prep_10_90_Train2.root","ML_Prep_10_90_Train4.root","ML_Prep_10_90_Train8.root","ML_Prep_40_60_Test-1.root"};
  string titles[5] = {"Truth p_{T} Distribution 10-90 GeV Train Bias 0 All Jets","Truth p_{T} Distribution 10-90 GeV Train Bias 2 All Jets","Truth p_{T} Distribution 10-90 GeV Train Bias 4 All Jets","Truth p_{T} Distribution 10-90 GeV Train Bias 8 All Jets","Truth p_{T} Distribution 40-60 GeV Train Bias 0 All Jets"};
  string types[5] = {"10-90 GeV, 0 Bias","10-90 GeV, p_{T}^{2} Bias","10-90 GeV, p_{T}^{4} Bias","10-90 GeV, p_{T}^{8} Bias","40-60 GeV, 0 Bias"};
  gStyle->SetOptStat(0);
  TH1F* pthist[5];
  auto legend = new TLegend(0.75,0.65,1.0,0.9);
  TCanvas *c1 = new TCanvas("c1","c1");
  for(int i=0; i<5; i++)
    {
      string name = names[i];
      c1->SetLogy();
      if(i<4)
	{
	  pthist[i] = new TH1F("temp","Truth p_{T} Distributions", 40, 10, 90);
	}
      else
	{
	  pthist[i] = new TH1F("temp","Truth p_{T} Distributions",10, 40,60);
	}
      TFile* file = TFile::Open(name.c_str());
      int jetn;
      float jetpt[100];
      TTree* tree = (TTree*) file->Get("Tree_Tree");
      tree->SetBranchAddress("jet_n", &jetn);
      tree->SetBranchAddress("jet_pt_true_pythia", jetpt);
      legend->AddEntry(pthist[i],types[i].c_str(),"p");
      for(int j=0; j<tree->GetEntries(); ++j)
	{
	  tree->GetEntry(j);
	  for(int k=0; k<jetn; ++k)
	    {
	      if(i<4)
		{
		  pthist[i]->SetMarkerColor(i+1);
		}
	      else
		{
		  pthist[i]->SetMarkerColor(i+2);
		}
	      pthist[i]->SetMarkerStyle(i+39);
	      if(i<4)
		{
		  if(jetpt[k] > 10 && jetpt[k] < 90) pthist[i]->Fill(jetpt[k]);
		}
	      else
		{
		  if(jetpt[k] > 40 && jetpt[k] < 60) pthist[i]->Fill(jetpt[k]);
		}
	    }
	}
      pthist[i]->GetYaxis()->SetRangeUser(1,1000000);
      pthist[i]->SetLineWidth(5);
      pthist[i]->GetYaxis()->SetTitle("Counts");
      pthist[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
    }
  float integral[5];
  integral[0] = pthist[0]->Integral(15,25);
  integral[1] = pthist[1]->Integral(15,25);
  integral[2] = pthist[2]->Integral(15,25);
  integral[3] = pthist[3]->Integral(15,25);
  integral[4] = pthist[4]->Integral(0,10);
  cout << integral[0] << endl;
  cout << integral[4] << endl;
  for(int i=0; i<5; ++i)
    {
      pthist[i]->Scale(integral[0]/integral[i]);
      pthist[i]->Draw("P0 HIST SAME");
    }
  legend->Draw();
  c1->SaveAs("Truthptdistalljetsmlprepped.pdf");
  return 0;
}
