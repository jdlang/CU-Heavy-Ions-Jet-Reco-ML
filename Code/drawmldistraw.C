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

int drawmldistraw()
{
  string names[1] = {"ML_Prep_40_60_Test8.root"};
  gStyle->SetOptStat(0);
  TH1F* pthist[5];
  auto legend = new TLegend(0.75,0.75,1.0,1.0);
  TCanvas *c1 = new TCanvas("c1","c1");
  for(int i=0; i<1; i++)
    {
      string name = names[i];
      TFile* file = TFile::Open(name.c_str());
      int jetn;
      float jetpt[100];
      float jetpttrue[100];
      TTree* tree = (TTree*) file->Get("Tree_Tree");
      ofstream outfile;
      ofstream outfile2;
      outfile2.open("trujets.csv");
      outfile.open("rawjets.csv");
      tree->SetBranchAddress("jet_n", &jetn);
      tree->SetBranchAddress("jet_pt_raw", jetpt);
      tree->SetBranchAddress("jet_pt_true_pythia", jetpttrue);
      for(int j=0; j<tree->GetEntries(); ++j)
	{
	  tree->GetEntry(j);
	  for(int k=0; k<jetn; ++k)
	    {
	      outfile << jetpt[k] << endl;
	      outfile2 << jetpttrue[k] << endl;
	    }
	}
      outfile.close();
    }
  return 0;
}
