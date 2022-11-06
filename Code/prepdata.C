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
int prepdata()
{
  int jetn;
  float features[100][10];
  float constnpt[100][10];
  vector<string> files{"ML_Prep_20_40_Test", "ML_Prep_40_60_Test", "ML_Prep_60_80_Test"};
  vector<string> trees{"Tree_20_40_Test", "Tree_40_60_Test", "Tree_60_80_Test"};
  vector<string> exten{".txt", ".root"};
  string dirname = "./compiled/Data/";

  TFile *file;
  TTree *tree;

  for(int i=0; i<files.size(); ++i)
    {

      for(int j=0; j<sizeof(features)/sizeof(features[0]); ++j)
	{
	  for(int k=0; k<sizeof(features[j])/sizeof(features[0][0]); ++k)
	    {
	      features[j][k] = -10000;
	    }
	}
      for(int j=0; j<sizeof(constnpt)/sizeof(constnpt[0]); ++j)
	{
	  for(int k=0; k<sizeof(constnpt[j])/sizeof(constnpt[0][0]); ++k)
	    {
	      constnpt[j][k] = -10000;
	    }
	}

      string concat = dirname + files[i] + exten[1];
      file = TFile::Open(concat.c_str());

      tree = file->Get<TTree>(trees[i].c_str());

      tree->SetBranchAddress("jet_n", &jetn);
      tree->SetBranchAddress("jet_pt_raw",features[0]);
      tree->SetBranchAddress("jet_pt_corr", features[1]);
      tree->SetBranchAddress("jet_mass", features[2]);
      tree->SetBranchAddress("jet_area", features[3]);
      tree->SetBranchAddress("jet_area_err",features[4]);
      //tree->SetBranchAddress("jet_const_n");
      tree->SetBranchAddress("const_pt_mean", features[5]);
      tree->SetBranchAddress("const_pt_median", features[6]);
      tree->SetBranchAddress("jet_y",features[7]);
      tree->SetBranchAddress("jet_phi", features[8]);
      tree->SetBranchAddress("jet_rho", features[9]);
      tree->SetBranchAddress("const_1_pt", &constnpt[0]);
      tree->SetBranchAddress("const_2_pt", &constnpt[1]);
      tree->SetBranchAddress("const_3_pt", &constnpt[2]);
      tree->SetBranchAddress("const_4_pt", &constnpt[3]);
      tree->SetBranchAddress("const_5_pt", &constnpt[4]);
      tree->SetBranchAddress("const_6_pt", &constnpt[5]);
      tree->SetBranchAddress("const_7_pt", &constnpt[6]);
      tree->SetBranchAddress("const_8_pt", &constnpt[7]);
      tree->SetBranchAddress("const_9_pt", &constnpt[8]);
      tree->SetBranchAddress("const_10_pt", &constnpt[9]);
      ofstream otherfile;
      concat = files[i] + exten[0];
      otherfile.open(concat);
      for(int j=0; j<tree->GetEntries(); ++i)
	{
	  tree->GetEntry(j);
	  otherfile << jetn << endl;
	  for(int l=0; l<sizeof(features)/sizeof(features[0]); ++l)
	    {
	      int m = 0;
	      while(features[l][m] != -10000)
		{
		  otherfile << features[l][m] << " ";
		  ++m;
		}
	      otherfile << endl;
	    }
	  otherfile << endl;
	  for(int l=0; l<sizeof(constnpt)/sizeof(constnpt[0]); ++l)
	    {
	      int m=0;
	      while(constnpt[l][m] != -10000)
		{
		  otherfile << constnpt[l][m] << " ";
		  ++m;
		}
	      otherfile << endl;
	    }
	  otherfile << endl << endl;
	}
      otherfile.close();
    }
  return 0;
}
