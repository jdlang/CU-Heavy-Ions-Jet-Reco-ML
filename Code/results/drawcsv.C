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

int drawcsv()
{
  gStyle->SetOptStat(0);
  int powers[1] = {-1};//,2,4,8};
  int realpow[1] = {0};
  int ptrange[4][2] = {18,22,28,32,38,42,48,52};
  string types[2] = {"corrected","LinearRegression"};
  int features[3] = {1,3,12};
  for(int i=0; i<sizeof(powers)/sizeof(powers[0]); ++i)
    {
      for(int l=0; l<sizeof(features)/sizeof(features[0]); ++l)
	{
	  for(int j=0; j<sizeof(ptrange)/sizeof(ptrange[0]); ++j)
	    {
	      TCanvas *c1 = new TCanvas("c1","c1");
	      auto legend = new TLegend(0.5,0.8,0.9,0.9);
	      TH1F * hists[2];
	      for(int k=0; k<2; ++k)
		{
		  string name = "results_10_90_train_ptbias" + to_string(powers[i]) + "_test_ptbias8_" + types[k] + "_" + to_string(features[l]) + "Features_set" + to_string(ptrange[j][0]) + "_" + to_string(ptrange[j][1]);
		  string title = "Train p_{T}^{" + to_string(realpow[i]) + "} 10-90 GeV, Test p_{T}^{8} " + to_string(ptrange[j][0]) + "-" + to_string(ptrange[j][1]) + " GeV, " + to_string(features[l]) + " Feature(s)";
		  ifstream infile;
		  infile.open(name + ".csv");
		  hists[k] = new TH1F("",title.c_str(), 150,-50,25);
		  string line;
		  while(getline(infile, line))
		    {
		      float value = stof(line);
		      hists[k]->Fill(value);
		    }
		  hists[k]->GetYaxis()->SetRangeUser(0,5000);
		  hists[k]->SetMarkerStyle(k+39);
		  hists[k]->SetMarkerColor(k+1);
		  hists[k]->GetYaxis()->SetTitle("Counts");
		  hists[k]->GetXaxis()->SetTitle("p_{T true}-p_{T pred}");
		  float mean = hists[k]->GetMean();
		  float sigma = hists[k]->GetStdDev();
		  stringstream sigstr;
		  sigstr << std::fixed << std::setprecision(2) << sigma;
		  stringstream menstr;
		  menstr << std::fixed << std::setprecision(2) << mean;
		  legend->AddEntry(hists[k], (types[k] + " #\sigma=" + sigstr.str() + " #\mu=" + menstr.str()).c_str(), "p");
		  hists[k]->Draw("P0 SAME HIST");
		}
	      legend->Draw();
	      c1->SaveAs(("results_10_90_train_ptbias" + to_string(powers[i]) + "_test_ptbias8_" + to_string(features[l]) + "Features_set" + to_string(ptrange[j][0]) + "_" + to_string(ptrange[j][1]) + ".pdf").c_str());//(name + ".pdf").c_str());
	    }
	}
    }
  return 0;
}
