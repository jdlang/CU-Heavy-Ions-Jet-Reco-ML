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

int drawcsvraw()
{
  gStyle->SetOptStat(0);
  int powers[1] = {-1};//,2,4,8};
  int realpow[1] = {0};
  int ptrange[4][2] = {18,22,28,32,38,42,48,52};
  string types[2] = {"corrected","LinearRegression"};
  int features[1] = {1};//,3,12};
  for(int i=0; i<sizeof(powers)/sizeof(powers[0]); ++i)
    {
      for(int l=0; l<sizeof(features)/sizeof(features[0]); ++l)
	{
	  TCanvas *c1 = new TCanvas("c1","c1");
	  auto legend = new TLegend(0.5,0.7,0.9,0.9);
	  TH1F * hists[4];
	  for(int j=0; j<sizeof(ptrange)/sizeof(ptrange[0]); ++j)
	    {

		  string name = "../rawjets.csv";
		  string title = "Different Truth p_{T} Regions over Raw p_{T} (Test Sample)";
		  ifstream infile;
		  ifstream infile2;
		  infile2.open("../trujets.csv");
		  infile.open(name);
		  hists[j] = new TH1F("",title.c_str(), 60,30,150);
		  string line;
		  string otherline;
		  while(getline(infile, line))
		    {
		      getline(infile2, otherline);
		      float value = stof(otherline);
		      float rawval = stof(line);
		      if((value < ptrange[j][1]) && (value > ptrange[j][0])) hists[j]->Fill(rawval);
		    }
		  hists[j]->GetYaxis()->SetRangeUser(0,2000);
		  hists[j]->SetMarkerStyle(j+39);
		  hists[j]->SetMarkerColor(j+1);
		  hists[j]->GetYaxis()->SetTitle("Counts");
		  hists[j]->GetXaxis()->SetTitle("Raw Jet p_{T}");
		  float mean = hists[j]->GetMean();
		  float sigma = hists[j]->GetStdDev();
		  stringstream sigstr;
		  sigstr << std::fixed << std::setprecision(2) << sigma;
		  stringstream menstr;
		  menstr << std::fixed << std::setprecision(2) << mean;
		  legend->AddEntry(hists[j], (to_string(ptrange[j][0]) + "-" + to_string(ptrange[j][1]) + " " + " #\sigma=" + sigstr.str() + " #\mu=" + menstr.str()).c_str(), "p");
		  hists[j]->Draw("P0 SAME HIST");
		}
	      legend->Draw();
	      c1->SaveAs("rawjetranges.pdf");//(name + ".pdf").c_str());
	}
    }
  return 0;
}
