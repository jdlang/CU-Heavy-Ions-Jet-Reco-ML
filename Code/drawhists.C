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

struct Jet
{
  float pt;
  float et;
  float ph;
  int n;
};

struct Particle
{
  float pt;
  float et;
  float ph;
  int jt;
  int id;
};

struct Result
{
  float pt;
  float et;
  float ph;
  float id;
};

bool comppt(const Jet &a, const Jet &b)
{
  return a.pt > b.pt;
}

float greaterval(float a, float b)
{
  if(a > b)
    {
      return a;
    }
  return b;
}

float lesserval(float a, float b)
{
  if(a < b)
    {
      return a;
    }
  return b;
}

float deltaR(Jet jet1, Jet jet2)
{
  float det = jet1.et - jet2.et;
  float dph1 = jet1.ph - jet2.ph;
  float dph2 = greaterval(jet1.ph, jet2.ph) - 2*M_PI - lesserval(jet1.ph, jet2.ph);
  float dph = lesserval(abs(dph1), abs(dph2));
  float dr = sqrt(det*det + dph*dph);
  return dr;
}

void swapjets(vector<Jet> &b, int index1, int index2)
{
  Jet temp;
  temp.pt = b[index1].pt;
  temp.et = b[index1].et;
  temp.ph = b[index1].ph;
  temp.n = b[index1].n;
  b[index1].pt = b[index2].pt;
  b[index1].et = b[index2].et;
  b[index1].ph = b[index2].ph;
  b[index1].n = b[index2].n;
  b[index2].pt = temp.pt;
  b[index2].et = temp.et;
  b[index2].ph = temp.ph;
  b[index2].n = temp.n;
}

float totalr2(vector<Jet> a, vector<Jet> b)
{
  float total = 0;
  for(int i=0; i<lesserval(a.size(), b.size()); i++)
    {
      total += deltaR(a[i], b[i]);
    }
  return total*total;
}
  

void matchjets(vector<Jet> &truth, vector<Jet> &reco)
{
  float tdR2;
  float ndR2;
  int mindex1;
  int mindex2;
  int counter = 0;
  float temp = 1;
  float r = 0.97;
  sort(truth.begin(), truth.end(), comppt);
  sort(reco.begin(), reco.end(), comppt);
  vector<Jet> recotemp;
  do
    {
      recotemp.clear();
      for(int i=0; i<reco.size(); i++)
	{
	  recotemp.push_back(reco[i]);
	}
      mindex1 = rand()%recotemp.size();
      mindex2 = rand()%recotemp.size();
      if(mindex1 >= recotemp.size() || mindex2 >= recotemp.size()) continue;
      float chance = ((float)rand())/RAND_MAX;
      if(chance < temp || counter < 100)
	{
	  swapjets(recotemp, mindex1, mindex2);
	  counter++;
	  temp *= r;
	  continue;
	}
      temp *= r;
      tdR2 = totalr2(recotemp, truth);
      swapjets(recotemp, mindex1, mindex2);
      ndR2 = totalr2(recotemp, truth);
      if(ndR2 < tdR2)
	{
	  swapjets(reco, mindex1, mindex2);
	}

    } while(temp > 0.00000001);
}

float eta_from_index(int index)
{
  return (0.0245*floor(index / 256) - 1.152);
}

float phi_from_index(int index)
{
  return (0.024*(index % 256) - M_PI);
}

float pt_from_pxpy(float px, float py)
{
  return sqrt(px*px+py*py);
}

float eta_from_p(float px, float py, float pz)
{
  float magp = sqrt(px*px+py*py+pz*pz);
  return atanh(pz/magp);
}

int drawhists(int display = 0)
{
  
  string feats[3] = {"1","3","12"};
  string trains[2] = {"10_90","30_70"};
  gStyle->SetOptStat(0);
  for(int k=0; k<2; k++)
    {
      string train = trains[k];
      for(int j=0; j<3; ++j)
	{
	  string feat = feats[j];
	  TFile *f4 = TFile::Open(("./compiled/Data/results_" + train +"corrected" + feat + "Features.root").c_str());
	  TFile *f3 = TFile::Open(("./compiled/Data/results_" + train +"MLPRegressor" + feat + "Features.root").c_str());
	  TFile *f2 = TFile::Open(("./compiled/Data/results_" + train +"LinearRegression" + feat + "Features.root").c_str());
	  TFile *f1 = TFile::Open(("./compiled/Data/results_" + train +"RandomForestRegressor" + feat + "Features.root").c_str());
	  
	  int ntypes = 4;
	  TH1D *h[ntypes];
	  TCanvas *c1 = new TCanvas("c1","c1");
	  auto legend = new TLegend(0.1,0.7,0.4,0.9);
	  h[0] = (TH1D*) f1->Get("hist");
	  h[1] = (TH1D*) f2->Get("hist");
	  h[2] = (TH1D*) f3->Get("hist");
	  h[3] = (TH1D*) f4->Get("hist");
	  
	  string names[4] = {"Random Forest","Linear Regressor","MLP Regressor","Area Correction"};
	  for(int i = 0; i<ntypes; i++)
	    {
	      string title(h[i]->GetTitle());
	      int pos = title.find("Falling");
	      pos = title.find("Falling", pos+1);
	      if(k==0) title.replace(pos,7,"Flat");
	      pos = title.find("Test ");
	      title.replace(pos,5,"");
	      h[i]->SetTitle(title.c_str());
	      h[i]->GetYaxis()->SetRangeUser(0,1000);
	      h[i]->SetLineColor(i+6);
	      legend->AddEntry(h[i],(names[i] + " #\sigma = " + to_string(h[i]->GetStdDev())).c_str(),"l");
	      h[i]->Draw("SAME");
	    }
	  legend->Draw();
	  c1->SaveAs((feat + "Features" + train + ".pdf").c_str());
	}
    }
  return 0;
}
