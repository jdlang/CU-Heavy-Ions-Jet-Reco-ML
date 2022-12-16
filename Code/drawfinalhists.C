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

int drawhists(int display = 0)
{
  string names[5] = {"Area Correction","Linear Regressor","MLP Regressor","Random Forest","ElasticNet Regressor"};
  float mean[3][3][3][5];
  float rmse[3][3][3][5];
  string feats[3] = {"1","3","12"};
  int powers[3] = {-1,4,8};
  string train = "10_90";
  gStyle->SetOptStat(0);
  TH1F* meanhist[3][3][5];
  TH1F* rmsehist[3][3][5];
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<5; j++)
	{
	  for(int k=0; k<3; ++k)
	    {
	      meanhist[i][k][j] = new TH1F((to_string(i)+to_string(j)).c_str(),(to_string(i)+to_string(j)).c_str(), 10, -1.5, 8.5);
	      rmsehist[i][k][j] = new TH1F((to_string(i)+to_string(j)).c_str(),(to_string(i)+to_string(j)).c_str(), 10, -1.5, 8.5);
	    }
	}
    }
  for(int l=0; l<3; l++)
    {
      for(int j=0; j<3; j++)
	{
	  for(int k=0; k<3; ++k)
	    {
	      string feat = feats[j];
	      TFile *f1 = TFile::Open(("./compiled/Data/results_" + train+ "train_ptbias" + to_string(powers[l])+ "test_ptbias" + to_string(powers[k]) +"corrected" + feat + "Features.root").c_str());
	      TFile *f3 = TFile::Open(("./compiled/Data/results_" + train+ "train_ptbias" + to_string(powers[l])+ "test_ptbias" + to_string(powers[k]) +"MLPRegressor" + feat + "Features.root").c_str());
	      TFile *f2 = TFile::Open(("./compiled/Data/results_" + train+ "train_ptbias" + to_string(powers[l])+ "test_ptbias" + to_string(powers[k]) +"LinearRegression" + feat + "Features.root").c_str());
	      TFile *f4 = TFile::Open(("./compiled/Data/results_" + train+ "train_ptbias" + to_string(powers[l])+ "test_ptbias" + to_string(powers[k]) +"RandomForestRegressor" + feat + "Features.root").c_str());
	      TFile *f5 = TFile::Open(("./compiled/Data/results_" + train+ "train_ptbias" + to_string(powers[l])+ "test_ptbias" + to_string(powers[k]) +"ElasticNet" + feat + "Features.root").c_str());
	      int ntypes = 5;
	      TH1D *h[ntypes];
	      TCanvas *c1 = new TCanvas("c1","c1");
	      auto legend = new TLegend(0.1,0.7,0.4,0.9);
	      h[0] = (TH1D*) f1->Get("hist");
	      h[1] = (TH1D*) f2->Get("hist");
	      h[2] = (TH1D*) f3->Get("hist");
	      h[3] = (TH1D*) f4->Get("hist");
	      h[4] = (TH1D*) f5->Get("hist");
	      
	      for(int i = 0; i<ntypes; i++)
		{
		  mean[l][k][j][i] = h[i]->GetMean();
		  rmse[l][k][j][i] = h[i]->GetStdDev();
		  cout << mean[l][k][j][i] << " " << rmse[l][k][j][i] << endl;
		  //if(abs(mean[l][k][j][i]) > 25) cout << mean[l][k][j][i] << endl; 
		  string title(feat + "Features, 10-90 GeV p_{T} bias " + to_string(powers[l]) + " Train, 40-60 GeV p_{T} bias " + to_string(powers[k]) + " Test");
		  h[i]->SetTitle(title.c_str());
		  h[i]->GetYaxis()->SetRangeUser(0,25000);
		  if(i == 4)
		    {
		      h[i]->SetMarkerColor(6);
		    }
		  else
		    {
		      h[i]->SetMarkerColor(i+1);
		    }
		  h[i]->SetMarkerStyle(i+39);
		  legend->AddEntry(h[i],(names[i] + " #\sigma = " + to_string(h[i]->GetStdDev())).c_str(),"p");
		  h[i]->Draw("P0 HIST SAME");
		}
	      legend->Draw();
	      c1->SaveAs((feat + "Features" + train + "train_ptbias" + to_string(powers[l])+ "test_ptbias" + to_string(powers[k]) + ".pdf").c_str());
	    }
	}
    }
  int nons[7] = {0,1,2,3,5,6,7};
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	{
	  TCanvas *c2 = new TCanvas("c2","c2");
	  TCanvas *c3 = new TCanvas("c3","c3");
	  auto legend2= new TLegend(0.35,0.7,0.65,0.9);
	  auto legend3= new TLegend(0.35,0.7,0.65,0.9);
	  for(int k=0; k<3; k++)
	    {
	      for(int l=0; l<5; l++)
		{
		  if(k==0)
		    {
		      legend2->AddEntry(meanhist[j][k][l], names[l].c_str(),"p");
		      legend3->AddEntry(rmsehist[j][k][l], names[l].c_str(),"p");
		    }
		  for(int m=0; m<7; m++)
		    {
		      meanhist[j][k][l]->Fill(nons[m], 10000);
		      rmsehist[j][k][l]->Fill(nons[m], 10000);
		    }
		  meanhist[j][k][l]->Fill(powers[k], mean[i][k][j][l]);
		  meanhist[j][k][l]->SetMarkerStyle(l+39);
		  if(l==4)
		    {
		      meanhist[j][k][l]->SetMarkerColor(6);
		    }
		  else
		    {
		      meanhist[j][k][l]->SetMarkerColor(l+1);
		    }
		  meanhist[j][k][l]->GetYaxis()->SetRangeUser(-30,30);
		  meanhist[j][k][l]->GetXaxis()->SetTitle("Test p_{T} Bias Power");
		  meanhist[j][k][l]->GetYaxis()->SetTitle("Mean of p_{T true} - p_{T corrected}");
		  meanhist[j][k][l]->SetTitle(("Mean as a Function of Test p_{T} Bias Power, " + feats[j] + " Features, " + powers[i] + " Train p_{T} Bias Power"));
		  rmsehist[j][k][l]->Fill(powers[k], rmse[i][k][j][l]);
		  rmsehist[j][k][l]->SetMarkerStyle(l+39);
		  if(l==4)
		    {
		      rmsehist[j][k][l]->SetMarkerColor(6);
		    }
		  else
		    {
		      rmsehist[j][k][l]->SetMarkerColor(l+1);
		    }
		  rmsehist[j][k][l]->GetYaxis()->SetRangeUser(0,15);
		  rmsehist[j][k][l]->GetXaxis()->SetTitle("Test p_{T} Bias Power");
		  rmsehist[j][k][l]->GetYaxis()->SetTitle("#\sigma of p_{T true} - p_{T corrected}");
		  rmsehist[j][k][l]->SetTitle(("RMSE as a Function of Test p_{T} Bias Power, " + feats[j] + " Features, " + powers[i] + " Train p_{T} Bias Power"));
		}
	      for(int l=0; l<5; ++l)
		{
		  c2->cd();
		  meanhist[j][k][l]->Draw("P0 HIST SAME");
		  c3->cd();
		  rmsehist[j][k][l]->Draw("P0 HIST SAME");
		}
	    }
	  c3->cd();
	  legend3->Draw();
	  c2->cd();
	  legend2->Draw();
	  c2->SaveAs(("meantrend_train" + to_string(powers[i]) + "_" + feats[j] + "features.pdf").c_str());
	  c3->SaveAs(("rmsetrend_train" + to_string(powers[i]) + "_" + feats[j] + "features.pdf").c_str());
	  c2->Clear();
	  c3->Clear();
	}
    }
  
  return 0;
}
