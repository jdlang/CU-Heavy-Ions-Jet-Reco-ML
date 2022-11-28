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

int drawhists2(int display = 0)
{
  string ranges[4] = {"40-60 GeV","40-60 GeV","48-52 GeV","48-52 GeV"};
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
	  for(int k=0; k<1; ++k)
	    {
	      string feat = feats[j];
	      TFile *f1 = TFile::Open(("./compiled/Data/results_" + train+ "train_ptbias" + to_string(powers[l])+ "test_ptbias" + to_string(powers[k]) +"corrected" + feat + "Features_setnumber" + to_string(display) + ".root").c_str());
	      TFile *f3 = TFile::Open(("./compiled/Data/results_" + train+ "train_ptbias" + to_string(powers[l])+ "test_ptbias" + to_string(powers[k]) +"MLPRegressor" + feat + "Features_setnumber" + to_string(display) + ".root").c_str());
	      TFile *f2 = TFile::Open(("./compiled/Data/results_" + train+ "train_ptbias" + to_string(powers[l])+ "test_ptbias" + to_string(powers[k]) +"LinearRegression" + feat + "Features_setnumber" + to_string(display) + ".root").c_str());
	      TFile *f4 = TFile::Open(("./compiled/Data/results_" + train+ "train_ptbias" + to_string(powers[l])+ "test_ptbias" + to_string(powers[k]) +"RandomForestRegressor" + feat + "Features_setnumber" + to_string(display) + ".root").c_str());
	      TFile *f5 = TFile::Open(("./compiled/Data/results_" + train+ "train_ptbias" + to_string(powers[l])+ "test_ptbias" + to_string(powers[k]) +"ElasticNet" + feat + "Features_setnumber" + to_string(display) + ".root").c_str());
	      
	      int ntypes = 5;
	      TH1D *h[ntypes];
	      TCanvas *c1 = new TCanvas("c1","c1");
	      auto legend = new TLegend(0.6,0.7,0.9,0.9);
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
		  string title(feat + " Feature(s), Train 10-90 GeV p_{T} bias " + to_string(powers[l]) + ", Test " + ranges[display] + " GeV p_{T} bias " + to_string(powers[k]));
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
	      c1->SaveAs((feat + "Features" + train + "train_ptbias_" + to_string(powers[l])+ "_test_ptbias_" + to_string(powers[k]) + "setnumber_" + to_string(display) + ".pdf").c_str());
	    }
	}
    }
  int nons[7] = {0,1,2,3,5,6,7};
  for(int k=0; k<1; k++)
    {
      for(int j=0; j<3; j++)
	{
	  TCanvas *c2 = new TCanvas("c2","c2");
	  TCanvas *c3 = new TCanvas("c3","c3");
	  auto legend2= new TLegend(0.2,0.7,0.5,0.9);
	  auto legend3= new TLegend(0.2,0.7,0.5,0.9);
	  for(int i=0; i<3; i++)
	    {
	      for(int l=0; l<5; l++)
		{
		  if(i==0)
		    {
		      legend2->AddEntry(meanhist[k][j][l], names[l].c_str(),"p");
		      legend3->AddEntry(rmsehist[k][j][l], names[l].c_str(),"p");
		    }
		  for(int m=0; m<7; m++)
		    {
		      meanhist[i][j][l]->Fill(nons[m], 10000);
		      rmsehist[i][j][l]->Fill(nons[m], 10000);
		    }
		  meanhist[k][j][l]->Fill(powers[i], mean[i][k][j][l]);
		  meanhist[k][j][l]->SetMarkerStyle(l+39);
		  if(l==4)
		    {
		      meanhist[k][j][l]->SetMarkerColor(6);
		    }
		  else
		    {
		      meanhist[k][j][l]->SetMarkerColor(l+1);
		    }
		  meanhist[k][j][l]->GetYaxis()->SetRangeUser(-30,10);
		  meanhist[k][j][l]->GetXaxis()->SetTitle("Train p_{T} Bias Power");
		  meanhist[k][j][l]->GetYaxis()->SetTitle("Mean of p_{T true} - p_{T corrected}");
		  string title1(("Mean as a Function of Train p_{T} Bias Power (10-90 GeV), " + feats[j] + " Feature(s), " + " Test p_{T} Bias Power " + powers[k] + " (" + ranges[display].c_str() + ")"));
		  meanhist[k][j][l]->SetTitle(title1.c_str());
		  rmsehist[k][j][l]->Fill(powers[i], rmse[i][k][j][l]);
		  rmsehist[k][j][l]->SetMarkerStyle(l+39);
		  if(l==4)
		    {
		      rmsehist[k][j][l]->SetMarkerColor(6);
		    }
		  else
		    {
		      rmsehist[k][j][l]->SetMarkerColor(l+1);
		    }
		  rmsehist[k][j][l]->GetYaxis()->SetRangeUser(0,15);
		  rmsehist[k][j][l]->GetXaxis()->SetTitle("Train p_{T} Bias Power");
		  rmsehist[k][j][l]->GetYaxis()->SetTitle("#\sigma of p_{T true} - p_{T corrected}");
		  string title2(("RMSE as a Function of Train p_{T} Bias Power (10-90 GeV), " + feats[j] + " Feature(s), " + " Test p_{T} Bias Power " + powers[k] + " (" + ranges[display].c_str() + ")"));
		  rmsehist[k][j][l]->SetTitle(title2.c_str());
		}
	      for(int l=0; l<5; ++l)
		{
		  c2->cd();
		  meanhist[i][j][l]->Draw("P0 HIST SAME");
		  c3->cd();
		  rmsehist[i][j][l]->Draw("P0 HIST SAME");
		}
	    }
	  c3->cd();
	  legend3->Draw();
	  c2->cd();
	  legend2->Draw();
	  c2->SaveAs(("meantrend_test" + to_string(powers[k]) + "_" + feats[j] + "features_" + "setnumber_" + to_string(display) + ".pdf").c_str());
	  c3->SaveAs(("rmsetrend_test" + to_string(powers[k]) + "_" + feats[j] + "features_" + "setnumber_" + to_string(display) + ".pdf").c_str());
	  c2->Clear();
	  c3->Clear();
	}
    }
  
  return 0;
}
