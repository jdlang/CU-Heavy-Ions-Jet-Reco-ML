#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include <cmath>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <fstream>
#include <map>
#include "_header.h"
using namespace Project_Constants;

const int label_arr_size = 2;
//char line2[100];
//int int1 = sprintf(line2, "2.76 TeV, N_{event} = %i", nEvent);
char* label_arr[label_arr_size] = {
  "FastJet, R = 0.4", "p_{T min, jet} > 5.0 GeV"};

const bool printOut = true;
const bool debug = false;

const int max_jets = 100;
const int max_parts = 400;
//const char* output_file_name = "Jet_ML_Prep.root";

void Jet_ML_Prep(const char* file_name, char* output_tree_description, float pt_min, float pt_max, int i, char* dir_data) {
    
  char output_file_name[200];
  sprintf(output_file_name, "./prepped/ML_Prep_%s%d.root", file_name, i);
  char csvfilename[200];
  sprintf(csvfilename, "./prepped/ML_Prep_%s%d.csv", file_name, i);
  string rawjetcsv("prepped/rawjets.csv");
  ofstream rawjetfile;
  rawjetfile.open(rawjetcsv);
  ofstream csvfile;
  csvfile.open(csvfilename);
  char output_tree_name[200];
  sprintf(output_tree_name, "Tree_%s", "Tree");
    
  std::cout << "----- Preparing " << output_tree_name << " -----" << std::endl;
    
  // Combined file of jets from thermal and PYTHIA data
  char combined_file_path[200];
  sprintf(combined_file_path, "%s/Pyth_%s_Trees%d.root", dir_data, file_name, i);
  TFile* combined_file = new TFile(combined_file_path);
  TTree* combJet_tree = (TTree*) combined_file->Get("Otherfj_tree");
    
  std::cout << "Reading combined file." << std::endl;
    
  int     c_jet_n;
  float   c_jet_pt[max_jets];
  float   c_jet_y[max_jets];
  float   c_jet_phi[max_jets];
  float   c_jet_E[max_jets];
  float   c_jet_mass[max_jets];
  float   c_jet_area[max_jets];
  float   c_jet_area_err[max_jets];
  int     c_jet_const_n[max_jets];
  float   c_jet_const_pt[max_jets][max_parts];
  float   c_jet_const_eta[max_jets][max_parts];
  float   c_jet_const_phi[max_jets][max_parts];
  float   c_jet_const_E[max_jets][max_parts];
    
  combJet_tree->SetBranchAddress("jet_n1",         &c_jet_n);
  combJet_tree->SetBranchAddress("jet_pt",        c_jet_pt);
  combJet_tree->SetBranchAddress("jet_y",         c_jet_y);
  combJet_tree->SetBranchAddress("jet_phi",       c_jet_phi);
  combJet_tree->SetBranchAddress("jet_E",         c_jet_E);
  combJet_tree->SetBranchAddress("jet_area",      c_jet_area);
  combJet_tree->SetBranchAddress("jet_area_err",  c_jet_area_err);
  combJet_tree->SetBranchAddress("jet_mass",      c_jet_mass);
  combJet_tree->SetBranchAddress("jet_const_n",   c_jet_const_n);
  combJet_tree->SetBranchAddress("jet_const_pt",  c_jet_const_pt);
  combJet_tree->SetBranchAddress("jet_const_eta", c_jet_const_eta);
  combJet_tree->SetBranchAddress("jet_const_phi", c_jet_const_phi);
  combJet_tree->SetBranchAddress("jet_const_E",   c_jet_const_E);
    
  TTree* pythJet_tree = (TTree*) combined_file->Get("FastJet_Tree");
    
  std::cout << "Reading PYTHIA file." << std::endl;
    
  int     p_jet_n;
  float   p_jet_pt[max_jets];
  float   p_jet_y[max_jets];
  float   p_jet_phi[max_jets];
  float   p_jet_E[max_jets];
  float   p_jet_mass[max_jets];
  float   p_jet_area[max_jets];
  float   p_jet_area_err[max_jets];
  int     p_jet_const_n[max_jets];
  float   p_jet_const_pt[max_jets][max_parts];
  float   p_jet_const_eta[max_jets][max_parts];
  float   p_jet_const_phi[max_jets][max_parts];
  float   p_jet_const_E[max_jets][max_parts];
  pythJet_tree->SetBranchAddress("jet_n",         &p_jet_n);
  pythJet_tree->SetBranchAddress("jet_pt",        p_jet_pt);
  pythJet_tree->SetBranchAddress("jet_y",         p_jet_y);
  pythJet_tree->SetBranchAddress("jet_phi",       p_jet_phi);
  pythJet_tree->SetBranchAddress("jet_E",         p_jet_E);
  pythJet_tree->SetBranchAddress("jet_area",      p_jet_area);
  pythJet_tree->SetBranchAddress("jet_area_err",  p_jet_area_err);
  pythJet_tree->SetBranchAddress("jet_mass",      p_jet_mass);
  pythJet_tree->SetBranchAddress("jet_const_n",   p_jet_const_n);
  pythJet_tree->SetBranchAddress("jet_const_pt",  p_jet_const_pt);
  pythJet_tree->SetBranchAddress("jet_const_eta", p_jet_const_eta);
  pythJet_tree->SetBranchAddress("jet_const_phi", p_jet_const_phi);
  pythJet_tree->SetBranchAddress("jet_const_E",   p_jet_const_E);
  // Sets output file and plot directory
  char output_file_path[200];
  sprintf(output_file_path, "%s/%s", "./", output_file_name);
  TFile* output_file = new TFile(output_file_path, "RECREATE");
  TTree* output_tree = new TTree(output_tree_name, output_tree_description);
  std::cout << "All files open." << std::endl;
  char subdir_plots[200];
  sprintf(subdir_plots, "%s/MachineLearning", dir_plots);
    
  // Plotted Data Output Tree
  // NOTE: This assumes jets have already been sorted from highest E to lowest E by FastJet!
    
  int nEvent = pythJet_tree->GetEntries();
    
  int jet_total_n = 0;
  int jet_hard_n = 0;
  int jet_soft_n = 0;
  float  jet_E_raw[max_jets];
  float  jet_E_reco[max_jets];
  float  jet_E_true[max_jets];
  //float  background_E_median[nEvent];
  //float  background_pt_median[nEvent];
  //float  thermal_E_median[nEvent];
  //float  jet_const_n_mean[nEvent];
  //float  jet_const_n_median[nEvent];
    
  // ML X value variables (features)
  int    jet_n;
  float  jet_pt_raw[max_jets];
  float  jet_pt_corr[max_jets];
  float  jet_pt_true_pythia[max_jets];
  float  jet_pt_true_paper[max_jets];
  float  jet_mass[max_jets];
  float  jet_area[max_jets];
  float  jet_area_err[max_jets];
  int    jet_const_n[max_jets];
  float  const_pt_mean[max_jets];
  float  const_pt_median[max_jets];
  float  const_1_pt[max_jets];
  float  const_2_pt[max_jets];
  float  const_3_pt[max_jets];
  float  const_4_pt[max_jets];
  float  const_5_pt[max_jets];
  float  const_6_pt[max_jets];
  float  const_7_pt[max_jets];
  float  const_8_pt[max_jets];
  float  const_9_pt[max_jets];
  float  const_10_pt[max_jets];
  float  jet_y[max_jets];
  float  jet_phi[max_jets];
  float  jet_rho[max_jets];
  float  background_rho;
  
  output_tree->Branch("jet_n",                &jet_n);
  output_tree->Branch("jet_pt_raw",           jet_pt_raw,         "jet_pt_raw[jet_n]/F");
  output_tree->Branch("jet_pt_corr",          jet_pt_corr,        "jet_pt_corr[jet_n]/F");
  output_tree->Branch("jet_pt_true_pythia",   jet_pt_true_pythia, "jet_pt_true_pythia[jet_n]/F");
  output_tree->Branch("jet_pt_true_paper",    jet_pt_true_paper,  "jet_pt_true_paper[jet_n]/F");
  output_tree->Branch("jet_mass",             jet_mass,           "jet_mass[jet_n]/F");
  output_tree->Branch("jet_area",             jet_area,           "jet_area[jet_n]/F");
  output_tree->Branch("jet_const_n",          jet_const_n,        "jet_const_n[jet_n]/I");
  output_tree->Branch("const_pt_mean",        const_pt_mean,      "const_pt_mean[jet_n]/F");
  output_tree->Branch("const_pt_median",      const_pt_median,    "const_pt_median[jet_n]/F");
  output_tree->Branch("const_1_pt",           const_1_pt,         "const_1_pt[jet_n]/F");
  output_tree->Branch("const_2_pt",           const_2_pt,         "const_2_pt[jet_n]/F");
  output_tree->Branch("const_3_pt",           const_3_pt,         "const_3_pt[jet_n]/F");
  output_tree->Branch("const_4_pt",           const_4_pt,         "const_4_pt[jet_n]/F");
  output_tree->Branch("const_5_pt",           const_5_pt,         "const_5_pt[jet_n]/F");
  output_tree->Branch("const_6_pt",           const_6_pt,         "const_6_pt[jet_n]/F");
  output_tree->Branch("const_7_pt",           const_7_pt,         "const_7_pt[jet_n]/F");
  output_tree->Branch("const_8_pt",           const_8_pt,         "const_8_pt[jet_n]/F");
  output_tree->Branch("const_9_pt",           const_9_pt,         "const_9_pt[jet_n]/F");
  output_tree->Branch("const_10_pt",          const_10_pt,        "const_10_pt[jet_n]/F");
  output_tree->Branch("jet_y",                jet_y,              "jet_y[jet_n]/F");
  output_tree->Branch("jet_phi",              jet_phi,            "jet_phi[jet_n]/F");
  output_tree->Branch("jet_rho",              jet_rho,            "jet_rho[jet_n]/F");
  output_tree->Branch("background_rho",       &background_rho);
    
  int jet_const_n_total = 0;
  int jet_true_pythia_counter = 0;
  int jet_true_paper_counter = 0;
  std::cout << "Setup complete; getting events." << std::endl;
  for ( int e = 0 ; e < combJet_tree->GetEntries() ; e++ ) {
    combJet_tree->GetEntry(e);
    pythJet_tree->GetEntry(e);
    // Finds the median pt of the background jets
    jet_n = 0;
    float  background_jet_pt[max_jets];
    float  background_jet_area[max_jets];
    float  background_rho_arr[max_jets];
    float  background_area_median;
    int    background_jet_n = 0;
    float  jet_const_n_arr[max_parts];
    int    jet_const_n_event = 0;
        
    for ( int j = 0 ; j < c_jet_n ; j++ )
      {
	jet_const_n_total++;
	jet_const_n_arr[j] = c_jet_const_n[j];
	jet_const_n_event++;
	if (j > 1) {
	  background_jet_pt[j-2]   = c_jet_pt[j];
	  background_jet_area[j-2] = c_jet_area[j];
	  background_rho_arr[j-2]  = c_jet_pt[j] / c_jet_area[j];
	  background_jet_n++;
	}
      }
    background_area_median  = TMath::Median(background_jet_n, background_jet_area);
    //background_pt_median[e] = TMath::Median(background_jet_n, background_jet_pt) * background_area_median;
    background_rho          = TMath::Median(background_jet_n, background_rho_arr);
    //        background_rho = background_pt_median[e] / background_area_median;
    //cout << "test0.5" << endl;
    //jet_const_n_mean[e] = TMath::Mean(jet_const_n_event, jet_const_n_arr);
    //jet_const_n_median[e] = TMath::Median(jet_const_n_event, jet_const_n_arr);
        
    // Loops through jets for machine learning predictors
        
    float const_total_pt;
    float const_pythia_pt[p_jet_n];
    int fill = 1;
    //cout << "test0" << endl;
    for(int q=0; q<p_jet_n; ++q)
      {
	jet_pt_true_paper[q]    = 0;
	jet_pt_true_pythia[q]   = 0;
	jet_pt_raw[q]       = 0;
	jet_pt_corr[q]      = 0;
	jet_mass[q]         = 0;
	jet_area[q]         = 0;
	jet_const_n[q]      = 0;
	const_pt_mean[q]    = 0;
	const_pt_median[q]  = 0;
	const_1_pt[q]       = 0;
	const_2_pt[q]       = 0;
	const_3_pt[q]       = 0;
	const_4_pt[q]       = 0;
	const_5_pt[q]       = 0;
	const_6_pt[q]       = 0;
	const_7_pt[q]       = 0;
	const_8_pt[q]       = 0;
	const_9_pt[q]       = 0;
	const_10_pt[q]      = 0;
	jet_y[q]            = 0;
	jet_phi[q]          = 0;
	jet_rho[q]          = 0;
      }
    for ( int j = 0 ; j < p_jet_n ; j++ ) {
      int c_jet_match_check[100] = {0};
      int c_jet_pot_pt[100] = {0};
      // Resets values for each iteration
      //int c_jet_max = 0;
      if(p_jet_pt[j] < pt_min || p_jet_pt[j] > pt_max)
	{
	  continue;
	}
      if(abs(p_jet_y[j]) > 0.5) continue;
      const_total_pt = 0.;
      for(int i=0; i<p_jet_n; ++i)
	{
	  const_pythia_pt[i] = 0.;
	}
      
      // Iterate through pythia jets to match jet pT_true
      int pythia_match = -1;
      int comb_match = -1;
      for ( int k = 0 ; k < c_jet_n ; k++ ) {
	if(c_jet_match_check[k]) continue;
	
	float dy = pow(p_jet_y[j]-c_jet_y[k],2);
	float dphi1;
	if(c_jet_phi[k] > p_jet_phi[j])
	  {
	    dphi1 = pow(c_jet_phi[k]-p_jet_phi[j]-2*3.1415,2);
	  }
	else
	  {
	    dphi1 = pow(p_jet_phi[j]-c_jet_phi[k]-2*3.1415,2);
	  }
	float dphi2 = pow(p_jet_phi[j]-c_jet_phi[k],2);
	float dphi;
	if(dphi2 < dphi1)
	  {
	    dphi = dphi2;
	  }
	else
	  {
	    dphi = dphi1;
	  }
	if (dphi + dy < fj_rSquared)
	  {
	  //if(c_jet_max < c_jet_pt[k])
	  //  {
	    c_jet_pot_pt[k] = 1;
	      //   c_jet_max = c_jet_pt[k];
	      //  }
	  }
      }
      int bestmatch = -1;
      float bestpt = 0;
      for(int k=0; k<c_jet_n; ++k)
	{
	  if(c_jet_pot_pt[k])
	    {
	      if(c_jet_pt[k] > bestpt)
		{
		  bestmatch = k;
		  pythia_match = j;
		  bestpt = c_jet_pt[k];
		}
	    }
	}
      comb_match = bestmatch;
      if(comb_match > 0)
	{
      c_jet_match_check[comb_match] = 1;
	}
      if(pythia_match < 0) {
	if(p_jet_pt[j] > pt_min && p_jet_pt[j] < pt_max)
	{
	  cout << "WARNING: truth jet "<< e << " " << j << " has no match!" << endl;
	}
	continue;
      }
      //cout << "test1" << endl;
      // Iterate through constituent particles to collect their pT for mean and median
      // This checks EVERY pythia particle in the jet, regardless of the jet match status
      
      float  jet_const_pt_arr[400];
      for(int k=0; k<p_jet_n; ++k)
	{
	  for(int q=0; q<p_jet_const_n[k]; ++q)
	    {
	      const_pythia_pt[k] += p_jet_const_pt[k][q];
	    }
	}
      //cout << "test1.0625" << endl;
      for ( int p = 0 ; p < c_jet_const_n[comb_match] ; p++ ) {
	jet_const_pt_arr[p] = c_jet_const_pt[comb_match][p];
	const_total_pt += c_jet_const_pt[comb_match][p];
      }
      //cout << "test1.125" << endl;
      //cout << jet_const_pt_arr[0] << endl;
      std::sort(jet_const_pt_arr, jet_const_pt_arr + c_jet_const_n[comb_match], greater<>());
      if((c_jet_pt[comb_match] > 400) || (c_jet_pt[comb_match] < 0.04))
	{
	  continue;
	}
      else
	{
	  jet_pt_raw[jet_n]       = c_jet_pt[comb_match];
	}
      /*
      if(c_jet_const_n[comb_match] < 2)
	{
	  continue;
	}
      */
      jet_pt_corr[jet_n]      = jet_pt_raw[jet_n] - (background_rho * c_jet_area[comb_match]);
      jet_mass[jet_n]         = c_jet_mass[comb_match];
      jet_area[jet_n]         = c_jet_area[comb_match];
      jet_const_n[jet_n]      = c_jet_const_n[comb_match];
      //cout << "test1.25" << endl;
      //YES, the below lines are correct. Look at the immediately preceding ones, dumbass!
      const_pt_mean[jet_n]    = TMath::Mean(jet_const_n[jet_n], jet_const_pt_arr);
      const_pt_median[jet_n]  = TMath::Median(jet_const_n[jet_n], jet_const_pt_arr);
      const_1_pt[jet_n]       = jet_const_pt_arr[0];
      const_2_pt[jet_n]       = jet_const_pt_arr[1];
      const_3_pt[jet_n]       = jet_const_pt_arr[2];
      const_4_pt[jet_n]       = jet_const_pt_arr[3];
      const_5_pt[jet_n]       = jet_const_pt_arr[4];
      const_6_pt[jet_n]       = jet_const_pt_arr[5];
      const_7_pt[jet_n]       = jet_const_pt_arr[6];
      const_8_pt[jet_n]       = jet_const_pt_arr[7];
      const_9_pt[jet_n]       = jet_const_pt_arr[8];
      const_10_pt[jet_n]      = jet_const_pt_arr[9];
      jet_y[jet_n]            = c_jet_y[comb_match];
      jet_phi[jet_n]          = c_jet_phi[comb_match];
      jet_rho[jet_n]          = background_rho;
      //cout << "test2" << endl;
      jet_pt_true_paper[jet_n]  = jet_pt_raw[comb_match] * (const_pythia_pt[pythia_match] / const_total_pt);
      if(jet_pt_true_paper[comb_match] > 200 || jet_pt_true_paper[comb_match] < 0.04)
	{
	  jet_pt_true_paper[jet_n] = 0;
	}
      if (e % 1000 == 0 ) std::cout << e << "-" << j << ": Raw Jet: " << jet_pt_raw[jet_n] << std::endl;
            
      if ( pythia_match < 0 ) continue;
            
      jet_pt_true_pythia[jet_n] = p_jet_pt[pythia_match];
      if(jet_pt_true_pythia[jet_n] < pt_max && jet_pt_true_pythia[jet_n] > pt_min)
	{
	  csvfile << jet_pt_true_pythia[jet_n] <<"," << jet_pt_corr[jet_n] <<","<< jet_pt_raw[jet_n] <<","<< jet_area[jet_n] <<","<< jet_mass[jet_n] <<","<< jet_const_n[jet_n] <<","<< const_pt_mean[jet_n] <<","<< const_pt_median[jet_n] <<","<< jet_rho[jet_n] <<","<< const_1_pt[jet_n] <<","<< const_2_pt[jet_n] <<","<< const_3_pt[jet_n] <<","<< const_4_pt[jet_n] << endl;
	  if((pt_min == 40) && (i == 8)) rawjetfile << jet_pt_raw[jet_n] << endl;
	}
      jet_true_pythia_counter++;
      if (e % 1000 == 0 ) std::cout << e << "-" << j << ": Truth_PYTHIA Jet: " << jet_pt_true_pythia[jet_n] << " ----- " << std::endl;
      if(pythia_match >= 0)
	{
	  jet_n++;
	}
    }
    if(fill) output_tree->Fill();
    
  }
  csvfile.close();
  rawjetfile.close();
  output_file->Write();
    
  std::cout << "Data written to file." << std::endl;
    
  combined_file->Close();
  output_file->Close();
    
  std::cout << "Files closed." << std::endl;
    
  std::cout << "----- Completed " << output_tree_name << " -----" << std::endl;
}

void joeymlprep() {
  vector<int> powers;
  vector<vector<int>> ranges;
  ifstream openfile;
  openfile.open("powers.txt");
  int a, b;
  while(openfile >> a)
    {
      powers.push_back(a);
    }
  openfile.close();
  openfile.open("ranges.txt");
  while(openfile >> a >> b)
    {
      vector<int> temp;
      temp.push_back(a);
      temp.push_back(b);
      ranges.push_back(temp);
    }
  for(int i=0; i<powers.size(); ++i)
    {
      Jet_ML_Prep(
		  (to_string(ranges[1][0]) + "_" + to_string(ranges[1][1]) + "_Test").c_str(),
		  "TESTING data for machine learning. 40-60 GeV with pT^n bias.",
		  ranges[1][0], ranges[1][1], powers[i],"/media/qq/MAIN");
      Jet_ML_Prep(
		  (to_string(ranges[0][0]) + "_" + to_string(ranges[0][1]) + "_Train").c_str(),
		  "TRAINING data for machine learning. 10-90 GeV with pT^n bias.",
		  ranges[0][0], ranges[0][1], powers[i],"/media/qq/MAIN");
    }
}
