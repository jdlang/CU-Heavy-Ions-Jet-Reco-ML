#include "Pythia8/Pythia.h"
using namespace Pythia8;
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/Selector.hh"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include <cmath>
#include <iostream>
#include "_header.h"
using namespace Project_Constants;

const bool print_out   = true;
const int  print_every_x = 100;
const bool debug       = false;
const bool use_voronoi = false; // Uses Voronoi for jet clustering



// ----- PYTHIA GENERATOR -----


TF1* n_func;
TF1* pt_func;
TF1* eta_func;
TF1* phi_func;

double Gaussian_Func(double *x,double *par) {
    double arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    double fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}

double Power_Func(double *x,double *par) {
    double arg = 0;
    if (par[2]!=0) arg = (x[0] - par[2]);
    double fitval = par[0] * pow(arg, par[1]);
    return fitval;
}

double FlatFunc(double *x,double *par) {
    double fitval = par[0] * 1;
    return fitval;
}

double ModifiedHagedorn_Func(double* x, double* par) {
    double arg1 = par[0] * pow(x[0], 2) / pow( (pow(x[0], 2) + pow(par[1], 2)), 0.5);
    double arg2 = pow( 1 + (x[0] /par[2]), par[3] );
    double arg3 = 1 / x[0];
    double fitval =  arg1 * arg2; // * arg3;
    return fitval;
}

int GenN() {
    int n = round( n_func->GetRandom(0, 2 * gaus_mean) );
    return n;
}

double GenPt() {
//    double pt = pt_func->GetRandom(0, 100);
    double pt = pt_func->GetRandom(0, 100);
    return pt;
}

double GenEta() {
    double eta = eta_func->GetRandom(-1 * eta_max, eta_max);
    return eta;
}

double GenPhi() {
    double phi = phi_func->GetRandom(-1 * math_pi, math_pi);
    return phi;
}

void Pythia_Generator(char* file_name, int func_eventCount, float func_beamPower, float func_ptHatMin,
                      float func_ptBiasPow, float func_slimMin, float func_slimMax, float func_slimRap) {
    
    char file_path[200];
    int temp_int = sprintf(file_path, "~/datatemp/Pyth_%s", file_name);
    TFile* output_file = new TFile(file_path, "RECREATE");
    TTree* pythia_tree = new TTree("Combined_Tree","Tree of particle jet events from PYTHIA p+p collisions");
    TTree* jet_tree    = new TTree("FastJet_Tree","Tree of jet clusters by event from PYTHIA");
    TTree* otherjet    = new TTree("Otherfj_tree","Tree of jet clusters from combined event.");
    
    if (print_out) std::cout << "Output file and trees successfully generated." << std::endl;

    n_func = new TF1("n_func", Gaussian_Func, 0, 2 * gaus_mean, 3);
    pt_func = new TF1("pt_func", ModifiedHagedorn_Func, 0, 100, 4);
    eta_func = new TF1("eta_func", FlatFunc, -1 * eta_max, eta_max, 1);
    phi_func = new TF1("phi_func", FlatFunc, -1 * math_pi, math_pi, 1);

    n_func->SetParameters(gaus_norm, gaus_mean, gaus_sigma);
    pt_func->SetParameters(64547, 3.076, 1.126, -8.491);
    eta_func->SetParameter(0, 1.);
    phi_func->SetParameter(0, 1.);
    
    // Output Variables
    
    // Initializes arrays for particle info, DO NOT USE DOUBLES!!!
    const int max_parts  = 6000;
    const int max_jets   = 100;
    const int max_consts = 200;
    
    int     particle_n;                 // particle index
    float   particle_pt[max_parts];     // transverse momentum
    float   particle_phi[max_parts];    // azimuthal angle
    float   particle_eta[max_parts];    // pseudorapidity
    int     particle_m[max_parts];      // particle mass
    
    int     jet_n;                      // jet index
    float   jet_pt[max_jets];           // transverse momentum
    float   jet_y[max_jets];            // rapidity (FastJet does not output eta)
    float   jet_phi[max_jets];          // azimuthal angle
    float   jet_E[max_jets];            // jet energy
    float   jet_mass[max_jets];
    float   jet_area[max_jets];
    float   jet_area_err[max_jets];
    int     jet_const_n[max_jets];
    float   jet_const_pt[max_jets][max_consts];
    float   jet_const_eta[max_jets][max_consts];
    float   jet_const_phi[max_jets][max_consts];
    float   jet_const_E[max_jets][max_consts];

    int     jet_n1;                      // jet index
    float   jet_pt1[max_jets];           // transverse momentum
    float   jet_y1[max_jets];            // rapidity (FastJet does not output eta)
    float   jet_phi1[max_jets];          // azimuthal angle
    float   jet_E1[max_jets];            // jet energy
    float   jet_mass1[max_jets];
    float   jet_area1[max_jets];
    float   jet_area_err1[max_jets];
    int     jet_const_n1[max_jets];
    float   jet_const_pt1[max_jets][max_consts];
    float   jet_const_eta1[max_jets][max_consts];
    float   jet_const_phi1[max_jets][max_consts];
    float   jet_const_E1[max_jets][max_consts];
    
    // Builds particle and jet branches
    // Note: Can only do variable size for first variable in array
    pythia_tree->Branch("particle_n",   &particle_n);
    pythia_tree->Branch("particle_pt",  particle_pt,  "particle_pt[particle_n]/F");
    pythia_tree->Branch("particle_eta", particle_eta, "particle_eta[particle_n]/F");
    pythia_tree->Branch("particle_phi", particle_phi, "particle_phi[particle_n]/F");
    pythia_tree->Branch("particle_m", particle_m, "particle_m[particle_n]/F");
    
    jet_tree->Branch("jet_n",        &jet_n);
    jet_tree->Branch("jet_pt",       jet_pt,       "jet_pt[jet_n]/F");
    jet_tree->Branch("jet_y",        jet_y,        "jet_y[jet_n]/F");
    jet_tree->Branch("jet_phi",      jet_phi,      "jet_phi[jet_n]/F");
    jet_tree->Branch("jet_E",        jet_E,        "jet_E[jet_n]/F");
    jet_tree->Branch("jet_mass",     jet_mass,     "jet_mass[jet_n]/F");
    jet_tree->Branch("jet_area",     jet_area,     "jet_area[jet_n]/F");
    jet_tree->Branch("jet_area_err", jet_area_err, "jet_area_err[jet_n]/F");
    jet_tree->Branch("jet_const_n",     jet_const_n,    "jet_const_n[jet_n]/I");
    jet_tree->Branch("jet_const_pt",    jet_const_pt,   "jet_const_pt[jet_n][400]/F");
    jet_tree->Branch("jet_const_eta",   jet_const_eta,  "jet_const_eta[jet_n][400]/F");
    jet_tree->Branch("jet_const_phi",   jet_const_phi,  "jet_const_phi[jet_n][400]/F");
    jet_tree->Branch("jet_const_E",     jet_const_E,    "jet_const_E[jet_n][400]/F");

    otherjet->Branch("jet_n1",        &jet_n1);
    otherjet->Branch("jet_pt",       jet_pt1,       "jet_pt1[jet_n1]/F");
    otherjet->Branch("jet_y",        jet_y1,        "jet_y[jet_n1]/F");
    otherjet->Branch("jet_phi",      jet_phi1,      "jet_phi[jet_n1]/F");
    otherjet->Branch("jet_E",        jet_E1,        "jet_E[jet_n1]/F");
    otherjet->Branch("jet_mass",     jet_mass1,     "jet_mass[jet_n1]/F");
    otherjet->Branch("jet_area",     jet_area1,     "jet_area[jet_n1]/F");
    otherjet->Branch("jet_area_err", jet_area_err1, "jet_area_err[jet_n1]/F");
    otherjet->Branch("jet_const_n",     jet_const_n1,    "jet_const_n[jet_n1]/I");
    otherjet->Branch("jet_const_pt",    jet_const_pt1,   "jet_const_pt[jet_n1][400]/F");
    otherjet->Branch("jet_const_eta",   jet_const_eta1,  "jet_const_eta[jet_n1][400]/F");
    otherjet->Branch("jet_const_phi",   jet_const_phi1,  "jet_const_phi[jet_n1][400]/F");
    otherjet->Branch("jet_const_E",     jet_const_E1,    "jet_const_E[jet_n1][400]/F");
    
    // --- PYTHIA Setup ---
    // --- LHC process and output selection. Initialization.
    Pythia pythia;
    Settings& pythia_settings = pythia.settings;
    pythia_settings.parm("Beams:eCM", func_beamPower);  // Uses units of GeV
    pythia_settings.parm("PhaseSpace:pTHatMin", func_ptHatMin );
    pythia_settings.parm("PhaseSpace:pTHatMax", 0.);
    pythia.readString("HardQCD:all = on");              // Turns on hard scattering
    if ( func_ptBiasPow >= 0. ) {
        pythia_settings.parm("PhaseSpace:bias2SelectionPow", func_ptBiasPow);
        pythia.readString("PhaseSpace:bias2Selection = on");
        pythia.readString("PhaseSpace:bias2SelectionRef = 100.");
    }
    pythia.init();
    // -- to try later
    // multiple ptHatMin
    // combine and reweight ptHatMins
    
    if ( print_out ) std::cout << "PYTHIA settings initialized." << std::endl;
    
    // --- FastJet Setup ---
    // --- Cluster each PYTHIA event
    fastjet::JetDefinition jet_definition(fastjet::antikt_algorithm, fj_jetR);
    fastjet::AreaDefinition area_definition;
    cout << fj_jetR << endl;
    if (!use_voronoi) {
        double ghost_etamax = 0.9;
        double ghost_area    = 0.05;
        int    active_area_repeats = 5;
        fastjet::GhostedAreaSpec ghost_spec(ghost_etamax, active_area_repeats, ghost_area);
        area_definition = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
    } else {
        double effective_Rfact = 1.0;
        area_definition = fastjet::VoronoiAreaSpec(effective_Rfact);
    }
    
    if ( print_out ) std::cout << "FastJet settings initialized." << std::endl;
    
    // Loop to generate PYTHIA events. Skip if error.
    int e = 0;
    while ( e < func_eventCount ) {
        // ---- PYTHIA -----
        if (!pythia.next()) continue;

        particle_n = 0;

        for ( int p = 0; p < pythia.event.size(); p++ ) {
            if (!pythia.event[p].isFinal()) continue;           // Ends if final particle
            if (!pythia.event[p].isCharged()) continue;         // Skips neutral particles
            if (pythia.event[p].pT() < 0.15) continue;          // Skips particles with pT < 0.15
            if (fabs(pythia.event[p].eta()) > 0.9) continue;    // Skips particles with eta outside +/-0.9

            particle_pt[particle_n]  = pythia.event[p].pT();
            particle_eta[particle_n] = pythia.event[p].eta();
            particle_phi[particle_n] = pythia.event[p].phi();
	    particle_m[particle_n] = -pythia.event[p].m();

            particle_n++;
        }

	int nParticles = GenN();
        
        // ---------- FASTJET ----------
        
        TLorentzVector input_vec;
        std::vector<fastjet::PseudoJet> input_particles;
        
        for ( int p = 0 ; p < particle_n ; p++) {
            input_vec.SetPtEtaPhiM(particle_pt[p], particle_eta[p], particle_phi[p], particle_m[p]);
            input_particles.push_back(fastjet::PseudoJet( input_vec.Px(), input_vec.Py(), input_vec.Pz(), input_vec.E() ));
        }
        
        fastjet::ClusterSequenceArea jet_clusters(input_particles, jet_definition, area_definition);
        std::vector<fastjet::PseudoJet> jet_list = sorted_by_pt(jet_clusters.inclusive_jets(8));
        // jet_clusters.inclusive_jets(fj_jetptmin)  // Use if NOT recording area
        
        jet_n = jet_list.size();
        if ( jet_n == 0 ) continue;
        if ( jet_n >= max_jets ) jet_n = max_jets;
        
        bool acceptJet = false;
        int flag = 0;
        for ( int j = 0 ; j < jet_n ; j++ ) {
            jet_pt[j]   = jet_list[j].pt();
            jet_y[j]    = jet_list[j].rap();
            jet_phi[j]  = jet_list[j].phi();
            jet_E[j]    = jet_list[j].E();
            jet_mass[j] = jet_list[j].m();
            jet_area[j] = jet_list[j].area();
            jet_area_err[j] = jet_list[j].area_error();
            
            if ( (jet_pt[j] >= func_slimMin) && (jet_pt[j] <= func_slimMax) &&
                (jet_y[j] >= -func_slimRap) && (jet_y[j] <= func_slimRap) ) acceptJet = true;
            
            std::vector<fastjet::PseudoJet> jet_constituents = jet_list[j].constituents();
            jet_const_n[j]  = jet_constituents.size();
            
            for ( int p = 0 ; p < jet_constituents.size() ; p++) {
                jet_const_pt[j][p]  = jet_constituents[p].pt();
                jet_const_eta[j][p] = jet_constituents[p].eta();
                jet_const_phi[j][p] = jet_constituents[p].phi();
                jet_const_E[j][p]   = jet_constituents[p].E();
            }
        }
        if ( acceptJet ) {
	  for (int p = 0 ; p < nParticles ; p++) {
            particle_pt[particle_n] = GenPt();
            particle_eta[particle_n] = GenEta();
            particle_phi[particle_n] = GenPhi();
	    particle_m[particle_n] = m_pion;
            
            // Increments particle_n for the next particle
            particle_n++;
	  }

	  input_particles.clear();
	  for ( int p = 0 ; p < particle_n ; p++) {
	    input_vec.SetPtEtaPhiM(particle_pt[p], particle_eta[p], particle_phi[p], particle_m[p]);
	    input_particles.push_back(fastjet::PseudoJet( input_vec.Px(), input_vec.Py(), input_vec.Pz(), input_vec.E() ));
	  }

	  fastjet::ClusterSequenceArea jet_clusters1(input_particles, jet_definition, area_definition);
	  jet_list = sorted_by_pt(jet_clusters1.inclusive_jets(8));


	  jet_n1 = jet_list.size();
	  if ( jet_n1 == 0 ) flag = 1;
	  if ( jet_n1 >= max_jets ) jet_n1 = max_jets;
	  
	  
	  for ( int j = 0 ; j < jet_n1; j++ ) {
            jet_pt1[j]   = jet_list[j].pt();
            jet_y1[j]    = jet_list[j].rap();
            jet_phi1[j]  = jet_list[j].phi();
            jet_E1[j]    = jet_list[j].E();
            jet_mass1[j] = jet_list[j].m();
            jet_area1[j] = jet_list[j].area();
            jet_area_err1[j] = jet_list[j].area_error();
            std::vector<fastjet::PseudoJet> jet_constituents = jet_list[j].constituents();
            jet_const_n1[j]  = jet_constituents.size();
            
            for ( int p = 0 ; p < jet_constituents.size() ; p++) {
	      jet_const_pt1[j][p]  = jet_constituents[p].pt();
	      jet_const_eta1[j][p] = jet_constituents[p].eta();
	      jet_const_phi1[j][p] = jet_constituents[p].phi();
	      jet_const_E1[j][p]   = jet_constituents[p].E();
            }
	  }
	  if((e % 100 == 0) && (e != 0))
	    {
	      cout << "Events " << e-100 <<"-"<< e-1 << " written to file." << endl;
	    }
	  if(!flag)
	    {
	      pythia_tree ->Fill();
	      jet_tree    ->Fill();
	      otherjet    ->Fill();
	      e++;
	    }
	}
    }
	
    if ( print_out ) pythia.stat(); // Prints PYTHIA stats at the end

    pythia_tree ->Write();
    jet_tree    ->Write();
    otherjet    ->Write();
    std::cout << "Output trees written to." << std::endl;
    output_file->Close();
    

    
    std::cout << "Files saved and closed." << std::endl;
}

void Event_Generator(char* file_name, int func_eventCount, float func_beamPower,
                     float func_ptBiasPow, float func_slimMin, float func_slimMax, float func_slimRap) {
    
    float func_ptHatMin = 0.75 * func_slimMin;
    
    std::cout << ">>> Generate PYTHIA Events <<<" << std::endl;
    Pythia_Generator(file_name, func_eventCount, func_beamPower, func_ptHatMin,
                     func_ptBiasPow, func_slimMin, func_slimMax, func_slimRap);

    
    std::cout << ">>> Event Generation and Jet Clustering Complete! <<<" << std::endl;
}


int main() {
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
  int ntrials = 500000;
  for(int i=0; i<powers.size(); ++i)
    {
      string s(to_string(ranges[0][0]) + "_" + to_string(ranges[0][1]) + "_Train_Trees");
      string t = to_string(powers[i]);
      string u(".root");
      string stu = s+t+u;
      string r(to_string(ranges[1][0]) + "_" + to_string(ranges[1][1]) + "_Train_Trees");
      string rtu = r+t+u;
      char* file = const_cast<char*>(stu.c_str());
      char* file2 = const_cast<char*>(rtu.c_str());
      Event_Generator(
		      file, // file name
		      ntrials,     // number of events to generate
		      beamPower,  // beam power
		      (float)powers[i],        // pt bias power (pt^x), set to -1. to disable bias
		      ranges[0][0],        // pt min for slimming
		      ranges[0][1],        // pt max for slimming
		      slim_rap);  // max rapidity for slimming
      Event_Generator(
		      file2, // file name
		      ntrials,     // number of events to generate
		      beamPower,  // beam power
		      (float)powers[i],        // pt bias power (pt^x), set to -1. to disable bias
		      ranges[1][0],        // pt min for slimming
		      ranges[1][1],        // pt max for slimming
		      slim_rap);  // max rapidity for slimming
    }
  return 0;
}
