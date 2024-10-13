#include "call_libraries.h"
//boolean to differentiate between MC and DATA
bool is_MC = false; // false for running in data and true for running in MonteCarlo_pythia

// jet cuts
float jet_radius = 0.4;
double jet_pt_min = 550.;  //defauly 550.0
double jet_eta_min = -1.6;
double jet_eta_max = 1.6;

// track cut
double newtrk_eta_min = -5.0;
double newtrk_eta_max = 5.0;
double newtrk_jt03_min = 0.;
double newtrk_jt03_max = 3.;
double newtrk_jt0p33_min = 0.3;
double newtrk_jt0p33_max = 3.;
double newtrk_jt0p53_min = 0.5;
double newtrk_jt0p53_max = 3.;

//number of events to mix
int bkgFactor = 5;

//event quantities
float vz_cut_min = -15.0; //-ve vz acceptance
float vz_cut_max = 15.0; //+ve vz acceptance
std::vector<TString> event_filter_str; // skimed event filters
float pthat_cut = 50.0;

TString jet_trigger;
TString pp_jet_trigger = "didHLTFire400";

TString jet_collection;
TString pp_jet_collection = "analyzer500"; // PF jets

TString JEC_file;
TString JEC_file_data;
TString pp_JEC_file = "Spring18_ppRef5TeV_V6_MC_L2Relative_AK4PF.txt"; // for pp ak4PF jets
TString pp_JEC_file_data = "Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK4PF.txt"; // residual for pp ak4PF jets for data

// Correction for pp

TString fmb="2017pp_TrkCorr_Sept25_Final.root";
