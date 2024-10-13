#include "call_libraries.h"  // call libraries from ROOT and C++
#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
// event quantities
float vertexz;
UInt_t event;
UInt_t lumi;
ULong64_t run;
// events quantities from gen
float pthatweight; // event weight --> pthat weight
float pthat;  // pthat (initial parton pT)


// trigger quantities
bool jet_trigger_bit; // jet HLT path trigger used for analysis (jet_trigger variable in input_variables.h)

//event filter
std::vector<int> event_filter_bool(10);

//~~~~~~~~~~~~~~~~~~~~~ jet and track quantities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// reco jets
const int nmax = 99999;
int nref;         // number of jets
float jetpt[nmax]; // jet pT
float jeteta[nmax]; // jet eta
float jetphi[nmax]; // jet phi
float rawpt[nmax]; // jet pT without JEC
float trackMax[nmax]; // track maximum pT in a jet
float WTAeta[nmax]; // WTA jet eta
float WTAphi[nmax]; // WTA jet phi

// reco tracks
int nTrk;                 // number of track
float trkPt[nmax];       // track pT
float trkEta[nmax];      // track eta
float trkPhi[nmax];      // track phi
float trkPtError[nmax];    // track pT error (uncertainty)
float trkDxy1[nmax];    // track dxy impact parameter (transverse distance between primary vertex and collision - distance of closest approuch - DCA)
float trkDz1[nmax];     // track dz impact parameter (longitudinal distance between primary vertex and collision - distance of closest approuch - DCA)
float trkDxyError1[nmax]; // track dxy error (uncertainty)
float trkDzError1[nmax];  // track dxy error (uncertainty)
float trkChi2[nmax];     // track reconstruction chi2 of the fitting
//float pfEcal[nmax];      // particle flow energy deposit in ECAL
//float pfHcal[nmax];      // particle flow energy deposit in HCAL
//float trkmva[nmax];      // track mva for each step
//int trkalgo[nmax];       // track algorithm/step
int trkNdof[nmax];       // track number of degrees of freedom in the fitting 
int trkCharge[nmax];     // track charge
int trkNHit[nmax];      // number of hits in the tracker
int trkNlayer[nmax];     // number of layers with measurement in the tracker
//bool highPurity[nmax];      // tracker steps MVA selection

// gen jets
int ngen;             // number of gen jets
float gen_jtpt[nmax];  // gen jet pT
float gen_jteta[nmax]; // gen jet eta
float gen_jtphi[nmax]; // gen jet phi
float gen_WTAeta[nmax]; // gen WTA jet eta
float gen_WTAphi[nmax]; // gen ETA jet phi

// matched jets
float refpt[nmax]; // jet pT matched with Gen pT
float refeta[nmax]; // jet eta matched with Gen eta
float refphi[nmax]; // jet phi matched with Gen phi
int refparton_flavor[nmax]; // jet flavor matched with Gen flavor
int refparton_flavorForB[nmax]; // jet flavor matched with Gen flavor

//gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
//gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

//higPurity

std::vector<bool> *highPurity = 0;



//jet_Vtx
std::vector<float> *xVtx = 0;
std::vector<float> *yVtx = 0;
std::vector<float> *zVtx = 0;
std::vector<float> *nTracksVtx = 0;


//jet
int jetN;
std::vector<float> *jetPt = 0;
std::vector<float> *jetEta = 0;
std::vector<float> *jetPhi = 0;
std::vector<float> *jetTheta = 0;

//daughter_tracks
std::vector<std::vector<float>> *dau_pt = 0;
std::vector<std::vector<float>> *dau_phi = 0;
std::vector<std::vector<float>> *dau_eta = 0;
std::vector<std::vector<int>> *dau_chg = 0;
std::vector<std::vector<int>> *dau_pid = 0;
std::vector<std::vector<float>> *dau_ptError = 0;
std::vector<std::vector<float>> *dau_theta = 0;
std::vector<std::vector<float>> *dau_pt_STAR = 0;
std::vector<std::vector<float>> *dau_eta_STAR = 0;
std::vector<std::vector<float>> *dau_phi_STAR = 0;
std::vector<std::vector<float>> *dau_theta_STAR = 0;
std::vector<std::vector<float>> *dau_vz = 0;
std::vector<std::vector<float>> *dau_vx = 0;
std::vector<std::vector<float>> *dau_vy = 0;
std::vector<std::vector<float>> *dau_vrefx = 0;
std::vector<std::vector<float>> *dau_vrefy = 0;
std::vector<std::vector<float>> *dau_vrefz = 0;
std::vector<std::vector<float>> *dau_XYDCAsig= 0;
std::vector<std::vector<float>> *dau_ZDCAsig= 0;
// gen tracks
std::vector<float> *gen_trkpt = 0;  // gen particle pT
std::vector<float> *gen_trketa = 0; // gen particle eta
std::vector<float> *gen_trkphi = 0; // gen particle phi
std::vector<int> *gen_trkchg = 0;   // gen particle charge
std::vector<int> *gen_trkpid  = 0;   // gen particle pid
std::vector<int> *gen_trksube = 0;   // gen particle pid


void read_tree(TTree *tree, TString jet_trigger, bool MC)
{
  tree->SetBranchStatus("*", 0); // disable all branches - this is important while reading big files
  
  // enable branches of interest -> see definition of each variables above

  // event quantities

  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");                                                                                                                                     gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
  
  tree->SetBranchStatus(Form("%s",jet_trigger.Data()), 1);
  tree->SetBranchAddress(Form("%s",jet_trigger.Data()), &jet_trigger_bit);
  /*
  tree->SetBranchStatus("vz", 1);
  tree->SetBranchAddress("vz", &vertexz);
  */
  
  tree->SetBranchStatus("nEv", 1);
  tree->SetBranchAddress("nEv", &event);

  tree->SetBranchStatus("nRun", 1);
  tree->SetBranchAddress("nRun", &run);

  tree->SetBranchStatus("nLumi", 1);
  tree->SetBranchAddress("nLumi", &lumi);
  /*
  for(int i = 0; i < (int) event_filterstr.size(); i++)
    {
      tree->SetBranchStatus(Form("%s",event_filterstr[i].Data()), 1);
      tree->SetBranchAddress(Form("%s",event_filterstr[i].Data()), &event_filter_bool[i]);
    }
  */
  
  if(MC){
  // events quantities from gen
  tree->SetBranchStatus("weight", 1);
  tree->SetBranchAddress("weight", &pthatweight);
  
  tree->SetBranchStatus("pthat", 1);
  tree->SetBranchAddress("pthat", &pthat);
  }


  // reco jet quantities                                                                                                    
  /*
  tree->SetBranchStatus("nref", 1);
  tree->SetBranchAddress("nref", &nref);
    
  tree->SetBranchStatus("rawpt", 1);
  tree->SetBranchAddress("rawpt", &rawpt);
  
  tree->SetBranchStatus("trackMax", 1);
  tree->SetBranchAddress("trackMax", &trackMax);
  */
  tree->SetBranchStatus("jetN", 1);
  tree->SetBranchAddress("jetN", &jetN);

  tree->SetBranchStatus("jetPt", 1);
  tree->SetBranchAddress("jetPt", &jetPt);

  tree->SetBranchStatus("jetEta", 1);
  tree->SetBranchAddress("jetEta", &jetEta);

  tree->SetBranchStatus("jetPhi", 1);
  tree->SetBranchAddress("jetPhi", &jetPhi);
  
  tree->SetBranchStatus("jetTheta", 1);
  tree->SetBranchAddress("jetTheta", &jetTheta);
  /*

  tree->SetBranchStatus("WTAeta", 1);
  tree->SetBranchAddress("WTAeta", &WTAeta);

  tree->SetBranchStatus("WTAphi", 1);
  tree->SetBranchAddress("WTAphi", &WTAphi);
  */
  
  if(MC){
  // gen jet quantities                                                                                                 
  tree->SetBranchStatus("ngen", 1);
  tree->SetBranchAddress("ngen", &ngen);
  
  tree->SetBranchStatus("genpt", 1);
  tree->SetBranchAddress("genpt", &gen_jtpt);
  
  tree->SetBranchStatus("geneta", 1);
  tree->SetBranchAddress("geneta", &gen_jteta);
  
  tree->SetBranchStatus("genphi", 1);
  tree->SetBranchAddress("genphi", &gen_jtphi);
  
  tree->SetBranchStatus("WTAgeneta", 1);
  tree->SetBranchAddress("WTAgeneta", &gen_WTAeta);
  
  tree->SetBranchStatus("WTAgenphi", 1);
  tree->SetBranchAddress("WTAgenphi", &gen_WTAphi);
  
  //matched gen jet quantities                                                                                          
  tree->SetBranchStatus("refpt", 1);
  tree->SetBranchAddress("refpt", &refpt);
  
  tree->SetBranchStatus("refeta", 1);
  tree->SetBranchAddress("refeta", &refeta);

  tree->SetBranchStatus("refphi", 1);
  tree->SetBranchAddress("refphi", &refphi);
  
  tree->SetBranchStatus("refparton_flavor", 1);
  tree->SetBranchAddress("refparton_flavor", &refparton_flavor);
  
  tree->SetBranchStatus("jtPartonFlavor", 1);
  tree->SetBranchAddress("jtPartonFlavor", &refparton_flavorForB);
  }

  //daughter quantities
  tree->SetBranchStatus("dau_chg", 1);                                                                                                                                                                  
  tree->SetBranchAddress("dau_chg", &dau_chg);

  tree->SetBranchStatus("dau_pt", 1);
  tree->SetBranchAddress("dau_pt", &dau_pt);

  tree->SetBranchStatus("dau_ptError", 1);
  tree->SetBranchAddress("dau_ptError", &dau_ptError);
  
  tree->SetBranchStatus("dau_eta", 1);
  tree->SetBranchAddress("dau_eta", &dau_eta);

  tree->SetBranchStatus("dau_phi", 1);
  tree->SetBranchAddress("dau_phi", &dau_phi);

  tree->SetBranchStatus("dau_theta", 1);
  tree->SetBranchAddress("dau_theta", &dau_theta);

  tree->SetBranchStatus("dau_pt_STAR", 1);
  tree->SetBranchAddress("dau_pt_STAR", &dau_pt_STAR);

  tree->SetBranchStatus("dau_eta_STAR", 1);
  tree->SetBranchAddress("dau_eta_STAR", &dau_eta_STAR);

  tree->SetBranchStatus("dau_phi_STAR", 1);
  tree->SetBranchAddress("dau_phi_STAR", &dau_phi_STAR);

  tree->SetBranchStatus("dau_theta_STAR", 1);
  tree->SetBranchAddress("dau_theta_STAR", &dau_theta_STAR);

  tree->SetBranchStatus("dau_pt_STAR", 1);
  tree->SetBranchAddress("dau_pt_STAR", &dau_pt_STAR);

  tree->SetBranchStatus("dau_vz", 1);
  tree->SetBranchAddress("dau_vz", &dau_vz);

  tree->SetBranchStatus("dau_vx", 1);
  tree->SetBranchAddress("dau_vx", &dau_vx);

  tree->SetBranchStatus("dau_vy", 1);
  tree->SetBranchAddress("dau_vy", &dau_vy);

  tree->SetBranchStatus("dau_vrefx", 1);
  tree->SetBranchAddress("dau_vrefx", &dau_vrefx);

  tree->SetBranchStatus("dau_vrefy", 1);
  tree->SetBranchAddress("dau_vrefy", &dau_vrefy);

  tree->SetBranchStatus("dau_vrefz", 1);
  tree->SetBranchAddress("dau_vrefz", &dau_vrefz);

  tree->SetBranchStatus("dau_XYDCAsig", 1);
  tree->SetBranchAddress("dau_XYDCAsig", &dau_XYDCAsig);

  tree->SetBranchStatus("dau_ZDCAsig", 1);
  tree->SetBranchAddress("dau_ZDCAsig", &dau_ZDCAsig);
  
  tree->SetBranchStatus("highPurity", 1);
  tree->SetBranchAddress("highPurity", &highPurity);

  
  // track quantities
  /*
  tree->SetBranchStatus("trkCharge", 1); 
  tree->SetBranchAddress("trkCharge", &trkCharge);

  tree->SetBranchStatus("trkPt", 1);
  tree->SetBranchAddress("trkPt", &trkPt);

  tree->SetBranchStatus("trkEta", 1); 
  tree->SetBranchAddress("trkEta", &trkEta);
  
  tree->SetBranchStatus("trkPhi", 1);
  tree->SetBranchAddress("trkPhi", &trkPhi);
  
  tree->SetBranchStatus("trkPtError", 1);
  tree->SetBranchAddress("trkPtError", &trkPtError);
  
  tree->SetBranchStatus("trkNHit", 1);
  tree->SetBranchAddress("trkNHit", &trkNHit);
  
  tree->SetBranchStatus("trkNlayer", 1);
  tree->SetBranchAddress("trkNlayer", &trkNlayer);
  
  tree->SetBranchStatus("trkDxy1", 1);
  tree->SetBranchAddress("trkDxy1", &trkDxy1); 

  tree->SetBranchStatus("trkDxyError1", 1);
  tree->SetBranchAddress("trkDxyError1", &trkDxyError1); 

  tree->SetBranchStatus("trkDz1", 1);
  tree->SetBranchAddress("trkDz1", &trkDz1); 

  tree->SetBranchStatus("trkDzError1", 1);
  tree->SetBranchAddress("trkDzError1", &trkDzError1); 

  tree->SetBranchStatus("highPurity", 1);
  tree->SetBranchAddress("highPurity", &highPurity);
  
  tree->SetBranchStatus("trkChi2", 1);
  tree->SetBranchAddress("trkChi2", &trkChi2);
  
  tree->SetBranchStatus("trkNdof", 1);
  tree->SetBranchAddress("trkNdof", &trkNdof);
  
  tree->SetBranchStatus("nTrk", 1);
  tree->SetBranchAddress("nTrk", &nTrk);
  */
  
  if(MC){
  // gen particle quantities
  tree->SetBranchStatus("pt", 1);
  tree->SetBranchAddress("pt", &gen_trkpt);
  
  tree->SetBranchStatus("eta", 1);
  tree->SetBranchAddress("eta", &gen_trketa);
  
  tree->SetBranchStatus("phi", 1);
  tree->SetBranchAddress("phi", &gen_trkphi);
  
  tree->SetBranchStatus("chg", 1);
  tree->SetBranchAddress("chg", &gen_trkchg);
  
  tree->SetBranchStatus("pdg", 1);
  tree->SetBranchAddress("pdg", &gen_trkpid);
  
  tree->SetBranchStatus("sube", 1);
  tree->SetBranchAddress("sube", &gen_trksube);
  }
}
