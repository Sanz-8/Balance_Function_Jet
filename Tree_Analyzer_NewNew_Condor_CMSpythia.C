#include "read_tree_New_CMSpythia.h"
#include "function_defination_CMSpythia.h" // call_libraries.h, input_variables.h and histogram_definition.h are called within function_defination.h
#include "coordinateTools.h" // for coordinate transformation w.r.t jet axis
#include "JetCorrector.h"
#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <iostream>

void Tree_Analyzer_NewNew_Condor_CMSpythia(TString input_file, int itxtoutFile, TString out_directory, TString name_of_roorFile)
{
  clock_t sec_start, sec_end;
  sec_start = clock();

  TDatime* date = new TDatime();

  print_start(); // start timing print

  // calling function to define array histograms
  array_hist_def();
  
  //added for CMSpythia: start
  jet_trigger = pp_jet_trigger;
  /*
  event_filter_str.resize(0);
  event_filter_str.push_back("pBeamScrapingFilter");
  event_filter_str.push_back("pPAprimaryVertexFilter");
  event_filter_str.push_back("HBHENoiseFilterResultRun2Loose");
  */
  jet_collection = pp_jet_collection;
  /*
  JEC_file = pp_JEC_file;
  
  vector<string> JECFiles;
  //JECFiles.push_back(Form("JEC_files/%s", JEC_file.c_str()));
  JECFiles.push_back(Form("JEC_files/%s", JEC_file.Data()));
  
  if(!is_MC)
    {
      JEC_file_data = pp_JEC_file_data;
      JECFiles.push_back(Form("JEC_files/%s", JEC_file_data.Data())); // for data only
      std::cout<<"It is runing for data (L2L3 residual file is attached)"<<std::endl;
      std::cout<<endl;
    }

  JetCorrector JEC(JECFiles);
  */
  //********* end ~~~~~~~~~~~~~~~~~

  TFile *f_mb;
  
  f_mb = TFile::Open(Form("TrkEffFiles/%s", fmb.Data()));
  th2f_eff=(TH2F*)f_mb->Get("rEff");
  th2f_fak=(TH2F*)f_mb->Get("rFak");
  th2f_sec=(TH2F*)f_mb->Get("rSec");
  
  fstream openInputFile;
  openInputFile.open(Form("%s",input_file.Data()), ios::in);
  if(!openInputFile.is_open())
    {
      cout << "List of input files not founded!" << endl;
      return;
    }
  
  std::vector<TString> file_name_vector;
  string file_chain;
  
  while(getline(openInputFile, file_chain))
    {
      //if(is_MC) file_name_vector.push_back(file_chain);
      //else file_name_vector.push_back(Form("root://cmsxrootd.fnal.gov/%s", file_chain.c_str()));
      file_name_vector.push_back(file_chain);
    }
  
  openInputFile.close();
  /*
  //added for CMSpythia: start
  TChain *hlt_tree = new TChain("hltanalysis/HltTree");
  TChain *hea_tree = new TChain("hiEvtAnalyzer/HiTree");
  TChain *ski_tree = new TChain("skimanalysis/HltTree");
  TChain *jet_tree = new TChain(Form("%s/t",jet_collection.Data()));
  TChain *trk_tree = new TChain("ppTrack/trackTree");
  TChain *gentrk_tree;
  if(is_MC) gentrk_tree = new TChain("HiGenParticleAna/hi");
  */


  TChain *jet_tree = new TChain(Form("%s/trackTree",jet_collection.Data()));
  //TTree *jet_tree = new TTree(Form("%s/trackTree",jet_collection.Data()));
  
  //********* end ~~~~~~~~~~~~~~~~~
  
  for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++)
    {
      TFile *testfile = TFile::Open(*listIterator);
      
      if(!testfile || testfile->IsZombie() || testfile->TestBit(TFile::kRecovered))
	{
	  cout << "File: " << *listIterator << " failed to open" << endl;
	  continue;
	}
      
      cout << "Adding file:--- " << *listIterator << "--- to the chains" << endl;
      jet_tree->Add(*listIterator);
      
      /*
      //added for CMSpythia: start 
      hlt_tree->Add(*listIterator);
      hea_tree->Add(*listIterator);
      ski_tree->Add(*listIterator);
      jet_tree->Add(*listIterator);
      trk_tree->Add(*listIterator);
      if(is_MC) gentrk_tree->Add(*listIterator);
      //********* end ~~~~~~~~~~~~~~~~~
      */
    }
  /*
  //added for CMSpythia: start
  hlt_tree->AddFriend(hea_tree);
  hlt_tree->AddFriend(ski_tree);
  hlt_tree->AddFriend(jet_tree);
  hlt_tree->AddFriend(trk_tree);
  if(is_MC) hlt_tree->AddFriend(gentrk_tree);
  */
  
  read_tree(jet_tree, jet_trigger.Data(), is_MC);
  //********* end ~~~~~~~~~~~~~~~~~

  int nevents = jet_tree->GetEntries(); // number of events
  
  cout<<endl;
  cout << "Total number of events in those files: "<< nevents << endl;
  cout<<endl;

  //gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  //gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

  //~~~~~~~~~~~~Define vectors use for Signal and Mixed event correlation
  // 1D event vector
  std::vector<float> pthatw_vec_1D;
  std::vector<int> processID_vec_1D;
  std::vector<float> vertexz_vec_1D;
  
  // 2D jet vector
  std::vector<std::vector<TVector3>> jet_vec_2D;
  std::vector<std::vector<double>> crtd_nch_trk_2D;
  // 3D trk vector
  std::vector<std::vector<std::vector<TVector3>>> jet_oldtrk_vec_3D;
  std::vector<std::vector<std::vector<TVector3>>> jet_newtrk_vec_3D;
  std::vector<std::vector<std::vector<TVector3>>> jet_newchtrk_vec_3D;
  std::vector<std::vector<std::vector<TVector3>>> jet_oldchtrk_vec_3D;
  std::vector<std::vector<std::vector<int>>> jet_trk_charge_3D;
  std::vector<std::vector<std::vector<double>>> weight_3D;
  int collected_events =0;

  for(int iev = 0; iev <nevents; iev++)
  //for(int iev = 0; iev < 10; iev++)
    {
      jet_tree->GetEntry(iev);
      //std::cout<<"ievt = "<<iev<<" & nref = "<<nref<<std::endl;
      
      
      //********************* All event cuts: start ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      
      //if(vertexz <= vz_cut_min || vertexz >= vz_cut_max) continue; // apply z vertex cut
      if(is_MC){
	if(pthat <= pthat_cut) continue; // apply pTHat cut
      }

      /*
      bool skimmed_evtfilter = false;
      for(int ii = 0; ii < (int) event_filter_str.size(); ii++)
        {
          if (event_filter_bool[ii] != 1) // condition for the skimmed event filters
	    {
              skimmed_evtfilter = true;
              break;
            }
        }
      */
      //std::cout<<"jet_trigger_bit:"<<jet_trigger_bit<<endl;
      //if(skimmed_evtfilter) continue; // apply the skimmed event filters                      
      //if(jet_trigger_bit!=1) continue; // apply jet trigger
      //if(jetN <= 0) continue; // if there is no jets in an event 
      //********************* All event cuts: end ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

      //~~~~~~~~~~~~Define 1D jet vectors use to fill 2D jet vector
      std::vector<TVector3> jet_vec_1D;
      std::vector<double> nch_weight_1D;
      //~~~~~~~~~~~~Define 2D trk vectors use to fill 3D trk vector
      std::vector<std::vector<TVector3>> oldtrk_vec_2D;
      std::vector<std::vector<TVector3>> newtrk_vec_2D;
      std::vector<std::vector<TVector3>> oldchtrk_vec_2D;
      std::vector<std::vector<TVector3>> newchtrk_vec_2D;
      std::vector<std::vector<int>> trk_charge_2D;
      std::vector<std::vector<double>> weight_2D;
       
      int nTotTrk_Jet = 0;
      
      //int jets_perevent = 0;
      //for (int ijet = 0; ijet < (int)njets; ijet++) // jet loop start
      for (int ijet = 0; ijet < dau_pt_STAR->size(); ijet++) // jet loop start
	{
	  //if(trackMax[ijet]/rawpt[ijet] < 0.01)continue; // Cut for jets for very low maxium pT track
	  //if(trackMax[ijet]/rawpt[ijet] > 0.98)continue; // Cut for jets where all the pT is taken by one track

	  float jet_pt = (*jetPt)[ijet];
	  float jet_eta = (*jetEta)[ijet];
	  float jet_phi = (*jetPhi)[ijet];
	  
	  //std::cout<<"Jet properties before JEC: "<<jet_pt_<<" "<<jet_eta<<" "<<jet_phi<<std::endl;
	  /*
	  JEC.SetJetPT(rawpt[ijet]);
	  JEC.SetJetEta(jet_eta);
	  JEC.SetJetPhi(jet_phi);
	  float jet_pt = JEC.GetCorrectedPT();
	  //float jet_pt = jet_pt_;
	  */
	  //std::cout<<"Jet properties after JEC: "<<jet_pt<<" "<<jet_eta<<" "<<jet_phi<<std::endl;
	  //std::cout<<" Automatic corrected jet pt = "<<jet_pt_<<"   &   "<<" Manually corrected jet pt = "<<jet_pt<<std::endl;

	  if (jet_pt < jet_pt_min) continue;  // jet pt cut
          if (fabs(jet_eta) >= jet_eta_max) continue;  // jet eta cut        
	  
	  nTotTrk_Jet++;

	  //push back to 1D jet vector
	  TVector3 saved_jet;
	  saved_jet.SetPtEtaPhi(jet_pt, jet_eta, jet_phi);
	  jet_vec_1D.push_back(saved_jet);
	  
	  //std::cout<<"Jet properties: "<<jet_pt<<" "<<jet_eta<<" "<<jet_phi<<std::endl;
	  //std::cout<<"Jet ncharge and Jet mult form tree: "<<JetnCh[ijet]<<"  "<<JetMult[ijet]<<std::endl;
	  
	  //~~~~~~~~~~~~Define 1D trk vectors use to fill 2D trk vector
	  std::vector<TVector3> oldtrk_vec_1D;
	  std::vector<TVector3> newtrk_vec_1D;
	  std::vector<TVector3> oldchtrk_vec_1D;
	  std::vector<TVector3> newchtrk_vec_1D;
	  std::vector<int> trk_charge_vec_1D;
	  std::vector<double> weight_1D;

	  //std::cout<<"High_purity Size:"<<(*highPurity)[ijet]<<endl;
	  //std::cout<<"Daughter size:"<<(*dau_pt_STAR)[ijet].size()<<endl;
	   
	  
	  int nchtrk_jet = 0, ntrk_jet = 0,nchtrk_pt03=0,nchtrk_pt0p33=0,nchtrk_pt0p53=0;
	  float dau_pt_sum = 0.0, dau_pt_STAR_sum =0.0,dau_pt_03_sum=0.0, dau_pt_0p33_sum=0.0,dau_pt_0p53_sum=0.0,dau_pt_03_STAR_sum =0.0,dau_pt_0p33_STAR_sum =0.0,dau_pt_0p53_STAR_sum =0.0;
	  double nch_weighted=0.0, weight=0.0;
	  //for (int itrk = 0; itrk < (int)JetMult[ijet]; itrk++) // jet track loop start
	  for (int itrk = 0; itrk < (*dau_pt_STAR)[ijet].size(); itrk++) // jet track loop start
	    {
	      float trk_pt = (*dau_pt)[ijet][itrk];
	      float trk_eta = (*dau_eta)[ijet][itrk];
	      float trk_phi = (*dau_phi)[ijet][itrk];
	      int trk_charge = (*dau_chg)[ijet][itrk];
	      float trk_pterr = (*dau_ptError)[ijet][itrk];
	      float new_trk_pt=(*dau_pt_STAR)[ijet][itrk];
	      float new_trk_eta=(*dau_eta_STAR)[ijet][itrk];
	      float new_trk_phi=(*dau_phi_STAR)[ijet][itrk];
	      float dauXYDCAsig=(*dau_XYDCAsig)[ijet][itrk];
	      float dauZDCAsig=(*dau_ZDCAsig)[ijet][itrk];
	      //bool trk_hp = (bool)(*highPurity)[itrk];
		/*
	      float trk_dxy = (float)trkDxy1[itrk];
	      float trk_dxyerr = (float)trkDxyError1[itrk];
	      float trk_dz = (float)trkDz1[itrk];
	      float trk_dzerr = (float)trkDzError1[itrk];
	      bool trk_hp = (bool)highPurity[itrk];
	      
	      if(trk_hp != 1) continue;
	      if(trk_charge == 0) continue;
	      if(fabs(trk_pterr) / trk_pt >= 0.1 ) continue;
		*/
	      if(fabs(dauXYDCAsig) >= 3.0 ) continue;
	      if(fabs(dauZDCAsig) >= 3.0 ) continue;
	      
	      if(trk_charge == 0) continue;
	      if(fabs(trk_pterr) / trk_pt >= 0.1 ) continue;
	      //std::cout<<" highPurity = "<<trk_hp<<std::endl;

	      //nTotTrk_Jet++;
	      
	      // to check whether the partcile is inside the the jet or not
	      
	      float DEta_jet_trk = (jet_eta - trk_eta);
	      float DPhi_jet_trk = TVector2::Phi_mpi_pi(jet_phi - trk_phi);
	      float DR_jet_trk = TMath::Sqrt(pow(DEta_jet_trk,2) + pow(DPhi_jet_trk,2));

	      h_pt_before->Fill(trk_pt);


	      if(DR_jet_trk > jet_radius)
		{
		  //std::cout<<"track is not inside jet cone: "<<trk_pt<<" "<<trk_eta<<" "<<trk_phi<<" "<<trk_charge<<std::endl;
		  continue;
		}

	       h_pt_after->Fill(trk_pt);

	      
	      
	      ntrk_jet++;
	      
	      if(fabs((int)trk_charge) > 0)
		{
		  nchtrk_jet++;
		  dau_pt_sum +=trk_pt;
		  dau_pt_STAR_sum +=new_trk_pt;

		  if(0<=trk_pt && trk_pt<3)
		    {
		      nchtrk_pt03++;
		      dau_pt_03_sum+=trk_pt;
		      dau_pt_03_STAR_sum+=new_trk_pt;
		      
		    }
		  if(0.3<=trk_pt && trk_pt<3)
		    {
		      nchtrk_pt0p33++;
                      dau_pt_0p33_sum+=trk_pt;
                      dau_pt_0p33_STAR_sum+=new_trk_pt;
		    }
		  if(0.5<=trk_pt && trk_pt<3)
                    {
                      nchtrk_pt0p53++;
                      dau_pt_0p53_sum+=trk_pt;
                      dau_pt_0p53_STAR_sum+=new_trk_pt;
                    }

		}
	      /*

	      if( th2f_eff->GetBinContent(th2f_eff->FindBin(new_trk_eta,new_trk_pt)) != 0. )
		{
		  weight = (((1.0 - th2f_fak->GetBinContent(th2f_fak->FindBin(new_trk_eta,new_trk_pt)))*(1.0 - th2f_sec->GetBinContent(th2f_sec->FindBin(new_trk_eta,new_trk_pt))))/(th2f_eff->GetBinContent(th2f_eff->FindBin(new_trk_eta,new_trk_pt))));
		}
	      else
		{
                 weight=1;
                }
	      */
	      if( th2f_eff->GetBinContent(th2f_eff->FindBin(trk_eta,trk_pt)) != 0. )
                {
                  weight = (((1.0 - th2f_fak->GetBinContent(th2f_fak->FindBin(trk_eta,trk_pt)))*(1.0 - th2f_sec->GetBinContent(th2f_sec->FindBin(trk_eta,trk_pt))))/(th2f_eff->GetBinContent(th2f_eff->FindBin(trk_eta,trk_pt))));
                }
              else
                {
		  weight=1;
		}


	      
	      /*
	      float new_trk_pt = ptWRTJet(jet_pt, jet_eta, jet_phi, trk_pt, trk_eta, trk_phi);
	      float new_trk_eta = etaWRTJet(jet_pt, jet_eta, jet_phi, trk_pt, trk_eta, trk_phi);
	      float new_trk_phi = phiWRTJet(jet_pt, jet_eta, jet_phi, trk_pt, trk_eta, trk_phi);
	      */
	      //std::cout<<trk_pt<<" "<<trk_eta<<" "<<trk_phi<<"    "<<new_trk_pt<<" "<<new_trk_eta<<" "<<new_trk_phi<<std::endl;
	      
	      //push back to 1D trk vector
	      trk_charge_vec_1D.push_back(trk_charge);
	     
              nch_weighted+=1*weight; 
	      weight_1D.push_back(1*weight);

                
	      TVector3 saved_oldtrk, saved_newtrk, saved_oldchtrk, saved_newchtrk;
	      saved_oldtrk.SetPtEtaPhi(trk_pt, trk_eta, trk_phi);
	      oldtrk_vec_1D.push_back(saved_oldtrk);	      
	      saved_newtrk.SetPtEtaPhi(new_trk_pt, new_trk_eta, new_trk_phi);
	      newtrk_vec_1D.push_back(saved_newtrk);
	     
	      if(fabs((int)trk_charge) > 0)
		{
		  saved_oldchtrk.SetPtEtaPhi(trk_pt, trk_eta, trk_phi);
		  oldchtrk_vec_1D.push_back(saved_oldchtrk);
		  
		  saved_newchtrk.SetPtEtaPhi(new_trk_pt, new_trk_eta, new_trk_phi);
		  newchtrk_vec_1D.push_back(saved_newchtrk);
		  
		}
	    } // jet track loop end
	  
	  //if(trk_charge_vec_1D.size() == 0) continue;

	  tp1d_dau_meanpT_nchjet->Fill(nchtrk_jet, dau_pt_sum/(float)nchtrk_jet);
	  tp1d_dau_STAR_meanpT_nchjet->Fill(nchtrk_jet, dau_pt_STAR_sum/(float)nchtrk_jet);
	  
	  tp1d_dau_meanpT_03_nchjet->Fill(nchtrk_pt03, dau_pt_03_sum/(float)nchtrk_pt03);
	  tp1d_dau_STAR_meanpT_03_nchjet->Fill(nchtrk_pt03, dau_pt_03_STAR_sum/(float)nchtrk_pt03);

	  tp1d_dau_meanpT_0p33_nchjet->Fill(nchtrk_pt0p33, dau_pt_0p33_sum/(float)nchtrk_pt0p33);
          tp1d_dau_STAR_meanpT_0p33_nchjet->Fill(nchtrk_pt0p33, dau_pt_0p33_STAR_sum/(float)nchtrk_pt0p33);

	  tp1d_dau_meanpT_0p53_nchjet->Fill(nchtrk_pt0p53, dau_pt_0p53_sum/(float)nchtrk_pt0p53);
          tp1d_dau_STAR_meanpT_0p53_nchjet->Fill(nchtrk_pt0p53, dau_pt_0p53_STAR_sum/(float)nchtrk_pt0p53);
         
          nch_weight_1D.push_back(nch_weighted);

	  

	  //push back to 2D trk vector
	  oldtrk_vec_2D.push_back(oldtrk_vec_1D);
	  newtrk_vec_2D.push_back(newtrk_vec_1D);
	  oldchtrk_vec_2D.push_back(oldchtrk_vec_1D);
	  newchtrk_vec_2D.push_back(newchtrk_vec_1D);
	  trk_charge_2D.push_back(trk_charge_vec_1D);
          weight_2D.push_back(weight_1D);

	  	  
	  //clear 1D trk vector
	  oldtrk_vec_1D.clear();
	  newtrk_vec_1D.clear();
	  oldchtrk_vec_1D.clear();
	  newchtrk_vec_1D.clear();
	  trk_charge_vec_1D.clear();
          weight_1D.clear();
	  
	} // jet loop end

      if(jet_vec_1D.size() == 0) continue;
      collected_events++;
      
      //std::cout<<"ievt = "<<iev<<" & Jets after cut  = "<<nTotTrk_Jet<<std::endl;

      //push back to 2D jet vector
      jet_vec_2D.push_back(jet_vec_1D);
      crtd_nch_trk_2D.push_back(nch_weight_1D);

      
      //clear 1D jet vector
      jet_vec_1D.clear();
      nch_weight_1D.clear();

      
      //push back to 3D trk vector
      jet_oldtrk_vec_3D.push_back(oldtrk_vec_2D);
      jet_newtrk_vec_3D.push_back(newtrk_vec_2D);
      jet_oldchtrk_vec_3D.push_back(oldchtrk_vec_2D);
      jet_newchtrk_vec_3D.push_back(newchtrk_vec_2D);
      jet_trk_charge_3D.push_back(trk_charge_2D);
      weight_3D.push_back(weight_2D);
      //std::cout<<" Size of jet  ="<< jet_newtrk_vec_3D.size()<<std::endl;
      
      //clear 2D trk vector
      oldtrk_vec_2D.clear();
      newtrk_vec_2D.clear();
      oldchtrk_vec_2D.clear();
      newchtrk_vec_2D.clear();
      trk_charge_2D.clear();
      weight_2D.clear();
      
      if(!is_MC) pthatweight = 1.0;
      int Subprocesscode = 1;
      //pthatweight = 1.; // for test
      
      //push back event quatities (pthatweight, Subprocesscode, vertexz, etc) before closing event loop
      pthatw_vec_1D.push_back(pthatweight);
      processID_vec_1D.push_back(Subprocesscode);
      //vertexz_vec_1D.push_back(vertexz);

      //fill event histograms
      //hntrkTotal->Fill(ntrkTotal, pthatweight);
      //hpthat->Fill(pthat, pthatweight);
      //hpthatw->Fill(pthatweight);
      //hsubprocessid->Fill(Subprocesscode, pthatweight);
      //hvertexz->Fill(vertexz);
      //h2D_run_lumi->Fill(run,lumi);
      //hnjets->Fill(njets, pthatweight);
    }
 // event loop end
  
  //std::cout<<" Total collected events after applying cuts is  = "<<std::endl;//collected_events<<std::endl;

  //return;

  //calling function for signal correlation
  signal_corr(pthatw_vec_1D, processID_vec_1D, jet_vec_2D, jet_oldtrk_vec_3D, jet_newtrk_vec_3D, jet_oldchtrk_vec_3D, jet_newchtrk_vec_3D, jet_trk_charge_3D,crtd_nch_trk_2D,weight_3D);
  //signal_corr(pthatw_vec_1D, processID_vec_1D, jet_vec_2D, jet_oldtrk_vec_3D, jet_oldchtrk_vec_3D, jet_trk_charge_3D);
  
  //calling function for mixing correlation  
  mixing_corr(pthatw_vec_1D, processID_vec_1D, jet_vec_2D, jet_oldtrk_vec_3D, jet_newtrk_vec_3D, jet_oldchtrk_vec_3D, jet_newchtrk_vec_3D, jet_trk_charge_3D,crtd_nch_trk_2D,weight_3D);
  // mixing_corr(pthatw_vec_1D, processID_vec_1D, jet_vec_2D, jet_oldtrk_vec_3D, jet_oldchtrk_vec_3D, jet_trk_charge_3D);
  //std::cout<<"   "<<jet_vec_2D

  //return;

  // clear event vector
  pthatw_vec_1D.clear();
  processID_vec_1D.clear();
  //vertexz_vec_1D.clear();

  //clear jet vector
  jet_vec_2D.clear();
  crtd_nch_trk_2D.clear();
  
  //clear trk vectorx
  jet_oldtrk_vec_3D.clear();
  jet_newtrk_vec_3D.clear();
  jet_oldchtrk_vec_3D.clear();
  jet_newchtrk_vec_3D.clear();
  jet_trk_charge_3D.clear();
  weight_3D.clear();


  
  
  /*
  // delete the tree to clean memory                                                                                    
  delete hlt_tree;
  delete hea_tree;
  delete ski_tree;
  delete jet_tree;
  delete trk_tree;
  if(is_MC)delete gentrk_tree;
  */
  delete jet_tree;
  // output file name
  std::string outfilename = Form("%s/%s_Outfile_%d", out_directory.Data(), name_of_roorFile.Data(), itxtoutFile);
  std::replace(outfilename.begin(), outfilename.end(), '.', 'p'); // replace . to p
  std::replace(outfilename.begin(), outfilename.end(), '-', 'N'); // replace - to N for negative
  
  TFile* fout = new TFile(Form("%s.root", outfilename.c_str()), "recreate");
  
  fout->mkdir("Event_Hist");
  fout->cd("Event_Hist");
  write_event_hist();
  
  fout->mkdir("Jet_Hist");
  fout->cd("Jet_Hist");
  write_jet_hist();
  
  fout->mkdir("Track_Hist");
  fout->cd("Track_Hist");
  write_track_hist();
  
  fout->mkdir("Corr_Hist");
  fout->cd("Corr_Hist");
  write_corr_hist();

  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~Finished, DONE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<<std::endl;
  
  sec_end = clock(); // stop time counting                                                                                       
  cout << "========================================" << endl;
  cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
  cout << "========================================" << endl;
  
  print_stop(); // Print time, date and hour when it stops
}

int main(int argc, char **argv)
{
  using namespace std;

  TString inputfile;
  int itxtout;
  TString outfile;
  TString outfile_name;
  if(argc == 5)
    {
      std::cout<<std::endl;
      std::cout<<"You have given "<< argc <<" arguments including the program name;  Your program will run"<<std::endl;
      std::cout<<std::endl;
      inputfile = argv[1];
      itxtout = atoi(argv[2]);
      outfile = argv[3];
      outfile_name  = argv[4];

      Tree_Analyzer_NewNew_Condor_CMSpythia(inputfile, itxtout, outfile, outfile_name);
    }
  return 0;
}
