#include "input_variables_CMSpythia.h"
#include "histogram_definition_CMSpythia.h"

void signal_corr(std::vector<float> pthatw_vec, std::vector<int> processID_vec, std::vector<std::vector<TVector3>> jet_vec, std::vector<std::vector<std::vector<TVector3>>> jet_oldtrk_vec, std::vector<std::vector<std::vector<TVector3>>> jet_newtrk_vec, std::vector<std::vector<std::vector<TVector3>>> jet_oldchtrk_vec, std::vector<std::vector<std::vector<TVector3>>> jet_newchtrk_vec, std::vector<std::vector<std::vector<int>>> jet_trk_charge,std::vector<std::vector<double>> nch_weight,std::vector<std::vector<std::vector<double>>>weight_3D)
{
  if(jet_newtrk_vec.size() != jet_vec.size() || jet_newtrk_vec.size() != processID_vec.size())
    {
      std::cout<<"event numbers are not same, pleaes check"<<std::endl;
    }
  
  else std::cout<<"Total events is: "<<jet_newtrk_vec.size()<<"  "<<jet_vec.size()<<"  "<<processID_vec.size()<<std::endl;

  std::cout<<endl;
  std::cout<<"~~~~~start signal correlation~~~~~~~~"<<std::endl;
  
  for(int ievt = 0; ievt < (int)jet_newtrk_vec.size(); ievt++) // event loop
    {
      if(ievt%10000 == 0) std::cout<<ievt<<"  events running for signal correlation of total events "<< jet_newtrk_vec.size() <<std::endl;
      //if(ievt%1 == 0) std::cout<<ievt<<"  events running for signal correlation of total events "<< jet_newtrk_vec.size() <<std::endl;
      
      double pthatw = pthatw_vec[ievt];

      //int subprocessID = processID_vec[ievt];

      int njets = 0;
      
      //std::cout<<" jet size =  "<<(int)jet_newtrk_vec[ievt].size()<<std::endl;

      for(int ijet = 0; ijet < (int)jet_newtrk_vec[ievt].size(); ijet++) // jet loop
	{
	  TVector3 jetVec = jet_vec[ievt][ijet];
	  double jet_pt = jetVec.Pt();
	  double jet_eta = jetVec.Eta();
	  double jet_phi = jetVec.Phi();
	  int ntrk_jet = jet_newtrk_vec[ievt][ijet].size();
	  //int nchtrk_jet = jet_newchtrk_vec[ievt][ijet].size(); // number of charge without correction
          double nchtrk_jet = nch_weight[ievt][ijet]; //// number of charge with correction

	  // extract the nchtrk bin

	  int nchbin  = hnch_bin_hist->FindBin(nchtrk_jet) - 1;

	  //apply jet pt, eta cut
	  //if(jet_pt < jet_pt_min || TMath::Abs(jet_eta) > jet_eta_max) continue;
	  //if(jet_pt < 120 || TMath::Abs(jet_eta) > jet_eta_max) continue;

	  //fill jet histograms
	  /*
	  hjet_pt->Fill(jet_pt, pthatw);
	  hjet_eta->Fill(jet_eta, pthatw);
	  hjet_phi->Fill(jet_phi, pthatw);
	  hJet_ntrk->Fill(ntrk_jet, pthatw);
	  hJet_nchtrk->Fill(nchtrk_jet, pthatw);
	  hJet_nchtrk_pT->Fill(nchtrk_jet, jet_pt, pthatw);
	  hJet_nchtrk_[nchbin]->Fill(nchtrk_jet, pthatw);
	  */
	  hjet_pt->Fill(jet_pt);                                                                                                                                                                  
          hjet_eta->Fill(jet_eta);  
          hjet_phi->Fill(jet_phi);  
          hJet_ntrk->Fill(ntrk_jet);                                                                                                                                                              
          hJet_nchtrk->Fill(nchtrk_jet);                                                           
          hJet_nchtrk_pT->Fill(nchtrk_jet, jet_pt);                                                                                                                                               
          hJet_nchtrk_[nchbin]->Fill(nchtrk_jet); 
	  nch_avg_[nchbin]->Fill(nchtrk_jet);

	  
	  //std::cout<<" trk size =  "<<(int)jet_newtrk_vec[ievt][ijet].size()<<std::endl;

	  for(int itrk = 0; itrk < (int)jet_newtrk_vec[ievt][ijet].size(); itrk++) // trigger trk loop
	    {
	      
	      TVector3 old_trgtrk_vec = jet_oldtrk_vec[ievt][ijet][itrk];
	      double old_trgtrk_pt = old_trgtrk_vec.Pt();
	      double old_trgtrk_eta = old_trgtrk_vec.Eta();
	      double old_trgtrk_phi = old_trgtrk_vec.Phi();
	      
	      
	      TVector3 new_trgtrk_vec = jet_newtrk_vec[ievt][ijet][itrk];
	      double new_trgtrk_pt = new_trgtrk_vec.Pt();
	      double new_trgtrk_eta = new_trgtrk_vec.Eta();
	      double new_trgtrk_phi = new_trgtrk_vec.Phi();
	      int trgtrk_charge = jet_trk_charge[ievt][ijet][itrk];
	      //double weight_trig=weight_3D[ievt][ijet][itrk];
	      //pthatw=1.0;
	      double weight_trig=1.0;
	      //if(new_trgtrk_pt > newtrk_jt03_max) continue; // trg pt cut

	      if(TMath::Abs(trgtrk_charge) > 0 )
		{
		  if(TMath::Abs(new_trgtrk_eta) <= newtrk_eta_max)
		    {
		      
		      hchtrk_pt->Fill(old_trgtrk_pt, pthatw*weight_trig);
		      hchtrk_eta->Fill(old_trgtrk_eta, pthatw*weight_trig);
		      hchtrk_phi->Fill(old_trgtrk_phi, pthatw*weight_trig);
		      
		      hchtrk_pt_jetaxis->Fill(new_trgtrk_pt, pthatw*weight_trig);
		      hchtrk_eta_jetaxis->Fill(new_trgtrk_eta, pthatw*weight_trig);
		      hchtrk_phi_jetaxis->Fill(new_trgtrk_phi, pthatw*weight_trig);
		      
		      if(trgtrk_charge > 0)
			{
			  
			  hchtrk_pt_p->Fill(old_trgtrk_pt, pthatw*weight_trig);
			  hchtrk_eta_p->Fill(old_trgtrk_eta, pthatw*weight_trig);
			  hchtrk_phi_p->Fill(old_trgtrk_phi, pthatw*weight_trig);
			  
			  hchtrk_pt_jetaxis_p->Fill(new_trgtrk_pt, pthatw*weight_trig);
			  hchtrk_eta_jetaxis_p->Fill(new_trgtrk_eta, pthatw*weight_trig);
			  hchtrk_phi_jetaxis_p->Fill(new_trgtrk_phi, pthatw*weight_trig);
			}
		      else if(trgtrk_charge < 0)
			{
			  
			  hchtrk_pt_m->Fill(old_trgtrk_pt, pthatw*weight_trig);
			  hchtrk_eta_m->Fill(old_trgtrk_eta, pthatw*weight_trig);
			  hchtrk_phi_m->Fill(old_trgtrk_phi, pthatw*weight_trig);
			  
			  hchtrk_pt_jetaxis_m->Fill(new_trgtrk_pt, pthatw*weight_trig);
			  hchtrk_eta_jetaxis_m->Fill(new_trgtrk_eta, pthatw*weight_trig);
			  hchtrk_phi_jetaxis_m->Fill(new_trgtrk_phi, pthatw*weight_trig);
			}
		      
		      if(new_trgtrk_pt >= newtrk_jt03_min && new_trgtrk_pt <= newtrk_jt03_max)
			{
			  hnchtrk_jt03->AddBinContent(1, pthatw*weight_trig);
			  hnchtrk_jt03_[nchbin]->AddBinContent(1, pthatw*weight_trig);
			  
			  if(trgtrk_charge > 0)
			    {
			      hnchtrk_jt03_pm->AddBinContent(1, pthatw*weight_trig);
			      hnchtrk_jt03_pm_[nchbin]->AddBinContent(1, pthatw*weight_trig);
			    }
			  else if(trgtrk_charge < 0)
			    {
			      hnchtrk_jt03_pm->AddBinContent(2, pthatw*weight_trig);
			      hnchtrk_jt03_pm_[nchbin]->AddBinContent(2, pthatw*weight_trig);
			    }
			}
		      if(new_trgtrk_pt >= newtrk_jt0p33_min && new_trgtrk_pt <= newtrk_jt0p33_max)
			{
			  hnchtrk_jt0p33->AddBinContent(1, pthatw*weight_trig);
			  hnchtrk_jt0p33_[nchbin]->AddBinContent(1, pthatw*weight_trig);
			  
			  if(trgtrk_charge > 0)
                            {
                              hnchtrk_jt0p33_pm->AddBinContent(1, pthatw*weight_trig);
                              hnchtrk_jt0p33_pm_[nchbin]->AddBinContent(1, pthatw*weight_trig);
                            }
                          else if(trgtrk_charge < 0)
                            {
                              hnchtrk_jt0p33_pm->AddBinContent(2, pthatw*weight_trig);
                              hnchtrk_jt0p33_pm_[nchbin]->AddBinContent(2, pthatw*weight_trig);
			    }
			}
		      if(new_trgtrk_pt >= newtrk_jt0p53_min && new_trgtrk_pt <= newtrk_jt0p53_max)
			{
			  hnchtrk_jt0p53->AddBinContent(1, pthatw*weight_trig);
			  hnchtrk_jt0p53_[nchbin]->AddBinContent(1, pthatw*weight_trig);

			  if(trgtrk_charge > 0)
                            {
                              hnchtrk_jt0p53_pm->AddBinContent(1, pthatw*weight_trig);
                              hnchtrk_jt0p53_pm_[nchbin]->AddBinContent(1, pthatw*weight_trig);
                            }
                          else if(trgtrk_charge < 0)
                            {
                              hnchtrk_jt0p53_pm->AddBinContent(2, pthatw*weight_trig);
                              hnchtrk_jt0p53_pm_[nchbin]->AddBinContent(2, pthatw*weight_trig);
                            }
			}
		    } // trk eta condition
		} // charge condition
	      
	      for(int jtrk = 0; jtrk < (int)jet_newtrk_vec[ievt][ijet].size(); jtrk++) // associate trk loop
		{
		  if(itrk == jtrk) continue; // to avoid auto correlation
		  
		  TVector3 old_assotrk_vec = jet_oldtrk_vec[ievt][ijet][jtrk];
		  double old_assotrk_pt = old_assotrk_vec.Pt();
		  double old_assotrk_eta = old_assotrk_vec.Eta();
		  double old_assotrk_phi = old_assotrk_vec.Phi();
		  
		  
		  TVector3 new_assotrk_vec = jet_newtrk_vec[ievt][ijet][jtrk];
		  double new_assotrk_pt = new_assotrk_vec.Pt();
		  double new_assotrk_eta = new_assotrk_vec.Eta();
		  double new_assotrk_phi = new_assotrk_vec.Phi();
		  int assotrk_charge = jet_trk_charge[ievt][ijet][jtrk];
				  
		  //if(new_assotrk_pt > newtrk_jt03_max) continue; // asso pt cut
		  
		  double deta_jetaxis = new_trgtrk_eta - new_assotrk_eta;
		  double dphi_jetaxis = new_trgtrk_phi - new_assotrk_phi;
		  double dphi2_jetaxis = new_assotrk_phi - new_trgtrk_phi;
		  double weight_ass=weight_3D[ievt][ijet][jtrk];
		  //pthatw=weight_trig*weight_ass;

		  if(dphi_jetaxis > 1.5*TMath::Pi())
		    {
		      dphi_jetaxis = dphi_jetaxis - 2.0*TMath::Pi();
		    }
		  else if(dphi_jetaxis < -0.5*TMath::Pi())
		    {
		      dphi_jetaxis = dphi_jetaxis + 2.0*TMath::Pi();
		    }
		  
		  if(dphi2_jetaxis > 1.5*TMath::Pi())
		    {
		      dphi2_jetaxis = dphi2_jetaxis - 2.0*TMath::Pi();
		    }
		  else if(dphi2_jetaxis < -0.5*TMath::Pi())
		    {
		      dphi2_jetaxis = dphi2_jetaxis + 2.0*TMath::Pi();
		    }

		  if(TMath::Abs(trgtrk_charge) > 0 && TMath::Abs(assotrk_charge) > 0)
		    {
		      if(TMath::Abs(new_trgtrk_eta) <= newtrk_eta_max && TMath::Abs(new_assotrk_eta) <= newtrk_eta_max)
			{
			  if((new_trgtrk_pt >= newtrk_jt03_min && new_trgtrk_pt <= newtrk_jt03_max) && (new_assotrk_pt >= newtrk_jt03_min && new_assotrk_pt <= newtrk_jt03_max))
			    {      
			      hsignal_newchtrk_jt03->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
			      hsignal_newchtrk_jt03->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
			      hsignal_newchtrk_jt03->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
			      hsignal_newchtrk_jt03->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
			      
			      hsignal_newchtrk_jt03_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt03_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt03_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt03_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

			      if(trgtrk_charge > 0 && assotrk_charge > 0)
				{
				  hsignal_newchtrk_jt03_pp->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				  hsignal_newchtrk_jt03_pp->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				  hsignal_newchtrk_jt03_pp->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				  hsignal_newchtrk_jt03_pp->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

				  hsignal_newchtrk_jt03_pp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_pp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_pp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_pp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				}
			      if(trgtrk_charge < 0 && assotrk_charge < 0)
				{
				  hsignal_newchtrk_jt03_mm->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_mm->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_mm->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_mm->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

				  hsignal_newchtrk_jt03_mm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_mm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_mm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_mm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				}
			      if(trgtrk_charge > 0 && assotrk_charge < 0)
				{
				  hsignal_newchtrk_jt03_pm->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_pm->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_pm->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_pm->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

				  hsignal_newchtrk_jt03_pm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_pm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_pm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_pm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				}
			      if(trgtrk_charge < 0 && assotrk_charge > 0)
				{
				  hsignal_newchtrk_jt03_mp->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_mp->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_mp->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_mp->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

				  hsignal_newchtrk_jt03_mp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_mp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_mp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt03_mp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				}
			    }
			  if((new_trgtrk_pt >= newtrk_jt0p33_min && new_trgtrk_pt <= newtrk_jt0p33_max) && (new_assotrk_pt >= newtrk_jt0p33_min && new_assotrk_pt <= newtrk_jt0p33_max))
			    {
			      hsignal_newchtrk_jt0p33->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt0p33->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt0p33->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt0p33->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

			      hsignal_newchtrk_jt0p33_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt0p33_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt0p33_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt0p33_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

			      if(trgtrk_charge > 0 && assotrk_charge > 0)
				{
				  hsignal_newchtrk_jt0p33_pp->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				  hsignal_newchtrk_jt0p33_pp->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				  hsignal_newchtrk_jt0p33_pp->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				  hsignal_newchtrk_jt0p33_pp->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

				  hsignal_newchtrk_jt0p33_pp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_pp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_pp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_pp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				}
			      if(trgtrk_charge < 0 && assotrk_charge < 0)
				{
				  hsignal_newchtrk_jt0p33_mm->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_mm->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_mm->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_mm->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

				  hsignal_newchtrk_jt0p33_mm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_mm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_mm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_mm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				}
			      if(trgtrk_charge > 0 && assotrk_charge < 0)
				{
				  hsignal_newchtrk_jt0p33_pm->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_pm->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_pm->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_pm->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

				  hsignal_newchtrk_jt0p33_pm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_pm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_pm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_pm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				}
			      if(trgtrk_charge < 0 && assotrk_charge > 0)
				{
				  hsignal_newchtrk_jt0p33_mp->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_mp->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_mp->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_mp->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

				  hsignal_newchtrk_jt0p33_mp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_mp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_mp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p33_mp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				}
			    }
			  if((new_trgtrk_pt >= newtrk_jt0p53_min && new_trgtrk_pt <= newtrk_jt0p53_max) && (new_assotrk_pt >= newtrk_jt0p53_min && new_assotrk_pt <= newtrk_jt0p53_max))
			    {
			      hsignal_newchtrk_jt0p53->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt0p53->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt0p53->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt0p53->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

			      hsignal_newchtrk_jt0p53_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt0p53_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt0p53_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                              hsignal_newchtrk_jt0p53_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

			      if(trgtrk_charge > 0 && assotrk_charge > 0)
				{
				  hsignal_newchtrk_jt0p53_pp->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				  hsignal_newchtrk_jt0p53_pp->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				  hsignal_newchtrk_jt0p53_pp->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				  hsignal_newchtrk_jt0p53_pp->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

				  hsignal_newchtrk_jt0p53_pp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_pp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_pp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_pp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				}
			      if(trgtrk_charge < 0 && assotrk_charge < 0)
				{
				  hsignal_newchtrk_jt0p53_mm->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_mm->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_mm->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_mm->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

				  hsignal_newchtrk_jt0p53_mm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_mm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_mm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_mm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				}
			      if(trgtrk_charge > 0 && assotrk_charge < 0)
				{
				  hsignal_newchtrk_jt0p53_pm->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_pm->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_pm->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_pm->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

				  hsignal_newchtrk_jt0p53_pm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_pm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_pm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_pm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				}
			      if(trgtrk_charge < 0 && assotrk_charge > 0)
				{
				  hsignal_newchtrk_jt0p53_mp->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_mp->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_mp->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_mp->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);

				  hsignal_newchtrk_jt0p53_mp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_mp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_mp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
                                  hsignal_newchtrk_jt0p53_mp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				}
			    }			      
			} // eta condition
		    } // charge condition

		} // asso trk loop
	    } // trg trk loop end
	  njets++;
	} // jet loop end
      if(njets > 0)
	{
	  hnjets_afterCut->Fill(njets);
	}
    } // event loop end
  std::cout<<endl;
  std::cout<<"~~~~~end signal correlation~~~~~~~~"<<std::endl;
} // function loop

void mixing_corr(std::vector<float> pthatw_vec, std::vector<int> processID_vec, std::vector<std::vector<TVector3>> jet_vec, std::vector<std::vector<std::vector<TVector3>>> jet_oldtrk_vec, std::vector<std::vector<std::vector<TVector3>>> jet_newtrk_vec, std::vector<std::vector<std::vector<TVector3>>> jet_oldchtrk_vec, std::vector<std::vector<std::vector<TVector3>>> jet_newchtrk_vec, std::vector<std::vector<std::vector<int>>> jet_trk_charge,std::vector<std::vector<double>> nch_weight,std::vector<std::vector<std::vector<double>>> weight_3D)
{
  std::cout<<endl;
  std::cout<<"~~~~~start mixed event correlation~~~~~~~~"<<std::endl;
  
  for(int ievt = 0; ievt < (int)jet_newtrk_vec.size(); ievt++) // 1st event loop
    {
      if(ievt%10000 == 0) std::cout<<ievt<<"  events running for mixing correlation of total events "<< jet_newtrk_vec.size() <<std::endl;
      //if(ievt%1 == 0) std::cout<<ievt<<"  events running for mixing correlation of total events "<< jet_newtrk_vec.size() <<std::endl;
      
      double pthatw = pthatw_vec[ievt];
      //int subprocessID = processID_vec[ievt];
      //float ivertexz = vertexz[ievt];

      // mixing algorithm
      int mixstart = ievt+1;
      int mixend   = (int)jet_newtrk_vec.size();
      
      if(mixstart > (0.7*(jet_newtrk_vec.size())))
	{
	  int mixstart = (0.5*(jet_newtrk_vec.size()));
	  int mixend   = (int)jet_newtrk_vec.size();
	}
      
      int nmix = 0;
      for(int jevt = mixstart; jevt < mixend; jevt++) // 2nd event loop
	{
	  if(ievt == jevt) continue;
	  //float jvertexz = vertexz[jevt];
	  
	  //if (fabs(jvertexz - ivertexz) > 2.0) continue;
	  
	  nmix++;
	  if(nmix > bkgFactor) break;
	  
	  for(int ijet = 0; ijet < (int)jet_newtrk_vec[ievt].size(); ijet++) // 1st jet loop
	    {
	      TVector3 jetVec = jet_vec[ievt][ijet];
	      double jet_pt = jetVec.Pt();
	      double jet_eta = jetVec.Eta();
	      double jet_phi = jetVec.Phi();
	      //int nchtrk_jet = jet_newchtrk_vec[ievt][ijet].size(); // Without Correction
	      int ntrk_jet = jet_newtrk_vec[ievt][ijet].size();
	      double nchtrk_jet = nch_weight[ievt][ijet];//With Correction
	      int nchbin  = hnch_bin_hist->FindBin(nchtrk_jet) - 1;
	      
	      //apply jet cut
	      //if(jet_pt < jet_pt_min || TMath::Abs(jet_eta) > jet_eta_max) continue;
	      
	      for(int jjet = 0; jjet < (int)jet_newtrk_vec[jevt].size(); jjet++) // 2nd jet loop
		{
		  TVector3 jetVec2 = jet_vec[jevt][jjet];
		  double jet_pt2 = jetVec.Pt();
		  double jet_eta2 = jetVec.Eta();
		  double jet_phi2 = jetVec.Phi();
		  int nchtrk_jet2 = jet_newchtrk_vec[jevt][jjet].size();
		  int ntrk_jet2 = jet_newtrk_vec[jevt][jjet].size();
		  
		  // apply jet cut
		  //if(jet_pt2 < jet_pt_min || TMath::Abs(jet_eta2) > jet_eta_max) continue;
		  
		  for(int itrk = 0; itrk < (int)jet_newtrk_vec[ievt][ijet].size(); itrk++) // trigger trk loop
		    {
		      TVector3 old_trgtrk_vec = jet_oldtrk_vec[ievt][ijet][itrk];
		      double old_trgtrk_pt = old_trgtrk_vec.Pt();
		      double old_trgtrk_eta = old_trgtrk_vec.Eta();
		      double old_trgtrk_phi = old_trgtrk_vec.Phi();
		      TVector3 new_trgtrk_vec = jet_newtrk_vec[ievt][ijet][itrk];
		      double new_trgtrk_pt = new_trgtrk_vec.Pt();
		      double new_trgtrk_eta = new_trgtrk_vec.Eta();
		      double new_trgtrk_phi = new_trgtrk_vec.Phi();
		      int trgtrk_charge = jet_trk_charge[ievt][ijet][itrk];
		      double weight_trig=weight_3D[ievt][ijet][itrk];		      
		      //if(new_trgtrk_pt > newtrk_jt03_max) continue; // trg pt cut
		      
		      for(int jtrk = 0; jtrk < (int)jet_newtrk_vec[jevt][jjet].size(); jtrk++) // associate trk loop
			{
			  TVector3 old_assotrk_vec = jet_oldtrk_vec[jevt][jjet][jtrk];
			  double old_assotrk_pt = old_assotrk_vec.Pt();
			  double old_assotrk_eta = old_assotrk_vec.Eta();
			  double old_assotrk_phi = old_assotrk_vec.Phi();
			  TVector3 new_assotrk_vec = jet_newtrk_vec[jevt][jjet][jtrk];
			  double new_assotrk_pt = new_assotrk_vec.Pt();
			  double new_assotrk_eta = new_assotrk_vec.Eta();
			  double new_assotrk_phi = new_assotrk_vec.Phi();
			  int assotrk_charge = jet_trk_charge[jevt][jjet][jtrk];
                          double weight_ass=weight_3D[jevt][jjet][jtrk];
			  //double pthatw=weight_trig*weight_ass;
		
			  //if(new_assotrk_pt > newtrk_jt03_max) continue; // asso pt cut
			  
			  double deta_jetaxis = new_trgtrk_eta - new_assotrk_eta;
			  double dphi_jetaxis = new_trgtrk_phi - new_assotrk_phi;
			  double dphi2_jetaxis = new_assotrk_phi - new_trgtrk_phi;
			  
			  if(dphi_jetaxis > 1.5*TMath::Pi())
			    {
			      dphi_jetaxis = dphi_jetaxis - 2.0*TMath::Pi();
			    }
			  else if(dphi_jetaxis < -0.5*TMath::Pi())
			    {
			      dphi_jetaxis = dphi_jetaxis + 2.0*TMath::Pi();
			    }
			  
			  if(dphi2_jetaxis > 1.5*TMath::Pi())
			    {
			      dphi2_jetaxis = dphi2_jetaxis - 2.0*TMath::Pi();
			    }
			  else if(dphi2_jetaxis < -0.5*TMath::Pi())
			    {
			      dphi2_jetaxis = dphi2_jetaxis + 2.0*TMath::Pi();
			    }
			  
			  if(TMath::Abs(trgtrk_charge) > 0 && TMath::Abs(assotrk_charge) > 0)
			    {
			      if(TMath::Abs(new_trgtrk_eta) <= newtrk_eta_max && TMath::Abs(new_assotrk_eta) <= newtrk_eta_max)
				{
				  if((new_trgtrk_pt >= newtrk_jt03_min && new_trgtrk_pt <= newtrk_jt03_max) && (new_assotrk_pt >= newtrk_jt03_min && new_assotrk_pt <= newtrk_jt03_max))
				    {
				      hmixing_newchtrk_jt03->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt03->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt03->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt03->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				      
				      hmixing_newchtrk_jt03_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt03_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt03_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt03_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				      
				      if(trgtrk_charge > 0 && assotrk_charge > 0)
					{
					  hmixing_newchtrk_jt03_pp->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_pp->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_pp->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_pp->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  
					  hmixing_newchtrk_jt03_pp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_pp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_pp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_pp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					}
				      if(trgtrk_charge < 0 && assotrk_charge < 0)
					{
					  hmixing_newchtrk_jt03_mm->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_mm->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_mm->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_mm->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  
					  hmixing_newchtrk_jt03_mm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_mm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_mm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_mm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					}
				      if(trgtrk_charge > 0 && assotrk_charge < 0)
					{
					  hmixing_newchtrk_jt03_pm->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_pm->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_pm->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_pm->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  
					  hmixing_newchtrk_jt03_pm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_pm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_pm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_pm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					}
				      if(trgtrk_charge < 0 && assotrk_charge > 0)
					{
					  hmixing_newchtrk_jt03_mp->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_mp->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_mp->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_mp->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  
					  hmixing_newchtrk_jt03_mp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_mp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_mp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt03_mp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					}
				    }
				  if((new_trgtrk_pt >= newtrk_jt0p33_min && new_trgtrk_pt <= newtrk_jt0p33_max) && (new_assotrk_pt >= newtrk_jt0p33_min && new_assotrk_pt <= newtrk_jt0p33_max))
				    {
				      hmixing_newchtrk_jt0p33->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt0p33->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt0p33->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt0p33->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				      
				      hmixing_newchtrk_jt0p33_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt0p33_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt0p33_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt0p33_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				      if(trgtrk_charge > 0 && assotrk_charge > 0)
					{
					  hmixing_newchtrk_jt0p33_pp->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_pp->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_pp->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_pp->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  
					  hmixing_newchtrk_jt0p33_pp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_pp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_pp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_pp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					}
				      if(trgtrk_charge < 0 && assotrk_charge < 0)
					{
					  hmixing_newchtrk_jt0p33_mm->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_mm->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_mm->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_mm->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  
					  hmixing_newchtrk_jt0p33_mm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_mm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_mm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_mm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					}
				      if(trgtrk_charge > 0 && assotrk_charge < 0)
					{
					  hmixing_newchtrk_jt0p33_pm->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_pm->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_pm->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_pm->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  
					  hmixing_newchtrk_jt0p33_pm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_pm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_pm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_pm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					}
				      if(trgtrk_charge < 0 && assotrk_charge > 0)
					{
					  hmixing_newchtrk_jt0p33_mp->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_mp->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_mp->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_mp->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  
					  hmixing_newchtrk_jt0p33_mp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_mp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_mp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p33_mp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					}
				    }
				  if((new_trgtrk_pt >= newtrk_jt0p53_min && new_trgtrk_pt <= newtrk_jt0p53_max) && (new_assotrk_pt >= newtrk_jt0p53_min && new_assotrk_pt <= newtrk_jt0p53_max))
				    {
				      hmixing_newchtrk_jt0p53->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt0p53->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt0p53->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt0p53->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				      
				      hmixing_newchtrk_jt0p53_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt0p53_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt0p53_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				      hmixing_newchtrk_jt0p53_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
				      
				      if(trgtrk_charge > 0 && assotrk_charge > 0)
					{
					  hmixing_newchtrk_jt0p53_pp->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_pp->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_pp->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_pp->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  
					  hmixing_newchtrk_jt0p53_pp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_pp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_pp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_pp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					}
				      if(trgtrk_charge < 0 && assotrk_charge < 0)
					{
					  hmixing_newchtrk_jt0p53_mm->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_mm->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_mm->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_mm->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  
					  hmixing_newchtrk_jt0p53_mm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_mm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_mm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_mm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					}
				      if(trgtrk_charge > 0 && assotrk_charge < 0)
					{
					  hmixing_newchtrk_jt0p53_pm->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_pm->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_pm->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_pm->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  
					  hmixing_newchtrk_jt0p53_pm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_pm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_pm_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_pm_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					}
				      if(trgtrk_charge < 0 && assotrk_charge > 0)
					{
					  hmixing_newchtrk_jt0p53_mp->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_mp->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_mp->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_mp->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  
					  hmixing_newchtrk_jt0p53_mp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_mp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_mp_[nchbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					  hmixing_newchtrk_jt0p53_mp_[nchbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, pthatw/4.);
					}
				    }
				} // eta condition
			    } // charge condition
			  
			} // asso trk loop
		    } // trg trk loop end
		} // 1st jet loop end
	    } // 2nd jet loop end
	} // 2nd event loop end
    } // 1st event loop end
  std::cout<<endl;
  std::cout<<"~~~~~end mixed event correlation~~~~~~~~"<<std::endl;
} // function loop
      
void print_start()
{
  cout << endl;
  time_t init = time(0);
  char* init_time = ctime(&init); // convert now to string form                                        
  cout << "Starting at : " << init_time << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
  cout << endl;
}

void print_stop()
{
  time_t end = time(0);
  char* end_time = ctime(&end); // convert now to string form                                                   
  cout << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << endl;
  cout << "Stopping at : " << end_time << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << endl;
}
