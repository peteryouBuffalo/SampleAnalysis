#ifndef Plots_cxx
#define Plots_cxx

#include "Obj.cxx"
#include "TString.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <TTreeReaderArray.h>

// These constants define the number of bins
// & pT ranges we want for different variables.
const unsigned NBIN_PT_JET = 2000;
const float X_PT_JET[2] = {0, 2000};
const unsigned NBIN_PT_Z = 2000;
const float X_PT_Z[2] = {0, 2000};
const unsigned NBIN_ETA = 480;
const float X_ETA[2] = {-5.0, 5.0};
const unsigned NBIN_PHI = 480;
const float X_PHI[2] = {-4.0, 4.0};

const float NBIN_M_H = 300;
const float X_M_H[2] = {0, 300};
const float NBIN_M_Z = 300;
const float X_M_Z[2] = {0, 300};

const float NBIN_DR = 90;
const float X_DR[2] = {-0.5, 8.5};

// Main Plots class
class VHPlots
{

  public:
    // Constructor
    VHPlots(TString name) : m_name(name) {
 
      h_pt_jet = new TH1D(name + "_pt_jet", "", NBIN_PT_JET, X_PT_JET[0], X_PT_JET[1]);
      h_eta_jet = new TH1D(name + "_eta_jet", "", NBIN_ETA, X_ETA[0], X_ETA[1]);
      h_phi_jet = new TH1D(name + "_phi_jet", "", NBIN_PHI, X_PHI[0], X_PHI[1]);
      h_mSV_jet = new TH1D(name + "_mSV_jet", "", 100, 0, 10);
      h_Njet = new TH1D(name + "_Njet", "", 15, -0.5, 14.5);   

      h_H_pt_jet0 = new TH1D(name + "_H_pt_jet0", "", NBIN_PT_JET, X_PT_JET[0], X_PT_JET[1]);
      h_H_pt_jet1 = new TH1D(name + "_H_pt_jet1", "", NBIN_PT_JET, X_PT_JET[0], X_PT_JET[1]);
      h_H_pt_jet2 = new TH1D(name + "_H_pt_jet2", "", NBIN_PT_JET, X_PT_JET[0], X_PT_JET[1]);
      h_H_eta_jet0 = new TH1D(name + "_H_eta_jet0", "", NBIN_ETA, X_ETA[0], X_ETA[1]);
      h_H_eta_jet1 = new TH1D(name + "_H_eta_jet1", "", NBIN_ETA, X_ETA[0], X_ETA[1]);
      h_H_eta_jet2 = new TH1D(name + "_H_eta_jet2", "", NBIN_ETA, X_ETA[0], X_ETA[1]);
      h_H_phi_jet0 = new TH1D(name + "_H_phi_jet0", "", NBIN_PHI, X_PHI[0], X_PHI[1]);
      h_H_phi_jet1 = new TH1D(name + "_H_phi_jet1", "", NBIN_PHI, X_PHI[0], X_PHI[1]);
      h_H_phi_jet2 = new TH1D(name + "_H_phi_jet2", "", NBIN_PHI, X_PHI[0], X_PHI[1]);

      h_Z_pt_jet0 = new TH1D(name + "_Z_pt_jet0", "", NBIN_PT_JET, X_PT_JET[0], X_PT_JET[1]);
      h_Z_pt_jet1 = new TH1D(name + "_Z_pt_jet1", "", NBIN_PT_JET, X_PT_JET[0], X_PT_JET[1]);
      h_Z_pt_jet2 = new TH1D(name + "_Z_pt_jet2", "", NBIN_PT_JET, X_PT_JET[0], X_PT_JET[1]);
      h_Z_eta_jet0 = new TH1D(name + "_Z_eta_jet0", "", NBIN_ETA, X_ETA[0], X_ETA[1]);
      h_Z_eta_jet1 = new TH1D(name + "_Z_eta_jet1", "", NBIN_ETA, X_ETA[0], X_ETA[1]);
      h_Z_eta_jet2 = new TH1D(name + "_Z_eta_jet2", "", NBIN_ETA, X_ETA[0], X_ETA[1]);
      h_Z_phi_jet0 = new TH1D(name + "_Z_phi_jet0", "", NBIN_PHI, X_PHI[0], X_PHI[1]);
      h_Z_phi_jet1 = new TH1D(name + "_Z_phi_jet1", "", NBIN_PHI, X_PHI[0], X_PHI[1]);
      h_Z_phi_jet2 = new TH1D(name + "_Z_phi_jet2", "", NBIN_PHI, X_PHI[0], X_PHI[1]);

      h_HMass = new TH1D(name + "_HMass", "", NBIN_M_H, X_M_H[0], X_M_H[1]);
      h_ZMass = new TH1D(name + "_ZMass", "", NBIN_M_Z, X_M_Z[0], X_M_Z[1]);

      h_dR_H = new TH1D(name + "_dR_H", "", 100, 0, 10);
      h_dR_Z = new TH1D(name + "_dR_Z", "", 100, 0, 10);
      h_dPhi_H = new TH1D(name + "_dPhi_H", "", 60, -TMath::Pi(), TMath::Pi());
      h_dPhi_Z = new TH1D(name + "_dPhi_Z", "", 60, -TMath::Pi(), TMath::Pi());

      h_pt_jet->Sumw2();
      h_eta_jet->Sumw2();
      h_phi_jet->Sumw2();
      h_mSV_jet->Sumw2();
      h_Njet->Sumw2();

      h_H_pt_jet0->Sumw2(); h_H_pt_jet1->Sumw2(); h_H_pt_jet2->Sumw2();
      h_H_eta_jet0->Sumw2(); h_H_eta_jet1->Sumw2(); h_H_eta_jet2->Sumw2();
      h_H_phi_jet0->Sumw2(); h_Z_phi_jet1->Sumw2(); h_H_phi_jet2->Sumw2();
      h_Z_pt_jet0->Sumw2(); h_Z_pt_jet1->Sumw2(); h_Z_pt_jet2->Sumw2();
      h_Z_eta_jet0->Sumw2(); h_Z_eta_jet1->Sumw2(); h_Z_eta_jet2->Sumw2();
      h_Z_phi_jet0->Sumw2(); h_Z_phi_jet1->Sumw2(); h_Z_phi_jet2->Sumw2();

      h_HMass->Sumw2();
      h_ZMass->Sumw2();
      h_dR_H->Sumw2();
      h_dR_Z->Sumw2();
      h_dPhi_H->Sumw2();
      h_dPhi_Z->Sumw2();
    } ;

    // Fill the general histograms.
    void Fill(HObj& H, ZObj& Z, float w=1) {

      h_HMass->Fill(H.m_lvec.M(), w);
      h_ZMass->Fill(Z.m_lvec.M(), w);
      
      if (H.nJets() >= 2) {
        h_dR_H->Fill(H.DeltaR(), w);
        h_dPhi_H->Fill(H.DPhi(), w);
  
        TLorentzVector jet0 = H.getJet(0).m_lvec;
        h_H_pt_jet0->Fill(jet0.Pt(), w);
        h_H_eta_jet0->Fill(jet0.Eta(), w);
        h_H_phi_jet0->Fill(jet0.Phi(), w);

        TLorentzVector jet1 = H.getJet(1).m_lvec;
        h_H_pt_jet1->Fill(jet1.Pt(), w);
        h_H_eta_jet1->Fill(jet1.Eta(), w);
        h_H_phi_jet1->Fill(jet1.Phi(), w);

        if (H.nJets() > 2) {
          TLorentzVector jet2 = H.getJet(2).m_lvec;
          h_H_pt_jet2->Fill(jet2.Pt(), w);
          h_H_eta_jet2->Fill(jet2.Eta(), w);
          h_H_phi_jet2->Fill(jet2.Phi(), w);
        }
      }
      if (Z.nJets() >= 2) {
        h_dR_Z->Fill(Z.DeltaR(), w);
        h_dPhi_Z->Fill(Z.DPhi(), w);

        TLorentzVector jet0 = Z.getJet(0).m_lvec;
        h_Z_pt_jet0->Fill(jet0.Pt(), w);
        h_Z_eta_jet0->Fill(jet0.Eta(), w);
        h_Z_phi_jet0->Fill(jet0.Phi(), w);

        TLorentzVector jet1 = Z.getJet(1).m_lvec;
        h_Z_pt_jet1->Fill(jet1.Pt(), w);
        h_Z_eta_jet1->Fill(jet1.Eta(), w);
        h_Z_phi_jet1->Fill(jet1.Phi(), w);

        if (Z.nJets() > 2) {
          TLorentzVector jet2 = Z.getJet(2).m_lvec;
          h_Z_pt_jet2->Fill(jet2.Pt(), w);
          h_Z_eta_jet2->Fill(jet2.Eta(), w);
          h_Z_phi_jet2->Fill(jet2.Phi(), w);
        }
      }
    };

    // Fill stuff about jets.
    void FillJets(std::vector<JetObj>& jets, float w=1) {
      for (auto it : jets) {
        h_pt_jet->Fill(it.m_lvec.Pt(), w);
        h_eta_jet->Fill(it.m_lvec.Eta(), w);
        h_phi_jet->Fill(it.m_lvec.Phi(), w);
        h_mSV_jet->Fill(it.m_mSV, w);
      }
    }

    void FillNjet(size_t nJet, float w=1) {
      h_Njet->Fill(nJet, w);
    }

    // Return a list of all the histograms.
    std::vector<TH1*> returnHisto() {
      std::vector<TH1*> histolist;
 
      // Jet-Related Plots
      histolist.push_back(h_pt_jet);
      histolist.push_back(h_eta_jet);
      histolist.push_back(h_phi_jet);
      histolist.push_back(h_mSV_jet);
      histolist.push_back(h_Njet);
 
      histolist.push_back(h_H_pt_jet0); histolist.push_back(h_H_pt_jet1); 
      histolist.push_back(h_H_pt_jet2); 
      histolist.push_back(h_H_eta_jet0); histolist.push_back(h_H_eta_jet1);
      histolist.push_back(h_H_eta_jet2);
      histolist.push_back(h_H_phi_jet0); histolist.push_back(h_H_phi_jet1);
      histolist.push_back(h_H_phi_jet2);

      histolist.push_back(h_Z_pt_jet0); histolist.push_back(h_Z_pt_jet1);
      histolist.push_back(h_Z_pt_jet2);
      histolist.push_back(h_Z_eta_jet0); histolist.push_back(h_Z_eta_jet1);
      histolist.push_back(h_Z_eta_jet2);
      histolist.push_back(h_Z_phi_jet0); histolist.push_back(h_Z_phi_jet1);
      histolist.push_back(h_Z_phi_jet2);

     
      // Mass-Related Plots
      histolist.push_back(h_HMass);
      histolist.push_back(h_ZMass);

      // VH-Related Plots
      histolist.push_back(h_dR_H);
      histolist.push_back(h_dR_Z);
      histolist.push_back(h_dPhi_H);
      histolist.push_back(h_dPhi_Z);

      return histolist;
    }

  protected:
    // Variables
    TString m_name;
    
    //Higgs Plots
    TH1D* h_HMass;    
    TH1D* h_dR_H;
    TH1D* h_dPhi_H;
    
    TH1D* h_H_pt_jet0;
    TH1D* h_H_eta_jet0;
    TH1D* h_H_phi_jet0;
    TH1D* h_H_pt_jet1;
    TH1D* h_H_eta_jet1;
    TH1D* h_H_phi_jet1;
    TH1D* h_H_pt_jet2;
    TH1D* h_H_eta_jet2;
    TH1D* h_H_phi_jet2;

    //Z Plots
    TH1D* h_ZMass;
    TH1D* h_dR_Z;
    TH1D* h_dPhi_Z;

    TH1D* h_Z_pt_jet0;
    TH1D* h_Z_eta_jet0;
    TH1D* h_Z_phi_jet0;
    TH1D* h_Z_pt_jet1;
    TH1D* h_Z_eta_jet1;
    TH1D* h_Z_phi_jet1;
    TH1D* h_Z_pt_jet2;
    TH1D* h_Z_eta_jet2;
    TH1D* h_Z_phi_jet2;
    
    //VH Plots
    TH1D* h_pt_jet;
    TH1D* h_eta_jet;
    TH1D* h_phi_jet;
    TH1D* h_mSV_jet;
    TH1D* h_Njet;
} ;

//=== Gen Plots ===
// For purposes of analyzing MC samples, we want to occasionally deal with the
// Gen-Level objects. Rather than continuously adding plots, it'd be nice to have
// a class that handles them all together.
class GenPlots
{
  public:
    // Constructor
    GenPlots(TString name) : m_name(name) {
      
      h_Higgs_mass = new TH1D(name + "_Higgs_mass", "", NBIN_M_H, X_M_H[0], X_M_H[1]);
      h_Higgs_pt = new TH1D(name + "_Higgs_pt", "", NBIN_PT_JET, X_PT_JET[0], X_PT_JET[1]);
      h_Higgs_phi = new TH1D(name + "_Higgs_phi", "", NBIN_PHI, X_PHI[0], X_PHI[1]);

      h_Z_mass = new TH1D(name + "_Z_mass", "", NBIN_M_H, X_M_H[0], X_M_H[1]);
      h_Z_pt = new TH1D(name + "_Z_pt", "", NBIN_PT_JET, X_PT_JET[0], X_PT_JET[1]);
      h_Z_phi = new TH1D(name + "_Z_phi", "", NBIN_PHI, X_PHI[0], X_PHI[1]);

      h_nHiggs = new TH1D(name + "_nHiggs", "", 30, -0.5, 29.5);
      h_nZ = new TH1D(name + "_nZ", "", 30, -0.5, 29.5);
      h_cjet_phi = new TH1D(name + "_cjet_phi", "", NBIN_PHI, X_PHI[0], X_PHI[1]);
      h_bjet_phi = new TH1D(name + "_bjet_phi", "", NBIN_PHI, X_PHI[0], X_PHI[1]);
      h_ljet_phi = new TH1D(name + "_ljet_phi", "", NBIN_PHI, X_PHI[0], X_PHI[1]);

      h_HZ_dPhi = new TH1D(name + "_HZ_dPhi", "", NBIN_PHI, X_PHI[0], X_PHI[1]);
      h_HZ_dR = new TH1D(name + "_HZ_dR", "", 40, -0.5, 12.5); 

      h_Hc_dR = new TH1D(name + "_Hc_dR", "", NBIN_DR, X_DR[0], X_DR[1]);
      h_Hb_dR = new TH1D(name + "_Hb_dR", "", NBIN_DR, X_DR[0], X_DR[1]);
      h_Zc_dR = new TH1D(name + "_Zc_dR", "", NBIN_DR, X_DR[0], X_DR[1]);
      h_Zb_dR = new TH1D(name + "_Zb_dR", "", NBIN_DR, X_DR[0], X_DR[1]);
      h_Zl_dR = new TH1D(name + "_Zl_dR", "", NBIN_DR, X_DR[0], X_DR[1]);
      h_cc_dR = new TH1D(name + "_cc_dR", "", NBIN_DR, X_DR[0], X_DR[1]);
      h_bb_dR = new TH1D(name + "_bb_dR", "", NBIN_DR, X_DR[0], X_DR[1]);
      h_ll_dR = new TH1D(name + "_ll_dR", "", NBIN_DR, X_DR[0], X_DR[1]);

      h_cc_dPhi = new TH1D(name + "_cc_dPhi", "", NBIN_PHI, X_PHI[0], X_PHI[1]);
      h_bb_dPhi = new TH1D(name + "_bb_dPhi", "", NBIN_PHI, X_PHI[0], X_PHI[1]);
    };

    void Fill(GenObj* H, GenObj* Z, float w=1) {
      h_Higgs_mass->Fill(H->m_lvec.M(), w);
      h_Higgs_pt->Fill(H->m_lvec.Pt(), w);
      h_Higgs_phi->Fill(H->m_lvec.Phi(), w);

      h_Z_mass->Fill(Z->m_lvec.M(), w);
      h_Z_pt->Fill(Z->m_lvec.Pt(), w);
      h_Z_phi->Fill(Z->m_lvec.Phi(), w);

      h_HZ_dPhi->Fill(Z->m_lvec.DeltaPhi(H->m_lvec), w);
      h_HZ_dR->Fill(Z->m_lvec.DeltaR(H->m_lvec), w);
    }

    void FillMult(Int_t nHiggs, Int_t nZ, float w=1) {
      h_nHiggs->Fill(nHiggs);
      h_nZ->Fill(nZ);
    }

    void FillJets(GenObj* obj, std::vector<JetObj>& jets, Int_t jetID, float w=1) {
      Int_t pdgID = obj->m_pdgID;
      if (abs(pdgID) == 25) { // Higgs
        for (auto it : jets) {
          float dR = it.m_lvec.DeltaR(obj->m_lvec);
          if (jetID == 4) h_Hc_dR->Fill(dR, w);
          else if (jetID == 5) h_Hb_dR->Fill(dR, w);
        }
 
        if (jets.size() < 2) return;
        for (Int_t i = 0; i < jets.size(); ++i) {
          for (Int_t j = i + 1; j < jets.size(); ++j) {
            float dR = jets.at(i).m_lvec.DeltaR(jets.at(j).m_lvec);
            float dPhi = jets.at(i).m_lvec.DeltaPhi(jets.at(j).m_lvec);
            if (jetID == 4) {
              h_cc_dPhi->Fill(dPhi, w);
              h_cc_dR->Fill(dR, w);
            }
          }
        }
      }
      else if (abs(pdgID) == 23) { // Z
        for (auto it : jets) {
          float dR = it.m_lvec.DeltaR(obj->m_lvec);
          float phi = it.m_lvec.Phi();
          if (jetID == 5) {
            h_Zb_dR->Fill(dR, w);
            h_cjet_phi->Fill(phi, w);
          }
          else if (jetID == 4){
            h_Zc_dR->Fill(dR, w);
            h_bjet_phi->Fill(phi, w);
          }
          else{
            h_Zl_dR->Fill(dR, w);
            h_ljet_phi->Fill(phi, w);
          }
        }
         
        if (jets.size() < 2) return;
        for (Int_t i = 0; i < jets.size(); ++i) {
          for (Int_t j = 0; j < jets.size(); ++j) {
            float dR = jets.at(i).m_lvec.DeltaR(jets.at(j).m_lvec);
            float dPhi = jets.at(i).m_lvec.DeltaPhi(jets.at(j).m_lvec);
            if (jetID == 5){
              h_bb_dR->Fill(dR, w);
              h_bb_dPhi->Fill(dPhi, w);
            }
            else if (jetID == 4) h_cc_dR->Fill(dR, w);
            else h_ll_dR->Fill(dR, w);
          }
        }
      }
    }//end-method 

    std::vector<TH1*> returnHisto() {
      std::vector<TH1*> histolist;
      
      histolist.push_back(h_Higgs_mass); histolist.push_back(h_Higgs_pt); histolist.push_back(h_Higgs_phi);
      histolist.push_back(h_Z_mass); histolist.push_back(h_Z_pt); histolist.push_back(h_Z_phi);
      histolist.push_back(h_cjet_phi); histolist.push_back(h_bjet_phi); histolist.push_back(h_ljet_phi);     

      histolist.push_back(h_nHiggs); histolist.push_back(h_nZ);
      histolist.push_back(h_HZ_dPhi); histolist.push_back(h_HZ_dR);
      histolist.push_back(h_Hc_dR); histolist.push_back(h_Hb_dR);
      histolist.push_back(h_Zc_dR); histolist.push_back(h_Zb_dR); histolist.push_back(h_Zl_dR);

      histolist.push_back(h_cc_dR); histolist.push_back(h_bb_dR); histolist.push_back(h_ll_dR);
      histolist.push_back(h_cc_dPhi); histolist.push_back(h_bb_dPhi);
      return histolist;
    }

  protected:
    // Variables
    TString m_name;

    // Plots Related to Gen Objects Themselves
    TH1D* h_nHiggs;
    TH1D* h_nZ;
    TH1D* h_Higgs_mass;
    TH1D* h_Higgs_pt;
    TH1D* h_Higgs_phi;
    TH1D* h_Z_mass;
    TH1D* h_Z_pt;
    TH1D* h_Z_phi;
    TH1D* h_cjet_phi;
    TH1D* h_bjet_phi;
    TH1D* h_ljet_phi;

    // Plots Related to Physics Analysis
    TH1D* h_HZ_dPhi;      //Phi separation of H & Z
    TH1D* h_HZ_dR;        //R separation of H & Z
    TH1D* h_Hc_dR;        //R separation of H & c-jets
    TH1D* h_Hb_dR;        //R separation of H & b-jets
    TH1D* h_Zc_dR;        //R separation of Z & c-jets
    TH1D* h_Zb_dR;        //R separation of Z & b-jets
    TH1D* h_Zl_dR;        //R separation of Z & l-jets
    TH1D* h_cc_dR;        //R separation of c-jets
    TH1D* h_bb_dR;        //R separation of b-jets
    TH1D* h_ll_dR;        //R separation of l-jets
    TH1D* h_cc_dPhi;      //Phi separation of c-jets
    TH1D* h_bb_dPhi;      //Phi separation of b-jets
} ;

#endif
