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
const float X_PT_JET[2] = {0,2000};
const unsigned NBIN_PT_Z = 2000;
const float X_PT_Z[2] = {0,2000};
const unsigned NBIN_ETA = 480;
const float X_ETA[2] = {-2.4,2.4};

class VHPlots 
{
  
  public:
    // Constructor
    VHPlots(TString name) : m_name(name) {

      // Initialize our histograms.
      h_pt_jet = new TH1D(name + "_pt_jet","",NBIN_PT_JET,X_PT_JET[0],X_PT_JET[1]) ;
      h_eta_jet = new TH1D(name + "_eta_jet","",NBIN_ETA,X_ETA[0],X_ETA[1]) ;
      h_Njet = new TH1D(name + "_Njet","",15,0,15) ;
      h_flavor_jet = new TH1D(name + "_flavor_jet", "", 10,-0.5, 9.5);

      h_Nbjet = new TH1D(name + "_Nbjet", "", 15, 0, 15);
      h_Ncjet = new TH1D(name + "_Ncjet", "", 15, 0, 15);
      h_Nljet = new TH1D(name + "_Nljet", "", 15, 0, 15);     

      h_nElec = new TH1D(name + "_Nelec", "", 15, 0, 15);
      h_nMuon = new TH1D(name + "_Nmuon", "", 15, 0, 15);
      h_nLep = new TH1D(name + "_Nlep", "", 15, 0, 15);
     
      h_HMass = new TH1D(name + "_HMass", "", 300, 0, 300);
      h_ZMass = new TH1D(name + "_ZMass", "", 300, 0, 300);

      h_pt_jet->Sumw2() ;
      h_eta_jet->Sumw2() ;
      h_Njet->Sumw2() ;
      h_flavor_jet->Sumw2() ;
      h_Nbjet->Sumw2();
      h_Ncjet->Sumw2();
      h_Nljet->Sumw2();

      h_HMass->Sumw2();
      h_ZMass->Sumw2();

      h_nElec->Sumw2() ;
      h_nMuon->Sumw2() ;
      h_nLep->Sumw2() ;
    } ;
 
    // Fill the general histograms.   
    void Fill(JetObj& J, float w=1) { 

      // Fill the basic information of the jet
      h_pt_jet->Fill(J.m_lvec.Pt(), w) ;
      h_eta_jet->Fill(J.m_lvec.Eta(), w) ;
      h_flavor_jet->Fill(J.m_flav, w);
    } ;

    // Fill the number of jets.
    void FillNjet(size_t nJet, float w=1) {
      h_Njet->Fill(nJet,w) ;
    } ;

    // Fill the number of jets by flavor.
    void FillNjet_flav(size_t nJet, int flav, float w=1) {
      if (flav == 5) h_Nbjet->Fill(nJet, w);
      else if (flav == 4) h_Ncjet->Fill(nJet, w);
      else if (flav > 0) h_Nljet->Fill(nJet, w);
    }

    // Fill the number of leptons
    void FillNlep(size_t nElec, size_t nMuon, float w=1) {
      h_nElec->Fill(nElec, w);
      h_nMuon->Fill(nMuon, w);
      h_nLep->Fill(nElec + nMuon, w);
    } ;

    // Fill Masses
    void FillHmass(Float_t mass, float w=1) {
      h_HMass->Fill(mass, w);
    }

    void FillZmass(Float_t mass, float w=1) {
      h_ZMass->Fill(mass, w);
    }

    // Return a list of all the histograms.
    std::vector<TH1*> returnHisto(){
      std::vector<TH1*> histolist ;
      
      histolist.push_back(h_pt_jet) ;
      histolist.push_back(h_eta_jet) ;   
      histolist.push_back(h_Njet) ;
      histolist.push_back(h_flavor_jet);
      histolist.push_back(h_Nbjet);
      histolist.push_back(h_Ncjet);
      histolist.push_back(h_Nljet);     

      histolist.push_back(h_HMass);
      histolist.push_back(h_ZMass);

      histolist.push_back(h_nElec);
      histolist.push_back(h_nMuon);
      histolist.push_back(h_nLep);      

      return histolist ; 
    }

  protected:
    // Variables
    TString m_name ;
    TH1D* h_pt_jet ;
    TH1D* h_eta_jet ;
    TH1D* h_Njet ;
    TH1D* h_flavor_jet;
    TH1D* h_Nbjet;
    TH1D* h_Ncjet;
    TH1D* h_Nljet;

    TH1D* h_HMass;
    TH1D* h_ZMass;

    TH1D* h_nElec;
    TH1D* h_nMuon;
    TH1D* h_nLep;
} ;


#endif
