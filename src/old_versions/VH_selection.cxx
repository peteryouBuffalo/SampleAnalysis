#define VH_selection_cxx

#include <math.h>

#include "TList.h"
#include "TParameter.h"

#include "VH_selection.h"
#include "Global.h"
#include "Obj.cxx"

//VH_selection::VH_selection(bool isData) : Selector(isData), h_zee_jet(0), h_zmm_jet(0) {}

VH_selection::~VH_selection() {
  if (h_VH) delete h_VH;
} 

void VH_selection::SlaveBegin(Reader* r) {
  
  h_evt = new TH1D("Nevt","",4,-1.5,2.5) ; //bin 1: total negative weight events, bin 2: total positive weight events, bin 3: total event weighted by genWeight (= bin2 - bin1, if genWeight are always -1,1

  h_VH = new VHPlots("VbbHcc") ;
  h_GenPlots = new GenPlots("GenObj") ;

  h_flavor_jet = new TH1D("flavor_jet", "", 25, -1.5, 23.5);
  h_Nbjet = new TH1D("Nbjet", "", 13, -0.5, 12.5);
  h_Ncjet = new TH1D("Ncjet", "", 13, -0.5, 12.5);
  h_Nljet = new TH1D("Nljet", "", 13, -0.5, 12.5);

  h_Higgs_nJet = new TH1D("Higgs_nJet", "", 13, -0.5, 12.5);
  h_evt_cutflow = new TH1D("evt_cutflow", "", 13, -0.5, 12.5);
  h_elec_cutflow = new TH1D("elec_cutflow", "", 10, -0.5, 9.5);
  h_muon_cutflow = new TH1D("muon_cutflow", "", 10, -0.5, 9.5);

  //Add histograms to fOutput so they can be saved in Processor::SlaveTerminate
  r->GetOutputList()->Add(h_evt) ;
  std::vector<TH1*> tmp = h_VH->returnHisto() ;
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_GenPlots->returnHisto() ;
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  
  r->GetOutputList()->Add(h_evt_cutflow);
  r->GetOutputList()->Add(h_elec_cutflow);
  r->GetOutputList()->Add(h_muon_cutflow);

  r->GetOutputList()->Add(h_flavor_jet);
  r->GetOutputList()->Add(h_Nbjet);
  r->GetOutputList()->Add(h_Ncjet);
  r->GetOutputList()->Add(h_Nljet);
  r->GetOutputList()->Add(h_Higgs_nJet); 
 
  const Int_t nevt = 12;
  const Int_t nx = 4, nx1 = 4;
  const char *evt_cut[nevt] = { "All Events", "Pass Lepton Selection",
    "Jet requirements", "p_{T}^{miss} filter", "Jet ID", "Pass b-tag selection", 
    "N_2^{DDT}", "Trigger", "Golden JSON", "CvB", "CvL"};
  const char *elec_cut[nx] = { "all", "ip", "kine", "ID" };
  const char *muon_cut[nx1] = { "all", "kine", "medium muon ID", "iso" };

  const Int_t nxBB = 4;
  const char *evt_cutBB[nxBB] = { "All Events", "Pass GenObj reconstruction", "Pass c-jet (H) requirement",
    "Pass b-jet (Z) requirement" };

  for (int i=1;i<=nxBB;i++) h_evt_cutflow->GetXaxis()->SetBinLabel(i+1.5, evt_cutBB[i-1]);

  for (int i=1;i<=nx;i++) h_elec_cutflow->GetXaxis()->SetBinLabel(i+1.5, elec_cut[i-1]);
  for (int i=1;i<=nx1;i++) h_muon_cutflow->GetXaxis()->SetBinLabel(i+1.5, muon_cut[i-1]);

}

void VH_selection::Process(Reader* r) {

  //Weights
  float genWeight = 1.;
  float puSF = 1.;
  float l1preW = 1.;
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
  if (*(r->genWeight) < 0) genWeight = -1. ;
  if (*(r->genWeight) == 0) {
    genWeight = 0;
    h_evt->Fill(0) ;
  }
  if (*(r->genWeight) < 0) h_evt->Fill(-1) ;
  if (*(r->genWeight) > 0) h_evt->Fill(1) ;
  if (m_centralGenWeight != 0)  genWeight = *(r->genWeight)/m_centralGenWeight; //use central gen weight to normalize gen weight
  puSF = PileupSF(*(r->Pileup_nTrueInt));
#endif

  h_evt->Fill(2,genWeight);

#if defined(MC_2016) || defined(MC_2017)
  l1preW = *(r->L1PreFiringWeight_Nom);
  //std::cout << "\nPrefiring: " << l1preW;
  if (m_l1prefiringType == "l1prefiringu") l1preW = *(r->L1PreFiringWeight_Up);
  if (m_l1prefiringType == "l1prefiringd") l1preW = *(r->L1PreFiringWeight_Dn);
#endif

#if defined(DATA_2016) || defined(DATA_2017) || defined(DATA_2018)
  h_evt->Fill(-1) ;
  if (!m_lumiFilter.Pass(*(r->run),*(r->luminosityBlock))) return;
  h_evt->Fill(1) ;
#endif

  float evtW = 1. ;

  if (!m_isData) evtW *= genWeight*puSF*l1preW;

  //=============Get objects============= 
  
  // Jets
  std::vector<JetObj> jets ;
  for (unsigned int i = 0 ; i < *(r->nJet) ; ++i) {
    int jetFlav = -1 ;
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
    jetFlav = (r->Jet_hadronFlavour)[i];
#endif
    JetObj jet((r->Jet_pt)[i],(r->Jet_eta)[i],(r->Jet_phi)[i],(r->Jet_mass)[i],jetFlav,(r->Jet_btagDeepB)[i],(r->Jet_btagDeepFlavB)[i],(r->Jet_puIdDisc)[i]) ;
    jets.push_back(jet) ;
  }

  //Make selection and fill histograms
  h_VH->FillJets(jets);
  h_VH->FillNjet(jets.size());
  
  //== Generator Object Reconstruction ========================================
  // For purposes of being able to analyze the MC samples, we want to locate 
  // the Z & Higgs boson data inside the generated particles. We want to store
  // them here. Remember that abs(pdgID) = 23 is Z, abs(pdgID) = 25 is Higgs.
  GenObj* genZ;
  GenObj* genHiggs;

  Int_t nHiggs, nZ;
  for (unsigned int i = 0; i < *(r->nGenPart); ++i) {
    
    // Get the PDG ID and keep this particle if it's what we want.
    int pdgID = (r->GenPart_pdgId)[i];
    if (abs(pdgID) != 23 && abs(pdgID) != 25) continue;

    // See if it's a Z boson
    if (abs(pdgID) == 23) {
      nZ++;
      genZ = new GenObj(pdgID, (r->GenPart_pt)[i], (r->GenPart_eta)[i],
                 (r->GenPart_phi)[i], (r->GenPart_mass)[i], i,
                 (r->GenPart_genPartIdxMother)[i], (r->GenPart_status)[i]);
    }

    // See if it's a Higgs boson
    if (abs(pdgID) == 25) {
      nHiggs++;
      genHiggs = new GenObj(pdgID, (r->GenPart_pt)[i], (r->GenPart_eta)[i],
                     (r->GenPart_phi)[i], (r->GenPart_mass)[i], i,
                     (r->GenPart_genPartIdxMother)[i], (r->GenPart_status)[i]);
    }
  }  

  h_GenPlots->FillMult(nHiggs, nZ, evtW);

  //== Event Selection ========================================================
  h_evt_cutflow->Fill(1);  // all events
  

  // Make sure we have the proper generator events. (They should exist in every
  // single event by definition of the file.)
  if (genZ != NULL and genHiggs != NULL) {
    h_evt_cutflow->Fill(2); // theoretically all events

    h_GenPlots->Fill(genHiggs, genZ, evtW);
  }  

  // Select the types of jets
  std::vector<JetObj> bjets, cjets, ljets;
  for (auto it : jets) {
    if (it.m_flav == 5) bjets.push_back(it);
    else if (it.m_flav == 4) cjets.push_back(it);
    else if (it.m_flav >= 0) ljets.push_back(it);

    h_flavor_jet->Fill(it.m_flav, evtW);
  } 

  h_Nbjet->Fill(bjets.size(), evtW);
  h_Ncjet->Fill(cjets.size(), evtW);
  h_Nljet->Fill(ljets.size(), evtW);

  //== Handle Stuff Related to the Jets & Gen Objects here ==
  h_GenPlots->FillJets(genHiggs, bjets, 5, evtW);
  h_GenPlots->FillJets(genHiggs, cjets, 4, evtW);
  h_GenPlots->FillJets(genZ, bjets, 5, evtW);
  h_GenPlots->FillJets(genZ, cjets, 4, evtW);
  h_GenPlots->FillJets(genZ, ljets, 0, evtW);

  float dRcut = TMath::Pi()/2;
  // ==== Start of Actual Selections ====
  // Note: In each case, we have Higgs forced to CC, so we never need
  // a secondary check for the Higgs jets.
  std::vector<JetObj> Hjets = GenObj::get_proper_jets(cjets, genHiggs, dRcut);
  h_Higgs_nJet->Fill(Hjets.size(), evtW);
  if (Hjets.size() >= 2) {

    h_evt_cutflow->Fill(3); // passed jet requirement (H)
    
    HObj H(Hjets); // Reconstruct the Higgs boson

    std::vector<JetObj> Zjets = GenObj::get_proper_jets(bjets, genZ, dRcut);
    if (Zjets.size() >= 2) {
      h_evt_cutflow->Fill(4); //passed jet requirement (Z)
      ZObj Z(Zjets); // Reconstruct the Z boson
      
      h_VH->Fill(H, Z, evtW);

    }//end-Z-Selection
  }//end-Higgs-Selection

} //end Process



void VH_selection::Terminate(TList* mergedList, std::string outFileName) {
  
  /*TList* aList = new TList() ;
  TParameter<double>* a = new TParameter<double>("lep_eta",CUTS.Get<float>("lep_eta")) ;
  aList->Add(a) ;
  a = new TParameter<double>("lep_pt",CUTS.Get<float>("lep_pt")) ;
  aList->Add(a) ;
  */

  //TFile* f = new TFile(outFileName.c_str(), "UPDATE") ;
  //aList->Write("VH_selectionCuts",TObject::kSingleKey) ;

  //f->Close() ;
}
