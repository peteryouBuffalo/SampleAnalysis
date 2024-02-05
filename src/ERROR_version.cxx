#define VH_selection_cxx

/*/////////////////////////////////////////////////////////////////////////////
IMPORT STATEMENTS
/////////////////////////////////////////////////////////////////////////////*/

#include <math.h>

#include "TList.h"
#include "TParameter.h"

#include "VH_selection.h"
#include "Global.h"
#include "Obj.cxx"

/*/////////////////////////////////////////////////////////////////////////////
MAIN SELECTOR METHODS
/////////////////////////////////////////////////////////////////////////////*/

VH_selection::~VH_selection() { 

}

/////////////////////////////////////////////
// Slave Begin
/////////////////////////////////////////////
void VH_selection::SlaveBegin(Reader *r) {

  // Set up the event plot
  h_evt = new TH1D("Nevt", "", 4, -1.5, 2.5);	// bin 1 - total negative w evt
      						// bin 2 - total positive w evt
                                        	// bin 3 - total event weight by		  						         // genWeight (= bin2 - bin1)    						  						  			
  // Set up the VHPlots instances
  h_VH_MC = new VHPlots("VbbHcc_MC");
  h_VH_all = new VHPlots("VbbHcc_all");

  h_VH_tags = new VHPlots("VbbHcc_tags");
  h_VH_algo = new VHPlots("VbbHcc_algo");
  h_VH_both = new VHPlots("VbbHcc_both");

  h_VH_tags_all = new VHPlots("VbbHcc_tags_all");
  h_VH_algo_all = new VHPlots("VbbHcc_algo_all");
  h_VH_both_all = new VHPlots("VbbHcc_both_all");

  h_VH_tags_afterTag = new VHPlots("VbbHcc_tags_afterTag");
  h_VH_algo_afterTag = new VHPlots("VbbHcc_algo_afterTag");
  h_VH_both_afterTag = new VHPlots("VbbHcc_both_afterTag");  

  // Set up the EffPlots instances
  h_eff_tags = new EffPlots("VbbHcc_tags");
  h_eff_algo = new EffPlots("VbbHcc_algo");
  h_eff_both = new EffPlots("VbbHcc_both");

  // Set up the CutFlows (for events)
  h_evt_MC_cutflow = new TH1D("VbbHcc_MC_CutFlow", "", 2, 0, 2);
  h_evt_MC_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_MC_cutflow->GetXaxis()->SetBinLabel(2, "Passed daughter selection");
  
  h_evt_tags_cutflow = new TH1D("VbbHcc_tags_CutFlow", "", 7, 0, 7);
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(2, "jet cuts");
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(3, "b-tags");
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(4, "c-tags");
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(5, "MET cut");
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(6, "pT(Z) cut");
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(7, "dPhi cut");
  
  h_evt_algo_cutflow = new TH1D("VbbHcc_algo_CutFlow", "", 7, 0, 7);
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(2, "jet cuts");
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(3, "b-tags");
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(4, "c-tags");
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(5, "MET cut");
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(6, "pT(Z) cut");
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(7, "dPhi cut");
  
  h_evt_both_cutflow = new TH1D("VbbHcc_both_CutFlow", "", 7, 0, 7);
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(2, "jet cuts");
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(3, "b-tags");
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(4, "c-tags");
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(5, "MET cut");
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(6, "pT(Z) cut");
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(7, "dPhi cut");

  // Set up the CutFlows (for obj selections)
  h_jet_cutflow = new TH1D("VbbHcc_CutFlow_jets", "", 4, 0, 4);
  h_jet_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_jet_cutflow->GetXaxis()->SetBinLabel(2, "pT cut");
  h_jet_cutflow->GetXaxis()->SetBinLabel(3, "eta cut");
  h_jet_cutflow->GetXaxis()->SetBinLabel(4, "iso req");
  
  h_elec_cutflow = new TH1D("VbbHcc_CutFlow_elec", "", 5, 0, 5);
  h_elec_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_elec_cutflow->GetXaxis()->SetBinLabel(2, "pT cut");
  h_elec_cutflow->GetXaxis()->SetBinLabel(3, "eta cut");
  h_elec_cutflow->GetXaxis()->SetBinLabel(4, "etaSC cut");
  h_elec_cutflow->GetXaxis()->SetBinLabel(5, "loose ID cut");
  
  h_muon_cutflow = new TH1D("VbbHcc_CutFlow_muon", "", 5, 0, 5);
  h_muon_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_muon_cutflow->GetXaxis()->SetBinLabel(2, "pT cut");
  h_muon_cutflow->GetXaxis()->SetBinLabel(3, "eta cut");
  h_muon_cutflow->GetXaxis()->SetBinLabel(4, "loose ID cut");
  h_muon_cutflow->GetXaxis()->SetBinLabel(5, "iso cut");

  // Set up the miscellaneous histograms.
  h_nJet = new TH1D("VbbHcc_nJet", "", 13, -0.5, 12.5);
  h_nBjet_loose = new TH1D("VbbHcc_nBjet_loose", "", 10, -0.5, 9.5);
  h_nCjet_loose = new TH1D("VbbHcc_nCjet_loose", "", 10, -0.5, 9.5);
  h_nBjet_medium = new TH1D("VbbHcc_nBjet_medium", "", 10, -0.5, 9.5);
  h_nCjet_medium = new TH1D("VbbHcc_nCjet_medium", "", 10, -0.5, 9.5);
  h_bScore = new TH1D("VbbHcc_bScore", "", 20, 0., 1.);
  h_cScore = new TH1D("VbbHcc_cScore", "", 20, 0., 1.);
  h_btag_v_ctag = new TH2D("VbbHcc_btag_v_ctag", "", 20, 0., 1., 20, 0., 1.);

  h_nJet_all = new TH1D("VbbHcc_nJet_all", "", 13, -0.5, 12.5);
  h_nBjet_loose_all = new TH1D("VbbHcc_nBjet_loose_all", "", 10, -0.5, 9.5);
  h_nCjet_loose_all = new TH1D("VbbHcc_nCjet_loose_all", "", 10, -0.5, 9.5);
  h_nBjet_medium_all = new TH1D("VbbHcc_nBjet_medium_all", "", 10, -0.5, 9.5);
  h_nCjet_medium_all = new TH1D("VbbHcc_nCjet_medium_all", "", 10, -0.5, 9.5);
  h_bScore_all = new TH1D("VbbHcc_bScore_all", "", 20, 0., 1.);
  h_cScore_all = new TH1D("VbbHcc_cScore_all", "", 20, 0., 1.);
  h_btag_v_ctag_all = new TH2D("VbbHcc_btag_v_ctag_all", "", 20, 0., 1., 20, 0., 1.);

  // Add histograms to fOutput so they can be saved in Processor::SlaveTermifor(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
nate
  r->GetOutputList()->Add(h_evt);

  std::vector<TH1*> tmp = h_VH_MC->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_all->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  tmp = h_VH_tags->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_algo->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_both->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  
  tmp = h_VH_tags_all->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_algo_all->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_both_all->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  tmp = h_VH_tags_afterTag->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_algo_afterTag->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_both_afterTag->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  tmp = h_eff_tags->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_eff_algo->returnHisto();for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_eff_both->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  
  r->GetOutputList()->Add(h_evt_MC_cutflow);
  r->GetOutputList()->Add(h_evt_tags_cutflow);
  r->GetOutputList()->Add(h_evt_algo_cutflow);
  r->GetOutputList()->Add(h_evt_both_cutflow);
  r->GetOutputList()->Add(h_jet_cutflow);
  r->GetOutputList()->Add(h_elec_cutflow);
  r->GetOutputList()->Add(h_muon_cutflow);

  r->GetOutputList()->Add(h_nJet);
  r->GetOutputList()->Add(h_nBjet_loose);
  r->GetOutputList()->Add(h_nCjet_loose);
  r->GetOutputList()->Add(h_nBjet_medium);
  r->GetOutputList()->Add(h_nCjet_medium);
  r->GetOutputList()->Add(h_bScore);
  r->GetOutputList()->Add(h_cScore);
  r->GetOutputList()->Add(h_btag_v_ctag);

  r->GetOutputList()->Add(h_nJet_all);
  r->GetOutputList()->Add(h_nBjet_loose_all);
  r->GetOutputList()->Add(h_nCjet_loose_all);
  r->GetOutputList()->Add(h_nBjet_medium_all);
  r->GetOutputList()->Add(h_nCjet_medium_all);
  r->GetOutputList()->Add(h_bScore_all);
  r->GetOutputList()->Add(h_cScore_all);
  r->GetOutputList()->Add(h_btag_v_ctag_all);

}// end SlaveBegin

/////////////////////////////////////////////
// Custom Methods
/////////////////////////////////////////////

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
std::vector<std::vector<int> > VH_selection::DauIdxs_ZH(Reader* r){

  // Store the indices of the H and Z daughters
  std::vector<std::vector<int> > dauIdxs;
  std::vector<int> dauIdxsZ;
  std::vector<int> dauIdxsH;

  // Loop over the gen particles
  for (unsigned i = 0; i < *(r->nGenPart); ++i) {
    
    int mIdx = (r->GenPart_genPartIdxMother)[i]; // mother index
    int mId  = (r->GenPart_pdgId)[mIdx];         // mother pdg ID
    int flav = (r->GenPart_pdgId)[i];            // ptcl pdg ID

    // b quark check (id = 5, and mother ID = 23 (Z))
    if (abs(flav) == 5 && mIdx > -1 && abs(mId) == 23) 
      dauIdxsZ.push_back(i);

    // c quark check (id = 4, and mother ID = 25 (H))
    if (abs(flav) == 4 && mIdx > -1 && abs(mId) == 25) 
      dauIdxsH.push_back(i);
  }

  // Push back and return the proper IDs
  dauIdxs.push_back(dauIdxsZ);
  dauIdxs.push_back(dauIdxsH);

  return dauIdxs;

} // end DauIdxs_H
#endif

/////////////////////////////////////////////
// Process
/////////////////////////////////////////////
void VH_selection::Process(Reader* r) {

  //===========================================================================
  // WEIGHTS & GEN INFORMATION
  //===========================================================================
  
  float genWeight = 1.;
  float puSF = 1.;
  float l1preW = 1.;

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
  if (*(r->genWeight) < 0) genWeight = -1.;
  if (*(r->genWeight) == 0) {
    genWeight = 0;
    h_evt->Fill(0);
  }
  if (*(r->genWeight) < 0) h_evt->Fill(-1);
  if (*(r->genWeight) > 0) h_evt->Fill(1);
  // use central gen weight to normalize gen weight
  if (m_centralGenWeight != 0) genWeight = *(r->genWeight)/m_centralGenWeight;
  puSF = PileupSF(*(r->Pileup_nTrueInt));
#endif

  h_evt->Fill(2, genWeight);

#if defined(MC_2016) || defined(MC_2017)
  l1preW = *(r->L1PreFiringWeight_Nom);
  if (m_l1prefiringType == "l1prefiringu") l1preW = *(r->L1PreFiringWeight_Up);
  if (m_l1prefiringType == "l1prefiringd") l1preW = *(r->L1PreFiringWeight_Dn);
#endif

#if defined(DATA_2016) || defined(DATA_2017) || defined(DATA_2018)
  h_evt->Fill(-1);
  if (!m_lumiFilter.Pass(*(r->run), *(r->luminosityBlock))) return;
  h_evt->Fill(1);
#endif

  float evtW = 1.;
  if (!m_isData) evtW *= genWeight * puSF * l1preW;

  //===========================================================================
  // GET OBJECTS
  //===========================================================================

  // ==== JETS ====
  std::vector<JetObj> jets;
  for (unsigned int i = 0; i < *(r->nJet); ++i) {

    int jetFlav = -1;
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
    jetFlav = (r->Jet_hadronFlavour)[i];
#endif

    JetObj jet((r->Jet_pt)[i], (r->Jet_eta)[i], (r->Jet_phi)[i], (r->Jet_mass)[i],
      jetFlav, (r->Jet_btagDeepFlavB)[i], (r->Jet_puIdDisc)[i]);

#if defined(NANOAODV9)
    jet.m_deepCvL = (r->Jet_btagDeepFlavCvL)[i];
#endif

#if defined(NANOAODV7)
    jet.m_deepCvL = (r->FatJet_btagDDCvL)[i];
#endif

    jets.push_back(jet);

  }//end-jet
  h_VH_all->FillJets(jets, evtW);

  // ==== ELECTRONS ====
  std::vector<LepObj> elecs;
  for (unsigned int i = 0; i < *(r->nElectron); ++i) {
    
    // Produce an electron from the information.
    h_elec_cutflow->Fill(0.5, genWeight); // All electrons

    float etaSC = (r->Electron_eta)[i] - (r->Electron_deltaEtaSC)[i];
    LepObj elec((r->Electron_pt)[i], (r->Electron_eta)[i], etaSC,
           (r->Electron_phi)[i], (r->Electron_mass)[i], i,
           (r->Electron_charge)[i], 0);

    // CUTS #1, 2 - pT and eta cuts ///////////////////////
    if (elec.m_lvec.Pt() < CUTS.Get<float>("lep_pt1")) continue;
    h_elec_cutflow->Fill(1.5, genWeight); // pass pT cut
    if (fabs(elec.m_lvec.Eta()) > CUTS.Get<float>("lep_eta")) continue;
    h_elec_cutflow->Fill(2.5, genWeight); // pass eta cut

    // CUT #3 - make a cut based on SC ////////////////////
    if (fabs(etaSC) < 1.566 && fabs(etaSC) > 1.442) continue;
    h_elec_cutflow->Fill(3.5, genWeight); // pass SC cut

    // CUT #4 - make sure to only keep proper electrons ///
    int elecID = r->Electron_cutBased[i];
    if (elecID < 2) continue; // loose electron ID
    h_elec_cutflow->Fill(4.5, genWeight); // pass ID cut

    // Add the electrons to the list
    elecs.push_back(elec);
  }//end-elec

  // ==== MUONS =====
  std::vector<LepObj> muons;
  for (unsigned int i = 0; i < *(r->nMuon); ++i) {
    
    // Produce a muon from the information
    h_muon_cutflow->Fill(0.5, genWeight); // all muons
    LepObj muon((r->Muon_pt)[i], (r->Muon_eta)[i], -1, (r->Muon_phi)[i],
           (r->Muon_mass)[i], i, (r->Muon_charge)[i], 
           (r->Muon_pfRelIso04_all)[i]);

    // CUTS #1,2 - pT and eta cuts ////////////////////////
    if (muon.m_lvec.Pt() < CUTS.Get<float>("lep_pt1")) continue;
    h_muon_cutflow->Fill(1.5, genWeight); // pass pT cut
    if (muon.m_lvec.Eta() > CUTS.Get<float>("lep_eta")) continue;
    h_muon_cutflow->Fill(2.5, genWeight); // pass eta cut
    
    // CUT #3 - check loose ID ////////////////////////////
    if (r->Muon_looseId[i] <= 0) continue;
    h_muon_cutflow->Fill(3.5, genWeight); // pass ID cut

    // CUT #4 - isolation cut /////////////////////////////
    if (muon.m_iso > CUTS.Get<float>("muon_iso")) continue;
    h_muon_cutflow->Fill(4.5, genWeight); // pass iso cut

    // Add the muon to the list
    muons.push_back(muon);
  }//end-muon

  //===========================================================================
  // START EVENT SELECTION
  //===========================================================================

  // All CutFlows need to show the total events.
  h_evt_MC_cutflow->Fill(0.5, genWeight);   // MC Truth
  h_evt_tags_cutflow->Fill(0.5, genWeight); // Tagging Only
  h_evt_algo_cutflow->Fill(0.5, genWeight); // Matching Prioritized
  h_evt_both_cutflow->Fill(0.5, genWeight); // Tagging Prioritized

  // Fill MET spectrum
  h_VH_all->FillMET(*(r->MET_pt), evtW);

  // Capture the MC truth values so taht we can use them for reference.
  bool canCompareToMC = false;
  std::vector<std::vector<int> > dauIdxs;
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)

  // Make sure we have two daughters for each particle.
  int idZ = 0, idH = 1;
  dauIdxs = DauIdxs_ZH(r);
  std::vector<JetObj> gen_cjets, gen_bjets;
  if (dauIdxs[idZ].size() == 2 && dauIdxs[idH].size() == 2) {
    
    canCompareToMC = true;
    h_evt_MC_cutflow->Fill(1.5, genWeight); // pass # daughters
    int idx1_Z = dauIdxs[idZ][0];
    int idx2_Z = dauIdxs[idZ][1];
    int idx1_H = dauIdxs[idH][0];
    int idx2_H = dauIdxs[idH][1];

    // Reconstruct from the GenPart data.
    // Two b quarks --> ZObj
    JetObj b0((r->GenPart_pt)[idx1_Z], (r->GenPart_eta)[idx1_Z],
      (r->GenPart_phi)[idx1_Z], (r->GenPart_mass)[idx1_Z], 5, 0., 0.);
    JetObj b1((r->GenPart_pt)[idx2_Z], (r->GenPart_eta)[idx2_Z],
      (r->GenPart_phi)[idx2_Z], (r->GenPart_mass)[idx2_Z], 5, 0., 0.);
    std::vector<JetObj> bjets{b0, b1};
    gen_bjets.push_back(b0);
    gen_bjets.push_back(b1);
    ZObj Z(bjets);

    // Two c quarks --> HObj
    JetObj c0((r->GenPart_pt)[idx1_H], (r->GenPart_eta)[idx1_H],
      (r->GenPart_phi)[idx1_H], (r->GenPart_mass)[idx1_H], 4, 0., 0.);
    JetObj c1((r->GenPart_pt)[idx2_H], (r->GenPart_eta)[idx2_H],
      (r->GenPart_phi)[idx2_H], (r->GenPart_mass)[idx2_H], 4, 0., 0.);
    std::vector<JetObj> cjets{c0, c1};
    gen_cjets.push_back(c0);
    gen_cjets.push_back(c1);
    HObj H(cjets);

    // Fill the histograms
    h_VH_MC->Fill(H, Z, evtW); 

  }//end-MC-truth

#endif

  // For the remainder of the selections, we only want to have jets
  // that meet our criteria. Let's make sure we have enough jets to
  // meet our requirements:
  // 1. pT(j) > 30 GeV
  // 2. |eta| < 2.5
  // 3. dR(small-R jet, lepton) < 0.4 = discard
  std::vector<JetObj> selected_jets;
  for (unsigned int i = 0; i < jets.size(); ++i) {

    h_jet_cutflow->Fill(0.5, genWeight); // all jets
    TLorentzVector vec = jets.at(i).m_lvec;
    
    if (vec.Pt() < 30) continue;
    h_jet_cutflow->Fill(1.5, genWeight); // passed pT cut

    if (fabs(vec.Eta()) > 2.5) continue;
    h_jet_cutflow->Fill(2.5, genWeight); // passed eta cut

    // Check the electrons & muons for overlap with the jets.
    // If dR(jet, lepton) < 0.4, discard the jet. 
    bool should_discard = false;
    for (unsigned int j = 0; j < elecs.size(); ++j) {
      float dR = vec.DeltaR(elecs.at(j).m_lvec);
      if (dR < 0.4) { should_discard = true; break; }
    }//end-elec
    if (should_discard) continue;

    should_discard = false;
    for (unsigned int j = 0; j < muons.size(); ++j) {
      float dR = vec.DeltaR(muons.at(j).m_lvec);
      if (dR < 0.4) { should_discard = true; break; }
    }
    if (should_discard) continue;
    h_jet_cutflow->Fill(3.5, genWeight); // passed iso

    // Add jet to our selected jet list
    selected_jets.push_back(jets.at(i));

  }//end-selected-jets

  h_VH_all->FillJets_selected(selected_jets, evtW);

  /****************************************************************************
  * JET ANALYSIS                                                              *
  ****************************************************************************/
  /*h_nJet->Fill(selected_jets.size(), evtW);
  h_nJet_all->Fill(jets.size(), evtW);  

  // Check our selected jets against the loose & medium WP.
  int nB_loose = 0, nB_medium = 0;
  int nC_loose = 0, nC_medium = 0;
  for (unsigned int i = 0; i < selected_jets.size(); ++i) {

    // Pull the scores from the jet
    float btag = selected_jets[i].m_deepCSV;
    float ctag = selected_jets[i].m_deepCvL;
    h_bScore->Fill(btag, evtW);
    h_cScore->Fill(ctag, evtW);
    h_btag_v_ctag->Fill(btag, ctag, evtW);

    // Check b-tagging against our WPs
    if (btag > 0.0508) nB_loose++;
    if (btag > 0.3) nB_medium++;

    // Check c-tagging against our WPs
    if (ctag > 0.327) nC_loose++;
    if (ctag > 0.370) nC_medium++;
  }
  h_nBjet_loose->Fill(nB_loose);
  h_nBjet_medium->Fill(nB_medium);
  h_nCjet_loose->Fill(nC_loose);
  h_nCjet_medium->Fill(nC_medium);

  // Check ALL jets against the loose & medium WP.
  nB_loose = 0; nB_medium = 0;
  nC_loose = 0; nC_medium = 0;
  for (unsigned int i = 0; i < jets.size(); ++i) {

    // Pull the scores from the jet
    float btag = jets[i].m_deepCSV;
    float ctag = jets[i].m_deepCvL;
    h_bScore_all->Fill(btag, evtW);
    h_cScore_all->Fill(ctag, evtW);
    h_btag_v_ctag_all->Fill(btag, ctag, evtW);
    
    // Check b-tagging against our WPs
    if (btag > 0.0508) nB_loose++;
    if (btag > 0.3) nB_medium++;

    // Check c-tagging against our WPs
    if (ctag > 0.327) nC_loose++;
    if (ctag > 0.370) nC_medium++;
  }
  h_nBjet_loose_all->Fill(nB_loose);
  h_nBjet_medium_all->Fill(nB_medium);
  h_nCjet_loose_all->Fill(nC_loose);
  h_nCjet_medium_all->Fill(nC_medium);*/

  if (selected_jets.size() >= 4) {

    h_evt_tags_cutflow->Fill(1.5, genWeight); // passed jet selection
    h_evt_algo_cutflow->Fill(1.5, genWeight);
    h_evt_both_cutflow->Fill(1.5, genWeight);

    /**************************************************************************
    * CASE #2 - TAGGING ONLY                                                  *
    **************************************************************************/
    
    // Select two jets with the largets btag value and then
    // make pass ~medium WP (Jet_btagDeepFlavB > 0.3)
    std::vector<JetObj> jets2 = selected_jets;
    std::vector<JetObj> bjets;
    float csv0 = -2000.0, csv1 = -2000.0; 
    int bIdx0 = -1, bIdx1 = -1;
 
    for (unsigned int i = 0; i < jets2.size(); ++i) {
      float jet_csv = jets2.at(i).m_deepCSV;
      if (jet_csv > csv0) { csv0 = jet_csv; bIdx0 = i; }
    } 
    bjets.push_back(jets2[bIdx0]);
    jets2.erase(jets2.begin() + bIdx0);

    for (unsigned int i = 0; i < jets2.size(); ++i) {
      float jet_csv = jets2.at(i).m_deepCSV;
      if (jet_csv > csv1) { csv1 = jet_csv; bIdx1 = i; }
    }
    bjets.push_back(jets2[bIdx1]);
    jets2.erase(jets2.begin() + bIdx1);

    if (csv0 > 0.05 && csv1 > 0.05) {
 
      h_evt_tags_cutflow->Fill(2.5, genWeight); // pass b-cuts
      ZObj Z(bjets);
      std::vector<JetObj> cjets;

      // Do the same matching process but with c-jets
      float cvl0 = -2000.0, cvl1 = -2000.0;
      int cIdx0 = -1, cIdx1 = -1;
       
      for (unsigned int i = 0; i < jets2.size(); ++i) {
        float jet_cvl = jets2.at(i).m_deepCvL;
        if (jet_cvl > cvl0) { cvl0 = jet_cvl; cIdx0 = i; }
      }
      cjets.push_back(jets2[cIdx0]);
      jets2.erase(jets2.begin() + cIdx0);

      for (unsigned int i = 0; i < jets2.size(); ++i) {
        float jet_cvl = jets2.at(i).m_deepCvL;
        if (jet_cvl > cvl1) { cvl1 = jet_cvl; cIdx1 = i; }
      }
      cjets.push_back(jets2[cIdx1]);
      jets2.erase(jets2.begin() + cIdx1);

      if (cvl0 > 0.33 && cvl1 > 0.33) {

        h_evt_tags_cutflow->Fill(3.5, genWeight); // pass c-cuts
        HObj H(cjets);

        // Plot all options of H & Z even if they don't survive
        // our cuts. We want to have them for reference of what
        // cuts we should make.
        h_VH_tags_all->Fill(H, Z, evtW);

        if(*(r->MET_pt) < 140) {
          h_evt_tags_cutflow->Fill(4.5, genWeight); // pass MET cut
          if (Z.m_lvec.Pt() > 50) {
            h_evt_tags_cutflow->Fill(5.5, genWeight); // pass pT(Z) cut
            float dPhi = Z.m_lvec.DeltaPhi(H.m_lvec);
            if (fabs(dPhi) > 2.5) {

              h_evt_tags_cutflow->Fill(6.5, genWeight); // pass dPhi cut
              h_VH_tags->Fill(H, Z, evtW);

              // If we're in an MC file, let's check to match our
              // selected objets to the MC truth. This should
              // help us see how accurate we are.            
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)

              if (canCompareToMC)
                h_eff_tags->Fill(H, Z, gen_cjets, gen_bjets, evtW);

#endif

            }//end-dPhi
          }//end-pT(Z)
        }//end-MET       

      }//end-c-tag

    }//end-b-tag

    /**************************************************************************
    * CASE #3 - MASS MATCHING PRIORITIZED                                     *
    **************************************************************************/
    std::vector<JetObj> jets3 = selected_jets;
    std::sort(jets3.begin(), jets3.end(), std::greater<JetObj>());

    // Create the objects for the distance calculations and make sure we sort
    // them in ascending order. (All necessary calculations for alogrithms 
    // are handled within the DObj class.)
    DObj d0(jets3, 0, 1, 2, 3);      // H(0,1) ; Z(2,3)
    DObj d1(jets3, 0, 2, 1, 3);      // H(0,2) ; Z(1,3)
    DObj d2(jets3, 0, 3, 1, 2);      // H(0,3) ; Z(1,2)
    std::vector<DObj> distances{ d0, d1, d2 };
    std::sort(distances.begin(), distances.end());

    // Determine the distance between the closest two.
    float deltaD = fabs(distances[0].m_d - distances[1].m_d);

    // Determine which pairing we want to use based on this deltaD.
    // If we are outside the resolution window (30 GeV), we can 
    // just choose the closest pair. Otherwise, we need to make a 
    // logical choice between d0 and d1.
    DObj chosenPair = distances[0];
    if (deltaD < 30) {
      float pt0 = distances[0].HPt();
      float pt1 = distances[1].HPt();
      int idx = 0;
      if (pt1 > pt0) idx = 1;
      chosenPair = distances[idx];
    }

    // Plot all the possible combinations for some references
    // before we go through the cuts.
    std::vector<JetObj> cjets0 { jets3[d0.m_hIdx0], jets3[d0.m_hIdx1] };
    std::vector<JetObj> bjets0 { jets3[d0.m_zIdx0], jets3[d0.m_zIdx1] };
    HObj H0(cjets0); ZObj Z0(bjets0);
    h_VH_algo_all->Fill(H0, Z0, evtW);

    std::vector<JetObj> cjets1 { jets3[d1.m_hIdx0], jets3[d1.m_hIdx1] };
    std::vector<JetObj> bjets1 { jets3[d1.m_zIdx0], jets3[d1.m_zIdx1] };
    HObj H1(cjets1); ZObj Z1(bjets1);
    h_VH_algo_all->Fill(H1, Z1, evtW);

    std::vector<JetObj> cjets2 { jets3[d2.m_hIdx0], jets3[d2.m_hIdx1] };
    std::vector<JetObj> bjets2 { jets3[d2.m_zIdx0], jets3[d2.m_zIdx1] };
    HObj H2(cjets2); ZObj Z2(bjets2);
    h_VH_algo_all->Fill(H2, Z2, evtW);

    // Reconstruct the objects here for various uses.
    std::vector<JetObj> bjets_algo;
    bjets_algo.push_back(jets3[chosenPair.m_zIdx0]);
    bjets_algo.push_back(jets3[chosenPair.m_zIdx1]);
    ZObj Z(bjets_algo);

    std::vector<JetObj> cjets_algo;
    cjets_algo.push_back(jets3[chosenPair.m_hIdx0]);
    cjets_algo.push_back(jets3[chosenPair.m_hIdx1]);
    HObj H(cjets_algo);
 
    // Now, check to see if we pass the tagging requirements
    // and our other cuts.
    if (chosenPair.Z_has_bjets()) {
      h_evt_algo_cutflow->Fill(2.5, genWeight); // pass b-tag
      if (chosenPair.H_has_cjets()) {
        h_evt_algo_cutflow->Fill(3.5, genWeight); // pass c-tag
        h_VH_algo_afterTag->Fill(H, Z, evtW);

        if (*(r->MET_pt) < 140) {
          h_evt_algo_cutflow->Fill(4.5, genWeight); // pass MET
          if (chosenPair.ZPt() > 50) {
            h_evt_algo_cutflow->Fill(5.5, genWeight); // pass pT(Z)
            float dPhi = fabs(chosenPair.DPhi());
            if (dPhi > 2.5) {

              h_evt_algo_cutflow->Fill(6.5, genWeight); // pass dPhi
              
              // Reconstruct the objects and fill our histograms.
              h_VH_algo->Fill(H, Z, evtW);

              // If we're in a MC file, let's check to match our
              // selected objects to teh MC truth. This should
              // help us see how accurate we are.
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
              if (canCompareToMC) 
                h_eff_algo->Fill(H, Z, gen_cjets, gen_bjets, evtW);
#endif

            }//end dPhi
          }//end pT(Z)
        }//end MET
      }//end-c-tag
    }//end-b-tag

    /**************************************************************************
    * CASE #4 - TAGGING PRIORITIZED                                           *
    **************************************************************************/
    std::vector<JetObj> jets4 = selected_jets;
    std::sort(jets4.begin(), jets4.end(), std::greater<JetObj>());

    // Check to see which jets pass which tagging
    std::vector<int> bIndices, cIndices;
    bool btags[4] = { false, false, false, false };
    bool ctags[4] = { false, false, false, false };
    for (int i = 0; i < 4; ++i) {
      if (jets4[i].m_deepCSV > 0.05) { btags[i] = true; bIndices.push_back(i); }
      if (jets4[i].m_deepCvL > 0.33) { ctags[i] = true; cIndices.push_back(i); }
    }    

    // Make sure we have no jets that pass both tags.
    bool properly_chosen = true;
    for (int i = 0; i < 4; ++i) {
      if (ctags[i] && btags[i]) {
        properly_chosen = false; break;
      }
    }

    // Make sure we have 2 b-jets and 2 c-jets
    if (properly_chosen && bIndices.size() == 2 && cIndices.size() == 2) {
    
      h_evt_both_cutflow->Fill(2.5, genWeight); // pass b-tag
      h_evt_both_cutflow->Fill(3.5, genWeight); // pass c-tag

      // Get the jets into proper lists. Reconstruct the bosons.
      std::vector<JetObj> bjets, cjets;
      bjets.push_back(jets4[bIndices[0]]);
      bjets.push_back(jets4[bIndices[1]]);
      cjets.push_back(jets4[cIndices[0]]);
      cjets.push_back(jets4[cIndices[1]]);
      ZObj Z(bjets); HObj H(cjets);

      // Go through the rest of the cuts.
      if (*(r->MET_pt) < 140) {
        h_evt_both_cutflow->Fill(4.5, genWeight); // pass MET
        if (Z.m_lvec.Pt() > 50) {
          h_evt_both_cutflow->Fill(5.5, genWeight); // pass pT(Z)
          float dPhi = fabs(Z.m_lvec.DeltaPhi(H.m_lvec));
          if (dPhi > 2.5) {

            h_evt_both_cutflow->Fill(6.5, genWeight); // pass dPhi
            h_VH_both->Fill(H, Z, evtW);             

            // If we're in a MC file, let's check to match our
            // selected objets to the MC truth. This should
            // help us see how accurate we are.
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
            if (canCompareToMC)
              h_eff_both->Fill(H, Z, gen_cjets, gen_bjets, evtW);
#endif

          }//end-dPhi
        }//end-pT(Z)
      }//end MET
    }//end-tagging
    // If we don't, we have to go through the mass matching algorithm.
    else {

      // Create the objects for the distance calculations & make sure we sort
      // them in ascending order. (All necessary calculations for algorithms
      // are handled within the DObj class.)
      DObj d0(jets4, 0, 1, 2, 3);    // H(0,1) ; Z(2,3)
      DObj d1(jets4, 0, 2, 1, 3);    // H(0,2) ; Z(1,3)
      DObj d2(jets4, 0, 3, 1, 2);    // H(0,3) ; Z(1,2)
      std::vector<DObj> distances { d0, d1, d2 };
      std::sort(distances.begin(), distances.end());

      // Determine the distance between the closest two.
      float deltaD = fabs(distances[0].m_d - distances[1].m_d);

      // Determine which pairing we want to use based on this deltaD.
      // If we are outside the resolution window (30 GeV), we can
      // just choose the closest pair. Otherwise, we need to make a
      // logical choice between d0 and d1.
      DObj chosenPair = distances[0];
      if (deltaD < 30) {
        float pt0 = distances[0].HPt();
        float pt1 = distances[1].HPt();
        int idx = 0;
        if (pt1 > pt0) idx = 1;
        chosenPair = distances[idx];
      }

      // Plot all the possible combinations for some references
      // before we go through the cuts.
      std::vector<JetObj> cjets0 { jets4[d0.m_hIdx0], jets4[d0.m_hIdx1] };
      std::vector<JetObj> bjets0 { jets4[d0.m_zIdx0], jets4[d0.m_zIdx1] };
      HObj H0(cjets0); ZObj Z0(bjets0);
      h_VH_both_all->Fill(H0, Z0, evtW);

      std::vector<JetObj> cjets1 { jets4[d1.m_hIdx0], jets4[d1.m_hIdx1] };
      std::vector<JetObj> bjets1 { jets4[d1.m_zIdx0], jets4[d1.m_zIdx1] };
      HObj H1(cjets1); ZObj Z1(bjets1);
      h_VH_both_all->Fill(H1, Z1, evtW);

      std::vector<JetObj> cjets2 { jets4[d2.m_hIdx0], jets4[d2.m_hIdx1] };
      std::vector<JetObj> bjets2 { jets4[d2.m_zIdx0], jets4[d2.m_zIdx1] };
      HObj H2(cjets2); ZObj Z2(bjets2);
      h_VH_both_all->Fill(H2, Z2, evtW);

      // Reconstruct our objects here for use.
      std::vector<JetObj> bjets_both;
      bjets_both.push_back(jets4[chosenPair.m_zIdx0]);
      bjets_both.push_back(jets4[chosenPair.m_zIdx1]);
      ZObj Z(bjets_both);

      std::vector<JetObj> cjets_both;
      cjets_both.push_back(jets4[chosenPair.m_hIdx0]);
      cjets_both.push_back(jets4[chosenPair.m_hIdx1]);
      HObj H(cjets_both);

      // Now, check to see if we pass the tagging requirements
      // and our other cuts.
      if (chosenPair.Z_has_bjets()) {
        h_evt_both_cutflow->Fill(2.5, genWeight); // pass b-tag
        if (chosenPair.H_has_cjets()) {
          h_evt_both_cutflow->Fill(3.5, genWeight); // pass c-tag
 
          // Fill here so we have all of them that are tagged properly.
          h_VH_both_afterTag->Fill(H, Z, evtW);          

          if(*(r->MET_pt) < 140) {
            h_evt_both_cutflow->Fill(4.5, genWeight); // pass MET
            if (chosenPair.ZPt() > 50) {
              h_evt_both_cutflow->Fill(5.5, genWeight); // pass pT(Z)
              float dPhi = fabs(chosenPair.DPhi());
              if (dPhi > 2.5) {

                h_evt_both_cutflow->Fill(6.5, genWeight); // pass dPhi

                // Fill our histograms.
                h_VH_both->Fill(H, Z, evtW);

                // If we're in a MC file, let's check to match our
                // selected objects to the MC truth. This should
                // help us see how accurate we are.
#if defined(MC_2016) || defined(MC_2017) || defined (MC_2018)
                if (canCompareToMC)
                  h_eff_both->Fill(H, Z, gen_cjets, gen_bjets, evtW);
#endif

              }//end dPhi
            }//end pT(Z)
          }//end MET
        }//end-c-tag
      }//end-b-tag

    }//end-mass

  }//end-jets-cut
}// end Process

/////////////////////////////////////////////
// Terminate
/////////////////////////////////////////////
void VH_selection::Terminate(TList* mergedList, std::string outFileName) {

}
