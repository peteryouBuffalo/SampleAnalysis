#define VH_selection_cxx

#include <math.h>

#include "TList.h"
#include "TParameter.h"

#include "VH_selection.h"
#include "Global.h"
#include "Obj.cxx"

//VH_selection::VH_selection(bool isData) : Selector(isData), h_zee_jet(0), h_zmm_jet(0) {}

VH_selection::~VH_selection() {
  //if (h_VH) delete h_VH;
} 

void VH_selection::SlaveBegin(Reader* r) {

  // Set up the plots we want
  h_evt = new TH1D("Nevt","",4,-1.5,2.5) ; //bin 1: total negative weight events
                                           // bin 2: total positive weight events
                                           // bin 3: total event weighted by genWeight
                                           // (= bin2 - bin1, if genWeight are always -1,1

  h_VH = new VHPlots("VbbHcc");
  h_VH_tags = new VHPlots("VbbHcc_tags");
  h_VH_algo = new VHPlots("VbbHcc_algo");
  h_VH_both = new VHPlots("VbbHcc_both");
  h_VH_bothTags = new VHPlots("VbbHcc_bothTags");
  h_VH_bothAlgo = new VHPlots("VbbHcc_bothAlgo");

  h_eff_tags = new EffPlots("VbbHcc_tags");
  h_eff_algo = new EffPlots("VbbHcc_algo");
  h_eff_both = new EffPlots("VbbHcc_both");

  // Cut flow to select events for analysis
  h_evt_cutflow = new TH1D("VbbHcc_CutFlow", "", 2, 0, 2);
  h_evt_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_cutflow->GetXaxis()->SetBinLabel(2, "Passed daughter selection");

  h_evt_tags_cutflow = new TH1D("VbbHcc_tags_CutFlow", "", 7, 0, 7);
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(2, "Passed jet requirements");
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(3, "Passed b-tagging for Z candidate");
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(4, "Passed c-tagging for H candidate");
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(5, "Passed MET cut");
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(6, "Passed V pT cut");
  h_evt_tags_cutflow->GetXaxis()->SetBinLabel(7, "Passed dPhi cut");

  h_evt_algo_cutflow = new TH1D("VbbHcc_algo_CutFlow", "", 7, 0, 7);
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(2, "Passed jet requirements");
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(3, "Passed b-tagging for Z candidate");
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(4, "Passed c-tagging for H candidate");
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(5, "Passed MET cut");
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(6, "Passed V pT cut");
  h_evt_algo_cutflow->GetXaxis()->SetBinLabel(7, "Passed dPhi cut");

  h_evt_both_cutflow = new TH1D("VbbHcc_both_CutFlow", "", 7, 0, 7);
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(2, "Passed jet requirements");
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(3, "Passed b-tagging for Z candidate");
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(4, "Passed c-tagging for H candidate");
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(5, "Passed MET cut");
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(6, "Passed V pT cut");
  h_evt_both_cutflow->GetXaxis()->SetBinLabel(7, "Passed dPhi cut");

  // Cut flow to select jets for event
  h_jet_cutflow = new TH1D("VbbHcc_CutFlow_jets", "", 4, 0, 4);
  h_jet_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_jet_cutflow->GetXaxis()->SetBinLabel(2, "pT cut");
  h_jet_cutflow->GetXaxis()->SetBinLabel(3, "eta cut");
  h_jet_cutflow->GetXaxis()->SetBinLabel(4, "iso req");

  // Cut flow to select electrons
  h_elec_cutflow = new TH1D("VbbHcc_CutFlow_elec", "", 4, 0, 4);
  h_elec_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_elec_cutflow->GetXaxis()->SetBinLabel(2, "Passed phase-space cut");
  h_elec_cutflow->GetXaxis()->SetBinLabel(3, "Passed SC cut");
  h_elec_cutflow->GetXaxis()->SetBinLabel(4, "Passed id cut");  

  // Cut flow to select muons
  h_muon_cutflow = new TH1D("VbbHcc_CutFlow_muon", "", 4, 0, 4);
  h_muon_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_muon_cutflow->GetXaxis()->SetBinLabel(2, "Passed phase-space cut");
  h_muon_cutflow->GetXaxis()->SetBinLabel(3, "Passed looseID cut");
  h_muon_cutflow->GetXaxis()->SetBinLabel(4, "Passed iso cut");

  // Other Histograms we want
  h_CSV = new TH1D("VbbHcc_CSV", "", 200, 0, 1);
  h_CvL = new TH1D("VbbHcc_CvL", "", 200, 0, 1);

  h_Nselected = new TH1D("VbbHcc_Nselected", "", 10, -0.5, 9.5);
  h_Nbjet = new TH1D("VbbHcc_Nbjet", "", 10, -0.5, 9.5);
  h_Ncjet = new TH1D("VbbHcc_Ncjet", "", 10, -0.5, 9.5);
  h_Z_dM = new TH1D("VbbHcc_Z_dM", "", 100, 0, 100);
  h_H_dM = new TH1D("VbbHcc_H_dM", "", 100, 0, 100);

  //h_HZ0 = new TH1D("VbbHcc_HZ0", "", 100, 0, 100);
  //h_HZ1 = new TH1D("VbbHcc_HZ1", "", 100, 0, 100);
  //h_HZ2 = new TH1D("VbbHcc_HZ2", "", 100, 0, 100);
  //h_dH = new TH1D("VbbHcc_dH", "", 100, 0, 100);

  h_tags_MH_v_MZ = new TH2D("VbbHcc_tags_MH_v_MZ", "", 200, 0, 200, 200, 0, 200);
  h_algo_MH_v_MZ = new TH2D("VbbHcc_algo_MH_v_MZ", "", 200, 0, 200, 200, 0, 200);
  h_both_MH_v_MZ = new TH2D("VbbHcc_both_MH_v_MZ", "", 200, 0, 200, 200, 0, 200);
  h_tags_MH_v_MZ_select = new TH2D("VbbHcc_tags_MH_v_MZ_select", "", 200, 0, 200, 200, 0, 200);
  h_algo_MH_v_MZ_select = new TH2D("VbbHcc_algo_MH_v_MZ_select", "", 200, 0, 200, 200, 0, 200);
  h_both_MH_v_MZ_select = new TH2D("VbbHcc_both_MH_v_MZ_select", "", 200, 0, 200, 200, 0, 200);
 
 
  //Add histograms to fOutput so they can be saved in Processor::SlaveTerminate
  r->GetOutputList()->Add(h_evt) ;
  std::vector<TH1*> tmp = h_VH->returnHisto() ;
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_tags->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_algo->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_both->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_bothTags->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_bothAlgo->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_eff_tags->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_eff_algo->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_eff_both->returnHisto();
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  r->GetOutputList()->Add(h_evt_cutflow);
  r->GetOutputList()->Add(h_evt_tags_cutflow);
  r->GetOutputList()->Add(h_evt_algo_cutflow);
  r->GetOutputList()->Add(h_evt_both_cutflow);
  r->GetOutputList()->Add(h_jet_cutflow);
  r->GetOutputList()->Add(h_elec_cutflow);
  r->GetOutputList()->Add(h_muon_cutflow);
  r->GetOutputList()->Add(h_CvL);
  r->GetOutputList()->Add(h_CSV);
  r->GetOutputList()->Add(h_Nselected);
  r->GetOutputList()->Add(h_Nbjet);
  r->GetOutputList()->Add(h_Ncjet);
  r->GetOutputList()->Add(h_Z_dM);
  r->GetOutputList()->Add(h_H_dM);
  r->GetOutputList()->Add(h_tags_MH_v_MZ);
  r->GetOutputList()->Add(h_algo_MH_v_MZ);
  r->GetOutputList()->Add(h_both_MH_v_MZ);
  r->GetOutputList()->Add(h_tags_MH_v_MZ_select);
  r->GetOutputList()->Add(h_algo_MH_v_MZ_select);
  r->GetOutputList()->Add(h_both_MH_v_MZ_select);

}//end-SlaveBegin

float D_HZ(float mH, float mZ) {
  float x0 = 125.0, y0 = 91.0;
  float k = x0 / y0;
  float m0 = std::max(mH, mZ);
  float m1 = std::min(mH, mZ);
  return fabs(m0 - k * m1) / sqrt(1 + pow(k,2));
}


float get_mass_from_indices(std::vector<JetObj>& jets, int idx0, int idx1) {
  TLorentzVector vec0 = jets.at(idx0).m_lvec;
  TLorentzVector vec1 = jets.at(idx1).m_lvec;
  TLorentzVector vec_total = vec0 + vec1;
  return vec_total.M(); 
}


#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
std::vector<std::vector<int> > VH_selection::DauIdxs_ZH(Reader* r) {
  // store the indices of the H and Z daughters
  std::vector<std::vector<int> > dauIdxs;
  std::vector<int> dauIdxsZ;
  std::vector<int> dauIdxsH;

  // loop over the gen parts
  for (unsigned i = 0; i < *(r->nGenPart); ++i) {
    // mother index
    int mIdx = (r->GenPart_genPartIdxMother)[i];
    // b quarks from Z (id = 5, has a mother, & mother id = 23)
    if (fabs((r->GenPart_pdgId)[i]) == 5 && mIdx > -1 && (r->GenPart_pdgId)[mIdx] == 23) dauIdxsZ.push_back(i);
    // c quarks from H (id = 4, has a mother, & mother id = 25)
    if (fabs((r->GenPart_pdgId)[i]) == 4 && mIdx > -1 && (r->GenPart_pdgId)[mIdx] == 25) dauIdxsH.push_back(i);
  }

  // push back & return the proper IDs
  dauIdxs.push_back(dauIdxsZ);
  dauIdxs.push_back(dauIdxsH);
  return dauIdxs;
}
#endif


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
    JetObj jet((r->Jet_pt)[i],(r->Jet_eta)[i],(r->Jet_phi)[i],(r->Jet_mass)[i], 
      jetFlav, (r->Jet_btagDeepFlavB)[i],(r->Jet_puIdDisc)[i]) ;

#if defined(NANOAODV9)
    jet.m_deepCvL = (r->Jet_btagDeepFlavCvL)[i];
#endif

#if defined(NANOAODV7)
    jet.m_deepCvL = (r->FatJet_btagDDCvL)[i];
#endif

    jets.push_back(jet) ;
  }

  // Electrons
  std::vector<LepObj> elecs;
  for (unsigned int i = 0; i < *(r->nElectron); ++i) {
    
    // Produce an electron from the information.
    h_elec_cutflow->Fill(0.5, genWeight); // All electrons

    float etaSC = (r->Electron_eta)[i] - (r->Electron_deltaEtaSC)[i];
    LepObj elec((r->Electron_pt)[i], (r->Electron_eta)[i], etaSC,
           (r->Electron_phi)[i], (r->Electron_mass)[i], i,
           (r->Electron_charge)[i], 0);

    // Cut #1 - check phase space /////////////////////
    // We want electrons with high enough pT and within
    // the eta limit of the tracker.
    if (elec.m_lvec.Pt() < CUTS.Get<float>("lep_pt1") || 
      fabs(elec.m_lvec.Eta()) > CUTS.Get<float>("lep_eta"))
       continue;

    h_elec_cutflow->Fill(1.5, genWeight); // passed phase cut

    // Cut #2 - make a cut based on SC /////////////////
    if (fabs(etaSC) < 1.566 && fabs(etaSC) > 1.442) continue;

    h_elec_cutflow->Fill(2.5, genWeight); // passed SC cut

    // Cut #3 - make sure to only keep proper electrons ////
    int elecID = r->Electron_cutBased[i];
    if (elecID < 2) continue; // loose electron ID
    h_elec_cutflow->Fill(3.5, genWeight); // passed ID cut

    elecs.push_back(elec);
  }

  // Muons
  std::vector<LepObj> muons;
  for (unsigned int i = 0; i < *(r->nMuon); ++i) {
  
    h_muon_cutflow->Fill(0.5, genWeight); // All muons
    LepObj muon((r->Muon_pt)[i], (r->Muon_eta)[i], -1, (r->Muon_phi)[i],
           (r->Muon_mass)[i], i, (r->Muon_charge)[i],
           (r->Muon_pfRelIso04_all)[i]);
 
    // Cut #1 - check phase space /////////////////////////
    if (muon.m_lvec.Pt() > CUTS.Get<float>("lep_pt1") ||
       fabs(muon.m_lvec.Eta()) > CUTS.Get<float>("lep_eta"))
       continue;

    h_muon_cutflow->Fill(1.5, genWeight); // passed phase cut
    
    // Cut #2 - check loose ID /////////////////////////////
    if (r->Muon_looseId[i] <= 0) continue;
    h_muon_cutflow->Fill(2.5, genWeight); // passed ID cut

    // Cut #3 - isolation cut //////////////////////////////
    if (muon.m_iso > CUTS.Get<float>("muon_iso")) continue;
    h_muon_cutflow->Fill(3.5, genWeight); // passed iso cut

    muons.push_back(muon);
  }
   

  // All cut flows need to show the total events //
  h_evt_cutflow->Fill(0.5, genWeight);
  h_evt_tags_cutflow->Fill(0.5, genWeight);
  h_evt_algo_cutflow->Fill(0.5, genWeight);
  h_evt_both_cutflow->Fill(0.5, genWeight);
 
  // So we have them for reference (in terms of making cuts),
  // make sure to place ALL jets into distributions.
  h_VH->FillJets(jets, evtW);
  h_VH->FillNjet(jets.size(), evtW);

  // Cut #1 - Proper Jets /////////////////////////
  // We want to keep jets that meet the following criteria:
  // pT(j) > 30 GeV, |eta| < 2.5
  // dR(small-R jet, lepton) < 0.4 = discard
  std::vector<JetObj> selected_jets;
  for (unsigned int i = 0; i < jets.size(); ++i) {
  
    TLorentzVector vec = jets.at(i).m_lvec;
    h_jet_cutflow->Fill(0.5, genWeight); // All jets
    if (vec.Pt() < 30) continue; 
    h_jet_cutflow->Fill(1.5, genWeight); // Passed pT cut
    if (fabs(vec.Eta()) > 2.5) continue;
    h_jet_cutflow->Fill(2.5, genWeight); // Passed eta cut
   
    // Check the electrons & muons for overlap with the jets.
    // If dR(jet, lepton) < 0.4, discard the jet.
    bool should_discard = false;
    for (unsigned int j = 0; j < elecs.size(); ++j) {
      float dR = vec.DeltaR(elecs.at(j).m_lvec);
      if (dR < 0.4) { should_discard = true; break; }
    };

    if (should_discard) continue;

    should_discard = false;
    for (unsigned int j = 0; j < muons.size(); ++j) {
      float dR = vec.DeltaR(muons.at(j).m_lvec);
      if (dR < 0.4) { should_discard = true; break; }
    };

    if (should_discard) continue;
    
    h_jet_cutflow->Fill(3.5, genWeight); // Passed iso

    selected_jets.push_back(jets.at(i));
  } 

  // We want to be able to plot the distributions of the
  // jets that survive our selections.
  h_Nselected->Fill(selected_jets.size(), genWeight);

// CASE #1 - MC TRUTH /////////////////////////////////////////////////////
// Remember, we only want to do this if we have a proper MC file.
// We ignore this part for data files.
bool canCompareToMC = false;
std::vector<std::vector<int> > dauIdxs;
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018) 

    // Make sure we have two daughters for each particle
    int Zid = 0, Hid = 1;
    dauIdxs = DauIdxs_ZH(r);
    std::vector<JetObj> genObjs;
    if (dauIdxs[Zid].size() == 2 && dauIdxs[Hid].size() == 2) {

      canCompareToMC = true;
      h_evt_cutflow->Fill(1.5, genWeight);
      int idx1_Z = dauIdxs[Zid][0];
      int idx2_Z = dauIdxs[Zid][1];
      int idx1_H = dauIdxs[Hid][0];
      int idx2_H = dauIdxs[Hid][1];
      
      // Reconstruct from the GenPart data.
      JetObj b0((r->GenPart_pt)[idx1_Z], (r->GenPart_eta)[idx1_Z],
         (r->GenPart_phi)[idx1_Z], (r->GenPart_mass)[idx1_Z], 5, 0., 0.);
      JetObj b1((r->GenPart_pt)[idx2_Z], (r->GenPart_eta)[idx2_Z],
         (r->GenPart_phi)[idx2_Z], (r->GenPart_mass)[idx2_Z], 5, 0., 0.);
      std::vector<JetObj> bjets{b0, b1};
      ZObj Z(bjets);
      genObjs.push_back(b0); genObjs.push_back(b1);

      JetObj c0((r->GenPart_pt)[idx1_H], (r->GenPart_eta)[idx1_H],
         (r->GenPart_phi)[idx1_H], (r->GenPart_mass)[idx1_H], 4, 0., 0.);
      JetObj c1((r->GenPart_pt)[idx2_H], (r->GenPart_eta)[idx2_H],
         (r->GenPart_phi)[idx2_H], (r->GenPart_mass)[idx2_H], 4, 0., 0.);
      std::vector<JetObj> cjets{c0, c1};
      HObj H(cjets);      
      genObjs.push_back(c0); genObjs.push_back(c1);

      // Fill the objects
      h_VH->Fill(H, Z, evtW);
    }

#endif 

  // Remember, we want at least 4 jets.
  if (selected_jets.size() >= 4) {

    h_evt_tags_cutflow->Fill(1.5, genWeight); // passed jet selection
    h_evt_algo_cutflow->Fill(1.5, genWeight);
    h_evt_both_cutflow->Fill(1.5, genWeight);   

    // CASE #2 - TAGGING //////////////////////////////////////////////////////
    
    // Select two jets with the largest btag value and then 
    // make pass ~medium WP (Jet_btagDeepFlavB > 0.3)
    std::vector<JetObj> jets2 = selected_jets;
    std::vector<JetObj> bjets;
    float csv0 = -2000.0, csv1 = -2000.0; int bIdx0 = -1, bIdx1 = -1;

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
    jets2.erase(jets2.begin() + bIdx0);

    if (csv0 > 0.3 && csv1 > 0.3)
    {
       h_evt_tags_cutflow->Fill(2.5, genWeight); // pass b-cuts
       ZObj Z(bjets);
       std::vector<JetObj> cjets;

       // Do the same matching process but with c-jets
       float cvl0 = -2000.0, cvl1 = -2000.0; int cIdx0 = -1, cIdx1 = -1;
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
       
       if (cvl0 > 0.37 && cvl1 > 0.37) {

         h_evt_tags_cutflow->Fill(3.5, genWeight); // pass c-cuts
         HObj H(cjets);
         if (*(r->MET_pt) < 140) {

           h_evt_tags_cutflow->Fill(4.5, genWeight); // pass MET cut
           if (Z.m_lvec.Pt() > 50) {

             h_evt_tags_cutflow->Fill(5.5, genWeight); // pass pT(Z) cut
             float dPhi = Z.m_lvec.DeltaPhi(H.m_lvec);
             if (fabs(dPhi) > 2.5) {

               h_evt_tags_cutflow->Fill(6.5, genWeight); // pass dPhi cut
               h_VH_tags->Fill(H, Z, evtW);
               h_tags_MH_v_MZ_select->Fill(H.m_lvec.M(), Z.m_lvec.M(), evtW);

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
               
             if (canCompareToMC) {
               // Get the gen-level objects. 
               int idx1_Z = dauIdxs[0][0], idx2_Z = dauIdxs[0][1];
               int idx1_H = dauIdxs[1][0], idx2_H = dauIdxs[1][1];
               JetObj b0((r->GenPart_pt)[idx1_Z], (r->GenPart_eta)[idx1_Z],
                 (r->GenPart_phi)[idx1_Z], (r->GenPart_mass)[idx1_Z], 5, 0., 0.);
               JetObj b1((r->GenPart_pt)[idx2_Z], (r->GenPart_eta)[idx2_Z],
                 (r->GenPart_phi)[idx2_Z], (r->GenPart_mass)[idx2_Z], 5, 0., 0.);
               std::vector<JetObj> bjets{b0, b1};

               JetObj c0((r->GenPart_pt)[idx1_H], (r->GenPart_eta)[idx1_H],
                 (r->GenPart_phi)[idx1_H], (r->GenPart_mass)[idx1_H], 4, 0., 0.);
               JetObj c1((r->GenPart_pt)[idx2_H], (r->GenPart_eta)[idx2_H],
                 (r->GenPart_phi)[idx2_H], (r->GenPart_mass)[idx2_H], 4, 0., 0.);
               std::vector<JetObj> cjets{c0, c1};

               h_eff_tags->Fill(H, Z, cjets, bjets, evtW);            
             }
#endif

             }//end-dPhi-cut
           }//end-pt(V)-cut
         }//end-met-cut          

       }//end-cjet-selection
    }//end-bjet-selection

    // CASE #3 - D_HZ ALGORITHM ///////////////////////////////////////////////

    // Make sure we have our jets properly sorted by pT(j)
    std::vector<JetObj> jets3 = selected_jets;
    std::sort(jets3.begin(), jets3.end(), std::greater<JetObj>());    
  
    // Create the objects for the distance calculations & make sure we sort
    // them in ascending order. NOTE: The algorithm for calculating the distance
    // is contained within the DObj class and is done immediately upon creation.
    // Additionally, the operators > and < within the class are designed to sort
    // the objects based on the distance calculated.
    DObj d0(jets3, 0, 1, 2, 3); // H(0,1) ; Z(2,3)
    DObj d1(jets3, 0, 2, 1, 3); // H(0,2) ; Z(1,3)
    DObj d2(jets3, 0, 3, 1, 2); // H(0,3) ; Z(1,2)
    std::vector<DObj> distances{ d0, d1, d2 };
    std::sort(distances.begin(), distances.end());

    //h_HZ0->Fill(d0.m_d, evtW); h_HZ1->Fill(d1.m_d, evtW); h_HZ2->Fill(d2.m_d, evtW);
    //h_VH_algo->FillAlgo(d0.m_d, d1.m_d, d2.m_d, evtW);
    h_VH_algo->FillAlgo(distances[0].m_d, distances[1].m_d, distances[2].m_d, evtW);
    h_algo_MH_v_MZ->Fill(d0.HM(), d0.ZM(), evtW);
    h_algo_MH_v_MZ->Fill(d1.HM(), d1.ZM(), evtW);
    h_algo_MH_v_MZ->Fill(d2.HM(), d2.ZM(), evtW); 

    // Determine the difference between the closest two.
    float deltaD = fabs(distances[0].m_d - distances[1].m_d);
    //h_dH->Fill(deltaD, evtW);    

    // Determine which pairing we want to use based on this deltaD.
    // If we are outside the resolution window (30 GeV), we can 
    // just choose the closest pair.
    DObj chosenPair = distances[0];
    if (deltaD >= 30) { chosenPair = distances[0];}
    // Otherwise, we want to choose the pairing of the lowest two
    // that has the largest pT(H).
    else {
      float pt0 = distances[0].HPt();
      float pt1 = distances[1].HPt();
      int idx;
      if (pt0 > pt1) idx = 0;
      else idx = 1;
      chosenPair = distances[idx];
    }

    // Now that we've passed the algorithm for matching, let's check the 
    // tagging for the b- and c-jets. The tagging will be taken care of
    // inside the DObj class.
    if (chosenPair.Z_has_bjets()) {

      h_evt_algo_cutflow->Fill(2.5, genWeight); // pass b-tagging
      if (chosenPair.H_has_cjets()) {

        h_evt_algo_cutflow->Fill(3.5, genWeight); // pass c-tagging
        if (*(r->MET_pt) < 140) {

          h_evt_algo_cutflow->Fill(4.5, genWeight); // pass MET cut
          if (chosenPair.ZPt() > 50) {

            h_evt_algo_cutflow->Fill(5.5, genWeight); // pass pT(Z) cut
            float dPhi = fabs(chosenPair.DPhi());
            if (dPhi > 2.5) {

              h_evt_algo_cutflow->Fill(6.5, genWeight); // pass dPhi cut
              // We need to build the objects to fill our histograms.
              std::vector<JetObj> bjets;
              bjets.push_back(jets3[chosenPair.m_zIdx0]);
              bjets.push_back(jets3[chosenPair.m_zIdx1]);
              ZObj Z(bjets);
              
              std::vector<JetObj> cjets;
              cjets.push_back(jets3[chosenPair.m_hIdx0]);
              cjets.push_back(jets3[chosenPair.m_hIdx1]);
              HObj H(cjets);

              h_VH_algo->Fill(H, Z, evtW);
              h_algo_MH_v_MZ_select->Fill(H.m_lvec.M(), Z.m_lvec.M(), evtW);

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)

            if (canCompareToMC) {
              // Get the gen-level objects. 
              int idx1_Z = dauIdxs[0][0], idx2_Z = dauIdxs[0][1];
              int idx1_H = dauIdxs[1][0], idx2_H = dauIdxs[1][1];
              JetObj b0((r->GenPart_pt)[idx1_Z], (r->GenPart_eta)[idx1_Z],
                (r->GenPart_phi)[idx1_Z], (r->GenPart_mass)[idx1_Z], 5, 0., 0.);
              JetObj b1((r->GenPart_pt)[idx2_Z], (r->GenPart_eta)[idx2_Z],
                (r->GenPart_phi)[idx2_Z], (r->GenPart_mass)[idx2_Z], 5, 0., 0.);
              std::vector<JetObj> bjets{b0, b1};

              JetObj c0((r->GenPart_pt)[idx1_H], (r->GenPart_eta)[idx1_H],
                (r->GenPart_phi)[idx1_H], (r->GenPart_mass)[idx1_H], 4, 0., 0.);
              JetObj c1((r->GenPart_pt)[idx2_H], (r->GenPart_eta)[idx2_H],
                (r->GenPart_phi)[idx2_H], (r->GenPart_mass)[idx2_H], 4, 0., 0.);
              std::vector<JetObj> cjets{c0, c1};

              h_eff_algo->Fill(H, Z, cjets, bjets, evtW);
             }
#endif

         
            }//end-dPhi-cut 
          }//end-pt(V)-cut
        }//end-MET-cut        

      }//end-b-tag

    }//end-c-tag

    // CASE #4 - TAGGING & ALGORITHM //////////////////////////////////////////

    // Make sure we have our jets properly sorted by pT(j)
    std::vector<JetObj> jets4 = selected_jets;
    std::sort(jets4.begin(), jets4.end(), std::greater<JetObj>());
    
    // Check to see which jets pass which tagging.
    std::vector<int> bIndices; std::vector<int> cIndices;
    bool ctags[4] = { false, false, false, false };
    bool btags[4] = { false, false, false, false };
    for (int i = 0; i < 4; ++i) {
      if (jets4[i].m_deepCSV > 0.3){ btags[i] = true; bIndices.push_back(i); }
      if (jets4[i].m_deepCvL > 0.37){ ctags[i] = true; cIndices.push_back(i); }
    }

    // Make sure we have no jets that pass both tags.
    bool properly_chosen = true;
    for (int i = 0; i < 4; ++i) {
      if (ctags[i] && btags[i]) {
        properly_chosen = false; break;
      }  
    }

    // Make sure we have 2 b-jets and 2 c-jets.
    if (properly_chosen && bIndices.size() == 2 && cIndices.size() == 2) {
      
      h_evt_both_cutflow->Fill(2.5, genWeight); //pass b-tagging
      h_evt_both_cutflow->Fill(3.5, genWeight); //pass c-tagging

      // Get the jets into proper lists. Reconstruct the bosons.
      std::vector<JetObj> bjets;
      bjets.push_back(jets4[bIndices[0]]);
      bjets.push_back(jets4[bIndices[1]]);
      ZObj Z(bjets);
 
      std::vector<JetObj> cjets;
      cjets.push_back(jets4[cIndices[0]]);
      cjets.push_back(jets4[cIndices[1]]);
      HObj H(cjets);

      // Go through the rest of our cuts.
      if (*(r->MET_pt) < 140) {
        h_evt_both_cutflow->Fill(4.5, genWeight); // pass MET cut
        if (Z.m_lvec.Pt() > 50) {

          h_evt_both_cutflow->Fill(5.5, genWeight); // pass pT(Z) cut
          float dPhi = fabs(Z.m_lvec.DeltaPhi(H.m_lvec));
 
          if (dPhi > 2.5) {
            h_evt_both_cutflow->Fill(6.5, genWeight); // pass dPhi cut
            h_VH_both->Fill(H, Z, evtW);
            h_both_MH_v_MZ_select->Fill(H.m_lvec.M(), Z.m_lvec.M(), evtW);
            h_VH_bothTags->Fill(H, Z, evtW);

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)

            if (canCompareToMC) {
               // Get the gen-level objects. 
               int idx1_Z = dauIdxs[0][0], idx2_Z = dauIdxs[0][1];
               int idx1_H = dauIdxs[1][0], idx2_H = dauIdxs[1][1];
               JetObj b0((r->GenPart_pt)[idx1_Z], (r->GenPart_eta)[idx1_Z],
                 (r->GenPart_phi)[idx1_Z], (r->GenPart_mass)[idx1_Z], 5, 0., 0.);
               JetObj b1((r->GenPart_pt)[idx2_Z], (r->GenPart_eta)[idx2_Z],
                 (r->GenPart_phi)[idx2_Z], (r->GenPart_mass)[idx2_Z], 5, 0., 0.);
               std::vector<JetObj> bjets{b0, b1};

               JetObj c0((r->GenPart_pt)[idx1_H], (r->GenPart_eta)[idx1_H],
                 (r->GenPart_phi)[idx1_H], (r->GenPart_mass)[idx1_H], 4, 0., 0.);
               JetObj c1((r->GenPart_pt)[idx2_H], (r->GenPart_eta)[idx2_H],
                 (r->GenPart_phi)[idx2_H], (r->GenPart_mass)[idx2_H], 4, 0., 0.);
               std::vector<JetObj> cjets{c0, c1};

               h_eff_both->Fill(H, Z, cjets, bjets, evtW);
             }
#endif
          }//end-dPhi-cut
        }//end-pt(V)-cut
      }//end-MET-cut   
    }//end-tagging

    // If we don't, we have to go through the mass matching algorithm.
    else {

      // Run the mass matching algorithm.
      DObj d0(jets4, 0, 1, 2, 3); // H(0,1) ; Z(2,3)
      DObj d1(jets4, 0, 2, 1, 3); // H(0,2) ; Z(1,3)
      DObj d2(jets4, 0, 3, 1, 2); // H(0,3) ; Z(1,2)
      std::vector<DObj> distances{ d0, d1, d2 };
      std::sort(distances.begin(), distances.end());    

      //h_HZ0->Fill(d0.m_d, evtW); h_HZ1->Fill(d1.m_d, evtW); h_HZ2->Fill(d2.m_d, evtW);
      //h_VH_both->FillAlgo(d0.m_d, d1.m_d, d2.m_d, evtW);
      h_VH_both->FillAlgo(distances[0].m_d, distances[1].m_d, distances[2].m_d, evtW);
      h_both_MH_v_MZ->Fill(d0.HM(), d0.ZM(), evtW);
      h_both_MH_v_MZ->Fill(d1.HM(), d1.ZM(), evtW);
      h_both_MH_v_MZ->Fill(d2.HM(), d2.ZM(), evtW);

      // Determine the difference between the closest two.
      float deltaD = fabs(distances[0].m_d - distances[1].m_d);
      //h_dH->Fill(deltaD, evtW);

      // Determine which pairing we want to use based on this deltaD.
      // If we are outside the resolution window (30 GeV), we can
      // just choose the closest pair.
      DObj chosenPair = distances[0];
      if (deltaD >= 30) { chosenPair = distances[0]; }
      // Otherwise, we want to choose the pairing of the lowest two
      // that has the largest pT(H).   
      else {
        float pt0 = distances[0].HPt();
        float pt1 = distances[0].HPt();
        int idx;
        if (pt0 > pt1) idx = 0;
        else idx = 1;
        chosenPair = distances[idx];
      }

      // Now that we've passed the algorithm, let's check the tagging
      // for b- and c-jets.
      if (chosenPair.Z_has_bjets()) {

        h_evt_both_cutflow->Fill(2.5, genWeight); // pass b-tagging
        if (chosenPair.H_has_cjets()) {

          h_evt_both_cutflow->Fill(3.5, genWeight); //  pass c-tagging
          if (*(r->MET_pt) < 140) {

            h_evt_both_cutflow->Fill(4.5, genWeight); // pass MET cut
            if (chosenPair.ZPt() > 50) {

              h_evt_both_cutflow->Fill(5.5, genWeight); // pass pT(Z) cut
              float dPhi = fabs(chosenPair.DPhi());
              if (dPhi > 2.5) {
                h_evt_both_cutflow->Fill(6.5, genWeight); // pass dPhi cut

                // We need to build the objects to fill our histograms.
                std::vector<JetObj> bjets;
	        bjets.push_back(jets4[chosenPair.m_zIdx0]);
	        bjets.push_back(jets4[chosenPair.m_zIdx1]);
	        ZObj Z(bjets);
                
                std::vector<JetObj> cjets;
	        cjets.push_back(jets4[chosenPair.m_hIdx0]);
	        cjets.push_back(jets4[chosenPair.m_hIdx1]);
	        HObj H(cjets);

	        h_VH_both->Fill(H, Z, evtW);
                h_VH_bothAlgo->Fill(H, Z, evtW);
                h_both_MH_v_MZ_select->Fill(H.m_lvec.M(), Z.m_lvec.M(), evtW);

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)

            if (canCompareToMC) {
               // Get the gen-level objects. 
               int idx1_Z = dauIdxs[0][0], idx2_Z = dauIdxs[0][1];
               int idx1_H = dauIdxs[1][0], idx2_H = dauIdxs[1][1];
               JetObj b0((r->GenPart_pt)[idx1_Z], (r->GenPart_eta)[idx1_Z],
                 (r->GenPart_phi)[idx1_Z], (r->GenPart_mass)[idx1_Z], 5, 0., 0.);
               JetObj b1((r->GenPart_pt)[idx2_Z], (r->GenPart_eta)[idx2_Z],
                 (r->GenPart_phi)[idx2_Z], (r->GenPart_mass)[idx2_Z], 5, 0., 0.);
               std::vector<JetObj> bjets{b0, b1};

               JetObj c0((r->GenPart_pt)[idx1_H], (r->GenPart_eta)[idx1_H],
                 (r->GenPart_phi)[idx1_H], (r->GenPart_mass)[idx1_H], 4, 0., 0.);
               JetObj c1((r->GenPart_pt)[idx2_H], (r->GenPart_eta)[idx2_H],
                 (r->GenPart_phi)[idx2_H], (r->GenPart_mass)[idx2_H], 4, 0., 0.);
               std::vector<JetObj> cjets{c0, c1};

               h_eff_both->Fill(H, Z, cjets, bjets, evtW);
             }
#endif

              }//end-dPhi-cut
            }//end-pT(V)-cut
          }//end-MET-cut
        }//end-c-tag
      }//end-b-tag
    }//end-mass-match


  }//end-jet-selection
 
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
