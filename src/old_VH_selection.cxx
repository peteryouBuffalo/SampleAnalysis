#define VH_selection_cxx

//== Include Statements =======================================================

#include <math.h>
#include "TList.h"
#include "TParameter.h"

#include "VH_selection.h"
#include "Global.h"
#include "Obj.cxx"

#include <bits/stdc++.h>

// ============================================================================
//  DECONSTRUCTOR
// ============================================================================

VH_selection::~VH_selection() { }

// ============================================================================
// CUSTOM METHODS
// ============================================================================


// == DauIdxs_ZH - get the indices of the MC truth particles == 
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
std::vector<std::vector<int> > VH_selection::DauIdxs_ZH(Reader *r) {

  // Store the indices of the Higgs and Z daughters
  std::vector<std::vector<int>> dauIdxs;
  std::vector<int> dauIdxsZ;
  std::vector<int> dauIdxsH;

  // Loop over the gen particles
  for (unsigned i = 0; i < *(r->nGenPart); ++i) {

    int mIdx = (r->GenPart_genPartIdxMother)[i]; // mother index     
    int mID  = (r->GenPart_pdgId)[mIdx];         // mother ID (PDG)
    int flav = (r->GenPart_pdgId)[i];            // ptcl ID (PDG)

    // b quark check (id = 5, mother ID = 23 [Z])
    if (abs(flav) == 5 && mIdx > -1 && abs(mID) == 23)
      dauIdxsZ.push_back(i);

    // c quark check (id = 4, mother ID = 25 [H])
    if (abs(flav) == 4 && mIdx > -1 && abs(mID) == 25)
      dauIdxsH.push_back(i);

    //if (mIdx > -1 && abs(mID) == 23)
    //  std::cout << "from Z, we get flav = " << flav << std::endl;
  }

  // Push back and return the proper IDs
  dauIdxs.push_back(dauIdxsZ);
  dauIdxs.push_back(dauIdxsH);

  return dauIdxs;
}
#endif


// == sort_by_second - Sort a pair by the second option == 
bool sort_by_second(const std::pair<int,float> &a,
                    const std::pair<int,float> &b)
{ return a.second < b.second; }

bool sort_by_second_descend(const std::pair<int,float> &a,
                            const std::pair<int,float> &b)
{ return a.second > b.second; } 


// == list_contains - This tells us if a list contains a given value ==
bool list_contains(std::vector<int> list, int value) {
  return std::count(list.begin(), list.end(), value);
}//end-list_contains


// == get_valid_ids - Check through a list of IDs and check against a list of ==
// == already used IDs to get the valid ones                                  ==
std::vector<int> get_valid_ids(std::vector<int> idxs, std::vector<int> invalid_idxs) {

  // Let's check for the top two valid idxs
  std::vector<int> valid_idxs;  
  for (size_t i = 0; i < idxs.size(); ++i) {

    // Check to see if idx #i is not already used.
    // If it's not used, let's add it to our list.
    if (!list_contains(invalid_idxs, idxs[i])) valid_idxs.push_back(idxs[i]);

    // If we've already found two valid idxs, let's break out
    if (valid_idxs.size() >= 2) break;

  }//end-for

  // If we have not found a valid index, return a fail-safe.
  return valid_idxs;

}


// == find_valid_combos - Check for all valid combinations of the given jets ==
std::vector<std::vector<int> > find_valid_combos(std::vector<int> bIdxs,
  std::vector<int> cIdxs) {

  std::vector<std::vector<int>> outs;

  // Check to make sure we even have enough points to check.
  if (bIdxs.size() >= 2 && cIdxs.size() >= 2) {

    // Option #1 - choose the top two bIndices
    int b0(bIdxs[0]), b1(bIdxs[1]);
    std::vector<int> bs{b0, b1};
    
    // Get the top two c indices that are not in the bList.
    std::vector<int> clist = get_valid_ids(cIdxs, bs);

    // If we have a proper number of cIDs and none of them are -1,
    // create the proper combination.
    if (clist.size() >= 2 && !list_contains(clist,-1)) {
      int c2(clist[0]); int c3(clist[1]);
      outs.push_back({b0, b1, c2, c3});
    }

    // Option #2 - choose the top two cIndices
    int c0(cIdxs[0]), c1(cIdxs[1]);
    std::vector<int> cs{c0, c1};   
 
    // Only check this way if there's overlap with the indices.
    if (list_contains(cs, b0) || list_contains(cs, b1)) { 

      // Get the top two b indices that are not in the cList.
      std::vector<int> blist = get_valid_ids(bIdxs, cs);
    
      // If we have a proper number of bIDs and none of them are -1,
      // create the proper combination.
      if (blist.size() >= 2 && !list_contains(blist,-1)) {
        int b2(blist[0]); int b3(blist[1]);
        outs.push_back({b2, b3, c0, c1});
      }
    }
  }//end-check

  // Return our list.
  return outs;

}//end-find_valid_combos


// == passes_btag - checks that a jet is b-tagged ==
bool passes_btag(JetObj& jet, float CSV_cut) {
  return jet.m_deepCSV > CSV_cut;
}


// == are_bjets - takes a list of jets and checks that they're all b-jets ==
// == (given a certain cut)                                               ==
bool are_bjets(std::vector<JetObj>& jets, float CSV_cut) {
  for (auto it : jets) {
    if (!passes_btag(it, CSV_cut)) return false;
  }
  return true;
}

// == passes_ctag - checks that a jet is c-tagged == 
bool passes_ctag(JetObj& jet, float CvL_cut, float CvB_cut) {
  bool passes_CvL = jet.m_deepCvL > CvL_cut;
  bool passes_CvB = jet.m_deepCvB > CvB_cut;
  return passes_CvL && passes_CvB;
}

// == are_cjets - takes a list of jets and checks that they're all c-jets ==
// == (given a certain cut)                                               == 
bool are_cjets(std::vector<JetObj>& jets, float CvL_cut, float CvB_cut) {
  for (auto it : jets) {
    if (!passes_ctag(it, CvL_cut, CvB_cut)) return false;
  }
  return true;
}

// == determine_trigger_SF - take in the list of jets, determine the HT of the ==
// == event, and then calculate return the appropriate value.                  ==
float determine_trigger_SF(std::vector<JetObj>& jets, int year) {
  float HT = 0.0;
  for (auto it : jets) HT += it.Pt();

  if (HT < 250 || HT >= 2000) return 1.0;
  else if (HT < 300) {
    if (year == 16) return 0.507;
    if (year == 17) return 0.488;
    if (year == 18) return 1.784;
  }
  else if (HT < 350) {
    if (year == 16) return 0.607;
    if (year == 17) return 0.563;
    if (year == 18) return 1.398;
  }
  else if (HT < 400) {
    if (year == 16) return 0.660;
    if (year == 17) return 0.716;
    if (year == 18) return 1.077;
  }
  else if (HT < 450) {
    if (year == 16) return 0.677;
    if (year == 17) return 0.701;
    if (year == 18) return 0.980;
  }
  else if (HT < 500) {
    if (year == 16) return 0.711;
    if (year == 17) return 0.827;
    if (year == 18) return 0.914;
  }
  else if (HT < 550) {
    if (year == 16) return 0.725;
    if (year == 17) return 0.845;
    if (year == 18) return 0.975;
  }
  else if (HT < 600) {
    if (year == 16) return 0.738;
    if (year == 17) return 0.802;
    if (year == 18) return 0.908;
  }
  else if (HT < 800) {
    if (year == 16) return 0.762;
    if (year == 17) return 0.867;
    if (year == 18) return 0.905;
  }
  else if (HT < 1000) {
    if (year == 16) return 0.765;
    if (year == 17) return 0.955;
    if (year == 18) return 0.914;
  }
  else if (HT < 2000) {
    if (year == 16) return 0.809;
    if (year == 17) return 0.942;
    if (year == 18) return 0.918;
  }

  return 1.0;
}

// ============================================================================
// SLAVE BEGIN
// ============================================================================

/*
 * This section is where we initialize all necessary plots
 * and then add all the plots to the output so we can use
 * them for analyses.
*/

void VH_selection::SlaveBegin(Reader *r) {

  // Set up all the necessary plots
  h_evt = new TH1D("Nevt", "", 4, -1.5, 2.5); // bin 1 = total negative w evt
                                              // bin 2 = total positive w evt
                                              // bin 3 = total event weight
                                              // genWeight (=bin2 - bin1)

  h_evt_VbbHcc = new TH1D("Nevt_VbbHcc", "", 4, 0, 4);
  h_evt_VbbHcc->GetXaxis()->SetBinLabel(1, "total");
  h_evt_VbbHcc->GetXaxis()->SetBinLabel(2, "VbbHcc MC");
  h_evt_VbbHcc->GetXaxis()->SetBinLabel(3, "found MC jets");

  // Set up the VHPlot instances
  h_VH_all = new VHPlots("VbbHcc_all");
  h_VH_select = new VHPlots("VbbHcc_select");
  h_VH_MC = new VHPlots("VbbHcc_MC");
  h_VH_MCjet_minDR = new VHPlots("VbbHcc_MCjet_minDR");
  h_VH_MCjet_minDR_noTag = new VHPlots("VbbHcc_MCjet_minDR_noTag");
  h_VH_MCjet_dRcollect = new VHPlots("VbbHcc_MCjet_dRcollect");
  h_VH_MCjet_dRcollect_noTag = new VHPlots("VbbHcc_MCjet_dRcollect_noTag");
  h_VH_MCjet_ideal = new VHPlots("VbbHcc_MCjet_ideal");
  h_VH_MCjet_DHZ = new VHPlots("VbbHcc_MCjet_DHZ");

  h_VH_tagOnly = new VHPlots("VbbHcc_tagOnly");
  h_VH_tagOnly_noMassCorr = new VHPlots("VbbHcc_tagOnly_noMassCorr");
  h_VH_tagOnly_noJEC = new VHPlots("VbbHcc_tagOnly_noJEC");
  h_VH_tagOnly_2b1c = new VHPlots("VbbHcc_tagOnly_2b1c");

  h_VH_algoFirst = new VHPlots("VbbHcc_algoFirst");
  h_VH_algoFirst_noMassCorr = new VHPlots("VbbHcc_algoFirst_noMassCorr");
  h_VH_algoFirst_noJEC = new VHPlots("VbbHcc_algoFirst_noJEC");
  h_VH_algoFirst_2b1c = new VHPlots("VbbHcc_algoFirst_2b1c");

  h_VH_tagFirst = new VHPlots("VbbHcc_tagFirst");
  h_VH_tagFirst_noMassCorr = new VHPlots("VbbHcc_tagFirst_noMassCorr");
  h_VH_tagFirst_noJEC = new VHPlots("VbbHcc_tagFirst_noJEC");
  h_VH_tagFirst_2b1c = new VHPlots("VbbHcc_tagFirst_2b1c");

  h_VH_alljet = new VHPlots("VbbHcc_alljet");
  h_VH_seljet = new VHPlots("VbbHcc_seljet");

  // Set up the EffPlot instances
  h_eff_tagOnly = new EffPlots("VbbHcc_eff_tagOnly");
  h_eff_algoFirst = new EffPlots("VbbHcc_eff_algoFirst");
  h_eff_tagFirst = new EffPlots("VbbHcc_eff_tagFirst");

  // Set up the JetPlots instances
  h_jets_selected = new JetPlots("VbbHcc_jets");
  h_jets_all = new JetPlots("VbbHcc_jets_all");
  h_jets_MC = new JetPlots("VbbHcc_MC_jets");

  h_jets_4b = new JetPlots("VbbHcc_4b");
  h_jets_2b2c = new JetPlots("VbbHcc_2b2c");

  // Set up the GenPlots instances
  h_genJet_all = new GenPlots("GenJet_all");
  h_genJet_cuts = new GenPlots("GenJet_cuts");
  h_genJet_VbbHcc = new GenPlots("GenJet_VbbHcc");
  
  // Set up the RecoPlots instances
  h_reco_minDR = new RecoPlots("Reco_minDR");
  h_reco_minDR_under30 = new RecoPlots("Reco_minDR_under30");
  h_reco_minDR_noTag = new RecoPlots("Reco_minDR_noTag");
  h_reco_minDR_noTag_under30 = new RecoPlots("Reco_minDR_noTag_under30");
  h_reco_dRcollect = new RecoPlots("Reco_dRcollect");
  h_reco_dRcollect_under30 = new RecoPlots("Reco_dRcollect_under30");
  h_reco_dRcollect_noTag = new RecoPlots("Reco_dRcollect_noTag");
  h_reco_dRcollect_noTag_under30 = new RecoPlots("Reco_dRcollect_noTag_under30");
  h_reco_ideal = new RecoPlots("Reco_ideal");
  h_reco_ideal_under30 = new RecoPlots("Reco_ideal_under30");
  h_reco_DHZ = new RecoPlots("Reco_DHZ");
  h_reco_DHZ_under30 = new RecoPlots("Reco_DHZ_under30");

  // Set up the Trigger efficiency plots
  h_trig_2016_QuadJet_TripleTag = new TriggerEffPlots("2016_QuadJet_TripleTag");
  h_trig_2016_QuadJet_DoubleTag = new TriggerEffPlots("2016_QuadJet_DoubleTag");
  h_trig_2016_DoubleJet_TripleTag = new TriggerEffPlots("2016_DoubleJet_TripleTag");
  h_trig_2016_DoubleJet_DoubleTag = new TriggerEffPlots("2016_DoubleJet_DoubleTag");
  h_trig_2017_QuadJet_TripleTag = new TriggerEffPlots("2017_QuadJet_TripleTag");
  h_trig_2017_QuadJet_noTag = new TriggerEffPlots("2017_QuadJet_noTag");
  h_trig_2018_QuadJet_TripleTag = new TriggerEffPlots("2018_QuadJet_TripleTag");
  h_trig_2018_QuadJet_noTag = new TriggerEffPlots("2018_QuadJet_noTag"); 
 
  h_trig_2016_QuadJet_TripleTag_tagged = new TriggerEffPlots("2016_QuadJet_TripleTag_tagged");
  h_trig_2016_QuadJet_DoubleTag_tagged = new TriggerEffPlots("2016_QuadJet_DoubleTag_tagged");
  h_trig_2016_DoubleJet_TripleTag_tagged = new TriggerEffPlots("2016_DoubleJet_TripleTag_tagged");
  h_trig_2016_DoubleJet_DoubleTag_tagged = new TriggerEffPlots("2016_DoubleJet_DoubleTag_tagged");
  h_trig_2017_QuadJet_TripleTag_tagged = new TriggerEffPlots("2017_QuadJet_TripleTag_tagged");
  h_trig_2017_QuadJet_noTag_tagged = new TriggerEffPlots("2017_QuadJet_noTag_tagged");
  h_trig_2018_QuadJet_TripleTag_tagged = new TriggerEffPlots("2018_QuadJet_TripleTag_tagged");
  h_trig_2018_QuadJet_noTag_tagged = new TriggerEffPlots("2018_QuadJet_noTag_tagged");

  h_trig_2016_QuadJet_TripleTag_ideal = new TriggerEffPlots("2016_QuadJet_TripleTag_ideal");
  h_trig_2016_QuadJet_DoubleTag_ideal = new TriggerEffPlots("2016_QuadJet_DoubleTag_ideal");
  h_trig_2016_DoubleJet_TripleTag_ideal = new TriggerEffPlots("2016_DoubleJet_TripleTag_ideal");
  h_trig_2016_DoubleJet_DoubleTag_ideal = new TriggerEffPlots("2016_DoubleJet_DoubleTag_ideal");
  h_trig_2017_QuadJet_TripleTag_ideal = new TriggerEffPlots("2017_QuadJet_TripleTag_ideal");
  h_trig_2017_QuadJet_noTag_ideal = new TriggerEffPlots("2017_QuadJet_noTag_ideal");
  h_trig_2018_QuadJet_TripleTag_ideal = new TriggerEffPlots("2018_QuadJet_TripleTag_ideal");
  h_trig_2018_QuadJet_noTag_ideal = new TriggerEffPlots("2018_QuadJet_noTag_ideal");

  h_trig_2016_QuadJet_TripleTag_3B = new TriggerEffPlots("2016_QuadJet_TripleTag_3B");
  h_trig_2016_DoubleJet_TripleTag_3B = new TriggerEffPlots("2016_DoubleJet_TripleTag_3B"); 
  h_trig_2017_QuadJet_TripleTag_3B = new TriggerEffPlots("2017_QuadJet_TripleTag_3B");
  h_trig_2018_QuadJet_TripleTag_3B = new TriggerEffPlots("2018_QuadJet_TripleTag_3B");
  h_trig_2016_QuadJet_TripleTag_2b2c = new TriggerEffPlots("2016_QuadJet_TripleTag_2b2c");
  h_trig_2016_DoubleJet_TripleTag_2b2c = new TriggerEffPlots("2016_DoubleJet_TripleTag_2b2c"); 
  h_trig_2017_QuadJet_TripleTag_2b2c = new TriggerEffPlots("2017_QuadJet_TripleTag_2b2c");
  h_trig_2018_QuadJet_TripleTag_2b2c = new TriggerEffPlots("2018_QuadJet_TripleTag_2b2c");

  h_trig_2017_QuadJet_TripleTag_RunB_ideal = new TriggerEffPlots("2017_QuadJet_TripleTag_RunB_ideal");
  h_trig_2017_QuadJet_TripleTag_RunB_3B = new TriggerEffPlots("2017_QuadJet_TripleTag_RunB_3B");
  h_trig_2017_QuadJet_TripleTag_RunB_2b2c = new TriggerEffPlots("2017_QuadJet_TripleTag_RunB_2b2c");

  h_trig_2017_QuadJet_noTagV2_ideal = new TriggerEffPlots("2017_QuadJet_noTagV2_ideal");
  h_trig_2017_QuadJet_noTagV3_ideal = new TriggerEffPlots("2017_QuadJet_noTagV3_ideal");
  h_trig_2017_QuadJet_noTagV4_ideal = new TriggerEffPlots("2017_QuadJet_noTagV4_ideal");
  h_trig_2017_QuadJet_noTagV5_ideal = new TriggerEffPlots("2017_QuadJet_noTagV5_ideal");

  h_trig_2018_QuadJet_noTagV2_ideal = new TriggerEffPlots("2018_QuadJet_noTagV2_ideal");
  h_trig_2018_QuadJet_noTagV3_ideal = new TriggerEffPlots("2018_QuadJet_noTagV3_ideal");
  h_trig_2018_QuadJet_noTagV4_ideal = new TriggerEffPlots("2018_QuadJet_noTagV4_ideal");
  h_trig_2018_QuadJet_noTagV5_ideal = new TriggerEffPlots("2018_QuadJet_noTagV5_ideal");

  // Set up the CutFlows (for events) 
  h_evt_MC_cutflow = new TH1D("VbbHcc_MC_CutFlow", "", 2, 0, 2);
  h_evt_MC_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_MC_cutflow->GetXaxis()->SetBinLabel(2, "Passed daughter selection");
  
  h_evt_tagOnly_cutflow = new TH1D("VbbHcc_tagOnly_CutFlow", "", 8, 0, 8);
  h_evt_tagOnly_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_tagOnly_cutflow->GetXaxis()->SetBinLabel(2, "MET cut");
  h_evt_tagOnly_cutflow->GetXaxis()->SetBinLabel(3, "jet cuts");
  h_evt_tagOnly_cutflow->GetXaxis()->SetBinLabel(4, "triggers");
  h_evt_tagOnly_cutflow->GetXaxis()->SetBinLabel(5, "b-tag #1");
  h_evt_tagOnly_cutflow->GetXaxis()->SetBinLabel(6, "b-tag #2");
  h_evt_tagOnly_cutflow->GetXaxis()->SetBinLabel(7, "c-tag #1");
  h_evt_tagOnly_cutflow->GetXaxis()->SetBinLabel(8, "c-tag #2");
  
  h_evt_algoFirst_cutflow = new TH1D("VbbHcc_algoFirst_CutFlow", "", 8, 0, 8);
  h_evt_algoFirst_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_algoFirst_cutflow->GetXaxis()->SetBinLabel(2, "MET cut");
  h_evt_algoFirst_cutflow->GetXaxis()->SetBinLabel(3, "jet cuts");
  h_evt_algoFirst_cutflow->GetXaxis()->SetBinLabel(4, "triggers");
  h_evt_algoFirst_cutflow->GetXaxis()->SetBinLabel(5, "b-tag #1");
  h_evt_algoFirst_cutflow->GetXaxis()->SetBinLabel(6, "b-tag #2");
  h_evt_algoFirst_cutflow->GetXaxis()->SetBinLabel(7, "c-tag #1");
  h_evt_algoFirst_cutflow->GetXaxis()->SetBinLabel(8, "c-tag #2");
  
  h_evt_tagFirst_cutflow = new TH1D("VbbHcc_tagFirst_CutFlow", "", 5, 0, 5);
  h_evt_tagFirst_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_tagFirst_cutflow->GetXaxis()->SetBinLabel(2, "MET cut");
  h_evt_tagFirst_cutflow->GetXaxis()->SetBinLabel(3, "jet cuts");
  h_evt_tagFirst_cutflow->GetXaxis()->SetBinLabel(4, "triggers");
  h_evt_tagFirst_cutflow->GetXaxis()->SetBinLabel(5, "tags cut");
  
  h_evt_MCjet_ideal_cutflow = new TH1D("VbbHcc_MCjet_ideal_CutFlow", "", 7, 0, 7);
  h_evt_MCjet_ideal_cutflow->GetXaxis()->SetBinLabel(1, "Total");
  h_evt_MCjet_ideal_cutflow->GetXaxis()->SetBinLabel(2, "Jet Multiplicity");
  h_evt_MCjet_ideal_cutflow->GetXaxis()->SetBinLabel(3, "b-jet 0 Matched");
  h_evt_MCjet_ideal_cutflow->GetXaxis()->SetBinLabel(4, "b-jet 1 Matched");
  h_evt_MCjet_ideal_cutflow->GetXaxis()->SetBinLabel(5, "c-jet 0 Matched");
  h_evt_MCjet_ideal_cutflow->GetXaxis()->SetBinLabel(6, "c-jet 1 Matched");
  h_evt_MCjet_ideal_cutflow->GetXaxis()->SetBinLabel(7, "Jets Tagged");

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

  // Set up miscellaneous histograms
  h_nCombos = new TH1D("nCombos", "", 5, -0.5, 4.5);
  h_dR_ccjet = new TH1D("ccjet_dR", "", 100, 0, 10);
  h_dPhi_ccjet = new TH1D("ccjet_dPhi", "", 240, 0, 4.0);
  h_dR_bbjet = new TH1D("bbjet_dR", "", 100, 0, 10);
  h_dPhi_bbjet = new TH1D("bbjet_dPhi", "", 240, 0, 4.0); 

  h_mistag_leading = new TH1D("mistag_leading", "", 3, 0, 3);
  h_mistag_leading->GetXaxis()->SetBinLabel(1, "Total");
  h_mistag_leading->GetXaxis()->SetBinLabel(2, "mistag");
  h_mistag_leading->GetXaxis()->SetBinLabel(3, "proper");
  h_mistag_all = new TH1D("mistag_all", "", 3, 0, 3);
  h_mistag_all->GetXaxis()->SetBinLabel(1, "Total");
  h_mistag_all->GetXaxis()->SetBinLabel(2, "mistag");
  h_mistag_all->GetXaxis()->SetBinLabel(3, "proper");

  h_bRegCorr = new TH1D("bRegCorr", "", 200, 0, 2);
  h_cRegCorr = new TH1D("cRegCorr", "", 200, 0, 2);
  h_JetMass = new TH1D("JetMass", "", 300, 0, 300);  

  h_genWeight = new TH1D("genWeight", "", 500, -5, 5);
  h_puSF = new TH1D("puSF", "", 500, -5, 5);
  h_l1preW = new TH1D("l1preW", "", 500, -5, 5);
  h_trigSF = new TH1D("trigSF", "", 500, -5, 5);
  h_btagW = new TH1D("btagW", "", 500, -5, 5);
  h_ctagW = new TH1D("ctagW", "", 500, -5, 5);
  h_evtW = new TH1D("evtW", "", 500, -5, 5);
  h_lheW = new TH1D("lheW", "", 500, -5, 5);
  h_genW = new TH1D("genW", "", 500, -5, 5);
 
  h_nJet = new TH1D("nJet", "", 20, -0.5, 19.5);
  h_nAnalysisJet = new TH1D("nAnalysisJet", "", 20, -0.5, 19.5);

  h_nMuon = new TH1D("nMuon", "", 10, -0.5, 9.5);
  h_nElec = new TH1D("nElec", "", 10, -0.5, 9.5);

  // Add them to the return list so we can use them in our analyses.
  r->GetOutputList()->Add(h_evt);

  r->GetOutputList()->Add(h_bRegCorr);
  r->GetOutputList()->Add(h_cRegCorr);
  r->GetOutputList()->Add(h_JetMass);

  r->GetOutputList()->Add(h_genWeight);
  r->GetOutputList()->Add(h_puSF);
  r->GetOutputList()->Add(h_l1preW);
  r->GetOutputList()->Add(h_trigSF);
  r->GetOutputList()->Add(h_btagW);
  r->GetOutputList()->Add(h_ctagW);
  r->GetOutputList()->Add(h_evtW);
  r->GetOutputList()->Add(h_lheW); 
  r->GetOutputList()->Add(h_genW); 

  r->GetOutputList()->Add(h_nJet);
  r->GetOutputList()->Add(h_nAnalysisJet);

  // Add each plot in the VHPlot instances
  std::vector<TH1*> tmp = h_VH_MC->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_MCjet_minDR->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_MCjet_minDR_noTag->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_MCjet_dRcollect->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_MCjet_dRcollect_noTag->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_MCjet_ideal->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_MCjet_DHZ->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  
  tmp = h_VH_tagOnly->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_tagOnly_noMassCorr->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_tagOnly_noJEC->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_tagOnly_2b1c->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  
  tmp = h_VH_algoFirst->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_algoFirst_noMassCorr->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_algoFirst_noJEC->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_algoFirst_2b1c->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  
  tmp = h_VH_tagFirst->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_tagFirst_noMassCorr->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_tagFirst_noJEC->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_tagFirst_2b1c->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  tmp = h_VH_all->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_select->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_alljet->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_VH_seljet->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  // Add each plot in the JetPlot instances
  tmp = h_jets_selected->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_jets_all->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_jets_MC->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_jets_4b->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_jets_2b2c->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  // Add each plot in the EffPlot instances
  tmp = h_eff_tagOnly->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_eff_algoFirst->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_eff_tagFirst->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  tmp = h_genJet_all->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_genJet_cuts->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_genJet_VbbHcc->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  // Add each plot in the TriggerEffPlot instances
  tmp = h_trig_2016_QuadJet_TripleTag->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2016_QuadJet_DoubleTag->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2016_DoubleJet_TripleTag->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2016_DoubleJet_DoubleTag->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2017_QuadJet_TripleTag->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2017_QuadJet_noTag->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2018_QuadJet_TripleTag->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2018_QuadJet_noTag->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);  

  tmp = h_trig_2016_QuadJet_TripleTag_tagged->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2016_QuadJet_DoubleTag_tagged->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2016_DoubleJet_TripleTag_tagged->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2016_DoubleJet_DoubleTag_tagged->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2017_QuadJet_TripleTag_tagged->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2017_QuadJet_noTag_tagged->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2018_QuadJet_TripleTag_tagged->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2018_QuadJet_noTag_tagged->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  tmp = h_trig_2016_QuadJet_TripleTag_ideal->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2016_QuadJet_DoubleTag_ideal->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2016_DoubleJet_TripleTag_ideal->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2016_DoubleJet_DoubleTag_ideal->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2017_QuadJet_TripleTag_ideal->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2017_QuadJet_noTag_ideal->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2018_QuadJet_TripleTag_ideal->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2018_QuadJet_noTag_ideal->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  tmp = h_trig_2016_QuadJet_TripleTag_3B->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2016_DoubleJet_TripleTag_3B->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2017_QuadJet_TripleTag_3B->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2018_QuadJet_TripleTag_3B->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2016_QuadJet_TripleTag_2b2c->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2016_DoubleJet_TripleTag_2b2c->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2017_QuadJet_TripleTag_2b2c->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2018_QuadJet_TripleTag_2b2c->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  tmp = h_trig_2017_QuadJet_TripleTag_RunB_ideal->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2017_QuadJet_TripleTag_RunB_3B->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_trig_2017_QuadJet_TripleTag_RunB_2b2c->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  
  // Add each plot in the RecoPlot instances
  tmp = h_reco_minDR->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_reco_minDR_under30->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_reco_minDR_noTag->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_reco_minDR_noTag_under30->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_reco_dRcollect->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_reco_dRcollect_under30->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_reco_dRcollect_noTag->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_reco_dRcollect_noTag_under30->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_reco_ideal->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_reco_ideal_under30->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_reco_DHZ->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
  tmp = h_reco_DHZ_under30->returnHisto();
  for(size_t i=0; i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  // Add all other plots we may desire.
  r->GetOutputList()->Add(h_evt_VbbHcc);
  r->GetOutputList()->Add(h_evt_MC_cutflow);
  r->GetOutputList()->Add(h_evt_tagOnly_cutflow);
  r->GetOutputList()->Add(h_evt_algoFirst_cutflow);
  r->GetOutputList()->Add(h_evt_tagFirst_cutflow);
  r->GetOutputList()->Add(h_evt_MCjet_ideal_cutflow);

  r->GetOutputList()->Add(h_jet_cutflow);
  r->GetOutputList()->Add(h_elec_cutflow);
  r->GetOutputList()->Add(h_muon_cutflow);  
  r->GetOutputList()->Add(h_nCombos);
  r->GetOutputList()->Add(h_dR_ccjet);
  r->GetOutputList()->Add(h_dPhi_ccjet);
  r->GetOutputList()->Add(h_dR_bbjet);
  r->GetOutputList()->Add(h_dPhi_bbjet);

  r->GetOutputList()->Add(h_mistag_leading);
  r->GetOutputList()->Add(h_mistag_all);

  r->GetOutputList()->Add(h_nMuon);
  r->GetOutputList()->Add(h_nElec);

}// end SlaveBegin

// ============================================================================
// PROCESS
// ============================================================================

/*
 * This is where we do the actual analysis. We need to do the 
 * following steps for each event:
 *
 * 1. determine the events for the weight
 * 2. reconstruct the Physics objects
 * 3. analyze them as desired
 * 4. pray that nothing fails
*/

void VH_selection::Process(Reader* r) { 

  //=================================================================
  // WEIGHTS & GEN INFORMATION
  //=================================================================

  float genWeight = 1.;
  float puSF = 1.;
  float l1preW = 1.;

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)

  // Determine the generator weights
  if (*(r->genWeight) < 0) genWeight = -1.;
  if (*(r->genWeight) == 0) {
    genWeight = 0;
    h_evt->Fill(0);
  }
  if (*(r->genWeight) < 0) h_evt->Fill(-1);
  if (*(r->genWeight) > 0) h_evt->Fill(1);

  h_genWeight->Fill(*(r->genWeight));

  // Use central general weights to normalize generator weights
  if (m_centralGenWeight != 0) genWeight *= *(r->genWeight)/m_centralGenWeight;
  puSF = PileupSF(*(r->Pileup_nTrueInt));
  h_puSF->Fill(puSF);

#endif

  h_evt->Fill(2, genWeight); 

#if defined(MC_2016) || defined(MC_2017)

  // Handle the L1 Pre-Firing weight for 2016-2017
  l1preW = *(r->L1PreFiringWeight_Nom);

  if (m_l1prefiringType == "l1prefiringu") l1preW = *(r->L1PreFiringWeight_Up);
  if (m_l1prefiringType == "l1prefiringd") l1preW = *(r->L1PreFiringWeight_Dn);

  h_l1preW->Fill(l1preW);

#endif

#if defined(DATA_2016) || defined(DATA_2017) || defined(DATA_2018)

  h_evt->Fill(-1);
  if (!m_lumiFilter.Pass(*(r->run), *(r->luminosityBlock))) return;
  h_evt->Fill(1);

#endif

  float evtW = 1.;

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
  //if (!m_isData) 
  // Only set the event weights with the scales if we're in MC, not data.
  evtW *= genWeight * puSF * l1preW;

  //h_lheW->Fill(*(r->LHEWeight));
  //h_genW->Fill(*(r->Generator_Weight));
#endif

  // We don't want to have to change the working point in several spots,
  // so we choose here and then we can use these variables wherever needed.
  float desired_BvL = CUTS.Get<float>("BvL_mediumWP_deepJet");
  float desired_CvL = CUTS.Get<float>("CvL_mediumWP_deepJet");
  float desired_CvB = CUTS.Get<float>("CvB_mediumWP_deepJet");

  //=================================================================
  // RECONSTRUCT PHYSICS OBJECTS
  //=================================================================
  
  // ==== JETS (NO KINEMATIC CUTS) ====
  std::vector<JetObj> jets;
  for (unsigned int i = 0; i < *(r->nJet); ++i) {

    // Get the flavor of the jets
    int jetFlav = -1;
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
    jetFlav = (r->Jet_hadronFlavour)[i];
#endif

    // Make sure to properly correct the jet energy
    Float_t pt = (r->Jet_pt)[i];
    Float_t jet_mass = (r->Jet_mass)[i];
    Float_t bRegCorr = (r->Jet_bRegCorr)[i];
    Float_t cRegCorr = (r->Jet_cRegCorr)[i];
    Float_t bRegRes  = (r->Jet_bRegRes)[i];
    Float_t cRegRes  = (r->Jet_cRegRes)[i];

    h_bRegCorr->Fill(bRegCorr, evtW);
    h_cRegCorr->Fill(cRegCorr, evtW);
    h_JetMass->Fill(jet_mass, evtW);

    // Reconstruct the jet
    JetObj jet(pt, (r->Jet_eta)[i], (r->Jet_phi)[i], jet_mass,
      jetFlav, (r->Jet_btagDeepFlavB)[i], (r->Jet_puId)[i]);

    jet.StoreRegInfo(bRegCorr, bRegRes, cRegCorr, cRegRes);   

    // Handle getting the c-tag value which is different
    // depending on the version of NanoAOD.
#if defined(NANOAODV9)
    jet.m_deepCvL = (r->Jet_btagDeepFlavCvL)[i];
    jet.m_deepCvB = (r->Jet_btagDeepFlavCvB)[i];
#endif

#if defined(NANOAODV7)
    jet.m_deepCvL = (r->Jet_btagDeepFlavC)[i];
    jet.m_deepCvB = (r->FatJet_btagDDCvB)[i];
#endif

    // If we're in a MC file, let's note which gen jet
    // is associated with thsi reco jet.
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
    jet.m_genJetIdx = (r->Jet_genJetIdx)[i];
#endif


    // Add the jet to our overall list.
    jet.SetIdxAll(i);
    jets.push_back(jet);

  }//end-jets

  int year = 0;
#if defined(MC_2016)
  year = 16;
#endif
#if defined(MC_2017)
  year = 17;
#endif
#if defined(MC_2018)
  year = 18;
#endif

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)

  // Modify the event weight account for the HT in the event
  float trig_SF = determine_trigger_SF(jets, year);
  evtW *= trig_SF;
  h_trigSF->Fill(trig_SF);

  // Modify the event weight to account for the btag & ctag weights
  float btagW = CalBtagWeight(jets, CUTS.GetStr("jet_main_btagWP"), m_btagUncType); 
  h_btagW->Fill(btagW);
  evtW *= btagW;
  
  float ctagW = CalCtagWeight(jets, CUTS.GetStr("jet_main_btagWP"), m_ctagUncType);
  h_ctagW->Fill(ctagW);
  evtW *= ctagW;

  h_evtW->Fill(evtW);
#endif
   
  // ==== ELECTRONS ====
  std::vector<LepObj> elecs;
  for (unsigned int i = 0; i < *(r->nElectron); ++i) {

    // Produce an electron from the information.
    h_elec_cutflow->Fill(0.5, genWeight); // all electrons

    float etaSC = (r->Electron_eta)[i] - (r->Electron_deltaEtaSC)[i];
    LepObj elec((r->Electron_pt)[i], (r->Electron_eta)[i], etaSC,
           (r->Electron_phi)[i], (r->Electron_mass)[i], i,
           (r->Electron_charge)[i], 0);

    // Cut based on the pT & eta values
    if (elec.Pt() < CUTS.Get<float>("lep_pt1")) continue;
    h_elec_cutflow->Fill(1.5, genWeight); // pass pT cut
    if (fabs(elec.Eta()) > CUTS.Get<float>("lep_eta")) continue;
    h_elec_cutflow->Fill(2.5, genWeight); // pass eta cut

    // Cut based on SC value
    if (fabs(etaSC) < 1.566 && fabs(etaSC) > 1.442) continue;
    h_elec_cutflow->Fill(3.5, genWeight); // pass SC cut

    // Cut based on the isolation
    int elecID = r->Electron_cutBased[i];
    if (elecID < 2) continue; // loose electron ID
    h_elec_cutflow->Fill(4.5, genWeight); // pass ID cut

    // Add the electrons to the list
    elecs.push_back(elec);

  }//end-elecs

  h_nElec->Fill(elecs.size(), evtW);

  // ==== MUONS ====
  std::vector<LepObj> muons;
  for (unsigned int i = 0; i < *(r->nMuon); ++i) {

    // Produce a muon from the information
    h_muon_cutflow->Fill(0.5, genWeight); // all muons
    LepObj muon((r->Muon_pt)[i], (r->Muon_eta)[i], -1, (r->Muon_phi)[i],
           (r->Muon_mass)[i], i, (r->Muon_charge)[i],
           (r->Muon_pfRelIso04_all)[i]);

    // Cut based on the pT & eta values
    if (muon.Pt() < CUTS.Get<float>("lep_pt0")) continue;
    h_muon_cutflow->Fill(1.5, genWeight); // pass pT cut
    if (fabs(muon.Eta()) > CUTS.Get<float>("lep_eta")) continue;
    h_muon_cutflow->Fill(2.5, genWeight); // pass eta cut

    // Cut based on the loose ID
    if (r->Muon_looseId[i] <= 0) continue;
    h_muon_cutflow->Fill(3.5, genWeight); // pass ID cut

    // Cut based on an isolation value
    if (muon.m_iso > CUTS.Get<float>("muon_iso")) continue;
    h_muon_cutflow->Fill(4.5, genWeight); // pass iso cut

    // Add the muon to the list.
    muons.push_back(muon);

  }//end-muons

  h_nMuon->Fill(muons.size(), evtW);

  //=================================================================
  // START EVENT SELECTION
  //=================================================================
 
  // All CutFlows need to show the total events.
  h_evt_MC_cutflow->Fill(0.5, genWeight);   // MC Truth
  h_evt_tagOnly_cutflow->Fill(0.5, genWeight); // Tagging Only
  h_evt_algoFirst_cutflow->Fill(0.5, genWeight); // Matching Prioritized
  h_evt_tagFirst_cutflow->Fill(0.5, genWeight); // Tagging Prioritized

  h_VH_all->FillMET(*(r->MET_pt), evtW);
  if (*(r->MET_pt) >= CUTS.Get<float>("MET")) return;

  // Show that we pass the MET cut.
  h_evt_MC_cutflow->Fill(1.5, genWeight);
  h_evt_tagOnly_cutflow->Fill(1.5, genWeight);
  h_evt_algoFirst_cutflow->Fill(1.5, genWeight);
  h_evt_tagFirst_cutflow->Fill(1.5, genWeight);

  /********************************************************
  * MONTE CARLO TRUTH - we want to have the MC truth      *
  * values to get an idea of how well we are selecting    *
  * events. We should have several delta-like functions   *
  * in this case (mass specifically).                     *
  ********************************************************/
  
  bool is_VbbHcc_event = false;
  std::vector<std::vector<int> > dauIdxs;
  std::vector<JetObj> gen_bs;
  std::vector<JetObj> gen_cs;

  h_evt_VbbHcc->Fill(0.5, evtW);

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)

  // Make sure we have two daughters for each particle of the proper type.
  int idZ = 0, idH = 1;
  dauIdxs = DauIdxs_ZH(r);
  
  if (dauIdxs[idZ].size() == 2 && dauIdxs[idH].size() == 2) {
    is_VbbHcc_event = true;
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

    //b0.ApplyRegression(5); b1.ApplyRegression(5);
    std::vector<JetObj> MC_bjets{b0, b1};
    gen_bs.push_back(b0); gen_bs.push_back(b1);
    ZObj MC_Z(MC_bjets);

    // Two c quarks --> HObj
    JetObj c0((r->GenPart_pt)[idx1_H], (r->GenPart_eta)[idx1_H],
      (r->GenPart_phi)[idx1_H], (r->GenPart_mass)[idx1_H], 4, 0., 0.);
    JetObj c1((r->GenPart_pt)[idx2_H], (r->GenPart_eta)[idx2_H],
      (r->GenPart_phi)[idx2_H], (r->GenPart_mass)[idx2_H], 4, 0., 0.);

    //c0.ApplyRegression(4); c1.ApplyRegression(4);
    std::vector<JetObj> MC_cjets{c0, c1};
    gen_cs.push_back(c0); gen_cs.push_back(c1);
    HObj MC_H(MC_cjets);

    // Fill the histograms
    h_VH_MC->FillVH(MC_Z, MC_H, evtW);

    is_VbbHcc_event = true;
    h_evt_VbbHcc->Fill(1.5, evtW);

  }//end

#endif

  /********************************************************
  * MONTE CARLO TRUTH JETS - In the previous MC truth, we *
  * only have the quarks. Here, we want to find the jets  *
  * that match the MC quarks.                             *    
  ********************************************************/

  std::vector<JetObj> gen_bjets;
  std::vector<JetObj> gen_cjets;
  bool found_MCjets = false;

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)

  // Go through all the gen jets
  std::vector<JetObj> genJet_all;
  std::vector<JetObj> genJet_cut;
  for (unsigned int i = 0; i < *(r->nGenJet); ++i) {
    
    // Build a jet from the information.
    Int_t flavor = (r->GenJet_partonFlavour)[i];
    Float_t eta = (r->GenJet_eta)[i], pt = (r->GenJet_pt)[i];
    JetObj gjet(pt, eta, (r->GenJet_phi)[i],
      (r->GenJet_mass)[i], flavor, 0, 0);

    genJet_all.push_back(gjet);
    
    // Check the cuts we're interested in for jets.
    if (pt > CUTS.Get<float>("jet_pt") && fabs(eta) < CUTS.Get<float>("jet_eta")) {
      genJet_cut.push_back(gjet);
    }
  }

  h_genJet_all->Fill(genJet_all, evtW);
  h_genJet_cuts->Fill(genJet_cut, evtW);
  
  // We only want to check MC events where we know that we have c-jets. 
  // We need this because of the overall file is H->cc and Z->QQ.
  if (is_VbbHcc_event) {
    
    // Make a copy of the jets
    std::vector<JetObj> jetlist; jetlist = jets;
    std::vector<JetObj> genJet_list;  
    std::vector<JetObj> genCjet_list;
    std::vector<JetObj> genBjet_list;

    // Make a list of the GenJets that we can pull from.
    int nL = 0, nC = 0, nB = 0, nGlu = 0;
    for (unsigned int i = 0; i < *(r->nGenJet); ++i) {
  
      // Build a jet from the information
      Int_t flavor = (r->GenJet_partonFlavour)[i];
      if ((r->GenJet_pt)[i] < 45.0 || abs((r->GenJet_eta)[i]) > 2.4) continue;     
      JetObj gjet((r->GenJet_pt)[i], (r->GenJet_eta)[i], (r->GenJet_phi)[i], 
        (r->GenJet_mass)[i], flavor, 0, 0);

      // Sort the gets into lists based on flavors
      genJet_list.push_back(gjet);
      if (abs(flavor) <= 3) nL++;
      else if (abs(flavor) == 4){ genCjet_list.push_back(gjet); nC++; }
      else if (abs(flavor) == 5){ genBjet_list.push_back(gjet); nB++; }
      else if (abs(flavor) == 21) nGlu++;
    }
    
    h_genJet_VbbHcc->Fill(genJet_list, evtW);
    
    // We now have 4 options we want to go through to see how we can properly match the jets.
    // Method 1 (original) - match partons to nearest properly-tagged jets.
    // Method 2 - match partons to nearest jets (without tagging)
    // Method 3 - combine any tagged jets within the desired range & match that combo-jet
    // Method 4 - combine ANY jets within the desired range & match that combo-jet
    
    // Methods 1 & 3 go here...
    float dR_cut = 0.4;
    if (genCjet_list.size() >= 2 && genBjet_list.size() >= 2) {
    
      // ========================
      // METHOD #1
      // ========================
      
      // Make a copy of the lists
      std::vector<JetObj> cjets = genCjet_list;
      std::vector<JetObj> bjets = genBjet_list;
    
      // Find the b-jets closest to the b-quarks
      std::vector<std::pair<int,float>> bpair_idx_dR;
      std::vector<int> chosenIdx;
      
      // Go through each b-parton...
      for (size_t i = 0; i < gen_bs.size(); ++i) {
        // ... and check the separation with each possible jet.
        for (size_t j = 0; j < bjets.size(); ++j) {
          float dR = fabs(gen_bs[i].m_lvec.DeltaR(bjets[j].m_lvec));
          bpair_idx_dR.push_back(std::make_pair(j,dR));
        }
        
        // Sort the matches via the dR values (second value in pair). We sort
        // in ascending order, so our proper choice will be the first option.
        std::sort(bpair_idx_dR.begin(), bpair_idx_dR.end(), sort_by_second);
        std::pair<int, float> proper_pair = bpair_idx_dR[0];
        
        // Move the proper jet to our list of chosen jets
        // and remove it as an option from the list.
        int idx = proper_pair.first;
        chosenIdx.push_back(idx);
        gen_bjets.push_back(bjets[idx]);
        bjets.erase(bjets.begin() + idx);
      }//end-bs
      
      // Find the c-jets closest to the c-quarks
      std::vector<std::pair<int,float>> cpair_idx_dR;
      std::vector<int> chosenIdx2;
      
      // Go through each c-parton...
      for (size_t i = 0; i < gen_cs.size(); ++i) {
        // ... and check the separation with each possible jet
        for (size_t j = 0; j < cjets.size(); ++j) {
          float dR = fabs(gen_cs[i].m_lvec.DeltaR(cjets[j].m_lvec));
          cpair_idx_dR.push_back(std::make_pair(j,dR));
        }
        
        // Sort the matches via the dR values (second value in pair). We sort
        // in ascending order, so our proper choice will be the first option.
        std::sort(bpair_idx_dR.begin(), bpair_idx_dR.end(), sort_by_second);
        std::pair<int, float> proper_pair = bpair_idx_dR[0];
        
        // Move the proper jet to our list of chosen jets
        // and remove it as an option from the list.
        int idx = proper_pair.first;
        chosenIdx2.push_back(idx);
        gen_cjets.push_back(cjets[idx]);
        cjets.erase(cjets.begin() + idx);
      }//end-cs
      
      // From these pairs, reconstruct the proper bosons
      // and fill our desired methods
      ZObj Z_MCjet(gen_bjets); HObj H_MCjet(gen_cjets);
      h_VH_MCjet_minDR->FillVH(Z_MCjet, H_MCjet, evtW);
      h_reco_minDR->Fill(gen_bs, gen_bjets, gen_cs, gen_cjets, evtW);
      if (Z_MCjet.M() < 30.0 || H_MCjet.M() < 30.0)
        h_reco_minDR_under30->Fill(gen_bs, gen_bjets, gen_cs, gen_cjets, evtW);
      
           
      // ==================
      // Method #3
      // ==================
      
      std::vector<JetObj> cjets2 = genCjet_list; 
      std::vector<JetObj> bjets2 = genBjet_list;
      std::vector<JetObj> genBjets; std::vector<JetObj> genCjets;
      
      // For each b-parton, take in any jet that's within our range
      //std::cout << "gen_bs.size() = " << gen_bs.size() << std::endl;
      for (size_t i = 0; i < gen_bs.size(); ++i) {
      
        //std::cout << ">> i = " << i << std::endl;
        JetObj comb_bjet(0, 0, 0, 0, 0, 0, 0);
        
        // For each jet in the b-list...
        //std::cout << ">> bjets2.size() = " << bjets2.size() << std::endl;
        for (Int_t j = bjets2.size() - 1; j >= 0; --j) {
          
          //std::cout << ">>>> j = " << j << std::endl; 
          // Check for the separation. If the separation is less than our cut...
          float dR = fabs(gen_bs[i].m_lvec.DeltaR(bjets2[j].m_lvec));
          if (dR < dR_cut) {
            // Combine the jet with our ones so far 
            // and remove it from the list
            comb_bjet.m_lvec += bjets2[j].m_lvec;
            bjets2.erase(bjets2.begin() + j);
          }
        }//end-j
        
        genBjets.push_back(comb_bjet);
        
      }//end-i
      
      // For each c-parton, take in any jet that's within our range
      for (size_t i = 0; i < gen_cs.size(); ++i) {
      
        JetObj comb_cjet(0, 0, 0, 0, 0, 0, 0);
        
        // For each jet in the c-list...
        for (Int_t j = cjets2.size() - 1; j >= 0; --j) {
        
          // Check for separation. If the separation is less than our cut...
          float dR = fabs(gen_cs[i].m_lvec.DeltaR(cjets2[j].m_lvec));
          if (dR < dR_cut) {
            // Combine the jet with our ones so far 
            // and remove it from the list
            comb_cjet.m_lvec += cjets2[j].m_lvec;
            cjets2.erase(cjets2.begin()+j);
          }
        
        }//end-j
        
        genCjets.push_back(comb_cjet);
      
      }//end-i
      
      // Check to make sure that we properly found each jet.
      // Each jet is required to be > 45.0 GeV, so if we see
      // values less than this, we did not find jets.
      bool found_jets = true;
      for (Int_t i = 0; i < genBjets.size(); ++i) {
        if (genBjets[i].Pt() < 45.0) found_jets = false; }
      for (Int_t j = 0; j < genCjets.size(); ++j) {
        if (genCjets[j].Pt() < 45.0) found_jets = false; }
      
      // From these pairs, reconstruct the proper bosons
      // and fill our desired methods. 
      if (found_jets) {
        ZObj Z_MCjet2(genBjets); HObj H_MCjet2(genCjets);
        h_VH_MCjet_dRcollect->FillVH(Z_MCjet2, H_MCjet2, evtW);
        h_reco_dRcollect->Fill(gen_bs, genBjets, gen_cs, genCjets, evtW);
        if (Z_MCjet.M() < 30.0 || H_MCjet.M() < 30.0)
          h_reco_dRcollect_under30->Fill(gen_bs, genBjets, gen_cs, genCjets, evtW);
      }
      
    }//end-method-1-3
    
    // Methods 2 & 4 go here...
    
    if (genJet_list.size() >= 4) {
    
      // Make a copy of the list of jets
      std::vector<JetObj> alljets = genJet_list;
      
      // ================
      // Method #2
      // ================
      
      std::vector<JetObj> bjets;
      std::vector<JetObj> cjets;
      
      // Compare each parton to each jet and see which jet gives us
      // the closest jet.
      float dR_min = 9999.0;
      int chosenIdx = -1;
      
      // Go through each b-parton
      for (size_t i = 0; i < gen_bs.size(); ++i) {
      
        for (size_t j = 0; j < alljets.size(); ++j) {
          float dR = fabs(gen_bs[i].m_lvec.DeltaR(alljets[j].m_lvec));
          if (dR < dR_min) { dR_min = dR; chosenIdx = j; }
        }
        
        // Take the jet we found that was closest and add it to our list
        bjets.push_back(alljets[chosenIdx]);
        alljets.erase(alljets.begin() + chosenIdx);
      
      }//end-i
      
      // Go through each c-parton with the remaining jets
      for (size_t i = 0; i < gen_cs.size(); ++i) {
      
        for (size_t j = 0; j < alljets.size(); ++j) {
          float dR = fabs(gen_cs[i].m_lvec.DeltaR(alljets[j].m_lvec));
          if (dR < dR_min) { dR_min = dR; chosenIdx = j; }
        }
        
        // Take the jet we found that was closest and add it to our list
        cjets.push_back(alljets[chosenIdx]);
        alljets.erase(alljets.begin() + chosenIdx);
      
      }//end-i
      
      // From these pairs, reconstruct the proper bosons
      // and fill our desired methods
      ZObj Z_MCjet(bjets); HObj H_MCjet(cjets);
      h_VH_MCjet_minDR_noTag->FillVH(Z_MCjet, H_MCjet, evtW);
      h_reco_minDR_noTag->Fill(gen_bs, bjets, gen_cs, cjets, evtW);
      if (Z_MCjet.M() < 30.0 || H_MCjet.M() < 30.0)
        h_reco_minDR_noTag_under30->Fill(gen_bs, bjets, gen_cs, cjets, evtW);
      
      // =====================
      // Method #4
      // =====================
      
      std::vector<JetObj> alljets2 = genJet_list;
      std::vector<JetObj> genBjets;
      std::vector<JetObj> genCjets;
      
      // For each parton, pull in any jets within our dR range
      for (size_t i = 0; i < gen_bs.size(); ++i) {
      
        JetObj comb_bjet(0, 0, 0, 0, 0, 0, 0);
        
        // For each jet in the b-list...
        for (Int_t j = alljets2.size() - 1; j >= 0; --j) {
        
          // Check for the separation. If the separation is less than our cut...
          float dR = fabs(gen_bs[i].m_lvec.DeltaR(alljets2[j].m_lvec));
          if (dR < dR_cut) {
            // Combine the jet with our ones so far 
            // and remove it from the list
            comb_bjet.m_lvec += alljets2[j].m_lvec;
            alljets2.erase(alljets2.begin() + j);
          }
        }//end-j
        
        genBjets.push_back(comb_bjet);
        
      }//end-i
      
      // For each c-parton, take in any jet that's within our range
      for (size_t i = 0; i < gen_cs.size(); ++i) {
      
        JetObj comb_cjet(0, 0, 0, 0, 0, 0, 0);
        
        // For each jet in the c-list...
        for (Int_t j = alljets2.size() - 1; j >= 0; --j) {
        
          // Check for separation. If the separation is less than our cut...
          float dR = fabs(gen_cs[i].m_lvec.DeltaR(alljets2[j].m_lvec));
          if (dR < dR_cut) {
            // Combine the jet with our ones so far 
            // and remove it from the list
            comb_cjet.m_lvec += alljets2[j].m_lvec;
            alljets2.erase(alljets2.begin() + j);
          }
        
        }//end-j
        
        genCjets.push_back(comb_cjet);
      
      }//end-i
      
      bool found_jets = true;
      for (Int_t i = 0; i < genBjets.size(); ++i) {
        if (genBjets[i].Pt() < 45.0) found_jets = false; }
      for (Int_t j = 0; j < genCjets.size(); ++j) {
        if (genCjets[j].Pt() < 45.0) found_jets = false; }
      
      // From these pairs, reconstruct the proper bosons
      // and fill our desired methods
      if (found_jets) {
        ZObj Z_MCjet2(genBjets); HObj H_MCjet2(genCjets);
        h_VH_MCjet_dRcollect_noTag->FillVH(Z_MCjet2, H_MCjet2, evtW);
        h_reco_dRcollect_noTag->Fill(gen_bs, genBjets, gen_cs, genCjets, evtW);
        if (Z_MCjet.M() < 30.0 || H_MCjet.M() < 30.0)
          h_reco_dRcollect_noTag_under30->Fill(gen_bs, genBjets, gen_cs, genCjets, evtW);
      }
    
    }//end-method-2-4
    
    // ===================
    // IDEAL METHOD
    // ===================
    h_evt_MCjet_ideal_cutflow->Fill(0.5, evtW);
    if (genJet_list.size() >= 4) {
    
      h_evt_MCjet_ideal_cutflow->Fill(1.5, evtW);
    
      // Make a copy of our list of jets
      std::vector<JetObj> alljets = genJet_list;
      std::vector<JetObj> genBjets;
      std::vector<JetObj> genCjets;
      
      // Check each parton versus our jets. Check to see how many jets
      // are within a dR < 0.4 cut of our parton. If we have more than
      // one jet in ANY of the cases, do not continue.
      bool can_continue = true;
      
      for (size_t i = 0; i < gen_bs.size(); ++i) {
      
        JetObj comb_bjet(0, 0, 0, 0, 0, 0, 0);
        int number_found = 0;
        
        if (!can_continue) continue;
        
        // For each jet in the list...
        for (Int_t j = alljets.size() - 1; j >= 0; --j) {
        
          // Check for the separation. If the separation is less than our cut...
          float dR = fabs(gen_bs[i].m_lvec.DeltaR(alljets[j].m_lvec));
          if (!can_continue) continue;
          
          if (dR < dR_cut) {
            // Combine the jet with our ones so far 
            // and remove it from the list
            number_found++;
            comb_bjet.m_lvec += alljets[j].m_lvec;
            comb_bjet.m_flav += alljets[j].m_flav;
            alljets.erase(alljets.begin() + j);
            if (number_found > 1) { can_continue = false; break; }
            if (i == 0) h_evt_MCjet_ideal_cutflow->Fill(2.5, evtW);
            if (i == 1) h_evt_MCjet_ideal_cutflow->Fill(3.5, evtW);
          }
        }//end-j
        
        genBjets.push_back(comb_bjet);
        
      }//end-i
      
      // Only continue to check the c-jets if we didn't have an issue
      // with the b-jets above.
      if (can_continue) {
      
        // For each c-parton...
        for (size_t i = 0; i < gen_cs.size(); ++i) {
        
          JetObj comb_cjet(0, 0, 0, 0, 0, 0, 0);
          int number_found = 0;
          
          if (!can_continue) continue;
          
          // For each jet in the list...
          for (Int_t j = alljets.size() - 1; j >= 0; --j) {
        
            // Check for the separation. If the separation is less than our cut...
            float dR = fabs(gen_cs[i].m_lvec.DeltaR(alljets[j].m_lvec));
            if (!can_continue) continue;
          
            if (dR < dR_cut) {
              // Combine the jet with our ones so far 
              // and remove it from the list
              number_found++;
              comb_cjet.m_lvec += alljets[j].m_lvec;
              comb_cjet.m_flav += alljets[j].m_flav;
              alljets.erase(alljets.begin() + j);
              if (number_found > 1) { can_continue = false; break; }
              if (j == 0) h_evt_MCjet_ideal_cutflow->Fill(4.5, evtW);
              if (j == 1) h_evt_MCjet_ideal_cutflow->Fill(5.5, evtW);
            }
          }//end-j
          
          genCjets.push_back(comb_cjet);
        
        }//end-i
        
        // Only continue on if we didn't have an issue with the c-jets.
        if (can_continue) {
        
          // If we've reached this point, we have one jet for each c-parton
          // and one jet for each b-parton. However, we do not know if the
          // jets are actually b/c-tagged or not. We need to check to see
          // if it's true.
          if (abs(genBjets[0].m_flav) == 5 && abs(genBjets[1].m_flav) == 5 &&
              abs(genCjets[0].m_flav) == 4 && abs(genCjets[1].m_flav) == 4) {
              
            h_evt_MCjet_ideal_cutflow->Fill(6.5, evtW);
          
            // If we've found the flavors of each jet, add them to our 
            // proper plots.
            ZObj Z_MCjet(genBjets); HObj H_MCjet(genCjets);
            h_VH_MCjet_ideal->FillVH(Z_MCjet, H_MCjet, evtW);
            h_reco_ideal->Fill(gen_bs, genBjets, gen_cs, genCjets, evtW);
            if (Z_MCjet.M() < 30.0 || H_MCjet.M() < 30.0)
              h_reco_ideal_under30->Fill(gen_bs, genBjets, gen_cs, genCjets, evtW);
          
          }
          
        }
      
      }
    
    }//end-ideal-method
    
    // ===================
    // DHZ METHOD
    // ===================
    if (genJet_list.size() >= 4) {
      
      //std::cout << "checking DHZ on gen jets..." << std::endl;    

      // Make a copy of the jets & sort them by pT
      std::vector<JetObj> jets = genJet_list;
      std::sort(jets.begin(), jets.end(), JetObj::JetCompPt());
      
      // Create the objects for the distance calculations and make sure
      // we sort them in ascending order. (All necessary calculations
      // for alogrithms are handled within the DObj class.)
      DHZObj d0(jets, 0, 1, 2, 3);          // H(0,1) ; Z(2,3)
      DHZObj d1(jets, 0, 2, 1, 3);          // H(0,2) ; Z(1,3)
      DHZObj d2(jets, 0, 3, 1, 2);          // H(0,3) ; Z(1,2)
      std::vector<DHZObj> distances {d0, d1, d2};
      std::sort(distances.begin(), distances.end());
      
      // Determine the distance between the closest two. Then, based on this
      // distance, determine which pair we want to use. If we are outside the
      // resolution window (30 GeV), we can just choose the closest pair. 
      // Otherwise, we need to make a logical choice between d0 and d1.
      float deltaD = fabs(distances[0].m_d - distances[1].m_d);
      DHZObj chosenPair = distances[0];
      
      if (deltaD < 30) {
        float pt0 = distances[0].HPt();
        float pt1 = distances[1].HPt();
        int idx = 0;
        if (pt1 > pt0) idx = 1;
        chosenPair = distances[idx];
      }
      
      // Reconstruct the objects here for various uses.
      std::vector<JetObj> bjets;
      bjets.push_back(jets[chosenPair.m_zIdx0]);
      bjets.push_back(jets[chosenPair.m_zIdx1]);
      ZObj Z(bjets);
      
      std::vector<JetObj> cjets;
      cjets.push_back(jets[chosenPair.m_hIdx0]);
      cjets.push_back(jets[chosenPair.m_hIdx1]);
      HObj H(cjets);
      
      // Check our tagging
      bool btag1 = abs(bjets[0].m_flav) == 5; //chosenPair.Z_has_bjet0(desired_BvL);
      bool btag2 = abs(bjets[1].m_flav) == 5; //chosenPair.Z_has_bjet1(desired_BvL);
      
      //std::cout << "btag1: " << btag1 << ", btag2: " << btag2 << std::endl;
      if (btag1 && btag2) {
        //std::cout << "passed btag" << std::endl;
        bool ctag1 = abs(cjets[0].m_flav) == 4; //chosenPair.H_has_cjet0(desired_CvL, desired_CvB);
        bool ctag2 = abs(cjets[1].m_flav) == 4; // chosenPair.H_has_cjet1(desired_CvL, desired_CvB);
        if (ctag1 && ctag2) {
          //std::cout << "passed ctag" << std::endl;
        
          h_VH_MCjet_DHZ->FillVH(Z, H, evtW);
          h_reco_DHZ->Fill(gen_bs, bjets, gen_cs, cjets, evtW);
          if (Z.M() < 30 || H.M() < 30)
            h_reco_DHZ_under30->Fill(gen_bs, bjets, gen_cs, cjets, evtW);
        
        }//end-ctag
      
      }//end-btag
      
    }//end-DHZ-method
    
    
  }//end-VbbHcc-check
#endif

  /********************************************************
  * JET ANALYSIS/SELECTION - we want to make sure we only *
  * select jets that are  useful for our analysis.        *
  ********************************************************/

  // For the remainder of the selections, we want to only have jets
  // that meet our criteria. Let's make sure we have enough jets to
  // meet our requirements.
  std::vector<JetObj> analysis_jets;
  int nJets = 0;
  for (unsigned int i = 0; i < jets.size(); ++i) {

    h_jet_cutflow->Fill(0.5, genWeight); // all jets
    TLorentzVector vec = jets[i].m_lvec;

    // Check to make sure that the jet pass our pT cuts.
    // NOTE: we have an extra cut for the leading jet.
    if (vec.Pt() < CUTS.Get<float>("jet_pt")) continue;
    if (i == 0 and vec.Pt() < CUTS.Get<float>("jet_pt0")) continue;
    h_jet_cutflow->Fill(1.5, genWeight); 

    // Check to make sure our jets are in the eta region we want.
    if (fabs(vec.Eta()) > CUTS.Get<float>("jet_eta")) continue;
    h_jet_cutflow->Fill(2.5, genWeight); 

    // NEW: Make sure they pass the pileup ID requirement
    // if the jet pT < 50 GeV. NOTE: There will be different
    // values for 2016 vs. 2017-2018. In either case, we're
    // checking if it passes the medium WP.
    if (vec.Pt() < 50.0) {
#if defined(MC_2016) || defined(DATA_2016)
      if (jets[i].m_puid < 1) continue;
#endif 

#if defined(MC_2017) || defined(MC_2018) || defined(DATA_2017) || defined(DATA_2018)
      if (jets[i].m_puid < 4) continue;
#endif
    }

    // Check if the electrons & muons overlap with the jets.
    // If dR(jet, lepton) < 0.4, discard the jet.
    bool isElectron = jets[i].IsLepton(elecs);
    bool isMuon     = jets[i].IsLepton(muons);
    if (isElectron || isMuon) continue;
    h_jet_cutflow->Fill(3.5, genWeight); // passed iso cut

    // Add jet to our selected jet list
    jets[i].SetIdx(nJets); nJets++;
    analysis_jets.push_back(jets.at(i));

  }//end-for-jets

  // Fill the proper jet plots with the jets we're interested in.
  h_jets_selected->Fill(analysis_jets, evtW);
  h_jets_all->Fill(jets, evtW);

  //===========================================================================
  // Any analysis after this point should be using only our selected jets which
  // pass our selection criteria. There should be at least four of them.
  //===========================================================================
  //std::cout << "no. jets = " << analysis_jets.size() << std::endl;
  h_nJet->Fill(jets.size());
  h_nAnalysisJet->Fill(analysis_jets.size());
  if (analysis_jets.size() >= 4) {

    //std::cout << "we have four jets!!!" << std::endl;
    h_evt_tagOnly_cutflow->Fill(2.5, genWeight); // passed jet selection
    h_evt_algoFirst_cutflow->Fill(2.5, genWeight); 
    h_evt_tagFirst_cutflow->Fill(2.5, genWeight); 

    // For our cases where we make it this far,
    // check what our MET is and check the HT
    // (both versions - sum pT and sum "all").
    float HT = 0.;
    for (size_t i = 0; i < analysis_jets.size(); ++i)
      HT += analysis_jets[i].m_lvec.Pt();
    float HTmod = *(r->MET_pt);
    
    h_VH_select->FillMET(*(r->MET_pt), evtW);
    h_VH_select->FillHT(HT, HTmod, evtW);

    // In these effficiency calculations, we want to know what the flavors are
    // as a check against our analysis requirements, so we need to tag each 
    // jet. We will assume that B takes priority over C in terms of double
    // tagging.
    bool has_bjets = false, has_cjets = false;

    // Check if the highest two b-tagOnly pass our requirement
    std::vector<JetObj> copyjets = analysis_jets;
    std::vector<JetObj> copyjets2 = analysis_jets;
    std::sort(copyjets.begin(), copyjets.end(), JetObj::JetCompBtag());
    std::sort(copyjets2.begin(), copyjets2.end(), JetObj::JetCompBtag());
    std::vector<JetObj> bTagged { copyjets[0], copyjets[1] };
    copyjets.erase(copyjets.begin() + 1); copyjets.erase(copyjets.begin() + 0);
    copyjets2.erase(copyjets2.begin()+1); copyjets2.erase(copyjets2.begin()+0);
    has_bjets = are_bjets(bTagged, desired_BvL);

    // Check if the highest two c-tags pass our requirement
    std::sort(copyjets2.begin(), copyjets2.end(), JetObj::JetCompCtag());
    std::vector<JetObj> bTagged2 { copyjets[0], copyjets[1] };
    std::vector<JetObj> cTagged { copyjets2[0], copyjets2[1] };
    copyjets.erase(copyjets.begin() + 1); copyjets.erase(copyjets.begin() + 0);
    has_cjets = are_cjets(cTagged, desired_CvL, desired_CvB);
    bool has_more_bjets = are_bjets(bTagged2, desired_BvL);

    bool properly_tagged_4B = has_bjets && has_more_bjets;
    bool properly_tagged_3B = has_bjets && passes_btag(bTagged2[0], desired_BvL);
    bool properly_tagged_2b2c = has_bjets && has_cjets;

    // For each year, we want to make sure that we pass all the criteria we're
    // interested in. We've already got a variable to handle tagging. We need
    // variables (per year) for pT cuts.

#if defined(MC_2016) || defined(DATA_2016)
    std::vector<JetObj> jets2016 = analysis_jets;
    std::sort(jets2016.begin(), jets2016.end(), JetObj::JetCompPt());

    // Trigger #1 - make sure all jets have pT > 45 GeV
    bool pT_criteria_met_2016v1 = true;
    for (size_t i = 0; i < jets2016.size(); ++i){
      if (jets2016[i].Pt() <= 45) pT_criteria_met_2016v1 = false;
    }

    // Trigger #2 - make sure two jets have pT > 90 GeV and two have pT > 30 GeV
    bool pT_criteria_met_2016v2 = true;
    for (size_t i = 0; i < 2; ++i) {
      if (jets2016[i].Pt() <= 90) pT_criteria_met_2016v2 = false;
      if (jets2016[2+i].Pt() <= 30) pT_criteria_met_2016v2 = false;
    }
#endif

#if defined(MC_2017) || defined(DATA_2017) || defined(MC_2018) || defined(DATA_2018)
    std::vector<JetObj> jets2017_18 = analysis_jets;
    std::sort(jets2017_18.begin(), jets2017_18.end(), JetObj::JetCompPt());

    // Both triggers for 2017-18 require the following pT cuts:
    // 1. > 75 GeV, 2. > 60 GeV, 3. > 45 GeV, 4. > 40 GeV
    float cuts[4] = { 75.0, 60.0, 45.0, 40.0 };
    bool pT_criteria_met_2017_18 = true;
    for (size_t i = 0; i < 4; ++i) {
      if (jets2017_18[i].Pt() <= cuts[i]) pT_criteria_met_2017_18 = false;
    }

    // These are the additional versions we want to check for 2017-18.
    bool pT_criteria_met_2017_18_v2 = true;
    float cuts_v2[4] = { 103.0, 88.0, 75.0, 15.0 };
    for (size_t i = 0; i < 4; ++i) {
      if (jets2017_18[i].Pt() <= cuts_v2[i]) pT_criteria_met_2017_18_v2 = false;
    }

    bool pT_criteria_met_2017_18_v3 = true;
    float cuts_v3[4] = { 105.0, 88.0, 76.0, 15.0 };
    for (size_t i = 0; i < 4; ++i) {
      if (jets2017_18[i].Pt() <= cuts_v3[i]) pT_criteria_met_2017_18_v3 = false;
    }

    bool pT_criteria_met_2017_18_v4 = true;
    float cuts_v4[4] = { 111.0, 90.0, 80.0, 15.0 };
    for (size_t i = 0; i < 4; ++i) {
      if (jets2017_18[i].Pt() <= cuts_v4[i]) pT_criteria_met_2017_18_v4 = false;
    }

    bool pT_criteria_met_2017_18_v5 = true;
    float cuts_v5[4] = { 98.0, 83.0, 71.0, 15.0 };
    for (size_t i = 0; i < 4; ++i) {
      if (jets2017_18[i].Pt() <= cuts_v5[i]) pT_criteria_met_2017_18_v5 = false;
    }

#endif

    // ========================================================================
    // Trigger Efficiency
    // ========================================================================

    bool is_ZH = true;  // Use this boolean for cases where we want to ignore
                        // our reference trigger & just use the base events.
    bool noTrigCheck = false; // Use this if we want to check how the events
                             // are selected without the trigger requirement.

    // Reference Trigger - HLT_IsoMu24 
    // Probe Triggers:
    // 2016 - HLT_QuadJet45_TripleBTagCSV_p087 OR
    //        HLT_DoubleJet90_Double30_TripleBTagCSV_p087
    // 2017 - HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0
    // 2018 - HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5

    // We need to make sure for sake of the sample that we're using that
    // we have at least one isolated muon. Luckily, in our code above, 
    // we check for muons, so all we have to do is make sure our list
    // contains at least one.
    if (muons.size() >= 1 || is_ZH) {

      // Make sure events pass our reference frame.
      if (*(r->HLT_IsoMu24) || is_ZH) {

        // We want the option to see if HT is being used as just sum(jet pT)
        // or if it's ALL objects. We need an HT mod that is the sum of the
        // MET and the muon pT
        float HTmod = *(r->MET_pt);
        for (size_t i = 0; i < muons.size(); ++i)
          HTmod += muons[i].m_lvec.Pt();

        float HT = 0.0;
        for (size_t i = 0; i < analysis_jets.size(); ++i)
          HT += analysis_jets[i].Pt();

        // Add the events to the reference plots (true = isReference)
        h_trig_2016_QuadJet_TripleTag->Fill(analysis_jets, true, HTmod, evtW);
        h_trig_2016_QuadJet_DoubleTag->Fill(analysis_jets, true, HTmod, evtW);
        h_trig_2016_DoubleJet_TripleTag->Fill(analysis_jets, true, HTmod, evtW);
        h_trig_2016_DoubleJet_DoubleTag->Fill(analysis_jets, true, HTmod, evtW);
        h_trig_2017_QuadJet_TripleTag->Fill(analysis_jets, true, HTmod, evtW);
        h_trig_2017_QuadJet_noTag->Fill(analysis_jets, true, HTmod, evtW);
        h_trig_2018_QuadJet_TripleTag->Fill(analysis_jets, true, HTmod, evtW);
        h_trig_2018_QuadJet_noTag->Fill(analysis_jets, true, HTmod, evtW);

        if (properly_tagged_4B) {
          h_trig_2016_QuadJet_TripleTag_tagged->Fill(analysis_jets, true, HTmod, evtW);
          h_trig_2016_QuadJet_DoubleTag_tagged->Fill(analysis_jets, true, HTmod, evtW);
          h_trig_2016_DoubleJet_TripleTag_tagged->Fill(analysis_jets, true, HTmod, evtW);
          h_trig_2016_DoubleJet_DoubleTag_tagged->Fill(analysis_jets, true, HTmod, evtW);
          h_trig_2017_QuadJet_TripleTag_tagged->Fill(analysis_jets, true, HTmod, evtW);
          h_trig_2017_QuadJet_noTag_tagged->Fill(analysis_jets, true, HTmod, evtW);
          h_trig_2018_QuadJet_TripleTag_tagged->Fill(analysis_jets, true, HTmod, evtW);
          h_trig_2018_QuadJet_noTag_tagged->Fill(analysis_jets, true, HTmod, evtW);
        }

#if defined(MC_2016) || defined(DATA_2016)

        /**********************************************************
        * First version we wanna try is when we pass proper
        * tagging and pT requirements associated with the triggers.
        **********************************************************/

        // 2016 - QuadJet (4B tagging)
        if (properly_tagged_4B && pT_criteria_met_2016v1) {

          // Fill reference 
          h_trig_2016_QuadJet_TripleTag_ideal->Fill(analysis_jets, true, HTmod, evtW);
          h_trig_2016_QuadJet_DoubleTag_ideal->Fill(analysis_jets, true, HTmod, evtW);        

          // Fill reference + probe   
          if (*(r->HLT_QuadJet45_TripleBTagCSV_p087) || noTrigCheck)
            h_trig_2016_QuadJet_TripleTag_ideal->Fill(analysis_jets, false, HTmod, evtW);
           if (*(r->HLT_QuadJet45_DoubleBTagCSV_p087) || noTrigCheck)
            h_trig_2016_QuadJet_DoubleTag_ideal->Fill(analysis_jets, false, HTmod, evtW);       
        }

        // 2016 - QuadJet (3B tagging)
        if (properly_tagged_3B && pT_criteria_met_2016v1) {

          // Fill reference
          h_trig_2016_QuadJet_TripleTag_3B->Fill(analysis_jets, true, HTmod, evtW);
          
          // Fill reference + probe   
          if (*(r->HLT_QuadJet45_TripleBTagCSV_p087) || noTrigCheck)
            h_trig_2016_QuadJet_TripleTag_3B->Fill(analysis_jets, false, HTmod, evtW);
        }

        // 2016 - QuadJet (2b2c tagging)
        if (properly_tagged_2b2c && pT_criteria_met_2016v1) {

          // Fill reference
          h_trig_2016_QuadJet_TripleTag_2b2c->Fill(analysis_jets, true, HTmod, evtW);
          
          // Fill reference + probe   
          if (*(r->HLT_QuadJet45_TripleBTagCSV_p087) || noTrigCheck)
            h_trig_2016_QuadJet_TripleTag_2b2c->Fill(analysis_jets, false, HTmod, evtW);

        }

        
        // 2016 - DoubleJet (4B tagging)
        if (properly_tagged_4B && pT_criteria_met_2016v2) {

          // Pass our reference            
          h_trig_2016_DoubleJet_TripleTag_ideal->Fill(analysis_jets, true, HTmod, evtW);
          
          // Pass our reference + probe
          if (*(r->HLT_DoubleJet90_Double30_TripleBTagCSV_p087) || noTrigCheck)
            h_trig_2016_DoubleJet_TripleTag_ideal->Fill(analysis_jets, false, HTmod, evtW); 
        }

        // 2016 - DoubleJet (3B tagging)
        if (properly_tagged_3B && pT_criteria_met_2016v2) {

          // Pass our reference            
          h_trig_2016_DoubleJet_TripleTag_3B->Fill(analysis_jets, true, HTmod, evtW);

          // Pass our reference + probe
          if (*(r->HLT_DoubleJet90_Double30_TripleBTagCSV_p087) || noTrigCheck)
            h_trig_2016_DoubleJet_TripleTag_3B->Fill(analysis_jets, false, HTmod, evtW);
        }
 
        // 2016 - DoubleJet (2b2c tagging)
        if (properly_tagged_2b2c && pT_criteria_met_2016v2) {

          // Pass our reference            
          h_trig_2016_DoubleJet_TripleTag_2b2c->Fill(analysis_jets, true, HTmod, evtW);
          
          // Pass our reference + probe
          if (*(r->HLT_DoubleJet90_Double30_TripleBTagCSV_p087) || noTrigCheck)
            h_trig_2016_DoubleJet_TripleTag_2b2c->Fill(analysis_jets, false, HTmod, evtW);
        }

        /************************************************************
        * We wanna look at the versions with none of the requirements
        * and then just the tagging requirement.
        ************************************************************/

        // 2016 - QuadJet (Triple Tag vs Double Tag)
        if (*(r->HLT_QuadJet45_TripleBTagCSV_p087)) {
          h_trig_2016_QuadJet_TripleTag->Fill(analysis_jets, false, HTmod, evtW);
          if (properly_tagged_4B) {
            h_trig_2016_QuadJet_TripleTag_tagged->Fill(analysis_jets,
              false, HTmod, evtW);
          }
        }
        if (*(r->HLT_QuadJet45_DoubleBTagCSV_p087)) {
          h_trig_2016_QuadJet_DoubleTag->Fill(analysis_jets, false, HTmod, evtW);
          if (properly_tagged_4B) {
            h_trig_2016_QuadJet_DoubleTag_tagged->Fill(analysis_jets,
              false, HTmod, evtW);
          }
        }

        // 2016 - DoubleJet (Triple Tag vs Double Tag)
        if (*(r->HLT_DoubleJet90_Double30_TripleBTagCSV_p087)) {
          h_trig_2016_DoubleJet_TripleTag->Fill(analysis_jets, false, HTmod, evtW);
          if (properly_tagged_4B) {
            h_trig_2016_DoubleJet_TripleTag_tagged->Fill(analysis_jets,
              false, HTmod, evtW);
          }
        }
        if (*(r->HLT_DoubleJet90_Double30_DoubleBTagCSV_p087)) {
          h_trig_2016_DoubleJet_DoubleTag->Fill(analysis_jets, false, HTmod, evtW);
          if (properly_tagged_4B) {
            h_trig_2016_DoubleJet_DoubleTag_tagged->Fill(analysis_jets,
              false, HTmod, evtW);
          }
        }

#endif

#if defined(MC_2017) || defined(DATA_2017)

        /**********************************************************
        * First version we wanna try is when we pass proper
        * tagging and pT requirements associated with the triggers.
        **********************************************************/
        
        // 4B tagging requirement
        if (properly_tagged_4B && pT_criteria_met_2017_18 && HT > 300) {

          // Pass our reference trigger.
          h_trig_2017_QuadJet_TripleTag_ideal->Fill(analysis_jets, true, HTmod, evtW);
          h_trig_2017_QuadJet_noTag_ideal->Fill(analysis_jets, true, HTmod, evtW);

          // Pass our probe triggers.
          #if defined(MC_2017) || !defined(DATA_2017B)
          if (*(r->HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0) || noTrigCheck)
            h_trig_2017_QuadJet_TripleTag_ideal->Fill(analysis_jets, false, HTmod, evtW);
          
          if (*(r->HLT_PFHT300PT30_QuadPFJet_75_60_45_40) || noTrigCheck)
            h_trig_2017_QuadJet_noTag_ideal->Fill(analysis_jets, false, HTmod, evtW);     
          #endif

          #if defined(DATA_2017B)
          if (*(r->HLT_HT300PT30_QuadJet_75_60_45_40_TripeCSV_p07) || noTrigCheck)
            h_trig_2017_QuadJet_TripleTag_RunB_ideal->Fill(analysis_jets, false, HTmod, evtW);
          #endif
        }

        // 3B tagging requirement
        if (properly_tagged_3B && pT_criteria_met_2017_18 && HT > 300) {

          // Pass our reference trigger.
          h_trig_2017_QuadJet_TripleTag_3B->Fill(analysis_jets, true, HTmod, evtW);

          // Pass our probe triggers.
          #if defined(MC_2017) || !defined(DATA_2017B)
          if (*(r->HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0) || noTrigCheck)
            h_trig_2017_QuadJet_TripleTag_3B->Fill(analysis_jets, false, HTmod, evtW);
          #endif

          #if defined(DATA_2017B)
          if (*(r->HLT_HT300PT30_QuadJet_75_60_45_40_TripeCSV_p07) || noTrigCheck)
            h_trig_2017_QuadJet_TripleTag_RunB_3B->Fill(analysis_jets, false, HTmod, evtW);
          #endif 
        }

        // 2b2c tagging requirement
        if (properly_tagged_2b2c && pT_criteria_met_2017_18 && HT > 300) {

          // Pass our reference trigger.
          h_trig_2017_QuadJet_TripleTag_2b2c->Fill(analysis_jets, true, HTmod, evtW);

          // Pass our probe triggers.
          #if defined(MC_2017) || !defined(DATA_2017B)
          if (*(r->HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0) || noTrigCheck)
            h_trig_2017_QuadJet_TripleTag_2b2c->Fill(analysis_jets, false, HTmod, evtW);
          #endif

          #if defined(DATA_2017B)
          if (*(r->HLT_HT300PT30_QuadJet_75_60_45_40_TripeCSV_p07) || noTrigCheck)
            h_trig_2017_QuadJet_TripleTag_RunB_2b2c->Fill(analysis_jets, false, HTmod, evtW);
          #endif
        }

        /************************************************************
        * We wanna look at the versions with none of the requirements
        * and then just the tagging requirement.
        ************************************************************/
        
        // 2017 (Triple Tag vs No Tag)
        /*if (*(r->HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0)) {
          h_trig_2017_QuadJet_TripleTag->Fill(analysis_jets, false, HTmod, evtW);
          if (properly_tagged_4B) {
            h_trig_2017_QuadJet_TripleTag_tagged->Fill(analysis_jets,
              false, HTmod, evtW);
          }
        }
        if (*(r->HLT_PFHT300PT30_QuadPFJet_75_60_45_40)) {
          h_trig_2017_QuadJet_noTag->Fill(analysis_jets, false, HTmod, evtW);
          if (properly_tagged_4B) {
            h_trig_2017_QuadJet_noTag_tagged->Fill(analysis_jets,
              false, HTmod, evtW);
          }
        }*/  
#endif

#if defined(MC_2018) || defined(DATA_2018)

        /**********************************************************
        * First version we wanna try is when we pass proper
        * tagging and pT requirements associated with the triggers.
        **********************************************************/

        // 4B tagging version
        if (properly_tagged_4B && pT_criteria_met_2017_18 && HT > 330) {

          // Pass our reference trigger.
          h_trig_2018_QuadJet_TripleTag_ideal->Fill(analysis_jets, true, HTmod, evtW);
          h_trig_2018_QuadJet_noTag_ideal->Fill(analysis_jets, true, HTmod, evtW);

          // Pass our probe triggers.
          if (*(r->HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5) || noTrigCheck)
            h_trig_2018_QuadJet_TripleTag_ideal->Fill(analysis_jets, false, HTmod, evtW);
          if (*(r->HLT_PFHT330PT30_QuadPFJet_75_60_45_40) || noTrigCheck)
            h_trig_2018_QuadJet_noTag_ideal->Fill(analysis_jets, false, HTmod, evtW);
        }

        // 3B tagging version
        if (properly_tagged_3B && pT_criteria_met_2017_18 && HT > 330) {

          // Pass our reference trigger.
          h_trig_2018_QuadJet_TripleTag_3B->Fill(analysis_jets, true, HTmod, evtW);

          // Pass our probe triggers.
          if (*(r->HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5) || noTrigCheck)
            h_trig_2018_QuadJet_TripleTag_3B->Fill(analysis_jets, false, HTmod, evtW);
        }

        // 2b2c tagging version
        if (properly_tagged_2b2c && pT_criteria_met_2017_18 && HT > 330) {

          // Pass our reference trigger.
          h_trig_2018_QuadJet_TripleTag_2b2c->Fill(analysis_jets, true, HTmod, evtW);

          // Pass our probe triggers.
          if (*(r->HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5) || noTrigCheck)
            h_trig_2018_QuadJet_TripleTag_2b2c->Fill(analysis_jets, false, HTmod, evtW);
        }

        /***********************************************************
        * We wanna look at the version with none of the requirements
        * and then just the tagging requirement.
        ***********************************************************/        

        // 2018 (Triple Tag vs No Tag)
        /*if (*(r->HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5)) {
          h_trig_2018_QuadJet_TripleTag->Fill(analysis_jets, false, HTmod, evtW);
          if (properly_tagged_4B) {
            h_trig_2018_QuadJet_TripleTag_tagged->Fill(analysis_jets, 
              false, HTmod, evtW);
         two-factor authentication   }
        }
        if (*(r->HLT_PFHT330PT30_QuadPFJet_75_60_45_40)) {
          h_trig_2018_QuadJet_noTag->Fill(analysis_jets, false, HTmod, evtW);
          if (properly_tagged_4B) {
            h_trig_2018_QuadJet_noTag_tagged->Fill(analysis_jets,
              false, HTmod, evtW);
          }
        }*/
#endif

      } // end ref-trigger
    }//end # muons

    /***************************************************************************
    * Check what we get for H and Z candidates by just using the leading pT    *
    * jets after any cuts (but before ANY tagging.                             *
    ***************************************************************************/
     
    // Make a copy of the jets and select/sort by pT.
    std::vector<JetObj> afterjets; afterjets = analysis_jets;
    std::sort(afterjets.begin(), afterjets.end(), JetObj::JetCompPt());
    
    // Reconstruct the Higgs boson with the highest two.
    std::vector<JetObj> acjets { afterjets[0], afterjets[1] };
    HObj Hafter(acjets);
    
    // Reconstruct the Z boson with the next highest two.
    std::vector<JetObj> abjets { afterjets[2], afterjets[3] };
    ZObj Zafter(abjets);
    
    // Fill the proper plots
    h_VH_seljet->FillVH(Zafter, Hafter, evtW); 

    /**************************************************************************
    * MISTAG RATE ANALYSIS                                                    *
    **************************************************************************/
    
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
    std::vector<JetObj> jetlistCopy; jetlistCopy = analysis_jets;
    std::sort(jetlistCopy.begin(), jetlistCopy.end(), JetObj::JetCompPt());
    
    bool found_leading = false;
    for (size_t i = 0; i < jetlistCopy.size(); ++i) {

      // check the leading light jet (if we haven't yet)
      if (!found_leading) {
        JetObj leadingJet = jetlistCopy[i];
        if (leadingJet.m_flav < 4) {
          found_leading = true;
          h_mistag_leading->Fill(0.5, evtW);
          if (leadingJet.m_deepCSV > desired_BvL)
            h_mistag_leading->Fill(1.5, evtW);
          else h_mistag_leading->Fill(2.5, evtW);
        }        
      }//end-leading

      // Check all the light jets we have
      if (jetlistCopy[i].m_flav < 4) {
        h_mistag_all->Fill(0.5, evtW);
        if (jetlistCopy[i].m_deepCSV > desired_BvL)
          h_mistag_all->Fill(1.5, evtW);
        else h_mistag_all->Fill(2.5, evtW);
      }
    }
#endif

    /**************************************************************************
    * JET ANALYSIS - 4B vs 2b2c                                               *
    **************************************************************************/
    
    // Go through the jets and get the distributions
    // for when we have 4b jets in the final state.
    std::vector<JetObj> jetsB; jetsB = analysis_jets;
    std::sort(jetsB.begin(), jetsB.end(), JetObj::JetCompBtag());
    std::vector<JetObj> bjetsB { jetsB[0], jetsB[1], jetsB[2], jetsB[3] };
    if (are_bjets(bjetsB, desired_BvL)) {
      h_jets_4b->Fill(bjetsB, evtW);
    }

    // Go through the jets and get the distirbutions
    // for when we have 2b 2c jets in the final state.
    std::vector<JetObj> jetsBC; jetsBC = analysis_jets;
    std::sort(jetsBC.begin(), jetsBC.end(), JetObj::JetCompBtag());
    std::vector<JetObj> bjetsBC { jetsBC[0], jetsBC[1] };
    jetsBC.erase(jetsBC.begin() + 1); jetsBC.erase(jetsBC.begin() + 0);
    if (are_bjets(bjetsBC, desired_BvL)) {

      std::sort(jetsBC.begin(), jetsBC.end(), JetObj::JetCompCtag());
      std::vector<JetObj> cjetsBC { jetsBC[0], jetsBC[1] };
      jetsBC.erase(jetsBC.begin() + 1); jetsBC.erase(jetsBC.begin() + 0);
      if (are_cjets(cjetsBC, desired_CvL, desired_CvB)) {
        h_jets_2b2c->Fill(bjetsBC, evtW);
        h_jets_2b2c->Fill(cjetsBC, evtW);
      }
 
    }

    // We want to only keep going with our selections if we trigger properly.
    // We need to check each year of the analysis properly.
    bool trigger = false;

#if defined(MC_2016) || defined(DATA_2016)
    trigger = (*(r->HLT_QuadJet45_TripleBTagCSV_p087) || 
      *(r->HLT_DoubleJet90_Double30_TripleBTagCSV_p087));
    //std::cout << "Checking the 2016 trigger!!!" << std::endl;
#endif

#if defined(MC_2017) || defined(DATA_2017)
  #if !defined(DATA_2017B)
    trigger = *(r->HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0);
  #endif
  #if defined(DATA_2017B)
    trigger = *(r->HLT_HT300PT30_QuadJet_75_60_45_40_TripeCSV_p07);
  #endif
  //std::cout << "Checking the 2017 trigger!!!" << std::endl;
#endif

#if defined(MC_2018) || defined(DATA_2018)
    trigger = *(r->HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5);
    //std::cout << "Checking the 2018 trigger!!!" << std::endl;
#endif

    // Do not continue on if we haven't triggered properly.
    if (!trigger) return;

    //std::cout << "Passed our triggers!!!!" << std::endl;
    h_evt_tagOnly_cutflow->Fill(3.5, genWeight);
    h_evt_algoFirst_cutflow->Fill(3.5, genWeight);
    h_evt_tagFirst_cutflow->Fill(3.5, genWeight);

    /**************************************************************************
    * Selection Method #1 - TAGGING ONLY                                      *
    **************************************************************************/

    bool force_pass_tagging = false;
    
    // Make an appropriate copy of the jets to use for this analysis.
    std::vector<JetObj> jets2; jets2 = analysis_jets;

    // Select two jets with the largest btag values and then
    // check against our working point of interest. (NOTE: when
    // we sort by btag (descending order), the largest two will
    // automatically be in the first two indices.)
    std::sort(jets2.begin(), jets2.end(), JetObj::JetCompBtag());
    std::vector<JetObj> bjets2 { jets2[0], jets2[1] };
    std::vector<JetObj> bjets22 { jets2[0], jets2[1] };
    std::vector<JetObj> bjets23 { jets2[0], jets2[1] };
    jets2.erase(jets2.begin() + 1); jets2.erase(jets2.begin() + 0);
    
    bool btag1 = passes_btag(bjets2[0], desired_BvL);
    bool btag2 = passes_btag(bjets2[1], desired_BvL);

    if (btag1 || force_pass_tagging) {

      h_evt_tagOnly_cutflow->Fill(4.5, genWeight); // pass b-cut #1
      if (btag2 || force_pass_tagging) {

        // Since we passed the cut, adjust our jets with the proper
        // JEC for b-tagged jets and reconstruct the Z boson
        h_evt_tagOnly_cutflow->Fill(5.5, genWeight); // pass b-cut #2

        bjets2[0].ApplyRegression(5); bjets2[1].ApplyRegression(5); // Full JEC version
        ZObj Z2(bjets2);
        bjets22[0].ApplyRegression(5, false); bjets22[1].ApplyRegression(5, false); // No mass correction
        ZObj Z22(bjets22);
        ZObj Z23(bjets23); // no JEC
      
        // Select two jets with the largest ctag values and then
        // check against our working point of interest.
        std::sort(jets2.begin(), jets2.end(), JetObj::JetCompCtag());
        std::vector<JetObj> cjets2 { jets2[0], jets2[1] };
        std::vector<JetObj> cjets22 { jets2[0], jets2[1] };
        std::vector<JetObj> cjets23 { jets2[0], jets2[1] };
        std::vector<JetObj> cjets_2b1c { jets2[0], jets2[1] };
        jets2.erase(jets2.begin() + 1); jets2.erase(jets2.begin() + 0);
        
        bool ctag1 = passes_ctag(cjets2[0], desired_CvL, desired_CvB);
        bool ctag2 = passes_ctag(cjets2[1], desired_CvL, desired_CvB);
        
        // Uncomment the following lines if you want to check 4b tagging
        // instead of 2b2c.
        //ctag1 = passes_btag(cjets2[0], desired_BvL);
        //ctag2 = passes_btag(cjets2[1], desired_BvL);

        if (ctag1 || force_pass_tagging) {

          h_evt_tagOnly_cutflow->Fill(6.5, genWeight); // pass c-cut #1
          if (ctag2 || force_pass_tagging) {
          
            // Since we passed the cuts, adjust our jets with the proper
            // JEC for c-tagged jets and reconstruct the Higgs boson
            h_evt_tagOnly_cutflow->Fill(7.5, genWeight); // pass c-cut #2

            cjets2[0].ApplyRegression(4); cjets2[1].ApplyRegression(4); // Full JEC version
            HObj H2(cjets2);
            cjets22[0].ApplyRegression(4,false); cjets22[1].ApplyRegression(4,false); // No mass correction
            HObj H22(cjets22);
            HObj H23(cjets23); // no JEC

            h_VH_tagOnly->FillVH(Z2, H2, evtW);
            h_VH_tagOnly_noMassCorr->FillVH(Z22, H22, evtW);
            h_VH_tagOnly_noJEC->FillVH(Z23, H23, evtW);         

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
            h_VH_tagOnly->FillGluCheck(Z2, evtW);
            
            // If we have a proper VbbHcc event and we've found the 
            // proper MC truth objects, check the efficiency.
            if (is_VbbHcc_event) {
              h_eff_tagOnly->Fill(Z2, H2, gen_bs, gen_cs, true, evtW);
            }
#endif 

          }//end-c-cut-2
        }//end-c-cut-1
        
        // Check the 2b1c version
        if (ctag1 || ctag2) {
        
          cjets_2b1c[0].ApplyRegression(4);
          cjets_2b1c[1].ApplyRegression(4);
          HObj H2(cjets_2b1c);
          
          h_VH_tagOnly_2b1c->FillVH(Z2, H2, evtW);
        
        }//end-2b1c
        
      }//end-b-cut-2
    }//end-b-cut-1
   
    /**************************************************************************
    * Selection Method #2 -  MASS-MATCHING PRIORITIZED                        *
    **************************************************************************/

    // Make an appropriate copy of the jets to use for this analysis.
    std::vector<JetObj> jets3; jets3 = analysis_jets;
    std::sort(jets3.begin(), jets3.end(), JetObj::JetCompPt());
    
    // Create the objects for the distance calculations and make sure we sort
    // them in ascending order. (All necessary calculations for algorithms
    // are handled within the DObj class.)
    DHZObj d0(jets3, 0, 1, 2, 3);          // H(0,1) ; Z(2,3)
    DHZObj d1(jets3, 0, 2, 1, 3);          // H(0,2) ; Z(1,3)
    DHZObj d2(jets3, 0, 3, 1, 2);          // H(0,3) ; Z(1,2)
    std::vector<DHZObj> distances {d0, d1, d2};
    std::sort(distances.begin(), distances.end());

    // Determine the distance between the closest two. Then, based on this
    // distance, determine which pair we want to use. If we are outside the
    // resolution window (30 GeV), we can just choose the closest pair. 
    // Otherwise, we need to make a logical choice between d0 and d1.
    float deltaD = fabs(distances[0].m_d - distances[1].m_d);
    DHZObj chosenPair = distances[0];

    if (deltaD < 30) {
      float pt0 = distances[0].HPt();
      float pt1 = distances[1].HPt();
      int idx = 0;
      if (pt1 > pt0) idx = 1;
      chosenPair = distances[idx];
    }

    // Reconstruct the objects here for various uses.
    h_VH_algoFirst->FillAlgo(distances, evtW);
    std::vector<JetObj> bjets3, bjets32, bjets33;

    bjets3.push_back(jets3[chosenPair.m_zIdx0]); // full JEC
    bjets3.push_back(jets3[chosenPair.m_zIdx1]);
    bjets3[0].ApplyRegression(5); bjets3[1].ApplyRegression(5);
    ZObj Z3(bjets3);

    bjets32.push_back(jets3[chosenPair.m_zIdx0]); // no mass correction
    bjets32.push_back(jets3[chosenPair.m_zIdx1]); 
    bjets32[0].ApplyRegression(5,false); bjets32[1].ApplyRegression(5,false);
    ZObj Z32(bjets32);

    bjets33.push_back(jets3[chosenPair.m_zIdx0]); // no JEC
    bjets33.push_back(jets3[chosenPair.m_zIdx1]);
    ZObj Z33(bjets33);

    std::vector<JetObj> cjets3, cjets32, cjets33; 
    cjets3.push_back(jets3[chosenPair.m_hIdx0]); // full JEC
    cjets3.push_back(jets3[chosenPair.m_hIdx1]);
    cjets3[0].ApplyRegression(4); cjets3[1].ApplyRegression(4);
    HObj H3(cjets3);

    cjets32.push_back(jets3[chosenPair.m_hIdx0]); // no mass correction
    cjets32.push_back(jets3[chosenPair.m_hIdx1]);
    cjets32[0].ApplyRegression(4, false); cjets32[1].ApplyRegression(4,false);
    HObj H32(cjets32);
    
    cjets33.push_back(jets3[chosenPair.m_hIdx0]); // no JEC
    cjets33.push_back(jets3[chosenPair.m_hIdx1]);
    HObj H33(cjets33);

    bool btag31 = chosenPair.Z_has_bjet0(desired_BvL);
    bool btag32 = chosenPair.Z_has_bjet1(desired_BvL);

    // Now check our tagging requirements and other cuts.
    if (btag31 || force_pass_tagging) {  
      h_evt_algoFirst_cutflow->Fill(4.5, genWeight); // pass b-tag #1
      if (btag32 || force_pass_tagging) {
        h_evt_algoFirst_cutflow->Fill(5.5, genWeight); // pass b-tag #2

        bool ctag1 = chosenPair.H_has_cjet0(desired_CvL, desired_CvB);
        bool ctag2 = chosenPair.H_has_cjet1(desired_CvL, desired_CvB);

        // Uncomment the following two lines if you want to
        // check 4b instead of 2b2c tagging.
        //ctag1 = chosenPair.H_has_bjet0(desired_BvL);
        //ctag2 = chosenPair.H_has_bjet1(desired_BvL);
 
        if (ctag1 || force_pass_tagging) {
          h_evt_algoFirst_cutflow->Fill(6.5, genWeight); // pass c-tag #1
          if (ctag2 || force_pass_tagging) {
            h_evt_algoFirst_cutflow->Fill(7.5, genWeight); // pass c-tag #2
            
            // Fill our histograms appropriately.
            h_VH_algoFirst->FillVH(Z3, H3, evtW);
            h_VH_algoFirst_noMassCorr->FillVH(Z32, H32, evtW);
            h_VH_algoFirst_noJEC->FillVH(Z33, H33, evtW);

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
            h_VH_algoFirst->FillGluCheck(Z3, evtW);

            // If we have a proper VbbHcc event and we've found the 
            // proper MC truth objects, check the efficiency.
            if (is_VbbHcc_event) {
              h_eff_algoFirst->Fill(Z3, H3, gen_bs, gen_cs, true, evtW);
            }
#endif

          }//end-c-cut-2
        }//end-c-cut-1
        
        // Check the 2b1c version
        if (ctag1 || ctag2) {
          h_VH_algoFirst_2b1c->FillVH(Z3, H3, evtW);
        }//end-2b1c
        
      }//end-b-cut-2
    }//end-b-cut-1 

    /**************************************************************************
    * Selection Method #3 - TAGGING PRIORITIZED                               *
    **************************************************************************/
  
    // Make an appropriate copy of the jets to use for this analysis.
    std::vector<JetObj> jets4; jets4 = analysis_jets;
    
    // Check our jets and push them to the appropriate lists by 
    // checking if the tagging scores meet our desired cuts.
    std::vector<std::pair<int,float> > jets_idx_BvL;
    std::vector<std::pair<int,float> > jets_idx_CvL;
    std::vector<std::pair<int,float> > jets_idx_CvL_2b1c;
    bool already_tagged_c = false;
    
    // For each jet...
    for (size_t i = 0; i < jets4.size(); ++i) {
      
      float csv = jets4[i].m_deepCSV;
      float cvl = jets4[i].m_deepCvL;

      bool btag1 = passes_btag(jets4[i], desired_BvL);
      bool ctag1 = passes_ctag(jets4[i], desired_CvL, desired_CvB); 
      //ctag1 = passes_btag(jets4[i], desired_BvL);

      if (btag1 || force_pass_tagging) {
        std::pair<int,float> pair0(i,csv);
        jets_idx_BvL.push_back(pair0); 
      }
      if (ctag1 || force_pass_tagging){
        std::pair<int,float> pair1(i,cvl);
        jets_idx_CvL.push_back(pair1);
        already_tagged_c = true;
      }
      
      if (ctag1 || already_tagged_c) {
        std::pair<int,float> pair2(i,cvl);
        jets_idx_CvL_2b1c.push_back(pair2);
      }
    }

    // Sort the jets so that the largest tagging-scores
    // are at the top of our lists.
    std::sort(jets_idx_BvL.begin(), jets_idx_BvL.end(), sort_by_second_descend);
    std::sort(jets_idx_CvL.begin(), jets_idx_CvL.end(), sort_by_second_descend);
    std::sort(jets_idx_CvL_2b1c.begin(), jets_idx_CvL_2b1c.end(), sort_by_second_descend);

    // Create the vectors containing of b- and c-jets indices
    std::vector<int> bIndices;
    std::vector<int> cIndices;
    std::vector<int> cIndices_2b1c;

    for (auto p : jets_idx_BvL){ bIndices.push_back(p.first); }
    for (auto p : jets_idx_CvL){ cIndices.push_back(p.first); }
    for (auto p : jets_idx_CvL_2b1c){ cIndices_2b1c.push_back(p.first); }
 
    // Find appropriate combinations of jets.
    std::vector<std::vector<int>> combos = find_valid_combos(bIndices, cIndices);
    h_nCombos->Fill(combos.size(), evtW);
    
    // If there are possible combos, let's check them.
    if (combos.size() > 0) {

      h_evt_tagFirst_cutflow->Fill(4.5, genWeight); // passed tags

      // From the combos we've found, run the DHZ algorithm and
      // determine which is the best combination.
      std::vector<DHZObj> DHZ_values;
      for (size_t i = 0; i < combos.size(); ++i) {

        std::vector<int> idxs = combos[i];
        DHZObj D(jets4, idxs[2], idxs[3], idxs[0], idxs[1]);
        DHZ_values.push_back(D);
      }

      // Sort the values by their distance. Then, determine the distance
      // between the closest two. Based on this distance, determine which
      // pair we want to use. (NOTE: we might only have one combination.)
      std::sort(DHZ_values.begin(), DHZ_values.end());
      DHZObj chosenPair = DHZ_values[0];

      if (DHZ_values.size() >= 2) {
        float deltaD = fabs(DHZ_values[0].m_d - DHZ_values[1].m_d);
        if (deltaD < 30) {
          float pt0 = DHZ_values[0].HPt();
          float pt1 = DHZ_values[1].HPt();
          int idx = 0;
          if (pt1 > pt0) idx = 1;
          chosenPair = DHZ_values[idx];
        }
      }

      // Reconstruct the objects for use.
      h_VH_tagFirst->FillAlgo(DHZ_values, evtW);
      std::vector<JetObj> bjets4, bjets42, bjets43;

      bjets4.push_back(jets4[chosenPair.m_zIdx0]); // full JEC
      bjets4.push_back(jets4[chosenPair.m_zIdx1]);
      bjets4[0].ApplyRegression(5); bjets4[1].ApplyRegression(5);
      ZObj Z4(bjets4);

      bjets42.push_back(jets4[chosenPair.m_zIdx0]); // no mass correction
      bjets42.push_back(jets4[chosenPair.m_zIdx1]);
      bjets42[0].ApplyRegression(5,false); bjets42[1].ApplyRegression(5,false);
      ZObj Z42(bjets42);

      bjets43.push_back(jets4[chosenPair.m_zIdx0]); // no JEC
      bjets43.push_back(jets4[chosenPair.m_zIdx1]);
      ZObj Z43(bjets43);

      std::vector<JetObj> cjets4, cjets42, cjets43;
      cjets4.push_back(jets4[chosenPair.m_hIdx0]); // full JEC
      cjets4.push_back(jets4[chosenPair.m_hIdx1]);
      HObj H4(cjets4);

      cjets42.push_back(jets4[chosenPair.m_hIdx0]); // no mass correction
      cjets42.push_back(jets4[chosenPair.m_hIdx1]);
      HObj H42(cjets42);

      cjets43.push_back(jets4[chosenPair.m_hIdx0]); // no JEC
      cjets43.push_back(jets4[chosenPair.m_hIdx1]);
      HObj H43(cjets43);

      // Fill our histograms appropriately.
      h_VH_tagFirst->FillVH(Z4, H4, evtW);
      h_VH_tagFirst_noMassCorr->FillVH(Z42, H42, evtW);
      h_VH_tagFirst_noJEC->FillVH(Z43, H43, evtW);

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
        h_VH_tagFirst->FillGluCheck(Z4, evtW);

        // If we have a proper VbbHcc event and we've found the 
        // proper MC truth objects, check the efficiency.
        if (is_VbbHcc_event) {
          h_eff_tagFirst->Fill(Z4, H4, gen_bs, gen_cs, true, evtW);
        }
#endif

    }//end-tagging
    
    // Now, do the 2b1c version
    std::vector<std::vector<int>> combos_2b1c = find_valid_combos(bIndices, cIndices_2b1c);
    if (combos_2b1c.size() > 0) {
    
      // From the combos we've found, run the DHZ algorithm and
      // determine which is the best combination.
      std::vector<DHZObj> DHZ_values;
      for (size_t i = 0; i < combos_2b1c.size(); ++i) {
        std::vector<int> idxs = combos_2b1c[i];
        DHZObj D(jets4, idxs[2], idxs[3], idxs[0], idxs[1]);
        DHZ_values.push_back(D);
      }
      
      // Sort the values by their distance. Then, determine the distance
      // between the closest two. Based on this distance, determine which
      // pair we want to use. (NOTE: we might only have one combination.)
      std::sort(DHZ_values.begin(), DHZ_values.end());
      DHZObj chosenPair = DHZ_values[0];

      if (DHZ_values.size() >= 2) {
        float deltaD = fabs(DHZ_values[0].m_d - DHZ_values[1].m_d);
        if (deltaD < 30) {
          float pt0 = DHZ_values[0].HPt();
          float pt1 = DHZ_values[1].HPt();
          int idx = 0;
          if (pt1 > pt0) idx = 1;
          chosenPair = DHZ_values[idx];
        }
      }
      
      std::vector<JetObj> bjets4;
      bjets4.push_back(jets4[chosenPair.m_zIdx0]); // full JEC
      bjets4.push_back(jets4[chosenPair.m_zIdx1]);
      ZObj Z4(bjets4);
      
      std::vector<JetObj> cjets4;
      cjets4.push_back(jets4[chosenPair.m_hIdx0]); // full JEC
      cjets4.push_back(jets4[chosenPair.m_hIdx1]);
      HObj H4(cjets4);
      
      h_VH_tagFirst_2b1c->FillVH(Z4, H4, evtW);
    
    }//end-combos-2b1c

  }//end-analysis-jets
}// end Process

//== TERMINATE ================================================================

void VH_selection::Terminate(TList* mergedList, std::string outFileName) { }

//== END OF FILE ==============================================================

