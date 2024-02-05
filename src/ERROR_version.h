#ifndef VH_selection_h
#define VH_selection_h

// Include Statements /////////////////////////////////////////////////////////

#include "Reader.h"
#include "Selector.h"
#include "Obj.cxx"
#include "Plots.cxx"

// Class Definition ///////////////////////////////////////////////////////////

class VH_selection : public Selector {

public:
  // Constructor & Deconstructor
  VH_selection() {};
  virtual ~VH_selection();

  // Methods
  virtual void SlaveBegin(Reader* r);
  virtual void Process(Reader* r);
  virtual void SlaveTerminate(Reader* r) {};
  virtual void Terminate(TList* mergedList, std::string outFileName);
  
  std::vector<std::vector<int> > DauIdxs_ZH(Reader* r);
  
private:
  
  // Histograms
  TH1D* h_evt;
  
  // VH Plots
  VHPlots* h_VH_MC;    // MC Truth events
  VHPlots* h_VH_all;   // Plots related to general events

  VHPlots* h_VH_tags;  // Tagging Only
  VHPlots* h_VH_algo;  // Mass Matching Prioritized
  VHPlots* h_VH_both;  // Tagging Prioritized

  VHPlots* h_VH_tags_all; // Tagging Only (before cuts)
  VHPlots* h_VH_algo_all; // Matching Prioritized (before cuts)
  VHPlots* h_VH_both_all; // Tagging Prioritized (before cuts)

  VHPlots* h_VH_tags_afterTag; // Tagging Only (after tagging cuts)
  VHPlots* h_VH_algo_afterTag; // Match Prioritized (after tagging cuts)
  VHPlots* h_VH_both_afterTag; // Tagging Prioritzied (after tagigng cuts) 
 
  // CutFlows for event selections
  TH1D* h_evt_MC_cutflow;
  TH1D* h_evt_tags_cutflow;
  TH1D* h_evt_algo_cutflow;
  TH1D* h_evt_both_cutflow;
  
  // CutFlows for reconstruction
  TH1D* h_jet_cutflow;
  TH1D* h_elec_cutflow;
  TH1D* h_muon_cutflow;
  
  // Efficiency Plots
  EffPlots* h_eff_tags;
  EffPlots* h_eff_algo;
  EffPlots* h_eff_both;
  
  // Miscellaneous Plots
  TH1D* h_nJet;
  TH1D* h_nBjet_loose;
  TH1D* h_nCjet_loose;
  TH1D* h_nBjet_medium;
  TH1D* h_nCjet_medium;
  TH1D* h_bScore;
  TH1D* h_cScore;
  TH2D* h_btag_v_ctag;

  TH1D* h_nJet_all;
  TH1D* h_nBjet_loose_all;
  TH1D* h_nCjet_loose_all;
  TH1D* h_nBjet_medium_all;
  TH1D* h_nCjet_medium_all;
  TH1D* h_bScore_all;
  TH1D* h_cScore_all;
  TH2D* h_btag_v_ctag_all;

};

#endif
