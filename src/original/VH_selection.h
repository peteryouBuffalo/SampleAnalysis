#ifndef VH_selection_h
#define VH_selection_h

//== Include Statements =======================================================

#include "Reader.h"
#include "Selector.h"
#include "Obj.cxx"
#include "Plots.cxx"

//== Class Definition =========================================================

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
    bool passes_btag(JetObj& jet, float CSV_cut);
    bool passes_ctag(JetObj& jet, float CvL_cut, float CvB_cut);

  private:
  
    // Histograms
    TH1D* h_evt;

    // Jet Plots
    JetPlots *h_VH_jets;     // Selected jets
    JetPlots *h_VH_jets_all; // All jets captured

    // VH Plots
    VHPlots *h_VH_MC;    // MC Truth events
    VHPlots *h_VH_tags;  // Tagging Only
    VHPlots *h_VH_algo;  // Mass-Matching Prioritized
    VHPlots *h_VH_both;  // Tagging Prioritized
    VHPlots *h_VH_duong; // Tagging Prioritized (Duong version)
    VHPlots *h_VH_all;   // Any plots related to ALL cut types

    // Eff Plots
    EffPlots *h_eff_tags; 
    EffPlots *h_eff_algo;
    EffPlots *h_eff_both;
    EffPlots *h_eff_duong;

    // CutFlows for event selections
    TH1D* h_evt_MC_cutflow;
    TH1D* h_evt_tags_cutflow;
    TH1D* h_evt_algo_cutflow;
    TH1D* h_evt_both_cutflow;
    TH1D* h_evt_duong_cutflow;

    // CutFlows for reconstruction
    TH1D* h_jet_cutflow;
    TH1D* h_elec_cutflow;
    TH1D* h_muon_cutflow;
};

#endif
