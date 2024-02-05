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
    //bool passes_btag(JetObj& jet, float CSV_cut);
    //bool passes_ctag(JetObj& jet, float CvL_cut, float CvB_cut);

  private:
  
    // Histograms
    TH1D* h_evt;

    // Jet Plots
    JetPlots *h_jets_selected; // Selected Jets
    JetPlots *h_jets_all;	  // All jets found
    JetPlots *h_jets_MC;       // MC jets

    JetPlots *h_jets_4b;   // Jets that pass 4 b-jet requirements
    JetPlots *h_jets_2b2c; // Jets that pass 2 b, 2-cjet requirements

    // VH Plots (for MC Truth & Jets)
    VHPlots *h_VH_MC;                    // MC Truth events (GenObj)
    VHPlots *h_VH_MCjet_minDR;           // MC jet events (selected by min dR)
    VHPlots *h_VH_MCjet_minDR_noTag;     // " " " (w/o tagging)
    VHPlots *h_VH_MCjet_dRcollect;       // " " " (selected within dR cone)
    VHPlots *h_VH_MCjet_dRcollect_noTag; // " " " (w/o tagging)
    VHPlots *h_VH_MCjet_ideal;           // " " " (ideal scenario)
    VHPlots *h_VH_MCjet_DHZ;             // " " " (selected using DHZ alogrithm)
    
    // Reco Jet Plots (corresponding to MC jet selection methods)
    // under30 is for cases where events have m_Z < 30 GeV or m_H < 30 GeV
    RecoPlots *h_reco_minDR;                     
    RecoPlots *h_reco_minDR_under30;
    RecoPlots *h_reco_minDR_noTag;
    RecoPlots *h_reco_minDR_noTag_under30;
    RecoPlots *h_reco_dRcollect;
    RecoPlots *h_reco_dRcollect_under30;
    RecoPlots *h_reco_dRcollect_noTag;
    RecoPlots *h_reco_dRcollect_noTag_under30;
    RecoPlots *h_reco_ideal;
    RecoPlots *h_reco_ideal_under30;
    RecoPlots *h_reco_DHZ;
    RecoPlots *h_reco_DHZ_under30;
    
    // VH Plots (for our selection methods)
    VHPlots *h_VH_tagOnly;              // Tagging Only
    VHPlots *h_VH_tagOnly_noMassCorr;   // w/o mass correction
    VHPlots *h_VH_tagOnly_noJEC;        // w/o Jet Energy correction
    VHPlots *h_VH_tagOnly_2b1c;         // one c-jet requirement removed

    VHPlots *h_VH_algoFirst;              // Mass-Matching Prioritized
    VHPlots *h_VH_algoFirst_noMassCorr;   // w/o mass correction
    VHPlots *h_VH_algoFirst_noJEC;        // w/o Jet Energy correction
    VHPlots *h_VH_algoFirst_2b1c;         // one c-jet requirement removed

    VHPlots *h_VH_tagFirst;              // Tagging Prioritized
    VHPlots *h_VH_tagFirst_noMassCorr;   // w/o mass correction
    VHPlots *h_VH_tagFirst_noJEC;        // w/o Jet Energy correction
    VHPlots *h_VH_tagFirst_2b1c;         // one c-jet requirement removed

    // I don't know if we use these 4 but keep them here as is for now...
    VHPlots *h_VH_all;   // Any plots related to ALL cut types
    VHPlots *h_VH_select; // Any plots related to ALL selection types (after our MET cut)
    VHPlots *h_VH_alljet; // Using just any jets
    VHPlots *h_VH_seljet; // Using selected jets
    
    // Trigger Efficiency Plots - original triggers
    TriggerEffPlots *h_trig_2016_QuadJet_TripleTag;    // HLT_QuadJet45_Triple...
    TriggerEffPlots *h_trig_2016_DoubleJet_TripleTag;  // HLT_DoubleJet90_Double30_Triple...
    TriggerEffPlots *h_trig_2017_QuadJet_TripleTag;    // HLT_PFHT300PT30_QuadJet...
    TriggerEffPlots *h_trig_2018_QuadJet_TripleTag;    // HLT_PFHT330PT30_QuadJet...

    TriggerEffPlots *h_trig_2016_QuadJet_DoubleTag;    // 2016 but w/ DoubleBTag
    TriggerEffPlots *h_trig_2016_DoubleJet_DoubleTag;  // 2016 but w/ DoubleBTag
    TriggerEffPlots *h_trig_2017_QuadJet_noTag;        // 2017 but w/ no tagging requirement
    TriggerEffPlots *h_trig_2018_QuadJet_noTag;        // 2018 but w/ no tagging requirement

    // Versions that pass tagging requirements
    TriggerEffPlots *h_trig_2016_QuadJet_TripleTag_tagged;
    TriggerEffPlots *h_trig_2016_DoubleJet_TripleTag_tagged;
    TriggerEffPlots *h_trig_2017_QuadJet_TripleTag_tagged;
    TriggerEffPlots *h_trig_2018_QuadJet_TripleTag_tagged;

    TriggerEffPlots *h_trig_2016_QuadJet_DoubleTag_tagged;
    TriggerEffPlots *h_trig_2016_DoubleJet_DoubleTag_tagged;
    TriggerEffPlots *h_trig_2017_QuadJet_noTag_tagged;
    TriggerEffPlots *h_trig_2018_QuadJet_noTag_tagged;

    // Versions that pass all requirements in tagging names
    TriggerEffPlots *h_trig_2016_QuadJet_TripleTag_ideal;
    TriggerEffPlots *h_trig_2016_DoubleJet_TripleTag_ideal;
    TriggerEffPlots *h_trig_2017_QuadJet_TripleTag_ideal;
    TriggerEffPlots *h_trig_2017_QuadJet_TripleTag_RunB_ideal;
    TriggerEffPlots *h_trig_2018_QuadJet_TripleTag_ideal;

    TriggerEffPlots *h_trig_2016_QuadJet_DoubleTag_ideal;
    TriggerEffPlots *h_trig_2016_DoubleJet_DoubleTag_ideal;
    TriggerEffPlots *h_trig_2017_QuadJet_noTag_ideal;
    TriggerEffPlots *h_trig_2018_QuadJet_noTag_ideal;

    TriggerEffPlots *h_trig_2017_QuadJet_noTagV2_ideal;      
    TriggerEffPlots *h_trig_2017_QuadJet_noTagV3_ideal;      
    TriggerEffPlots *h_trig_2017_QuadJet_noTagV4_ideal;      
    TriggerEffPlots *h_trig_2017_QuadJet_noTagV5_ideal;      

    TriggerEffPlots *h_trig_2018_QuadJet_noTagV2_ideal;     
    TriggerEffPlots *h_trig_2018_QuadJet_noTagV3_ideal;
    TriggerEffPlots *h_trig_2018_QuadJet_noTagV4_ideal;
    TriggerEffPlots *h_trig_2018_QuadJet_noTagV5_ideal;

    // Versions where the tagging requirement is 3b
    TriggerEffPlots *h_trig_2016_QuadJet_TripleTag_3B;
    TriggerEffPlots *h_trig_2016_DoubleJet_TripleTag_3B;
    TriggerEffPlots *h_trig_2017_QuadJet_TripleTag_3B;
    TriggerEffPlots *h_trig_2017_QuadJet_TripleTag_RunB_3B;
    TriggerEffPlots *h_trig_2018_QuadJet_TripleTag_3B;

    // Versions where the tagging requirement is 2b2c
    TriggerEffPlots *h_trig_2016_QuadJet_TripleTag_2b2c;
    TriggerEffPlots *h_trig_2016_DoubleJet_TripleTag_2b2c;
    TriggerEffPlots *h_trig_2017_QuadJet_TripleTag_2b2c;
    TriggerEffPlots *h_trig_2017_QuadJet_TripleTag_RunB_2b2c;
    TriggerEffPlots *h_trig_2018_QuadJet_TripleTag_2b2c;

    // Efficiency Plots
    EffPlots *h_eff_tagOnly; // Efficiencies for Tagging Only
    EffPlots *h_eff_algoFirst; // Efficiencies for Mass-Matching Prioritized
    EffPlots *h_eff_tagFirst; // Efficiencies for Tagging Prioritized

    // CutFlows for event selections
    TH1D* h_evt_VbbHcc;
    TH1D* h_evt_MC_cutflow;
    TH1D* h_evt_tagOnly_cutflow;
    TH1D* h_evt_algoFirst_cutflow;
    TH1D* h_evt_tagFirst_cutflow;
    TH1D* h_evt_MCjet_ideal_cutflow;

    // CutFlows for reconstruction
    TH1D* h_jet_cutflow;
    TH1D* h_elec_cutflow;
    TH1D* h_muon_cutflow;

    // Miscellaneous
    TH1D* h_nCombos;
    TH1D* h_dR_ccjet;
    TH1D* h_dR_bbjet;
    TH1D* h_dPhi_ccjet;
    TH1D* h_dPhi_bbjet;

    TH1D* h_mistag_leading;
    TH1D* h_mistag_all;

    TH1D* h_bRegCorr;
    TH1D* h_cRegCorr;
    TH1D* h_JetMass;

    TH1D* h_JetMassBefore;
    TH1D* h_JetMassAfter;

    GenPlots *h_genJet_all;
    GenPlots *h_genJet_cuts;
    GenPlots *h_genJet_VbbHcc;

    TH1D* h_genWeight;
    TH1D* h_puSF;
    TH1D* h_l1preW;
    TH1D* h_trigSF;
    TH1D* h_btagW;
    TH1D* h_ctagW;
    TH1D* h_evtW;

    TH1D* h_lheW;
    TH1D* h_genW;

    TH1D* h_nJet;
    TH1D* h_nAnalysisJet;

    TH1D* h_nMuon;
    TH1D* h_nElec;
};

#endif
