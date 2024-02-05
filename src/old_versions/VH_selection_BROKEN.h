#ifndef VH_selection_h
#define VH_selection_h

// Template for creating new selection
#include "Reader.h"
#include "Selector.h"
#include "Obj.cxx"
#include "Plots.cxx"

//The Selection does not have Begin, since we do not have anything to do at the begining (overall) at client
//The histograms, ..., are booked and added to output list at SlaveBegin
//We need to have terminate since we might want to do overall tasks related to this selection only. At termination, all inforamtion from slaves is added. Example task is cutflow for this selection, which need information from all slaves.
//SlaveTerminate and Terminate are here for reservation now. histograms are written back in SlaveTerminate of Processor class

class VH_selection : public Selector {
public:

  // Constructor & Deconstructor
  VH_selection() {};
  virtual ~VH_selection() ;
  
  // These methods are called at the corresponding stage of processing of TSelector
  virtual void SlaveBegin(Reader* r) ;
  virtual void Process(Reader* r) ;
  virtual void SlaveTerminate(Reader* r) {} ;
  virtual void Terminate(TList* mergedList, std::string outFileName) ;

  std::vector<std::vector<int> > DauIdxs_ZH(Reader *r);

private:

  // Events
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
  
  // CutFlows for event selection
  TH1D* h_evt_MC_cutflow;
  TH1D* h_evt_tags_cutflow;
  TH1D* h_evt_algo_cutflow;
  TH1D* h_evt_both_cutflow;
  
  // CutFlows for reconstruction
  TH1D* h_jet_cutflow;
  TH1D* h_elec_cutflow;
  TH1D* h_muon_cutflow;
  
  // Tagging Scores
  TH1D* h_CSV;	// b-tag scores
  TH1D* h_CvL;	// c-tag scores
  
  // Efficiency Plots
  EffPlots* h_eff_tags;
  EffPlots* h_eff_algo;
  EffPlots* h_eff_both;
  
};
#endif
