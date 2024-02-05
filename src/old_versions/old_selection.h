#ifndef VH_selection_h
#define VH_selection_h

//Template for creating new selection

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
  VH_selection() {};
  virtual ~VH_selection() ;
  // These methods are called at the corresponding stage of processing of TSelector
  virtual void SlaveBegin(Reader* r) ;
  virtual void Process(Reader* r) ;
  virtual void SlaveTerminate(Reader* r) {} ;
  virtual void Terminate(TList* mergedList, std::string outFileName) ;

  std::vector<std::vector<int> > DauIdxs_ZH(Reader* r);
  //float check_efficiency(HObj& H, ZObj& Z, std::vector<JetObj> genObjs);

private:

  float D_HZ(float mH, float mZ);
  float get_mass_from_indices(std::vector<JetObj>& jets, int idx0, int idx1);

  //histograms
  TH1D* h_evt;
  
  //VH Plots & histograms
  VHPlots* h_VH;
  VHPlots* h_VH_tags;
  VHPlots* h_VH_algo;
  VHPlots* h_VH_both;
  VHPlots* h_VH_bothTags;
  VHPlots* h_VH_bothAlgo; 
 
  // CutFlows
  TH1D* h_evt_cutflow;
  TH1D* h_evt_tags_cutflow;
  TH1D* h_evt_algo_cutflow;
  TH1D* h_evt_both_cutflow;
  TH1D* h_jet_cutflow;
  TH1D* h_elec_cutflow;
  TH1D* h_muon_cutflow;
  
  // Tagging Scores
  TH1D* h_CSV;
  TH1D* h_CvL;

  // Misc.
  TH1D* h_Nselected;
  TH1D* h_Nbjet;
  TH1D* h_Ncjet;
  TH1D* h_Z_dM;
  TH1D* h_H_dM;

  // Mass-Matching Plots
  TH2D* h_tags_MH_v_MZ; TH2D* h_tags_MH_v_MZ_select;
  TH2D* h_algo_MH_v_MZ; TH2D* h_algo_MH_v_MZ_select;
  TH2D* h_both_MH_v_MZ; TH2D* h_both_MH_v_MZ_select;

  // Efficiency Plots
  EffPlots* h_eff_tags;
  EffPlots* h_eff_algo;
  EffPlots* h_eff_both;
} ;

#endif
