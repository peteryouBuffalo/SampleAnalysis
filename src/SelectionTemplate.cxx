#define SelectionTemplate_cxx

#include <math.h>

#include "TList.h"
#include "TParameter.h"

#include "SelectionTemplate.h"
#include "Global.h"
#include "Obj.cxx"

SelectionTemplate::SelectionTemplate(bool isData) : Selector(isData), h_zee_jet(0), h_zmm_jet(0) {}

SelectionTemplate::~SelectionTemplate() {
  if (h_zee_jet) delete h_zee_jet;
  if (h_zmm_jet) delete h_zmm_jet;
} 

void SelectionTemplate::SlaveBegin(Reader* r) {
  h_zee_jet = new ZbPlots("Zee_jet") ;
  h_zmm_jet = new ZbPlots("Zmm_jet") ;

  //Add histograms to fOutput so they can be saved in Processor::SlaveTerminate
  std::vector<TH1*> tmp = h_zee_jet->returnHisto() ;
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);

  tmp = h_zmm_jet->returnHisto() ;
  for(size_t i=0;i<tmp.size();i++) r->GetOutputList()->Add(tmp[i]);
}

void SelectionTemplate::Process(Reader* r) {
  //std::cout << "\n Data type: " << m_isData ;
  //std::cout << "\n Number of jets: " << *(r->nJet) ;
  //if (*(r->nJet)>0) std::cout << "\n First jet: " << (r->Jet_eta)[0] ;

  //=============Get objects============= 
  std::vector<LepObj> eles ;
  for (unsigned int i = 0 ; i < *(r->nElectron) ; ++i) {
    LepObj ele((r->Electron_pt)[i],(r->Electron_eta)[i],(r->Electron_phi)[i],(r->Electron_mass)[i],0) ;
    if (ele.m_lvec.Pt() > CUTS.Get<float>("lep_pt") && fabs(ele.m_lvec.Eta()) < CUTS.Get<float>("lep_eta")) eles.push_back(ele) ;
  }

  std::vector<LepObj> muons ;
  for (unsigned int i = 0 ; i < *(r->nMuon) ; ++i) {
    LepObj muon((r->Muon_pt)[i],(r->Muon_eta)[i],(r->Muon_phi)[i],(r->Muon_mass)[i],(r->Muon_pfRelIso03_all)[i]) ;
    muons.push_back(muon) ;
  }

  std::vector<JetObj> jets ;
  for (unsigned int i = 0 ; i < *(r->nJet) ; ++i) {
    JetObj jet((r->Jet_pt)[i],(r->Jet_eta)[i],(r->Jet_phi)[i],(r->Jet_mass)[i],(r->Jet_btagDeepB)[i]) ;
    jets.push_back(jet) ;
  }
  
  //Zee + jets
  if (eles.size() >= 2) {
    if (jets.size() >= 1) {
      ZObj Z(eles[0],eles[1]) ;
      h_zee_jet->Fill(Z) ;
    }

  }

  //Zmm + jets
  if (muons.size() >= 2) {
    if (jets.size() >= 1) {
      ZObj Z(muons[0],muons[1]) ;
      h_zmm_jet->Fill(Z) ;
    }

  }


} //end Process



void SelectionTemplate::Terminate(TList* mergedList, std::string outFileName) {
  
  TList* aList = new TList() ;
  TParameter<double>* a = new TParameter<double>("lep_eta",CUTS.Get<float>("lep_eta")) ;
  aList->Add(a) ;
  a = new TParameter<double>("lep_pt",CUTS.Get<float>("lep_pt")) ;
  aList->Add(a) ;
  
  TFile* f = new TFile(outFileName.c_str(), "UPDATE") ;
  aList->Write("SelectionTemplateCuts",TObject::kSingleKey) ;

  f->Close() ;
}
