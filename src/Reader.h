//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 14 16:02:21 2020 by ROOT version 6.14/09
// from TTree Events/Events
// found on file: 8F93A522-A364-BE4A-8A5D-26591CD081F1.root
//////////////////////////////////////////////////////////

#ifndef Reader_h
#define Reader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class Reader : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
   
   // Readers to access the data (delete the ones you do not need).
#if defined(DATA_2016) || defined(DATA_2017) || defined(DATA_2018)
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<UInt_t> luminosityBlock = {fReader, "luminosityBlock"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
#endif

#if defined(MC_2016) || defined(MC_2017)
   TTreeReaderValue<Float_t> L1PreFiringWeight_Dn = {fReader, "L1PreFiringWeight_Dn"};
   TTreeReaderValue<Float_t> L1PreFiringWeight_Nom = {fReader, "L1PreFiringWeight_Nom"};
   TTreeReaderValue<Float_t> L1PreFiringWeight_Up = {fReader, "L1PreFiringWeight_Up"};
#endif
   
   //Jet
   TTreeReaderValue<UInt_t> nJet = {fReader, "nJet"};
   TTreeReaderArray<Float_t> Jet_pt = {fReader, "Jet_pt"};
   TTreeReaderArray<Float_t> Jet_eta = {fReader, "Jet_eta"};
   TTreeReaderArray<Float_t> Jet_phi = {fReader, "Jet_phi"};
   TTreeReaderArray<Float_t> Jet_mass = {fReader, "Jet_mass"};
   TTreeReaderArray<Float_t> Jet_btagDeepB = {fReader, "Jet_btagDeepB"};
   TTreeReaderArray<Float_t> Jet_btagDeepFlavB = {fReader, "Jet_btagDeepFlavB"};
   TTreeReaderArray<Int_t> Jet_puId = {fReader, "Jet_puId"};
   TTreeReaderArray<Float_t> Jet_puIdDisc = {fReader, "Jet_puIdDisc"};
   TTreeReaderArray<Float_t> Jet_bRegCorr = {fReader, "Jet_bRegCorr"};
   TTreeReaderArray<Float_t> Jet_cRegCorr = {fReader, "Jet_cRegCorr"};
   TTreeReaderArray<Float_t> Jet_bRegRes = {fReader, "Jet_bRegRes"};
   TTreeReaderArray<Float_t> Jet_cRegRes = {fReader, "Jet_cRegRes"};

#if defined(NANOAODV9) || defined(DATA_2016) || defined(DATA_2017) || defined(DATA_2018)
   TTreeReaderArray<Float_t> Jet_btagDeepFlavCvL = {fReader, "Jet_btagDeepFlavCvL"};
   TTreeReaderArray<Float_t> Jet_btagDeepFlavCvB = {fReader, "Jet_btagDeepFlavCvB"};
#endif

#if defined(NANOAODV7)
   TTreeReaderArray<Float_t> Jet_btagDeepFlavC = {fReader, "Jet_btagDeepFlavC"};
   TTreeReaderArray<Float_t> FatJet_btagDDCvB = {fReader, "FatJet_btagDDCvB"};
   TTreeReaderArray<Float_t> FatJet_btagDDCvL = {fReader, "FatJet_btagDDCvL"};
   TTreeReaderArray<Float_t> FatJet_btagDDBvL = {fReader, "FatJet_btagDDBvL"};
#endif
#if defined(NANOAODV9)
   TTreeReaderArray<Float_t> FatJet_btagDDCvB = {fReader, "FatJet_btagDDCvBV2"};
   TTreeReaderArray<Float_t> FatJet_btagDDCvL = {fReader, "FatJet_btagDDCvLV2"};
   TTreeReaderArray<Float_t> FatJet_btagDDBvL = {fReader, "FatJet_btagDDBvLV2"};
#endif
 
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
   TTreeReaderValue<Float_t> Pileup_nTrueInt = {fReader, "Pileup_nTrueInt"};
   TTreeReaderValue<Float_t> genWeight = {fReader, "genWeight"};
   TTreeReaderArray<Int_t> Jet_hadronFlavour = {fReader, "Jet_hadronFlavour"};
#endif

   // Electron
   TTreeReaderValue<UInt_t> nElectron = {fReader, "nElectron"};
   TTreeReaderArray<Float_t> Electron_pt = {fReader, "Electron_pt"};
   TTreeReaderArray<Float_t> Electron_eta = {fReader, "Electron_eta"};
   TTreeReaderArray<Float_t> Electron_phi = {fReader, "Electron_phi"};
   TTreeReaderArray<Float_t> Electron_mass = {fReader, "Electron_mass"};
   TTreeReaderArray<Int_t> Electron_charge = {fReader, "Electron_charge"};
   TTreeReaderArray<Float_t> Electron_deltaEtaSC = {fReader, "Electron_deltaEtaSC"};
   TTreeReaderArray<Int_t> Electron_cutBased = {fReader, "Electron_cutBased"};

   // Muon
   TTreeReaderValue<UInt_t> nMuon = {fReader, "nMuon"};
   TTreeReaderArray<Float_t> Muon_pt = {fReader, "Muon_pt"};
   TTreeReaderArray<Float_t> Muon_eta = {fReader, "Muon_eta"};
   TTreeReaderArray<Float_t> Muon_phi = {fReader, "Muon_phi"};
   TTreeReaderArray<Float_t> Muon_mass = {fReader, "Muon_mass"};
   TTreeReaderArray<Int_t> Muon_charge = {fReader, "Muon_charge"};
   TTreeReaderArray<Float_t> Muon_pfRelIso04_all = {fReader, "Muon_pfRelIso04_all"};
   TTreeReaderArray<Bool_t> Muon_mediumId = {fReader, "Muon_mediumId"};
   TTreeReaderArray<Bool_t> Muon_looseId = {fReader, "Muon_looseId"};

   TTreeReaderValue<Float_t> MET_pt = {fReader, "MET_pt"};

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
TTreeReaderValue<UInt_t> nGenPart = {fReader, "nGenPart"};
TTreeReaderArray<Float_t> GenPart_eta = {fReader, "GenPart_eta"};
TTreeReaderArray<Float_t> GenPart_mass = {fReader, "GenPart_mass"};
TTreeReaderArray<Float_t> GenPart_phi = {fReader, "GenPart_phi"};
TTreeReaderArray<Float_t> GenPart_pt = {fReader, "GenPart_pt"};
TTreeReaderArray<Int_t> GenPart_genPartIdxMother = {fReader, "GenPart_genPartIdxMother"};
TTreeReaderArray<Int_t> GenPart_pdgId = {fReader, "GenPart_pdgId"};
TTreeReaderArray<Int_t> GenPart_status = {fReader, "GenPart_status"};

TTreeReaderValue<UInt_t> nGenJet = {fReader, "nGenJet"};
TTreeReaderArray<Float_t> GenJet_eta = {fReader, "GenJet_eta"};
TTreeReaderArray<Float_t> GenJet_mass = {fReader, "GenJet_mass"};
TTreeReaderArray<Float_t> GenJet_phi = {fReader, "GenJet_phi"};
TTreeReaderArray<Float_t> GenJet_pt = {fReader, "GenJet_pt"};
TTreeReaderArray<Int_t> GenJet_partonFlavour = {fReader, "GenJet_partonFlavour"};
TTreeReaderArray<UChar_t> GenJet_hadronFlavour = {fReader, "GenJet_hadronFlavour"};
TTreeReaderArray<Int_t> Jet_genJetIdx = {fReader, "Jet_genJetIdx"};
#endif

   // Triggers
   TTreeReaderValue<Bool_t> HLT_IsoMu24 = {fReader, "HLT_IsoMu24"};

#if defined(MC_2016) || defined(DATA_2016)
   TTreeReaderValue<Bool_t> HLT_QuadJet45_TripleBTagCSV_p087 = {fReader, "HLT_QuadJet45_TripleBTagCSV_p087"};
   TTreeReaderValue<Bool_t> HLT_QuadJet45_DoubleBTagCSV_p087 = {fReader, "HLT_QuadJet45_DoubleBTagCSV_p087"};
   TTreeReaderValue<Bool_t> HLT_DoubleJet90_Double30_TripleBTagCSV_p087 = {fReader, "HLT_DoubleJet90_Double30_TripleBTagCSV_p087"};
   TTreeReaderValue<Bool_t> HLT_DoubleJet90_Double30_DoubleBTagCSV_p087 = {fReader, "HLT_DoubleJet90_Double30_DoubleBTagCSV_p087"};
#endif

// NOTE: There is an issue with the desired trigger is that it only exists for
// 2017C-F. We need a different trigger for 2017B.
#if defined(MC_2017) || defined(DATA_2017)
   #if defined(MC_2017) || !defined(DATA_2017B)
   TTreeReaderValue<Bool_t> HLT_PFHT300PT30_QuadPFJet_75_60_45_40 = {fReader, "HLT_PFHT300PT30_QuadPFJet_75_60_45_40"};
   TTreeReaderValue<Bool_t> HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0 = {fReader, "HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0"};
   #endif
#endif

#if defined(DATA_2017B)
   TTreeReaderValue<Bool_t> HLT_HT300PT30_QuadJet_75_60_45_40_TripeCSV_p07 = {fReader, "HLT_HT300PT30_QuadJet_75_60_45_40_TripeCSV_p07"};
#endif


#if defined(MC_2018) || defined(DATA_2018)
   TTreeReaderValue<Bool_t> HLT_PFHT330PT30_QuadPFJet_75_60_45_40 = {fReader, "HLT_PFHT330PT30_QuadPFJet_75_60_45_40"};
   TTreeReaderValue<Bool_t> HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5 = {fReader, "HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5"};
#endif

#if defined(MC_2017) || defined(MC_2018) || defined(DATA_2017) || defined(DATA_2018)

   /*#if !defined(MC_2017B)
   TTreeReaderValue<Bool_t> HLT_QuadPFJet98_83_71_15 = {fReader, "HLT_QuadPFJet98_83_71_15"};
   TTreeReaderValue<Bool_t> HLT_QuadPFJet103_88_75_15 = {fReader, "HLT_QuadPFJet103_88_75_15"};
   TTreeReaderValue<Bool_t> HLT_QuadPFJet105_88_76_15 = {fReader, "HLT_QuadPFJet105_88_76_15"};
   TTreeReaderValue<Bool_t> HLT_QuadPFJet111_90_80_15 = {fReader, "HLT_QuadPFJet111_90_80_15"};
   #endif*/
#endif

   Reader(TTree * /*tree*/ =0) {}

   virtual ~Reader() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(Reader,0);

};

#endif

#ifdef Reader_cxx
void Reader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t Reader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef Reader_cxx
