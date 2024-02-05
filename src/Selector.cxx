#include "Selector.h"
#include "Global.h"
#include "math.h"

void Selector::Process(Reader* r) { 
} 

void Selector::SetRandom() {
  m_rand = new TRandom3() ;
}

void Selector::SetLumiMaskFilter(std::string fName_lumiMaskFilter) {
  m_lumiFilter.Set(fName_lumiMaskFilter) ;
}

// This is where we set up the Btag calibration. We have four parameters with it:
// 1. csvFileName - this is the BTV csv file that contains the values we're interested in
// 2. 
void Selector::SetBtagCalib(std::string csvFileName, std::string taggerName, std::string effFileName, std::string btagUncType) {
  m_btagCal = BTagCalibration(taggerName, csvFileName) ;
  m_btagReader = BTagCalibrationReader(BTagEntry::OP_MEDIUM,  // operating point
                                       "central",            //central sys type
                                       {"up","down"});       //other sys type
  
  m_btagReader.load(m_btagCal,     // calibration instance
            BTagEntry::FLAV_B,    // btag flavour
            "comb") ;             // measurement type
  m_btagReader.load(m_btagCal,  
            BTagEntry::FLAV_C,    
            "comb") ;            
  m_btagReader.load(m_btagCal, 
            BTagEntry::FLAV_UDSG,
            "incl") ;           
  m_btagUncType = btagUncType;
  m_btagEffFile = new TFile(effFileName.c_str(),"READ") ;

}

// This is where we setup the Ctag calibration. We have four parameters with it.
void Selector::SetCtagCalib(std::string csvFileName, std::string taggerName, std::string effFileName, std::string ctagUncType) {

  std::cout << "Setting ctagCal & ctagReader..." << std::endl;
  std::cout << ">>> ctagCal..." << std::endl;
  m_ctagCal = BTagCalibration(taggerName, csvFileName);
  std::cout << ">>> ctagReader..." << std::endl;
  m_ctagReader = BTagCalibrationReader(BTagEntry::OP_MEDIUM, // operating point
	                               "central",            // central sys type
				       {"up","down"});       // other sys type
  
  std::cout << "Loading ctagReader..." << std::endl;
  m_ctagReader.load(m_ctagCal, BTagEntry::FLAV_B, "comb");
  m_ctagReader.load(m_ctagCal, BTagEntry::FLAV_C, "comb");
  m_ctagReader.load(m_ctagCal, BTagEntry::FLAV_UDSG, "incl");

  std::cout << "Loading ctagEff file..." << std::endl;
  m_ctagUncType = ctagUncType;
  m_ctagEffFile = new TFile(effFileName.c_str(), "READ");
  std::cout << "CtagCalib completed." << std::endl;
}

void Selector::SetEleEffCorr(std::vector<std::string> fName_trig,std::string fName_recSF, std::string fName_IDSF, std::vector<float> w_trig, std::string eleUncType) {
  std::string trigN("EGamma_SF2D");
  TFile* fRec = new TFile(fName_recSF.c_str(),"READ") ;
  TFile* fID = new TFile(fName_IDSF.c_str(),"READ") ;
  m_hSF_eleRec = (TH2F*)fRec->Get("EGamma_SF2D") ;
  m_hSF_eleID = (TH2F*)fID->Get("EGamma_SF2D") ;
  for(std::string fN : fName_trig) {
    TFile* f = new TFile(fN.c_str(),"READ");
    m_hSF_eleTrig.push_back((TH2F*)f->Get(trigN.c_str()));
    m_hSF_eleTrig.back()->SetDirectory(0);
  }

  for(float w : w_trig) m_eleTrig_w.push_back(w) ;
  m_eleUncType = eleUncType;
}

//multiple inputs to deal with different SFs for different run periods 
void Selector::SetMuonEffCorr(std::vector<std::string> fName_trig, std::vector<std::string> fName_ID, std::vector<std::string> fName_iso, std::vector<float> w_trig, std::vector<float> w_ID, std::vector<float> w_iso, std::string muonUncType) {
  std::string trigN("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
  //std::string idN("NUM_MediumID_DEN_genTracks_eta_pt_syst");
  //std::string isoN("NUM_TightRelIso_DEN_MediumID_eta_pt_syst");//FIXME tight iso?
  std::string idN("NUM_MediumID_DEN_genTracks_eta_pt");//use this to get total unc = stat + syst
  std::string isoN("NUM_TightRelIso_DEN_MediumID_eta_pt");//FIXME tight iso?
  std::string isoN1("NUM_TightRelIso_DEN_MediumID_eta_pt");//FIXME: only stat available for GH for 2016 legacy
#if defined(MC_2017)
  trigN = "IsoMu27_PtEtaBins/abseta_pt_ratio";
  //idN = "NUM_MediumID_DEN_genTracks_pt_abseta_syst";
  //isoN = "NUM_TightRelIso_DEN_MediumID_pt_abseta_syst";
  idN = "NUM_MediumID_DEN_genTracks_pt_abseta"; //use this to get total unc = stat + syst
  isoN = "NUM_TightRelIso_DEN_MediumID_pt_abseta";
#endif
#if defined(MC_2018)
  trigN = "IsoMu24_PtEtaBins/abseta_pt_ratio";
  //idN = "NUM_MediumID_DEN_TrackerMuons_pt_abseta_syst";
  //isoN = "NUM_TightRelIso_DEN_MediumID_pt_abseta_syst";
  idN = "NUM_MediumID_DEN_TrackerMuons_pt_abseta";
  isoN = "NUM_TightRelIso_DEN_MediumID_pt_abseta";
#endif
  for(std::string fN : fName_trig) {
    TFile* f = new TFile(fN.c_str(),"READ");
    m_hSF_muonTrig.push_back((TH2F*)f->Get(trigN.c_str()));
    m_hSF_muonTrig.back()->SetDirectory(0);
  }
  for(std::string fN : fName_ID) {
    TFile* f = new TFile(fN.c_str(),"READ");
    m_hSF_muonID.push_back((TH2F*)f->Get(idN.c_str()));
    m_hSF_muonID.back()->SetDirectory(0);
  }
  for(std::string fN : fName_iso) {
    TFile* f = new TFile(fN.c_str(),"READ");
    if (fN.find("RunGH") != std::string::npos) {
      m_hSF_muonIso.push_back((TH2F*)f->Get(isoN1.c_str()));
    }
    else {
      m_hSF_muonIso.push_back((TH2F*)f->Get(isoN.c_str()));
    }
    m_hSF_muonIso.back()->SetDirectory(0);
  }
  for(float w : w_trig) m_muonTrig_w.push_back(w) ;
  for(float w : w_iso) m_muonIso_w.push_back(w) ;
  for(float w : w_ID) m_muonID_w.push_back(w) ;

  m_muonUncType = muonUncType;
}

void Selector::SetPileupSF(std::string fName_puSF) {
  TFile* f = new TFile(fName_puSF.c_str(),"READ") ;
  m_hSF_pu = (TH1D*)f->Get("pileup_ratio") ;
  m_hSF_pu->SetDirectory(0);
}

float Selector::PileupSF(int nTrueInt) {
  int iBin = m_hSF_pu->FindFixBin(nTrueInt) ;
  return m_hSF_pu->GetBinContent(iBin); 
}

void Selector::SetPUjetidCalib(std::string fName_PUjetID_SF,std::string fName_PUjetID_eff,std::string jetPUidUncType){
  //get files
  PUjetID_SF = new TFile(fName_PUjetID_SF.c_str(),"READ") ;
  PUjetID_eff = new TFile(fName_PUjetID_eff.c_str(),"READ") ;
  m_jetPUidUncType = jetPUidUncType ;
}

float Selector::PUjetWeight(std::vector<JetObj>& jets){
  //get SF
  std::string effSF = "h2_eff_sf"+m_year+"_M";
  TH2D* heffSF = (TH2D*)PUjetID_SF->Get(effSF.c_str()) ;

  //h2_eff_sf2016_M_Systuncty
  std::string effSFunc = "h2_eff_sf"+m_year+"_M_Systuncty";
  TH2D* heffSFunc = (TH2D*)PUjetID_SF->Get(effSFunc.c_str()) ;
  
  //get eff
  std::string effMC = "h2_eff_mc"+m_year+"_M";
  TH2D* heffMC = (TH2D*)PUjetID_eff->Get(effMC.c_str()) ;

  float pMC(1.0);
  float pData(1.0);
  float weight(1.) ;
  for (std::vector<JetObj>::iterator jetIt = jets.begin() ; jetIt != jets.end() ; ++jetIt) {
  float jetPt = (jetIt->m_lvec).Pt() ;
  float jetEta = (jetIt->m_lvec).Eta() ;
  float PUjet = jetIt->m_puid;
  if (fabs(jetEta) > CUTS.Get<float>("jet_eta") || jetPt < CUTS.Get<float>("jet_pt"))continue;

  int iBin = heffSF->FindFixBin(jetPt,jetEta) ; 
  float sf = heffSF->GetBinContent(iBin); 
  float sf_unc = heffSFunc->GetBinContent(iBin);
  if (m_jetPUidUncType == "up") sf = sf + sf_unc;
  if (m_jetPUidUncType == "down") sf = sf - sf_unc;
  if (sf < 0.5) sf = 0.5;
  if (sf > 1.5) sf = 1.5;

  int bin = heffMC->FindFixBin(jetPt,jetEta) ; 
  float eff = heffMC->GetBinContent(bin); 

   if (jetPt < 50){
      if (PUjet > 0.61){
      pData = pData*sf*eff ;
      pMC = pMC*eff ;
    }
    else {
      pData = pData*(1-sf*eff) ;
      pMC = pMC*(1-eff) ;
     }
   }//end jet pt > 50 GeV
  }// end jet loop
  if (pMC > 0) weight = pData/pMC ;
   return weight ;
}

// This method is what we use to determine the scale factor / event weight for b-tagging.
// These annotations are added by P.W. Young so we know what's going on.
float Selector::CalBtagWeight(std::vector<JetObj>& jets, std::string jet_main_btagWP, std::string uncType) {

  // Get the calibration files (suggested by BTV group)
  std::string bN = "b_pt_eff_"+m_year;
  std::string cN = "c_pt_eff_"+m_year;
  std::string lN = "l_pt_eff_"+m_year;
  if (jet_main_btagWP.find("deepJet") != std::string::npos) {
    bN = "bdj_pt_eff_"+m_year;
    cN = "cdj_pt_eff_"+m_year;
    lN = "ldj_pt_eff_"+m_year;
  }
  TH1D* hEff_b = (TH1D*)m_btagEffFile->Get(bN.c_str());
  TH1D* hEff_c = (TH1D*)m_btagEffFile->Get(cN.c_str());
  TH1D* hEff_l = (TH1D*)m_btagEffFile->Get(lN.c_str());
  float pMC(1.);
  float pData(1.);

  // Go through each jet in the given event (from function parameter)
  for (std::vector<JetObj>::iterator jetIt = jets.begin() ; jetIt != jets.end() ; ++jetIt) {

    // Get the pT of the jet & its flavor
    float jetPt = (jetIt->m_lvec).Pt() ;
    int iBin = hEff_b->FindFixBin(jetPt) ; //return overflow bin if jetPt > max pt range
    unsigned flav = jetIt->m_flav ;

    // Determine what type of uncertainty we want.
    std::string uncTypeInput = "central";
    float eff = hEff_l->GetBinContent(iBin); //jet with pt > max pt range of efficinecy histogram will get the eff of overflow bins
    //if (eff <= 0) std::cout << "\n Warning: Efficiency <=0, " << eff ; //we do not want eff = 0 
   
    // Check the flavors of the jets & get the appropriate efficiency.
    BTagEntry::JetFlavor flavCode(BTagEntry::FLAV_UDSG) ;
    if (flav == 5) {
      eff = hEff_b->GetBinContent(iBin);
      flavCode = BTagEntry::FLAV_B;
    }
    if (flav == 4) {
      eff = hEff_c->GetBinContent(iBin);
      flavCode = BTagEntry::FLAV_C;
    }

    if (uncType == "up" || uncType == "down") uncTypeInput = uncType; //all bc and light are up together
    if (uncType == "light_up" && flav != 4 && flav != 5) uncTypeInput = "up";
    if (uncType == "light_down" && flav != 4 && flav != 5) uncTypeInput = "down";
    if (uncType == "bc_up" && (flav == 4 || flav == 5)) uncTypeInput = "up";
    if (uncType == "bc_down" && (flav == 4 || flav == 5)) uncTypeInput = "down";

    // Determine the scale factor for this jet
    float sf = m_btagReader.eval_auto_bounds(
                 uncTypeInput, 
                 flavCode, 
                 fabs((jetIt->m_lvec).Eta()), // absolute value of eta
                 jetPt
    );

    //pass b-tagging requirement
    if (jetIt->m_deepCSV > CUTS.Get<float>("jet_"+jet_main_btagWP+"_" + m_year)) {
      pData = pData*sf*eff ;
      pMC = pMC*eff ;
    }
    else {
      pData = pData*(1-sf*eff) ;
      pMC = pMC*(1-eff) ;
    }
  } //end loop over jet 

  // Determine the sf from the Data & MC values
  float sf(1.) ;
  if (pMC > 0) sf = pData/pMC ;
  return sf ;
}

//This method is what we use to determine the scale factor / event weight for c-tagging.
float Selector::CalCtagWeight(std::vector<JetObj>& jets, std::string jet_main_ctagWP, std::string uncType) {

  // Get the calibration files (suggested by BTV group)
  std::string bN = "b_pt_eff_" + m_year;
  std::string cN = "c_pt_eff_" + m_year;
  std::string lN = "l_pt_eff_" + m_year;
  if (jet_main_ctagWP.find("deepJet") != std::string::npos) {
    bN = "bdj_pt_eff_" + m_year;
    cN = "cdj_pt_eff_" + m_year;
    lN = "ldj_pt_eff_" + m_year;
  }
  
  TH1D* hEff_b = (TH1D*)m_ctagEffFile->Get(bN.c_str());
  TH1D* hEff_c = (TH1D*)m_ctagEffFile->Get(cN.c_str());
  TH1D* hEff_l = (TH1D*)m_ctagEffFile->Get(lN.c_str());
  float pMC(1.);
  float pData(1.);
  
  // Go through each jet in the given event (from function parameter)
  for (std::vector<JetObj>::iterator jetIt = jets.begin(); jetIt != jets.end(); ++jetIt) {
    
    // Get the pT of the jet & its flavor
    float jetPt = (jetIt->m_lvec).Pt();
    int iBin = hEff_c->FindFixBin(jetPt) ; //return overflow bin if jetPt > max pt range
    unsigned flav = jetIt->m_flav ;
    
    // Determine what type of uncertainty we want.
    std::string uncTypeInput = "central";
    float eff = hEff_l->GetBinContent(iBin); //jet with pt > max pt range of efficinecy histogram will get the eff of overflow bins
    //if (eff <= 0) std::cout << "\n Warning: Efficiency <=0, " << eff ; //we do not want eff = 0 
    
    // Check the flavors of the jets & get the appropriate efficiency.
    BTagEntry::JetFlavor flavCode(BTagEntry::FLAV_UDSG) ;
    if (flav == 5) {
      eff = hEff_b->GetBinContent(iBin);
      flavCode = BTagEntry::FLAV_B;
    }
    if (flav == 4) {
      eff = hEff_c->GetBinContent(iBin);
      flavCode = BTagEntry::FLAV_C;
    }
    
    if (uncType == "up" || uncType == "down") uncTypeInput = uncType; //all bc and light are up together
    if (uncType == "light_up" && flav != 4 && flav != 5) uncTypeInput = "up";
    if (uncType == "light_down" && flav != 4 && flav != 5) uncTypeInput = "down";
    if (uncType == "bc_up" && (flav == 4 || flav == 5)) uncTypeInput = "up";
    if (uncType == "bc_down" && (flav == 4 || flav == 5)) uncTypeInput = "down";
    
    // Determine the scale factor for this jet
    float sf = m_ctagReader.eval_auto_bounds(
                 uncTypeInput,
                 flavCode,
                 fabs((jetIt->m_lvec).Eta()),
                 jetPt
    );
    
    // Pass c-tagging requirement
    if (jetIt->m_deepCvL > CUTS.Get<float>("jet_"+jet_main_ctagWP+"_CvL_" + m_year) && 
        jetIt->m_deepCvB > CUTS.Get<float>("jet_"+jet_main_ctagWP+"_CvB_" + m_year)) {
      pData = pData*sf*eff;
      pMC = pMC*eff;
    }
    else {
      pData = pData * (1-sf*eff);
      pMC = pMC * (1-eff);
    }
  }//end-loop over jet
  
  // Determine the sf from the Data & MC values
  float sf(1.) ;
  if (pMC > 0) sf = pData / pMC;
  return sf;
}

//Get scale factors from a list of calibration histograms h (each histo corresponds to a run periods, for example muon in 2016 has scale factors for B->F and G->H sets. w are weights for each sets 
std::vector<float> Selector::GetSF_2DHist(float x, float y, std::vector<TH2F*> h, std::vector<float> w) {
  std::vector<float> o{1,0,0}; //value, absolute error, relative error
  unsigned nX = h[0]->GetNbinsX() ;
  unsigned nY = h[0]->GetNbinsY() ;
  unsigned iX = h[0]->GetXaxis()->FindFixBin(x) ;
  unsigned iY = h[0]->GetYaxis()->FindFixBin(y) ;
  if (iX == 0 || iY == 0 || iX > nX || iY > nY) { //underflow and overflow bins
    return o ;
  }
  float sf(0) ;
  float e_sf(0) ;
  for (unsigned i = 0 ; i < h.size() ; ++i) {
    sf += w[i]*h[i]->GetBinContent(iX,iY);
    e_sf += w[i]*w[i]*h[i]->GetBinError(iX,iY)*h[i]->GetBinError(iX,iY);
  }
  e_sf = sqrt(e_sf) ;
  o[0] = sf ;
  o[1] = e_sf ;
  if(sf!=0) o[2] = e_sf/sf; //relative error
  else o[2] = 0.;

  return o ;
}

float Selector::CalEleSF(LepObj e1, LepObj e2) {
  std::vector<float> w{1.0};
  float sf = 1;
  float err = 0; //relative error treated as uncorrelated = (dy/y)^2 = sum[(dxi/xi)^2]
  float errRec = 0.0;
  float errID = 0.0;
  //reconstruction
  std::vector<TH2F*> h{m_hSF_eleRec};
  //first lepton
  std::vector<float> sfs = GetSF_2DHist(e1.m_lvec.Eta(),e1.m_lvec.Pt(),h,w) ;
  sf = sfs[0]; 
  err += sfs[2]*sfs[2]; //2=relative error
  errRec += sfs[2]*sfs[2]; //2=relative error
  //second lepton
  sfs = GetSF_2DHist(e2.m_lvec.Eta(),e2.m_lvec.Pt(),h,w);
  sf *= sfs[0] ;
  err += sfs[2]*sfs[2]; //2=relative error
  errRec += sfs[2]*sfs[2]; //2=relative error
  
  //ID
  h[0] = m_hSF_eleID;
  //first lepton
  sfs = GetSF_2DHist(e1.m_lvec.Eta(),e1.m_lvec.Pt(),h,w);
  sf *= sfs[0] ;
  err += sfs[2]*sfs[2];
  errID += sfs[2]*sfs[2];
  //second lepton
  sfs = GetSF_2DHist(e2.m_lvec.Eta(),e2.m_lvec.Pt(),h,w);
  sf *= sfs[0] ; 
  err += sfs[2]*sfs[2];
  errID += sfs[2]*sfs[2];

  err = sqrt(err);
  errRec = sqrt(errRec);
  errID = sqrt(errID);

  if (m_eleUncType == "up") sf = sf*(1+err);
  if (m_eleUncType == "down") sf = sf*(1-err);
  if (m_eleUncType == "idup") sf = sf*(1+errID);
  if (m_eleUncType == "iddown") sf = sf*(1-errID);
  if (m_eleUncType == "recup") sf = sf*(1+errRec);
  if (m_eleUncType == "recdown") sf = sf*(1-errRec);

  return sf; 
}

float Selector::CalSingleEleSF(LepObj e1) {
  std::vector<float> w{1.0};
  float sf = 1;
  float err = 0; //relative error treated as uncorrelated = (dy/y)^2 = sum[(dxi/xi)^2]
  float errRec = 0.0;
  float errID = 0.0;
  //reconstruction
  std::vector<TH2F*> h{m_hSF_eleRec};
  //first lepton
  std::vector<float> sfs = GetSF_2DHist(e1.m_lvec.Eta(),e1.m_lvec.Pt(),h,w) ;
  sf = sfs[0]; 
  err += sfs[2]*sfs[2]; //2=relative error
  errRec += sfs[2]*sfs[2]; //2=relative error
  
  //ID
  h[0] = m_hSF_eleID;
  //first lepton
  sfs = GetSF_2DHist(e1.m_lvec.Eta(),e1.m_lvec.Pt(),h,w);
  sf *= sfs[0] ;
  err += sfs[2]*sfs[2];
  errID += sfs[2]*sfs[2];

  err = sqrt(err);
  errRec = sqrt(errRec);
  errID = sqrt(errID);
  if (m_eleUncType == "up") sf = sf*(1+err);
  if (m_eleUncType == "down") sf = sf*(1-err);
  if (m_eleUncType == "recup") sf = sf*(1+errRec);
  if (m_eleUncType == "recdown") sf = sf*(1-errRec);
  if (m_eleUncType == "idup") sf = sf*(1+errID);
  if (m_eleUncType == "iddown") sf = sf*(1-errID);

  return sf; 
}

float Selector::CalMuonSF_id_iso(LepObj e1, LepObj e2) {
  float sf(1.0);
  float err(0.0);
  float errID(0.0);
  float errIso(0.0);
#ifdef MC_2016
  //////////////
  //ID
  //////////////
  //first lepton
  std::vector<float> sfs = GetSF_2DHist(e1.m_lvec.Eta(),e1.m_lvec.Pt(), m_hSF_muonID, m_muonID_w);
  sf *= sfs[0];
  err += sfs[2]*sfs[2]; //relative error
  errID += sfs[2]*sfs[2]; //relative error
  //second lepton
  sfs = GetSF_2DHist(e2.m_lvec.Eta(),e2.m_lvec.Pt(), m_hSF_muonID, m_muonID_w);
  sf *= sfs[0];
  err += sfs[2]*sfs[2];
  errID += sfs[2]*sfs[2];
  
  /////////////////
  //Iso
  /////////////////
  //first muon
  sfs = GetSF_2DHist(e1.m_lvec.Eta(),e1.m_lvec.Pt(), m_hSF_muonIso, m_muonIso_w);
  sf *= sfs[0];
  err += sfs[2]*sfs[2];
  errIso += sfs[2]*sfs[2];
  //second lepton
  sfs = GetSF_2DHist(e2.m_lvec.Eta(),e2.m_lvec.Pt(), m_hSF_muonIso, m_muonIso_w);
  sf *= sfs[0];
  err += sfs[2]*sfs[2];
  errIso += sfs[2]*sfs[2];
#endif
#if defined(MC_2017) || defined(MC_2018)
  ///////////
  //ID
  ///////////
  //first muon
  std::vector<float> sfs = GetSF_2DHist(e1.m_lvec.Pt(),fabs(e1.m_lvec.Eta()), m_hSF_muonID, m_muonID_w);
  sf *= sfs[0];
  err += sfs[2]*sfs[2];
  errID += sfs[2]*sfs[2];
  //second muon
  sfs = GetSF_2DHist(e2.m_lvec.Pt(),fabs(e2.m_lvec.Eta()), m_hSF_muonID, m_muonID_w);
  sf *= sfs[0];
  err += sfs[2]*sfs[2];
  errID += sfs[2]*sfs[2];
  
  ///////////
  //Iso
  ///////////
  //first muon
  sfs = GetSF_2DHist(e1.m_lvec.Pt(),fabs(e1.m_lvec.Eta()), m_hSF_muonIso, m_muonIso_w);
  sf *= sfs[0];
  err += sfs[2]*sfs[2];
  errIso += sfs[2]*sfs[2];
  //second muon
  sfs = GetSF_2DHist(e2.m_lvec.Pt(),fabs(e2.m_lvec.Eta()), m_hSF_muonIso, m_muonIso_w);
  sf *= sfs[0];
  err += sfs[2]*sfs[2];
  errIso += sfs[2]*sfs[2];
#endif
  err = sqrt(err);
  errID = sqrt(errID);
  errIso = sqrt(errIso);
  if (m_muonUncType == "up") sf = sf*(1+err);
  if (m_muonUncType == "down") sf = sf*(1-err);
  if (m_muonUncType == "idup") sf = sf*(1+errID);
  if (m_muonUncType == "iddown") sf = sf*(1-errID);
  if (m_muonUncType == "isoup") sf = sf*(1+errIso);
  if (m_muonUncType == "isodown") sf = sf*(1-errIso);

  return sf ;
}

float Selector::CalSingleMuonSF_id_iso(LepObj e1) {
  float sf(1.0);
  float err(0.0);
  float errID(0.0);
  float errIso(0.0);
#ifdef MC_2016
  //////////////
  //ID
  //////////////
  //first lepton
  std::vector<float> sfs = GetSF_2DHist(e1.m_lvec.Eta(),e1.m_lvec.Pt(), m_hSF_muonID, m_muonID_w);
  sf *= sfs[0];
  err += sfs[2]*sfs[2]; //relative error
  errID += sfs[2]*sfs[2]; //relative error
  
  /////////////////
  //Iso
  /////////////////
  //first muon
  sfs = GetSF_2DHist(e1.m_lvec.Eta(),e1.m_lvec.Pt(), m_hSF_muonIso, m_muonIso_w);
  sf *= sfs[0];
  err += sfs[2]*sfs[2];
  errIso += sfs[2]*sfs[2];
#endif
#if defined(MC_2017) || defined(MC_2018)
  ///////////
  //ID
  ///////////
  //first muon
  std::vector<float> sfs = GetSF_2DHist(e1.m_lvec.Pt(),fabs(e1.m_lvec.Eta()), m_hSF_muonID, m_muonID_w);
  sf *= sfs[0];
  err += sfs[2]*sfs[2];
  errID += sfs[2]*sfs[2];
  
  ///////////
  //Iso
  ///////////
  //first muon
  sfs = GetSF_2DHist(e1.m_lvec.Pt(),fabs(e1.m_lvec.Eta()), m_hSF_muonIso, m_muonIso_w);
  sf *= sfs[0];
  err += sfs[2]*sfs[2];
  errIso += sfs[2]*sfs[2];
#endif
  err = sqrt(err);
  errID = sqrt(errID);
  errIso = sqrt(errIso);
  if (m_muonUncType == "up") sf = sf*(1+err);
  if (m_muonUncType == "down") sf = sf*(1-err);
  if (m_muonUncType == "idup") sf = sf*(1+errID);
  if (m_muonUncType == "iddown") sf = sf*(1-errID);
  if (m_muonUncType == "isoup") sf = sf*(1+errIso);
  if (m_muonUncType == "isodown") sf = sf*(1-errIso);

  return sf ;
}
/*
TLorentzVector Selector::GetTrigObj(Reader* r, unsigned id, unsigned bits, float ptThresh) { 

  int id_trigObj = -1 ;
  float maxPt = ptThresh ;

  for (unsigned i = 0 ; i < *(r->nTrigObj) ; ++i) {
    if ((abs(r->TrigObj_id[i]) == id) && ((r->TrigObj_filterBits)[i] & bits)) { //has correct id, bits
      //std::cout << "\n" << (r->TrigObj_id)[i] << " " << (r->TrigObj_filterBits)[i] << "  " << bits ;
      if ((r->TrigObj_pt)[i] > maxPt) { //choose maximum pt trigger object > threshold
        id_trigObj = i ;
        maxPt = (r->TrigObj_pt)[i] ;
      }
    }
  }
  
  TLorentzVector tl_out(0,0,0,0) ;
  float mass(0.1) ;
  if (id == 11) mass = 0.0005 ;
  if (id_trigObj>=0) tl_out.SetPtEtaPhiM((r->TrigObj_pt)[id_trigObj],(r->TrigObj_eta)[id_trigObj],(r->TrigObj_phi)[id_trigObj],mass) ;
  return tl_out ;
}
*/

float Selector::CalTrigSF(int id, LepObj lep1, LepObj lep2, TLorentzVector trigObj, TH1D* h_dR1_trig, TH1D* h_dR2_trig, TH1D* h_pt1_trig, TH1D* h_pt2_trig) {
  
  float trigSF = 1.0 ;
  if (trigObj.Pt() < 0.01) return trigSF ; //empty trigger object
  float dR1 = lep1.m_lvec.DeltaR(trigObj) ;
  float dR2 = lep2.m_lvec.DeltaR(trigObj) ;
  h_dR1_trig->Fill(dR1) ;
  h_dR2_trig->Fill(dR2) ;
  if ((dR1 < dR2) && (dR1 < 0.2)) {
    h_pt1_trig->Fill(lep1.m_lvec.Pt()) ;
    if (id == 13) {
      std::vector<float> sfTmp = GetSF_2DHist(fabs(lep1.m_lvec.Eta()),lep1.m_lvec.Pt(), m_hSF_muonTrig, m_muonTrig_w) ;
      trigSF = sfTmp[0];
      if (m_muonUncType == "trigup") trigSF = sfTmp[0]*(1+sfTmp[2]);
      if (m_muonUncType == "trigdown") trigSF = sfTmp[0]*(1-sfTmp[2]);
    }
    if (id == 11) {
      //SC eta
      std::vector<float> sfTmp = GetSF_2DHist(lep1.m_scEta,lep1.m_lvec.Pt(), m_hSF_eleTrig, m_eleTrig_w);
      trigSF = sfTmp[0];
      if (m_eleUncType == "trigup") trigSF = sfTmp[0]*(1+sfTmp[2]);
      if (m_eleUncType == "trigdown") trigSF = sfTmp[0]*(1-sfTmp[2]);
    }
  }    
  if ((dR2 < dR1) && (dR2 < 0.2)) {
    h_pt2_trig->Fill(lep2.m_lvec.Pt()) ;
    if (id == 13) {
      std::vector<float> sfTmp = GetSF_2DHist(fabs(lep2.m_lvec.Eta()),lep2.m_lvec.Pt(), m_hSF_muonTrig, m_muonTrig_w) ;
      trigSF = sfTmp[0];
      if (m_muonUncType == "trigup") trigSF = sfTmp[0]*(1+sfTmp[2]);
      if (m_muonUncType == "trigdown") trigSF = sfTmp[0]*(1-sfTmp[2]);
    }
    if (id == 11) {
      std::vector<float> sfTmp = GetSF_2DHist(lep2.m_scEta,lep2.m_lvec.Pt(), m_hSF_eleTrig, m_eleTrig_w) ;
      trigSF = sfTmp[0];
      if (m_eleUncType == "trigup") trigSF = sfTmp[0]*(1+sfTmp[2]);
      if (m_eleUncType == "trigdown") trigSF = sfTmp[0]*(1-sfTmp[2]);
    }
  }
  
  return trigSF ;
}

float Selector::CalTrigSF_singleLepton(int id, LepObj lep1, TLorentzVector trigObj) {
  
  float trigSF = 1.0 ;
  if (trigObj.Pt() < 0.01) return trigSF ; //empty trigger object
  float dR1 = lep1.m_lvec.DeltaR(trigObj) ;
  if ((dR1 < 0.2)) {
    if (id == 13) {
      std::vector<float> sfTmp = GetSF_2DHist(fabs(lep1.m_lvec.Eta()),lep1.m_lvec.Pt(), m_hSF_muonTrig, m_muonTrig_w) ;
      trigSF = sfTmp[0];
      if (m_muonUncType == "trigup") trigSF = sfTmp[0]*(1+sfTmp[2]);
      if (m_muonUncType == "trigdown") trigSF = sfTmp[0]*(1-sfTmp[2]);
    }
    if (id == 11) {
      //SC eta
      std::vector<float> sfTmp = GetSF_2DHist(lep1.m_scEta,lep1.m_lvec.Pt(), m_hSF_eleTrig, m_eleTrig_w);
      trigSF = sfTmp[0];
      if (m_eleUncType == "trigup") trigSF = sfTmp[0]*(1+sfTmp[2]);
      if (m_eleUncType == "trigdown") trigSF = sfTmp[0]*(1-sfTmp[2]);
    }
  }    
  
  return trigSF ;
}
/*
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
unsigned Selector::MatchGenLep(Reader* r, LepObj lep, int pdgId) {
  float dRmin(1000) ;
  int indO(-1) ;
  for (unsigned i = 0 ; i < *(r->nGenDressedLepton) ; ++i) {
    if (pdgId == fabs((r->GenDressedLepton_pdgId)[i])) {
      TLorentzVector gLep_lvec ;
      gLep_lvec.SetPtEtaPhiM((r->GenDressedLepton_pt)[i], (r->GenDressedLepton_eta)[i], (r->GenDressedLepton_phi)[i], (r->GenDressedLepton_mass)[i]) ;
      float dRtmp = lep.m_lvec.DeltaR(gLep_lvec) ;
      if (dRtmp < dRmin) {
        dRmin = dRmin ;
        indO = i ;
      }
    }
  }
  return indO ;
}

#endif

float Selector::MuonRcSF(Reader* r, LepObj lep, int pdgId) {
  float sf(1.) ;
  sf = m_rc.kScaleDT(lep.m_charge, lep.m_lvec.Pt(), lep.m_lvec.Eta(), lep.m_lvec.Phi());
#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
  sf = 1 ;
  //int gLepInd = MatchGenLep(r, lep, 13) ;
  int gLepInd = r->Muon_genPartIdx[lep.m_idx] ;
  if(gLepInd >= 0) {
    sf = m_rc.kSpreadMC(lep.m_charge, lep.m_lvec.Pt(), lep.m_lvec.Eta(), lep.m_lvec.Phi(),(r->GenPart_pt[gLepInd]));
  }
  else {
    sf = m_rc.kSmearMC(lep.m_charge, lep.m_lvec.Pt(), lep.m_lvec.Eta(), lep.m_lvec.Phi(), (r->Muon_nTrackerLayers)[lep.m_idx], m_rand->Rndm());
  }
#endif
  return sf ;
}*/

