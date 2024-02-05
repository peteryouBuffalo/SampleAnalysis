#ifndef Obj_cxx
#define Obj_cxx

//The objects are wrapper so that the analysis is independent on variable name on the ntuple, only when the objects are created such dependences occur 

#include "TLorentzVector.h"

//-------------------------------objs------------------------------
class LepObj { // Lepton

  public:
    
    // Constructor & Deconstructor
    LepObj(float pt, float eta, float scEta, float phi, float mass, unsigned idx, int charge, 
           float iso) : m_scEta(scEta), m_idx(idx), m_charge(charge), m_iso(iso) { 
             m_lvec.SetPtEtaPhiM(pt, eta, phi, mass) ; 
    } ;

    virtual ~LepObj() {} ;

    // Variables
    TLorentzVector m_lvec ;  //4-vector
    
    float m_scEta ;  // eta value
    unsigned m_idx ; // index
    int m_charge ;   // charge
    float m_iso ;    // isolation
    
} ;


class JetObj { // Jets
  
  public:
    
    // Constructor & Deconstructor
    JetObj(float pt, float eta, float phi, float mass, unsigned flav, 
           float deepCSV, float deepJet, float PUjetID) : m_flav(flav), m_deepCSV(deepCSV), m_deepJet(deepJet), m_puid(PUjetID)  {
              m_lvec.SetPtEtaPhiM(pt, eta, phi, mass) ; 
    } ;

    virtual ~JetObj() {} ;

    // IsLepton - checks to make sure our jet isn't actually a lepton
    bool IsLepton(std::vector<LepObj>& leps) {
      float minDr = 1000 ;
      for (std::vector<LepObj>::iterator it = leps.begin() ; it != leps.end() ; ++it) {
        float dRtmp = m_lvec.DeltaR(it->m_lvec) ;
        if (dRtmp < minDr) minDr = dRtmp ;
      }
      return minDr <= 0.4 ;
    }

    // SetSV - determines the best secondary vertex for the jet
    void SetSV(std::vector<TLorentzVector>& sv) {
      float maxPt = -1;
      m_svIdx = -1;
      m_mSV = -1;
      for (unsigned isv = 0 ; isv < sv.size() ; ++isv) {
        float dRtmp = m_lvec.DeltaR(sv[isv]);
        if (dRtmp <= 0.4 && sv[isv].Pt() > maxPt) {
          m_svIdx = isv; 
          m_mSV = sv[isv].M();
        }
      }
    } ;

    // Variables
    TLorentzVector m_lvec ; // 4-vector
    unsigned m_flav ;       // jet flavor
    float m_deepCSV ;       // DeepCSV value
    float m_deepJet ;       // DeepJet value
    unsigned m_svIdx ;      // SV index
    float m_mSV;            // SV mass
    float m_puid;           // PU ID

} ;


class ZObj { // Z Boson

  public:

    // Constructor & Deconstructor
    ZObj(std::vector<JetObj> jetlist) : m_jets(jetlist) {
      for (size_t idx = 0; idx < jetlist.size(); ++idx) {
        m_lvec += jetlist.at(idx).m_lvec;
      }
    }
  
    virtual ~ZObj() {} ;

    // Methods
    size_t nJets() {
      return m_jets.size();
    }

    JetObj getJet(Int_t idx) {
      return m_jets.at(idx);
    }

    float DeltaR() {
      return m_jets.at(0).m_lvec.DeltaR(m_jets.at(1).m_lvec);
    }
 
    float Dphi(float x, float y) { return (x-y); }
    float DphiC(double x) {
      if ((x <= TMath::Pi() && x >= 0) or (x<0 && x > -TMath::Pi())) return x;
      else if (x >= TMath::Pi()) return DphiC(x-2*TMath::Pi());
      else if (x < -TMath::Pi()) return DphiC(x+2*TMath::Pi());
      return 0; // make sure we return a value
    }
    
    float DPhi() {
      float dphi = Dphi(m_jets.at(0).m_lvec.Phi(), m_jets.at(1).m_lvec.Phi());
      return DphiC(dphi);
    }

    // Variables
    TLorentzVector m_lvec;      // 4-vector
    std::vector<JetObj> m_jets; // list of jets that compose Z boson

} ;


class HObj { // Higgs boson

  public:

    // Constructor & Deconstructor
    HObj(std::vector<JetObj> jetlist) : m_jets(jetlist) {
      for (size_t idx = 0; idx < jetlist.size(); ++idx) {
        m_lvec += jetlist.at(idx).m_lvec;
      }
    }

    virtual ~HObj() {} ;

    // Methods
    size_t nJets() {
      return m_jets.size();
    }

    JetObj getJet(Int_t idx) {
      return m_jets.at(idx);
    }

    float DeltaR() {
      return m_jets.at(0).m_lvec.DeltaR(m_jets.at(1).m_lvec);
    }

    float Dphi(float x, float y) { return (x-y); }
    float DphiC(double x) {
      if ((x <= TMath::Pi() && x >= 0) or (x<0 && x > -TMath::Pi())) return x;
      else if (x >= TMath::Pi()) return DphiC(x-2*TMath::Pi());
      else if (x < -TMath::Pi()) return DphiC(x+2*TMath::Pi());
      return 0; // make sure we return a value
    }

    float DPhi() {
      float dphi = Dphi(m_jets.at(0).m_lvec.Phi(), m_jets.at(1).m_lvec.Phi());
      return DphiC(dphi);
    }


    // Variables
    TLorentzVector m_lvec;      // 4-vector
    std::vector<JetObj> m_jets; // list of jets that compose Higgs boson
};

class GenObj { // Generator Objects

  public:
    
    // Constructor & Deconstructor
    GenObj(int pdgID, float pt, float eta, float phi, float mass, unsigned idx,
           int mother, int status) : m_pdgID(pdgID), m_idx(idx), m_motherID(mother),
           m_status(status) {
             m_lvec.SetPtEtaPhiM(pt, eta, phi, mass);
    };
    virtual ~GenObj() {};

    // Variables
    TLorentzVector m_lvec; //4-vector
    int m_pdgID;
    unsigned m_idx;
    int m_motherID;
    int m_status;
    
    //Methods
    static std::vector<JetObj> get_proper_jets(std::vector<JetObj>& arr, GenObj* gen_obj,
    float dRcut) {
      std::vector<JetObj> constituents;
      for (size_t i = 0; i < arr.size(); ++i) {
        float dR = gen_obj->m_lvec.DeltaR(arr.at(i).m_lvec);
        if (dR < dRcut) constituents.push_back(arr.at(i));
      }
 
      return constituents;
    }
};


#endif
