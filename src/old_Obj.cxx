#ifndef Obj_cxx
#define Obj_cxx

//== Obj class ================================================================
// This file handles multiple classes that allow us to reconstruct objects. 
// In the general form, each was its own object, but we have some common items
// between each, so I'm rewriting this to have a base Obj class.
//=============================================================================

//== Include statements =======================================================

#include "TLorentzVector.h" // ROOT object for four-vectors
#include <math.h>           // General math functions

/******************************************************************************
* General Object                                                              *
******************************************************************************/
class GenObj {

  public:
  
    // Constructor & Deconstructor
    GenObj(float pt, float eta, float phi, float mass) {
      m_lvec.SetPtEtaPhiM(pt, eta, phi, mass);
      //std::cout << "GenObj created." << std::endl;
    };
    virtual ~GenObj() { };
    
    // Methods
    float Pt() { return m_lvec.Pt(); }
    float Eta() { return m_lvec.Eta(); }
    float M() { return m_lvec.M(); }
    float Phi() { return m_lvec.Phi(); } 
 
    // Variables
    TLorentzVector m_lvec; //4-vector

};

/******************************************************************************
* Lepton Object                                                               *
******************************************************************************/
class LepObj : public GenObj {

  public:
  
    // Constructor & Deconstructor
    LepObj(float pt, float eta, float scEta, float phi, float mass, unsigned idx,
      int charge, float iso) : GenObj(pt, eta, phi, mass) {
      m_scEta = scEta;
      m_idx = idx;
      m_charge = charge;
      m_iso = iso;  
      //std::cout << ">>> Lep Obj created." << std::endl;
    };
    virtual ~LepObj() {};
 
    // Methods
    float SCEta() { return m_scEta; }
    float Idx() { return m_idx; }
    float Charge() { return m_charge; }
    float Iso() { return m_iso; }
  
    // Variables
    float m_scEta;  // SC eta value
    unsigned m_idx; // index
    int m_charge;   // charge
    float m_iso;    // isolation
    int m_objIdx = 1;
  

};

/******************************************************************************
* Jet Object                                                                  *
******************************************************************************/
class JetObj : public GenObj {

  public:

    // Constructor & Deconstructor
    JetObj(float pt, float eta, float phi, float mass, unsigned flav, 
      float deepCSV, float PUjetID) : GenObj(pt, eta, phi, mass) {
      m_flav = flav;
      m_deepCSV = deepCSV;
      m_puid = PUjetID;

      m_pt = pt;
      m_mass = mass;
    };
    virtual ~JetObj() {};

    // This is a struct that allows us to compare jets by a number
    // of different attributes.
    struct JetCompBtag {
      JetCompBtag() {}
      bool operator()(const JetObj& j0, const JetObj& j1) const 
      { return j0.m_deepCSV > j1.m_deepCSV; }
    };

    struct JetCompCtag {
      JetCompCtag() {}
      bool operator() (const JetObj& j0, const JetObj& j1) const
      { return j0.m_deepCvL > j1.m_deepCvL; }
    };

    struct JetCompPt {
      JetCompPt() {}
      bool operator()(const JetObj& j0, const JetObj& j1) const
      { return j0.m_lvec.Pt() > j1.m_lvec.Pt(); }
    };

    // Comparison Methods - override the > operator so we can sort
    // the jets by pT (in descending order).
    //bool operator>(const JetObj &r) const
    //{ return m_lvec.Pt() > r.m_lvec.Pt(); }
    
    // IsLepton - checks to make sure our jet isn't actually a lepton
    bool IsLepton(std::vector<LepObj>& leps) {
      float minDr = 1000;
      for (std::vector<LepObj>::iterator it = leps.begin();
      it != leps.end(); ++it) {
        float dRtmp = m_lvec.DeltaR(it->m_lvec);
        if (dRtmp < minDr) minDr = dRtmp;
      }
      return minDr <= 0.4;
    };

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
    };

    // GetPComp - get the i'th component of the 3-vector momentum
    float GetPComp(int idx) {
      if (idx == 0) return m_lvec.Px();
      if (idx == 1) return m_lvec.Py();
      if (idx == 2) return m_lvec.Pz();
      if (idx == 3) return m_lvec.P();
      return m_lvec.Pt();
    }

    // StoreRegInfo - store the regression variables
    void StoreRegInfo(float bRC, float bRR, float cRC, float cRR) {
      m_bRegCorr = bRC; m_bRegRes = bRR;
      m_cRegCorr = cRC; m_cRegRes = cRR;
    }

    // ApplyRegression - applies the corrects given a jet type
    // abs(flav) = 4 -> c, abs(flav) = 5 -> b
    void ApplyRegression(float flav, bool applyToMass=true) {

      float pt = m_lvec.Pt();    
      m_pt = pt;
      if (abs(flav) == 4) pt *= m_cRegCorr;
      else if (abs(flav) == 5) pt *= m_bRegCorr;
      m_ptJEC = pt;

      // == apply RegRes here ==
      
      float phi = m_lvec.Phi();
      float eta = m_lvec.Eta();
      float mass = m_lvec.M();

      if (applyToMass) {
        m_mass = mass;
        if (abs(flav) == 4) mass *= m_cRegCorr;
        else if (abs(flav) == 5) mass *= m_bRegCorr;
        m_massJEC = mass;
      }

      m_lvec.SetPtEtaPhiM(pt, phi, eta, mass);

    }

    // SetIdx - set the index within the list
    void SetIdx(int idx) { m_Idx = idx; }
    void SetIdxAll(int idx) { m_IdxAll = idx; }
    int IdxAll() { return m_IdxAll; }
    int Idx() { return m_Idx; }

    // Variables
    int m_flav;  // jet flavor
    float m_deepCSV;  // b-tagging
    float m_deepCvL;  // c-tagging
    float m_deepCvB;  // C vs B score
    unsigned m_svIdx; // SV index
    float m_mSV;      // SV mass
    float m_puid;     // PU ID
    int m_Idx;        // index within list
    int m_IdxAll;     // index within event list
    int m_genJetIdx;  // index for gen jet associated with this jet

    float m_pt;       // pt before correction
    float m_ptJEC;    // pt after correction
    float m_mass;     // mass before correction
    float m_massJEC;  // mass after correction

    float m_bRegCorr; // pt correction for b-jet energy regression
    float m_cRegCorr; // pt correction for c-jet energy regression
    float m_bRegRes;  // res on pt corrected with b-jet regression
    float m_cRegRes;  // res on pt corrected with c-jet regression
};

/******************************************************************************
* Boson Object                                                                *
******************************************************************************/
class BosonObj : public GenObj {

  public:

    // Constructor & Deconstructor
    BosonObj(std::vector<JetObj> jetlist) : GenObj(0, 0, 0, 0) {
      for (unsigned int idx = 0; idx < jetlist.size(); ++idx) {
        m_lvec += jetlist.at(idx).m_lvec;
        m_jets.push_back(jetlist.at(idx));
      }
    }

    BosonObj(JetObj jet) : GenObj(0, 0, 0, 0) {
      m_jets.push_back(jet);
      m_lvec = jet.m_lvec;
    }

    virtual ~BosonObj() {};

    // Methods
    JetObj getJet(int idx) { return m_jets[idx]; }
    size_t nJets() { return m_jets.size(); }
    float DeltaR() {
      return m_jets[0].m_lvec.DeltaR(m_jets[1].m_lvec);
    }
    float Dphi(float x, float y) { return (x-y); }
    float DphiC(double x) {
      if ((x <= TMath::Pi() && x >= 0) or (x<0 && x > -TMath::Pi())) return x;
      else if (x >= TMath::Pi()) return DphiC(x-2*TMath::Pi());
      else if (x < -TMath::Pi()) return DphiC(x+2*TMath::Pi());
      return x;
    }

    float DPhi() { 
      float dphi = Dphi(m_jets[0].m_lvec.Phi(), m_jets[1].m_lvec.Phi());
      return DphiC(dphi);
    }

    // Variables
    std::vector<JetObj> m_jets;
    int bosonID = 0;
};

/******************************************************************************
* Z boson Object                                                              *
******************************************************************************/
class ZObj : public BosonObj {
  
  public:

    // Constructor & Deconstructor
    ZObj(std::vector<JetObj> jetlist) : BosonObj(jetlist) {
      bosonID = 23;
    }
    virtual ~ZObj() {};
};

/******************************************************************************
* Higgs boson Object                                                          *
******************************************************************************/
class HObj : public BosonObj {

  public:

    // Constructor & Deconstructor
    HObj(std::vector<JetObj> jetlist) : BosonObj(jetlist) {
      bosonID = 25;
    }
    virtual ~HObj() {};

};

/******************************************************************************
* D_HZ Algorithm Object                                                       *
******************************************************************************/
class DHZObj {

  public:

    // Constructor & Deconstructor
    DHZObj(std::vector<JetObj> jets, int h0, int h1, int z0, int z1) {
      m_jets = jets;
      m_hIdx0 = h0;  m_hIdx1 = h1;
      m_zIdx0 = z0;  m_zIdx1 = z1;
      calculate_d();
    };
    virtual ~DHZObj() {}

    // Overload the < and > operators. This allows us to sort the objects
    // by their distance values in ascending and descending order.
    bool operator>(const DHZObj &r) const
    { return m_d > r.m_d; }

    bool operator<(const DHZObj &r) const
    { return m_d < r.m_d; }

    // Methods - calculate the distance value
    void calculate_d() {

      // Reconstruct the 4-vectors of the bosons
      m_Hvec = m_jets[m_hIdx0].m_lvec + m_jets[m_hIdx1].m_lvec;
      m_Zvec = m_jets[m_zIdx0].m_lvec + m_jets[m_zIdx1].m_lvec;

      // Reconstruct the masses & calcualte the distance values.
      float m0 = m_Hvec.M(); float m1 = m_Zvec.M();
      float numerator = fabs(m0 - k*m1);
      float denominator = sqrt(1 + pow(k,2));
      m_d = numerator / denominator;

    };

    // Methods - Setters & Getters
    float HPt() { return m_Hvec.Pt(); }
    float HM() { return m_Hvec.M(); }
    float ZPt() { return m_Zvec.Pt(); }
    float ZM() { return m_Zvec.M(); }
    float DPhi() { return m_Zvec.DeltaPhi(m_Hvec); }

    // Methods - check if the jets pass our criteria
    bool H_has_cjet0(float desired_CvL, float desired_CvB) {
      float cvl0 = m_jets[m_hIdx0].m_deepCvL;
      float cvb0 = m_jets[m_hIdx0].m_deepCvB;
      return (cvl0 > desired_CvL && cvb0 > desired_CvB);
    }

    bool H_has_cjet1(float desired_CvL, float desired_CvB) {
      float cvl1 = m_jets[m_hIdx1].m_deepCvL;
      float cvb1 = m_jets[m_hIdx1].m_deepCvB;
      return (cvl1 > desired_CvL && cvb1 > desired_CvB);
    }
 
    bool H_has_cjets(float desired_CvL, float desired_CvB) {
      bool pass0 = H_has_cjet0(desired_CvL, desired_CvB);
      bool pass1 = H_has_cjet1(desired_CvL, desired_CvB);
      return pass0 && pass1;
    };
 
    bool H_has_bjet0(float desired_BvL) {
      float csv0 = m_jets[m_hIdx0].m_deepCSV;
      return (csv0 > desired_BvL);  
    };

    bool H_has_bjet1(float desired_BvL) {
      float csv1 = m_jets[m_hIdx0].m_deepCSV;
      return (csv1 > desired_BvL);
    };

    bool Z_has_bjet0(float desired_BvL) {
      float csv0 = m_jets[m_zIdx0].m_deepCSV;
      return (csv0 > desired_BvL);
    }
    
    bool Z_has_bjet1(float desired_BvL) {
      float csv1 = m_jets[m_zIdx1].m_deepCSV;
      return (csv1 > desired_BvL);
    }

    bool Z_has_bjets(float desired_BvL) {
      bool pass0 = Z_has_bjet0(desired_BvL);
      bool pass1 = Z_has_bjet1(desired_BvL);
      return pass0 && pass1;
    };

    // Variables
    float k = 125.0 / 91.0;
    int m_hIdx0, m_hIdx1;          // indices for jets selected to Higgs
    int m_zIdx0, m_zIdx1;          // indices for jets selected to Z boson
    std::vector<JetObj> m_jets;    // list of jets that we've selected
    float m_d;                     // distance calculated
    TLorentzVector m_Hvec;         // 4-vector for Higgs
    TLorentzVector m_Zvec;         // 4-vector for Z boson
};

#endif
