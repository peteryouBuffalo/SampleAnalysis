#ifndef Obj_cxx
#define Obj_cxx

/******************************************************************************
 ####### ######        # #######  #####  #######  #####  
 #     # #     #       # #       #     #    #    #     # 
 #     # #     #       # #       #          #    #       
 #     # ######        # #####   #          #     #####  
 #     # #     # #     # #       #          #          # 
 #     # #     # #     # #       #     #    #    #     # 
 ####### ######   #####  #######  #####     #     #####  
*******************************************************************************

Description: This is a series of classes that allow us to reconstruct objects.
             In the general form, each was its own object, but we have some
             common items between them, so all objects should extend from the
             base 'GenObj' class.

******************************************************************************/

// == [10] Included Libraries =================================================

// Standard C++ libraries/classes
#include <math.h>

// ROOT classes
#include "TLorentzVector.h"

// == [20] General Object class ===============================================

class GenObj
{

  /* PUBLIC: elements we want easy access to */
  public:
    
    // ============================
    // Constructor & Deconstructor
    // ============================
    GenObj(float pt, float eta, float phi, float mass) 
    {
      SetPtEtaPhiM(pt, eta, phi, mass);
    };
    virtual ~GenObj() { };

    // ========
    // Methods
    // ========
    float Pt() { return m_lvec.Pt(); }
    float Eta() { return m_lvec.Eta(); }
    float M() { return m_lvec.M(); }
    float Phi() { return m_lvec.Phi(); }
    void AddVec4(TLorentzVector vec) { m_lvec += vec; }
    void SetPtEtaPhiM(float pt, float phi, float eta, float mass)
    { m_lvec.SetPtEtaPhiM(pt, phi, eta, mass); }
    TLorentzVector Vec4() { return m_lvec; } 
 
  /* PRIVATE: elements we want to keep isolated to the class*/
  //private:
    
    // ==========
    // Variables
    // ==========
    TLorentzVector m_lvec;  // particle 4-vector

};//end-GenObj

// == [30] Lepton Object class ================================================

class LepObj : public GenObj
{

  /* PUBLIC: elements we want easy access to */
  public:
    
    // ============================
    // Constructor & Deconstructor
    // ============================
    LepObj(float pt, float eta, float scEta, float phi, float mass,
      unsigned idx, int charge, float iso) : GenObj(pt, eta, phi, mass)
    {
      m_scEta = scEta;
      m_idx = idx;
      m_charge = charge;
      m_iso = iso;
    };
    virtual ~LepObj() {};
    
    // ========
    // Methods
    // ========
    float SCEta() { return m_scEta; }
    float Idx() { return m_idx; }
    float Charge() { return m_charge; }
    float Iso() { return m_iso; }
    
  /* PRIVATE: variables we want to keep access within the class */
  //private:
    
    float m_scEta;       // SC eta value
    float m_idx;         // index
    float m_charge;      // charge
    float m_iso;         // isolation
    float m_objIdx = 1;

};//end-LepObj

// == [40] Jet Object class ===================================================

class JetObj : public GenObj 
{

  /* PUBLIC: elements that we want easy access to */
  public:
    
    // ============================
    // Constructor & Deconstructor
    // ============================
    JetObj(float pt, float eta, float phi, float mass, unsigned flav,
      float deepCSV, float PUjetID) : GenObj(pt, eta, phi, mass)
    {
      m_flav = flav;
      m_deepCSV = deepCSV;
      m_puid = PUjetID;
      
      m_pt = pt;
      m_mass = mass;
    };
    virtual ~JetObj() {};

    // The following structs allow us to compare jets 
    // by a number of different attributes.
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

    // ========
    // Methods
    // ========

    /* IsLepton - checks to make sure our jet isn't actually a lepton */
    bool IsLepton(std::vector<LepObj>& leps) 
    {
      float minDr = 1000;
      for (std::vector<LepObj>::iterator it = leps.begin();
      it != leps.end(); ++it) {
        float dRtmp = Vec4().DeltaR(it->Vec4());
        if (dRtmp < minDr) minDr = dRtmp;
      }
      return minDr <= 0.4;
    };

    /* SetSV - determines the best secondary vertex for the jet */
    void SetSV(std::vector<TLorentzVector>& sv) 
    {
      float maxPt = -1;
      m_svIdx = -1;
      m_mSV = -1;
      for (unsigned isv = 0; isv < sv.size(); ++isv) 
      {
        float dRtmp = Vec4().DeltaR(sv[isv]);
        if (dRtmp <= 0.4 && sv[isv].Pt() > maxPt)
        {
          m_svIdx = isv;
          m_mSV = sv[isv].M();
        }
      }
    };
    
    /* GetPComp - get the i'th component of the 3-vector momentum */
    float GetPComp(int idx)
    {
      if (idx == 0) return Vec4().Px();
      if (idx == 1) return Vec4().Py();
      if (idx == 2) return Vec4().Pz();
      if (idx == 3) return Vec4().P();
      return Pt();
    };

    /* StoreRegInfo - store the regression variables */
    void StoreRegInfo(float bRC, float bRR, float cRC, float cRR)
    {
      m_bRegCorr = bRC; m_bRegRes = bRR;
      m_cRegCorr = cRC; m_cRegRes = cRR;
    };

    /* ApplyRegression - applies the corrects given a jet type
      |abs(flav) = 4 -> c, abs(flav) = 5 -> b| */
    void ApplyRegression(float flav, bool applyToMass=true)
    {
      float pt = Pt();
      m_pt = pt;
      if (abs(flav) == 4) pt *= m_cRegCorr;
      else if (abs(flav) == 5) pt *= m_bRegCorr;
      m_ptJEC = pt;

      float phi = Phi();
      float eta = Eta();
      float mass = M();
      
      if (applyToMass) 
      {
        m_mass = mass;
        if (abs(flav) == 4) mass *= m_cRegCorr;
        else if (abs(flav) == 5) mass *= m_bRegCorr;
        m_massJEC = mass;
      }

      SetPtEtaPhiM(pt, phi, eta, mass);
    };

    /* Setters & Getters */
    void SetIdx(int idx) { m_Idx = idx; }
    void SetIdxAll(int idx) { m_IdxAll = idx; }
    int IdxAll() { return m_IdxAll; }
    int Idx() { return m_Idx; }

  /* PRIVATE: variables we want less access to */
  //private: 
    int m_flav;
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

};//end-JetObj

// == [50] Boson Object class =================================================

class BosonObj : public GenObj
{

  /* PUBLIC */
  public:
    
    // ============================
    // Constructor & Deconstructor
    // ============================
    BosonObj(std::vector<JetObj> jetlist) : GenObj(0, 0, 0, 0) 
    {
      for (unsigned int idx = 0; idx < jetlist.size(); ++idx) 
      {
        AddVec4(jetlist.at(idx).Vec4());
        m_jets.push_back(jetlist.at(idx));
      }
    };
    virtual ~BosonObj() {};

    // ========
    // Methods
    // ========
    JetObj getJet(int idx) { return m_jets[idx]; }
    size_t nJets() { return m_jets.size(); }
    float DeltaR()
    {
      return m_jets[0].Vec4().DeltaR(m_jets[1].Vec4());
    }
    float Dphi(float x, float y) { return (x-y); }
    float DphiC(double x)
    {
      if ((x <= TMath::Pi() && x >= 0) or (x<0 && x > -TMath::Pi())) return x;
      else if (x >= TMath::Pi()) return DphiC(x-2*TMath::Pi());
      else if (x < -TMath::Pi()) return DphiC(x+2*TMath::Pi());
      return x;
    }
    
    float DPhi() { 
      float dphi = Dphi(m_jets[0].Phi(), m_jets[1].Phi());
      return DphiC(dphi);
    }
    
    void SetBosonID(int id) { bosonID = id; }

  /* PRIVATE */
  private:
    // ==========
    // Variables
    // ==========
    std::vector<JetObj> m_jets;
    int bosonID = 0;

};//end-BosonObj

// == [61] Z Boson Object class ===============================================

class ZObj : public BosonObj
{
  public:
    ZObj(std::vector<JetObj> jetlist) : BosonObj(jetlist)
    { SetBosonID(23); }
    virtual ~ZObj() {};
};

// == [62] Higgs Boson Object class ===========================================

class HiggsObj : public BosonObj
{
  public:
    HiggsObj(std::vector<JetObj> jetlist) : BosonObj(jetlist)
    { SetBosonID(25); }
    virtual ~HiggsObj() {};
};

// == [70] DHZ Object class ===================================================

class DHZObj
{

  public:
    // ============================
    // Constructor & Deconstructor
    // ============================
    DHZObj(std::vector<JetObj> jets, int p0id0, int p0id1, int p1id0, int p1id1)
    {
      m_jets = jets;
      std::vector<int> idx0 {p0id0, p0id1};
      m_indices.push_back(idx0);
      std::vector<int> idx1 {p1id0, p1id1};
      m_indices.push_back(idx1);
      calculate_d(); 
    };
    virtual ~DHZObj() {}

    // Overload the < and > operators. This allows us to sort the objects
    // by their distance values in ascending and descending order.
    bool operator>(const DHZObj &r) const
    { return m_distance > r.m_distance; }

    bool operator<(const DHZObj &r) const
    { return m_distance < r.m_distance; }

    // ========
    // Methods 
    // ========
    void calculate_d() 
    {
      // Create the two possible candidates
      TLorentzVector v0 = (m_jets[m_indices[0][0]].m_lvec +
        m_jets[m_indices[0][1]].m_lvec);
      TLorentzVector v1 = (m_jets[m_indices[1][0]].m_lvec +
        m_jets[m_indices[1][1]].m_lvec);
 
      // Make sure that H is the higher pT pairing
      if (v0.Pt() > v1.Pt()){ 
        m_Hvec = v0; m_Zvec = v1; 
        m_Hidx = 0; m_Zidx = 1;
      }
      else { 
        m_Hvec = v1; m_Zvec = v0;
        m_Hidx = 1; m_Zidx = 0;
      }

      // Reconstruct the masses & calcualte their distances
      float m0 = m_Hvec.M(); float m1 = m_Zvec.M();
      float numerator = fabs(m0 - k*m1);
      float denominator = sqrt(1 + pow(k,2));
      m_distance = numerator / denominator;
    };

    // Setters & Getters
    float HPt() { return m_Hvec.Pt(); }
    float HM() { return m_Hvec.M(); }
    float ZPt() { return m_Zvec.Pt(); }
    float ZM() { return m_Zvec.M(); }
    float DPhi() { return m_Zvec.DeltaPhi(m_Hvec); }
    float Hidx(int i) { return m_indices[m_Hidx][i]; }
    float Zidx(int i) { return m_indices[m_Zidx][i]; }

    bool H_has_cjet0(float desired_CvL, float desired_CvB) {
      float cvl0 = m_jets[Hidx(0)].m_deepCvL;
      float cvb0 = m_jets[Hidx(0)].m_deepCvB;
      return (cvl0 > desired_CvL && cvb0 > desired_CvB);
    }

    bool H_has_cjet1(float desired_CvL, float desired_CvB) {
      float cvl1 = m_jets[Hidx(1)].m_deepCvL;
      float cvb1 = m_jets[Hidx(1)].m_deepCvB;
      return (cvl1 > desired_CvL && cvb1 > desired_CvB);
    }
 
    bool H_has_cjets(float desired_CvL, float desired_CvB) {
      bool pass0 = H_has_cjet0(desired_CvL, desired_CvB);
      bool pass1 = H_has_cjet1(desired_CvL, desired_CvB);
      return pass0 && pass1;
    };
 
    bool H_has_bjet0(float desired_BvL) {
      float csv0 = m_jets[Hidx(0)].m_deepCSV;
      return (csv0 > desired_BvL);  
    };

    bool H_has_bjet1(float desired_BvL) {
      float csv1 = m_jets[Hidx(1)].m_deepCSV;
      return (csv1 > desired_BvL);
    };

    bool Z_has_bjet0(float desired_BvL) {
      float csv0 = m_jets[Zidx(0)].m_deepCSV;
      return (csv0 > desired_BvL);
    }
    
    bool Z_has_bjet1(float desired_BvL) {
      float csv1 = m_jets[Zidx(1)].m_deepCSV;
      return (csv1 > desired_BvL);
    }

    bool Z_has_bjets(float desired_BvL) {
      bool pass0 = Z_has_bjet0(desired_BvL);
      bool pass1 = Z_has_bjet1(desired_BvL);
      return pass0 && pass1;
    };

    // Variables
    float k = 125.0 / 91.0;
    std::vector<std::vector<int>> m_indices;
    int m_Zidx, m_Hidx;
    std::vector<JetObj> m_jets;
    float m_distance;
    TLorentzVector m_Hvec;
    TLorentzVector m_Zvec;
};

#endif

// == [NN] END OF FILE ========================================================
