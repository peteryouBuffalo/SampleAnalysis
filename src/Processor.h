//////////////////////////////////////////////////////////
//Processor is essentially a TSelector that will run through the
//  data and execute the registered Selection in SlaveBegin, Process
//  and SlaveTerminate and thus allows to run several selection in one
//  pass over the data.
//////////////////////////////////////////////////////////

#ifndef Processor_h
#define Processor_h

#include "Reader.h"
#include "Selector.h"


//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondTools/BTau/interface/BTagCalibrationReader.h"

// Headers needed by this particular selector
//For processing sequence, please see
//https://root.cern.ch/developing-tselector

class Processor : public Reader {
public :

   Processor(TTree * /*tree*/ t =0) ;
   virtual ~Processor() { } ;
   virtual void    Begin(TTree *tree) ;
   virtual void    SlaveBegin(TTree *tree) ;
   virtual void    Init(TTree *tree) ;
   virtual Bool_t  Process(Long64_t entry) ;
   virtual void    SlaveTerminate() ;
   virtual void    Terminate() ;
   
   void addSelector(Selector* sel) {m_selectors.push_back(sel) ; } ;
   void setOutPutFileName(std::string name) {m_outputfilename=name; } ;
   void SetDataInfo(bool isData, std::string year) {m_isData = isData; m_year = year ; };

   ClassDef(Processor,0);

private:

  std::vector<Selector*> m_selectors ;
  std::string m_outputfilename ;
  bool m_isData ;
  std::string m_year ;

};

#endif
