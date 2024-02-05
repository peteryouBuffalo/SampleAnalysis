#ifndef SelectionTemplate_h
#define SelectionTemplate_h

//Template for creating new selection

#include "Reader.h"
#include "Selector.h"
#include "Plots.cxx"

//The Selection does not have Begin, since we do not have anything to do at the begining (overall) at client
//The histograms, ..., are booked and added to output list at SlaveBegin
//We need to have terminate since we might want to do overall tasks related to this selection only. At termination, all inforamtion from slaves is added. Example task is cutflow for this selection, which need information from all slaves.
//SlaveTerminate and Terminate are here for reservation now. histograms are written back in SlaveTerminate of Processor class

class SelectionTemplate : public Selector {
public:
  SelectionTemplate(bool isData=false) ;
  virtual ~SelectionTemplate() ;
    // These methods are called at the corresponding stage of processing of TSelector
  virtual void SlaveBegin(Reader* r) ;
  virtual void Process(Reader* r) ;
  virtual void SlaveTerminate(Reader* r) {} ;
  virtual void Terminate(TList* mergedList, std::string outFileName) ;
private:
  //histograms
  ZbPlots* h_zee_jet ;
  ZbPlots* h_zmm_jet ;
} ;

#endif
