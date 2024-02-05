#define Processor_cxx

#include "Processor.h"

Processor::Processor(TTree * /*tree*/ t) { 
  m_outputfilename="ProcOutput.root"; 
} ;

//sequence is Begin, SlaveBegin, Init, Process, SlaveTerminate, Terminate
void Processor::Begin(TTree *tree) {
  Reader::Begin(tree);
}

void Processor::SlaveBegin(TTree *tree) {
  Reader::SlaveBegin(tree);
  for (std::vector<Selector*>::iterator it = m_selectors.begin(); it != m_selectors.end(); it++) {
    (*it)->SlaveBegin(this);
    (*it)->SetDataInfo(m_isData, m_year) ;
  }
}

void Processor::Init(TTree* tree) {
  Reader::Init(tree) ;
}


Bool_t Processor::Process(Long64_t entry) {
 
  //Process this entry, basically set all values for this events
  Reader::Process(entry) ;
  
  int pow10 = 100;
  while(entry > pow10*10) pow10 *= 10;
  if(entry>=pow10 && entry % pow10*((int) entry / pow10) == 0) std::cout << "Entry " << entry << std::endl;
  
  //Process all the selection
  for ( std::vector<Selector*>::iterator it = m_selectors.begin(); 
	it != m_selectors.end(); it++  )
  {
    (*it)->Process(this);
  }

  return kTRUE ; 
}

void Processor::SlaveTerminate()
{
  //SlaveTerminate all selectors
  for ( std::vector<Selector*>::iterator it = m_selectors.begin(); 
	it != m_selectors.end(); it++  )
  {
    (*it)->SlaveTerminate(this);
  }

  //SlaveTerminate reader
  Reader::SlaveTerminate();
  
  //write objects (plots, ...). These will be added automatically by after SlaveTerminate
  char* h=const_cast<char*>(m_outputfilename.c_str());
  TFile* outF = new TFile(h,"recreate");
  TIter next(fOutput);
  while(TObject* obj = next() ) {
    obj->Write();
  }
  outF->Close();

}

void Processor::Terminate()
{
  //currently there is nothing done in selector Terminate. This is reserved for printing cutflow for example, since all cutflow histograms in slaves were added. This will print out overall cut flow
  for ( std::vector<Selector*>::iterator it = m_selectors.begin(); 
	it != m_selectors.end(); it++  )
  {
    (*it)->Terminate(fOutput,m_outputfilename);
  }
}

