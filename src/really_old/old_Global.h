#ifndef Global_h
#define Global_h

#include <iostream>
#include <algorithm>
#include <vector>
#include "boost/variant.hpp"

namespace glob {
  
  class Parameters {
    
    public:
      
      Parameters() {
        parameterNames.push_back("jet_pt") ;
        parameterNames.push_back("jet_eta") ;
        parameterNames.push_back("jet_main_btagWP") ;
        parameterNames.push_back("lep_eta") ;
        parameterNames.push_back("lep_pt0") ;
        parameterNames.push_back("lep_pt1") ;
        parameterNames.push_back("lep_jetOverlap_pt");
        parameterNames.push_back("lep_jetOverlap_eta");
        parameterNames.push_back("muon_iso");
        parameterNames.push_back("ZMassL");
        parameterNames.push_back("ZMassH");
      } ;
      
      //passing a constant string https://stackoverflow.com/questions/4475634/c-pass-a-string
      //template<class T> T Get(const std::string& name) ;
      //template<class T> void Set(const std::string& name, T val) ;
      template<class T> T Get(const std::string& name) {
      //float Get(const std::string& name) {
        if (std::count(parameterNames.begin(),parameterNames.end(),name)) {
          if (name == "jet_pt") return jet_pt ;
          if (name == "jet_eta") return jet_eta ;
        }
        else {
          std::cout << "\n There is no parameter " << name << ". Will terminate" << std::endl ;
          exit(1) ; 
        }
        return 0 ;
      } ;
      
      std::string GetStr(const std::string& name) {
        if (std::count(parameterNames.begin(),parameterNames.end(),name)) {
          if (name == "jet_main_btagWP") return jet_main_btagWP;
        }
        else {
          std::cout << "\n There is no parameter " << name << ". Will terminate" << std::endl ;
          exit(1) ; 
        }
        return "" ;
      } ;

      void SetStr(const std::string& name, std::string val) {
        if (std::count(parameterNames.begin(),parameterNames.end(),name)) {
          if (name == "jet_main_btagWP") jet_main_btagWP = val;
        }
        else std::cout << "\n Can not set value for parameter named: " << name << ". Not exist in list of parameters " << std::endl; 
      };

      template<class T> void Set(const std::string& name, T val) {
        if (std::count(parameterNames.begin(),parameterNames.end(),name)) {
          if (name == "jet_pt") jet_pt = val;
          if (name == "jet_eta") jet_eta = val;
        }
        else std::cout << "\n Can not set value for parameter named: " << name << ". Not exist in list of parameters " << std::endl; 
      } ;


    private:
      float jet_pt ;
      float jet_eta ;
      std::string jet_main_btagWP ;
      std::vector<std::string> initializedVars ;
      std::vector<std::string> parameterNames ;

  } ;
  
}


//http://www.cplusplus.com/forum/unices/5784/
extern glob::Parameters CUTS ;

#endif
