#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>

#include "StdArg.hpp"
#include "src/Reader.h" 
#include "src/Processor.h"
#include "src/Selector.h"
#include "src/VH_selection.h"

#include "src/Global.h"

std::vector<std::string> splitNames(const std::string& files, std::string sep = ",")
{
  std::vector<std::string> fileList;
  for (size_t i=0,n; i <= files.length(); i=n+1)
  {
    n = files.find_first_of(sep,i);
    if (n == string::npos)
      n = files.length();
    std::string tmp = files.substr(i,n-i);
    std::string ttmp;
    for(unsigned int j=0; j<tmp.size(); j++)
    {
      if(tmp[j]==' ' || tmp[j]=='\n') continue;
      ttmp+=tmp[j];
    }
    fileList.push_back(ttmp);
  }
  return fileList;

}

void SetParameters(std::string fName, glob::Parameters& para) {
  std::string line;
  std::ifstream myfile (fName);
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      
      //skip comment line start with "//"
      if (line.find("//") != std::string::npos) continue ;
      
      std::vector<std::string> cont ;
      std::vector<std::string> cont_no_space ;
      std::string delim(" ") ;
      boost::split(cont, line, boost::is_any_of(delim));
      for (size_t i = 0 ; i < cont.size() ; ++i) {
        //std::cout << "\n" << cont[i] << cont[i].find(" ") << " " << std::string::npos;
        if (cont[i].find_first_not_of(' ') != std::string::npos) {
          cont_no_space.push_back(cont[i]) ;
        }
      }

      if (cont_no_space.size() != 3) {
        std::cout << "\n Need name, value, value type. Will skip this \"" << line << "\"" << std::endl ;
      }
      else {
        std::string name = cont_no_space[0] ;
        if (cont_no_space[2] == "int") para.Set(name, std::stoi(cont[1])) ;
        if (cont_no_space[2] == "float") para.Set(name, std::stof(cont_no_space[1])) ;
        if (cont_no_space[2] == "string") para.SetStr(name, cont_no_space[1]) ;
      }
    }

    myfile.close();
  }
  
  else cout << "Unable to open file"; 

}



int main(int argc, char *argv[]) {
  
  StdArg arg(argc,argv);
  std::cout << "======================================================================" << std::endl;
  std::cout << "ZH analysis call:" << std::endl;
  for (int i=0;i<argc;i++) std::cout << argv[i] << " ";
  std::cout << endl;
  std::cout << "======================================================================" << std::endl;
  arg.keys << "-in" << "-filelist" << "-cfg" << "-out" << "-data" << "-year" << "-syst" << "-centralGenWeight" 
           << "-firstentry" << "-lastentry" ; 
  arg.Process();

  std::vector<std::string> filenames;
  if ( arg.Key("-in")) 
    { 
      string inputfilenames;
      inputfilenames  = arg.Get<string>("-in"); 
      filenames = splitNames(inputfilenames);
    }
  else 
    { 
      string inputfilename;
      if (arg.Key("-filelist"))
	{
	  inputfilename  = arg.Get<string>("-filelist");
	  std::ifstream fList(inputfilename.c_str());
	  if (!fList){
	    cerr << ">>Can't open file " << inputfilename << endl;
	    return 1;
	  }
	  char lineFromFile[10000];
	  while (fList.getline(lineFromFile,10000)){
	    if ( strlen(lineFromFile)>0 ) filenames.push_back(string(lineFromFile));
	  }
	}
    }

  std::string cfgfilename =  "Configs/inputParameters.txt";
  if (arg.Key("-cfg")) cfgfilename = arg.Get<string>("-cfg");
  else std::cout << "\nUsing default configuration " << cfgfilename << std::endl; 
  
  std::string outputfilename =  "output.root";
  if (arg.Key("-out")) outputfilename = arg.Get<string>("-out");
  
  bool isData(false) ;
  
  int intisdata=0;
  int intfirstentry=0;
  long int intlastentry=-1;
  string syst = "NONE";
  string year = "";
  double centralGenWeight = 0;

  if (arg.Key("-data")) intisdata = arg.Get<int>("-data");
  if (arg.Key("-year")) year = arg.Get<string>("-year");
  if (arg.Key("-centralGenWeight")) centralGenWeight = arg.Get<double>("-centralGenWeight");
  if (arg.Key("-firstentry")) intfirstentry = arg.Get<int>("-firstentry");
  if (arg.Key("-lastentry")) intlastentry = arg.Get<int>("-lastentry");
  if (arg.Key("-syst")) syst = arg.Get<string>("-syst");

  if(intisdata!=0) isData=true;
  
  std::cout << "\n=================================" << std::endl ;
  std::cout << "\nIs data:              " << isData ;
  std::cout << "\nYear:                 " << year ;
  std::cout << "\nCentral gen weight    " << centralGenWeight;
  std::cout << "\nFirst and last entry: " << intfirstentry << " " << intlastentry ;
  std::cout << "\nSystematic:           " << syst ;
  
  std::cout << std::endl ;

  SetParameters(cfgfilename,CUTS) ;

#if defined(TFILE)
  TFile* f = TFile::Open(filenames[0].c_str());
  TTree* chain = (TTree*)f->Get("Events");
#endif
#if defined(TCHAIN)
  std::string chainName("Events") ;
  TChain* chain = new TChain(chainName.c_str()) ;
  for ( std::vector<std::string>::const_iterator it = filenames.begin();it != filenames.end(); it++) {
    cout << "reading file " << it->c_str() << endl;
    int retcode = chain->Add(it->c_str(),-1);
    if ( retcode == 0 ) throw std::invalid_argument("the file "+*it+" does not exist of does not contain the tree named "+chainName);
  }
#endif 
  
  std::cout << "\n Number of events: " << chain->GetEntries() ;
  if (intlastentry == -1) intlastentry = chain->GetEntries() ;
  
  Processor ana ;
  ana.setOutPutFileName(outputfilename) ;
  ana.SetDataInfo(isData,year) ;
  
  //collection of all selectors
  std::vector<Selector*> sels;
  
  //Selection for VH 
  VH_selection sel ;
  
  std::string fName_btagSF;
  std::string fName_ctagSF;
  std::string fName_puSF;
  std::string fName_PUjetID_SF;
  std::string fName_PUjetID_eff;
 
  //Syst
  if (syst == "L1PREFIRINGU") sel.SetL1prefiring("l1prefiringu");
  if (syst == "L1PREFIRINGD") sel.SetL1prefiring("l1prefiringd");
  std::string jetPUidUncType = "central";
  if (syst == "JETPUIDU") jetPUidUncType = "up";
  if (syst == "JETPUIDD") jetPUidUncType = "down";

#ifdef MC_2016
  fName_btagSF = "CalibData/DeepCSV_2016LegacySF_WP_V1.csv";
  fName_btagSF = "CalibData/DeepJet_2016LegacySF_WP_V1.csv";
  //fName_btagSF = "CalibData/wp_deepJet_2016.csv";
  fName_ctagSF = "CalibData/ctagging_wp_deepJet_106XUL16postVFP_v1.csv";
  fName_ctagSF = fName_btagSF;
  fName_puSF = "CalibData/2016_pileup_ratio.root";
  if (syst == "PUU") fName_puSF = "CalibData/2016_pileup_ratio_up.root";
  if (syst == "PUD") fName_puSF = "CalibData/2016_pileup_ratio_down.root";
#endif

#ifdef MC_2017
  fName_btagSF = "CalibData/DeepCSV_94XSF_WP_V4_B_F.csv";
  fName_btagSF = "CalibData/DeepFlavour_94XSF_WP_V3_B_F.csv";
  //fName_btagSF = "CalibData/wp_deepJet_2017.csv";
  fName_ctagSF = "CalibData/ctagging_wp_deepJet_106XUL17_v1.csv";
  fName_ctagSF = fName_btagSF;
  fName_puSF = "CalibData/2017_pileup_ratio.root";
  if (syst == "PUU") fName_puSF = "CalibData/2017_pileup_ratio_up.root";
  if (syst == "PUD") fName_puSF = "CalibData/2017_pileup_ratio_down.root";
#endif

#ifdef MC_2018
  fName_btagSF = "CalibData/DeepCSV_102XSF_WP_V1.csv";
  fName_btagSF = "CalibData/DeepJet_102XSF_WP_V1.csv";
  //fName_btagSF = "CalibData/wp_deepJet_2018.csv";
  fName_ctagSF = "CalibData/ctagging_wp_deepJet_106XUL18_v2.csv";
  fName_ctagSF = fName_btagSF;
  fName_puSF = "CalibData/2018_pileup_ratio.root";
  if (syst == "PUU") fName_puSF = "CalibData/2018_pileup_ratio_up.root";
  if (syst == "PUD") fName_puSF = "CalibData/2018_pileup_ratio_down.root";
#endif

#if defined(DATA_2016)
  //TODO: update to UL
  std::string fName_lumiMaskFilter("CalibData/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt");
#endif
#if defined(DATA_2017)
  //TODO: update to UL
  std::string fName_lumiMaskFilter("CalibData/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt");
#endif
#if defined(DATA_2018)
  //TODO: update to UL
  std::string fName_lumiMaskFilter("CalibData/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt");
#endif

  std::string btagUncType = "central";
  if (syst == "BTAGU") btagUncType = "up";
  if (syst == "BTAGD") btagUncType = "down";
  if (syst == "BTAGLU") btagUncType = "light_up";
  if (syst == "BTAGLD") btagUncType = "light_down";
  if (syst == "BTAGBCU") btagUncType = "bc_up";
  if (syst == "BTAGBCD") btagUncType = "bc_down";

  std::string ctagUncType = "central";

#if defined(MC_2016) || defined(MC_2017) || defined(MC_2018)
  sel.SetCentralGenWeight(centralGenWeight);
  sel.SetPileupSF(fName_puSF);
  fName_PUjetID_SF = "CalibData/scalefactorsPUID_81Xtraining.root";
  fName_PUjetID_eff = "CalibData/effcyPUID_81Xtraining.root";
  sel.SetPUjetidCalib(fName_PUjetID_SF,fName_PUjetID_eff,jetPUidUncType); //pileup jet ID SF
  
  // Add what allows us to calibrate based on b-tags
  std::cout << "Setting b-tagging calibration..." << std::endl;
  if (CUTS.GetStr("jet_main_btagWP")=="deepCSVT") sel.SetBtagCalib(fName_btagSF,"DeepCSV","CalibData/effT.root",btagUncType);
  if (CUTS.GetStr("jet_main_btagWP")=="deepJetT") sel.SetBtagCalib(fName_btagSF,"DeepJet","CalibData/effT.root",btagUncType);
  if (CUTS.GetStr("jet_main_btagWP")=="deepCSVM") sel.SetBtagCalib(fName_btagSF,"DeepCSV","CalibData/effM.root",btagUncType);
  if (CUTS.GetStr("jet_main_btagWP")=="deepJetM") sel.SetBtagCalib(fName_btagSF,"DeepJet","CalibData/effM.root",btagUncType);

  // Add what allows us to calibrate based on c-tags
  std::cout << "Setting c-tagging calibration..." << std::endl;
  if (CUTS.GetStr("jet_main_btagWP")=="deepCSVT") sel.SetCtagCalib(fName_ctagSF,"DeepCSV","CalibData/effT.root",ctagUncType);
  if (CUTS.GetStr("jet_main_btagWP")=="deepJetT") sel.SetCtagCalib(fName_ctagSF,"DeepJet","CalibData/effT.root",ctagUncType);
  if (CUTS.GetStr("jet_main_btagWP")=="deepCSVM") sel.SetCtagCalib(fName_ctagSF,"DeepCSV","CalibData/effM.root",ctagUncType);
  if (CUTS.GetStr("jet_main_btagWP")=="deepJetM") sel.SetCtagCalib(fName_ctagSF,"DeepJet","CalibData/effM.root",ctagUncType);
  std::cout << "jet_main_btagWP = " << CUTS.GetStr("jet_main_btagWP") << std::endl;
#endif
#if defined(DATA_2016) || defined(DATA_2017) || defined(DATA_2018)
  sel.SetLumiMaskFilter(fName_lumiMaskFilter);
#endif


  sels.push_back(&sel) ;
  
  //add all selectors to processors
  for (std::vector<Selector*>::iterator it = sels.begin() ; it != sels.end() ; it++) ana.addSelector(*it) ;
  std::cout << "\n Selections added" << std::endl; 
  chain->Process(&ana,"",intlastentry,intfirstentry) ;
  std::cout << "\n End processing" << std::endl;
  return 0 ;
}
