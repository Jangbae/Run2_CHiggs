#include "BTagCalibForLJMet.cpp"
#include "step1.cc"
#include<vector>
using namespace std;

void makeStep1(TString macroDir, TString inputFile, TString outputFile){

  gROOT->SetMacroPath(macroDir);

  gSystem->AddIncludePath("-I$CMSSW_BASE/src/");

  TString incl("-I");
  incl+=macroDir;
  gSystem->AddIncludePath(incl);

  if ( inputFile.Contains("Run2017") || inputFile.Contains("Run2018") || inputFile.Contains("Single") || inputFile.Contains("Double") || inputFile.Contains("MuonEG") || inputFile.Contains("EGamma") ) { 
    step1 t(inputFile,outputFile.ReplaceAll(".root","nominal.root"));
    t.Loop("ljmet", "ljmet"); 
    }

  else {
  	vector<TString> shifts = { "nominal" };//, "JECup", "JECdown", "JERup", "JERdown" };
  	for (size_t i =0; i<shifts.size(); ++i) {
      cout << endl << "Running shift " << shifts[i] << endl;
      TString tName = "ljmet";
      if ( !shifts[i].Contains("nominal") ) { tName.Append("_"); tName.Append(shifts[i]); }
      step1 t(inputFile,outputFile.ReplaceAll(".root",shifts[i].Append(".root"))); //"shifts[i]" is now changed to "shifts[i].root"
      t.saveHistograms();
      t.Loop(tName, "ljmet");
      outputFile.ReplaceAll(shifts[i],".root"); //Change outputFile back to its original name.
      }
    }
}

// void makeStep1(TString macroDir, TString inputFile, TString outputFile){
// 
//   gROOT->SetMacroPath(macroDir);
// 
//   gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
// 
//   TString incl("-I");
//   incl+=macroDir;
//   gSystem->AddIncludePath(incl);
// 
//   step1 t(inputFile,outputFile);
//   if ( inputFile.Contains("Run2017") || inputFile.Contains("Run2018") || inputFile.Contains("Single") || inputFile.Contains("Double") || inputFile.Contains("MuonEG") || inputFile.Contains("EGamma") )
//     { t.Loop("ljmet", "ljmet"); }
// 
//   else {
//     t.saveHistograms();
//     vector<TString> shifts = { "ljmet", "ljmet_JECup", "ljmet_JECdown", "ljmet_JERup", "ljmet_JERdown" };
//     for (size_t i =0; i<shifts.size(); ++i) {
//       cout << endl << "Running shift " << shifts[i] << endl;
//       t.Loop(shifts[i], shifts[i]);
//     }
//   }
// 
// }
