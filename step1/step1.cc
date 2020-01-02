#define step1_cxx
#include "step1.h"
#include <fstream>
#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <sstream>
#include <string>
#include <vector>
#include <TH3.h>

#include "BTagCalibForLJMet.h"
using namespace std;

// ----------------------------------------------------------------------------
// Define functions
// ----------------------------------------------------------------------------

bool comparepair( const std::pair<double,int> a, const std::pair<double,int> b) { return a.first > b.first; }

TRandom3 Rand;

const double MTOP  = 173.5;
const double MW    = 80.4; 

bool step1::applySF(bool& isTagged, float tag_SF, float tag_eff){
  
  bool newTag = isTagged;
  if (tag_SF == 1) return newTag; //no correction needed 

  //throw die
  float coin = Rand.Uniform(1.);    

  if(tag_SF > 1){  // use this if SF>1

    if( !isTagged ) {

      //fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - tag_SF) / (1.0 - (tag_SF/tag_eff) );

      //upgrade to tagged
      if( coin < mistagPercent ) {newTag = true;}
    }

  }else{  // use this if SF<1
      
    //downgrade tagged to untagged
    if( isTagged && coin > tag_SF ) {newTag = false;}

  }

  return newTag;
}

// -------------------------------------
// Function- SAVE Histograms from ljmet
//--------------------------------------
void step1::saveHistograms() 
{ 
TH1D* numhist = (TH1D*)inputFile->Get("mcweightanalyzer/NumTrueHist");
TH1D* wgthist = (TH1D*)inputFile->Get("mcweightanalyzer/weightHist");
outputFile->cd();
numhist->Write();
wgthist->Write();

}


// ----------------------------------------------------------------------------
// MAIN EVENT LOOP
// ----------------------------------------------------------------------------

void step1::Loop(TString inTreeName, TString outTreeName ) 
{

  // setup calibration + reader
  BTagCalibrationForLJMet calib("deepcsv", "DeepCSV_94XSF_V4_B_F.csv");
  BTagCalibrationForLJMetReader reader(BTagEntryForLJMet::OP_RESHAPING,  // operating point
			       "central",             // central sys type
			       {
                    "down_cferr1", "up_cferr1",
                    "down_cferr2", "up_cferr2",
                    "down_hf", "up_hf",
                    "down_hfstats1", "up_hfstats1",
                    "down_hfstats2", "up_hfstats2",
                    "down_jes", "up_jes",
                    "down_jesAbsoluteMPFBias", "up_jesAbsoluteMPFBias",
                    "down_jesAbsoluteScale", "up_jesAbsoluteScale",
                    "down_jesAbsoluteStat", "up_jesAbsoluteStat",
                    "down_jesFlavorQCD", "up_jesFlavorQCD",
                    "down_jesFragmentation", "up_jesFragmentation",
                    "down_jesPileUpDataMC", "up_jesPileUpDataMC",
                    "down_jesPileUpPtBB", "up_jesPileUpPtBB",
                    "down_jesPileUpPtEC1", "up_jesPileUpPtEC1",
                    "down_jesPileUpPtEC2", "up_jesPileUpPtEC2",
                    "down_jesPileUpPtHF", "up_jesPileUpPtHF",
                    "down_jesPileUpPtRef", "up_jesPileUpPtRef",
                    "down_jesRelativeBal", "up_jesRelativeBal",
                    "down_jesRelativeFSR", "up_jesRelativeFSR",
                    "down_jesRelativeJEREC1", "up_jesRelativeJEREC1",
                    "down_jesRelativeJEREC2", "up_jesRelativeJEREC2",
                    "down_jesRelativeJERHF", "up_jesRelativeJERHF",
                    "down_jesRelativePtBB", "up_jesRelativePtBB",
                    "down_jesRelativePtEC1", "up_jesRelativePtEC1",
                    "down_jesRelativePtEC2", "up_jesRelativePtEC2",
                    "down_jesRelativePtHF", "up_jesRelativePtHF",
                    "down_jesRelativeStatEC", "up_jesRelativeStatEC",
                    "down_jesRelativeStatFSR", "up_jesRelativeStatFSR",
                    "down_jesRelativeStatHF", "up_jesRelativeStatHF",
                    "down_jesSinglePionECAL", "up_jesSinglePionECAL",
                    "down_jesSinglePionHCAL", "up_jesSinglePionHCAL",
                    "down_jesTimePtEta", "up_jesTimePtEta",
                    "down_lf", "up_lf",
                    "down_lfstats1", "up_lfstats1",
                    "down_lfstats2", "up_lfstats2"				   
				   }
				   );      // other sys types
  //C-only: cferr*. L-only: hf, lfstats*. B-only: lf, hfstats*
  
  reader.load(calib,                 // calibration instance
	      BTagEntryForLJMet::FLAV_B,     // btag flavour
	      "iterativefit");       // measurement type
  reader.load(calib, BTagEntryForLJMet::FLAV_C, "iterativefit");     // for FLAV_C
  reader.load(calib, BTagEntryForLJMet::FLAV_UDSG, "iterativefit");     // for FLAV_UDSG
  
  // ----------------------------------------------------------------------------
  // Turn on input tree branches
  // ----------------------------------------------------------------------------
  inputTree=(TTree*)inputFile->Get(inTreeName+"/"+inTreeName);
  if(inputTree->GetEntries()==0) {
    std::cout<<"WARNING! Found 0 events in the tree!!!!"<<std::endl;;
    exit(1);
  }
  Init(inputTree);

   if (inputTree == 0) return;
   
   inputTree->SetBranchStatus("*",0);

   //Event info
   inputTree->SetBranchStatus("event_CommonCalc",1);
   inputTree->SetBranchStatus("run_CommonCalc",1);
   inputTree->SetBranchStatus("lumi_CommonCalc",1);
   //   inputTree->SetBranchStatus("nPV_MultiLepCalc",1);
   inputTree->SetBranchStatus("nTrueInteractions_MultiLepCalc",1);
   inputTree->SetBranchStatus("MCWeight_MultiLepCalc",1);
   inputTree->SetBranchStatus("LHEweightids_MultiLepCalc",1);
   inputTree->SetBranchStatus("LHEweights_MultiLepCalc",1);
   inputTree->SetBranchStatus("NewPDFweights_MultiLepCalc",1);
   inputTree->SetBranchStatus("HTfromHEPUEP_MultiLepCalc",1);
   inputTree->SetBranchStatus("L1NonPrefiringProb_CommonCalc",1);
   inputTree->SetBranchStatus("L1NonPrefiringProbUp_CommonCalc",1);
   inputTree->SetBranchStatus("L1NonPrefiringProbDown_CommonCalc",1);
   
   //triggers
   inputTree->SetBranchStatus("vsSelMCTriggersEl_MultiLepCalc",1);
   inputTree->SetBranchStatus("viSelMCTriggersEl_MultiLepCalc",1);
   inputTree->SetBranchStatus("vsSelMCTriggersMu_MultiLepCalc",1);
   inputTree->SetBranchStatus("viSelMCTriggersMu_MultiLepCalc",1);
   inputTree->SetBranchStatus("vsSelTriggersEl_MultiLepCalc",1);
   inputTree->SetBranchStatus("viSelTriggersEl_MultiLepCalc",1);
   inputTree->SetBranchStatus("vsSelTriggersMu_MultiLepCalc",1);
   inputTree->SetBranchStatus("viSelTriggersMu_MultiLepCalc",1);
   
   //electrons
   inputTree->SetBranchStatus("elPt_MultiLepCalc",1);
   inputTree->SetBranchStatus("elEta_MultiLepCalc",1);
   inputTree->SetBranchStatus("elPhi_MultiLepCalc",1);
   inputTree->SetBranchStatus("elEnergy_MultiLepCalc",1);
   inputTree->SetBranchStatus("elMiniIso_MultiLepCalc",1);  
   
   //muons
   inputTree->SetBranchStatus("muPt_MultiLepCalc",1);
   inputTree->SetBranchStatus("muEta_MultiLepCalc",1);
   inputTree->SetBranchStatus("muPhi_MultiLepCalc",1);
   inputTree->SetBranchStatus("muEnergy_MultiLepCalc",1);
   inputTree->SetBranchStatus("muMiniIso_MultiLepCalc",1);  
   
   //missing et
   inputTree->SetBranchStatus("corr_met_MultiLepCalc",1);
   inputTree->SetBranchStatus("corr_met_phi_MultiLepCalc",1);
   inputTree->SetBranchStatus("corr_metmod_MultiLepCalc",1);
   inputTree->SetBranchStatus("corr_metmod_phi_MultiLepCalc",1);

   //boosted truth
   inputTree->SetBranchStatus("HadronicVHtPt_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtEta_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtPhi_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtEnergy_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtID_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtD0Pt_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtD0Eta_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtD0Phi_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtD0E_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtD1Pt_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtD1Eta_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtD1Phi_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtD1E_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtD2Pt_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtD2Eta_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtD2Phi_JetSubCalc",1);
   inputTree->SetBranchStatus("HadronicVHtD2E_JetSubCalc",1);

   //genParticles
   inputTree->SetBranchStatus("genPt_MultiLepCalc",1);
   inputTree->SetBranchStatus("genEta_MultiLepCalc",1);
   inputTree->SetBranchStatus("genPhi_MultiLepCalc",1);
   inputTree->SetBranchStatus("genEnergy_MultiLepCalc",1);
   inputTree->SetBranchStatus("genStatus_MultiLepCalc",1);
   inputTree->SetBranchStatus("genID_MultiLepCalc",1);
   inputTree->SetBranchStatus("genIndex_MultiLepCalc",1);
   inputTree->SetBranchStatus("genMotherID_MultiLepCalc",1);
   inputTree->SetBranchStatus("genMotherIndex_MultiLepCalc",1);
   inputTree->SetBranchStatus("genJetPt_MultiLepCalc",1);
   inputTree->SetBranchStatus("genJetEta_MultiLepCalc",1);
   inputTree->SetBranchStatus("genJetPhi_MultiLepCalc",1);
   inputTree->SetBranchStatus("genJetEnergy_MultiLepCalc",1);
   if(isTT){
   inputTree->SetBranchStatus("genJetPtNoClean_MultiLepCalc",1);
   inputTree->SetBranchStatus("genJetEtaNoClean_MultiLepCalc",1);
   inputTree->SetBranchStatus("genJetPhiNoClean_MultiLepCalc",1);
   inputTree->SetBranchStatus("genJetEnergyNoClean_MultiLepCalc",1);
   }                

   //JetSubCalc
   inputTree->SetBranchStatus("theJetHFlav_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetPFlav_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetPt_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetEta_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetPhi_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetEnergy_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetDeepCSVb_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetDeepCSVbb_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetDeepCSVc_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetDeepCSVudsg_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetBTag_bSFup_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetBTag_bSFdn_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetBTag_lSFup_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetBTag_lSFdn_JetSubCalc",1);
   inputTree->SetBranchStatus("theJetBTag_JetSubCalc",1);
   
   //top
   inputTree->SetBranchStatus("ttbarMass_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topEnergy_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topEta_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topMass_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topPhi_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topPt_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topID_TTbarMassCalc",1);

   inputTree->SetBranchStatus("allTopsEnergy_TTbarMassCalc",1);
   inputTree->SetBranchStatus("allTopsEta_TTbarMassCalc",1);
   inputTree->SetBranchStatus("allTopsPhi_TTbarMassCalc",1);
   inputTree->SetBranchStatus("allTopsPt_TTbarMassCalc",1);
   inputTree->SetBranchStatus("allTopsID_TTbarMassCalc",1);
   inputTree->SetBranchStatus("allTopsStatus_TTbarMassCalc",1);

   inputTree->SetBranchStatus("genTtbarIdCategory_TTbarMassCalc",1);

   //top W
   inputTree->SetBranchStatus("topWEnergy_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topWEta_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topWPhi_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topWPt_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topWID_TTbarMassCalc",1);

   //top b
   inputTree->SetBranchStatus("topbEnergy_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topbEta_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topbPhi_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topbPt_TTbarMassCalc",1);
   inputTree->SetBranchStatus("topbID_TTbarMassCalc",1);

   //HOTTaggerCalc
   inputTree->SetBranchStatus("topJet1Index_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topJet2Index_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topJet3Index_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topNAK4_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topNtops_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topBestGenEnergy_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topBestGenEta_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topBestGenPhi_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topBestGenPt_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topDRmax_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topDThetaMax_HOTTaggerCalc", 1);
   inputTree->SetBranchStatus("topDThetaMin_HOTTaggerCalc", 1);
   inputTree->SetBranchStatus("topDiscriminator_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topEta_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topMass_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topNconstituents_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topPhi_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topPt_HOTTaggerCalc",1);
   inputTree->SetBranchStatus("topType_HOTTaggerCalc",1);

   
  // ----------------------------------------------------------------------------
  // Create output tree and define branches
  // ----------------------------------------------------------------------------
   
   // OUTPUT FILE
   outputFile->cd();
   TTree *outputTree = new TTree(outTreeName, outTreeName);


   //Common
   outputTree->Branch("event_CommonCalc",&event_CommonCalc,"event_CommonCalc/L");
   outputTree->Branch("run_CommonCalc",&run_CommonCalc,"run_CommonCalc/I");
   outputTree->Branch("lumi_CommonCalc",&lumi_CommonCalc,"lumi_CommonCalc/I");
   outputTree->Branch("nTrueInteractions_MultiLepCalc",&nTrueInteractions_MultiLepCalc,"nTrueInteractions_MultiLepCalc/I");
   outputTree->Branch("isElectron",&isElectron,"isElectron/I");
   outputTree->Branch("isMuon",&isMuon,"isMuon/I");
   outputTree->Branch("MCPastTrigger",&MCPastTrigger,"MCPastTrigger/I");
   outputTree->Branch("DataPastTrigger",&DataPastTrigger,"DataPastTrigger/I");
   outputTree->Branch("L1NonPrefiringProb_CommonCalc",&L1NonPrefiringProb_CommonCalc,"L1NonPrefiringProb_CommonCalc/D");
   outputTree->Branch("L1NonPrefiringProbUp_CommonCalc",&L1NonPrefiringProbUp_CommonCalc,"L1NonPrefiringProbUp_CommonCalc/D");
   outputTree->Branch("L1NonPrefiringProbDown_CommonCalc",&L1NonPrefiringProbDown_CommonCalc,"L1NonPrefiringProbDown_CommonCalc/D");
   
   //weights
   outputTree->Branch("MCWeight_MultiLepCalc",&MCWeight_MultiLepCalc,"MCWeight_MultiLepCalc/D");
   outputTree->Branch("renormWeights",&renormWeights);
   outputTree->Branch("pdfWeights",&pdfWeights);
   outputTree->Branch("pdfNewWeights",&pdfNewWeights);
   outputTree->Branch("pdfNewNominalWeight",&pdfNewNominalWeight,"pdfNewNominalWeight/F");
   outputTree->Branch("pileupWeight",&pileupWeight,"pileupWeight/F");
   outputTree->Branch("pileupWeightUp",&pileupWeightUp,"pileupWeightUp/F");
   outputTree->Branch("pileupWeightDown",&pileupWeightDown,"pileupWeightDown/F");
   outputTree->Branch("HTSF_Pol",&HTSF_Pol,"HTSF_Pol/F");
   outputTree->Branch("HTSF_PolUp",&HTSF_PolUp,"HTSF_PolUp/F");
   outputTree->Branch("HTSF_PolDn",&HTSF_PolDn,"HTSF_PolDn/F");
   outputTree->Branch("topPtWeight13TeV",&topPtWeight13TeV,"topPtWeight13TeV/F");          
   outputTree->Branch("EGammaGsfSF",&EGammaGsfSF,"EGammaGsfSF/F");
   outputTree->Branch("lepIdSF",&lepIdSF,"lepIdSF/F");
   outputTree->Branch("triggerSF",&triggerSF,"triggerSF/F");
   outputTree->Branch("isoSF",&isoSF,"isoSF/F");
   
   //ttbar generator
   outputTree->Branch("ttbarMass_TTbarMassCalc",&ttbarMass_TTbarMassCalc,"ttbarMass_TTbarMassCalc/D");
   outputTree->Branch("genTopPt",&genTopPt,"genTopPt/F");
   outputTree->Branch("genAntiTopPt",&genAntiTopPt,"genAntiTopPt/F");
   outputTree->Branch("topEnergy_TTbarMassCalc",&topEnergy_TTbarMassCalc);
   outputTree->Branch("topEta_TTbarMassCalc",&topEta_TTbarMassCalc);
   outputTree->Branch("topMass_TTbarMassCalc",&topMass_TTbarMassCalc);
   outputTree->Branch("topPhi_TTbarMassCalc",&topPhi_TTbarMassCalc);   
   outputTree->Branch("topPt_TTbarMassCalc",&topPt_TTbarMassCalc);      
   outputTree->Branch("topID_TTbarMassCalc",&topID_TTbarMassCalc);
   
   //leptons
   outputTree->Branch("corr_met_MultiLepCalc",&corr_met_MultiLepCalc,"corr_met_MultiLepCalc/D");
   outputTree->Branch("corr_met_phi_MultiLepCalc",&corr_met_phi_MultiLepCalc,"corr_met_phi_MultiLepCalc/D");
   outputTree->Branch("corr_metmod_MultiLepCalc",&corr_metmod_MultiLepCalc,"corr_metmod_MultiLepCalc/D");
   outputTree->Branch("corr_metmod_phi_MultiLepCalc",&corr_metmod_phi_MultiLepCalc,"corr_metmod_phi_MultiLepCalc/D");
   outputTree->Branch("leptonPt_MultiLepCalc",&leptonPt_MultiLepCalc,"leptonPt_MultiLepCalc/F");
   outputTree->Branch("leptonEta_MultiLepCalc",&leptonEta_MultiLepCalc,"leptonEta_MultiLepCalc/F");
   outputTree->Branch("leptonPhi_MultiLepCalc",&leptonPhi_MultiLepCalc,"leptonPhi_MultiLepCalc/F");
   outputTree->Branch("leptonEnergy_MultiLepCalc",&leptonEnergy_MultiLepCalc,"leptonEnergy_MultiLepCalc/F");
   outputTree->Branch("leptonMVAValue_MultiLepCalc",&leptonMVAValue_MultiLepCalc,"leptonMVAValue_MultiLepCalc/F");
   outputTree->Branch("leptonMiniIso_MultiLepCalc",&leptonMiniIso_MultiLepCalc,"leptonMiniIso_MultiLepCalc/F");
   outputTree->Branch("MT_lepMet",&MT_lepMet,"MT_lepMet/F");
   outputTree->Branch("MT_lepMetmod",&MT_lepMetmod,"MT_lepMetmod/F");
   outputTree->Branch("minDPhi_MetJet",&minDPhi_MetJet,"minDPhi_MetJet/F");

   // AK4
   outputTree->Branch("theJetPt_JetSubCalc_PtOrdered",&theJetPt_JetSubCalc_PtOrdered);
   outputTree->Branch("theJetEta_JetSubCalc_PtOrdered",&theJetEta_JetSubCalc_PtOrdered);
   outputTree->Branch("theJetPhi_JetSubCalc_PtOrdered",&theJetPhi_JetSubCalc_PtOrdered);
   outputTree->Branch("theJetEnergy_JetSubCalc_PtOrdered",&theJetEnergy_JetSubCalc_PtOrdered);
   outputTree->Branch("theJetDeepCSVb_JetSubCalc_PtOrdered",&theJetDeepCSVb_JetSubCalc_PtOrdered);
   outputTree->Branch("theJetDeepCSVbb_JetSubCalc_PtOrdered",&theJetDeepCSVbb_JetSubCalc_PtOrdered);
   outputTree->Branch("theJetDeepCSVc_JetSubCalc_PtOrdered",&theJetDeepCSVc_JetSubCalc_PtOrdered);
   outputTree->Branch("theJetDeepCSVudsg_JetSubCalc_PtOrdered",&theJetDeepCSVudsg_JetSubCalc_PtOrdered);
   outputTree->Branch("theJetHFlav_JetSubCalc_PtOrdered",&theJetHFlav_JetSubCalc_PtOrdered);
   outputTree->Branch("theJetPFlav_JetSubCalc_PtOrdered",&theJetPFlav_JetSubCalc_PtOrdered);
   outputTree->Branch("theJetBTag_JetSubCalc_PtOrdered",&theJetBTag_JetSubCalc_PtOrdered);
   outputTree->Branch("AK4HTpMETpLepPt",&AK4HTpMETpLepPt,"AK4HTpMETpLepPt/F");
   outputTree->Branch("AK4HT",&AK4HT,"AK4HT/F");
   outputTree->Branch("NJets_JetSubCalc",&NJets_JetSubCalc,"NJets_JetSubCalc/I");
   outputTree->Branch("NJetsCSV_JetSubCalc",&NJetsCSV_JetSubCalc,"NJetsCSV_JetSubCalc/I");
   outputTree->Branch("NJetsCSVwithSF_JetSubCalc",&NJetsCSVwithSF_JetSubCalc,"NJetsCSVwithSF_JetSubCalc/I");
   outputTree->Branch("NJetsCSVwithSF_JetSubCalc_bSFup",&NJetsCSVwithSF_JetSubCalc_bSFup,"NJetsCSVwithSF_JetSubCalc_bSFup/I");
   outputTree->Branch("NJetsCSVwithSF_JetSubCalc_bSFdn",&NJetsCSVwithSF_JetSubCalc_bSFdn,"NJetsCSVwithSF_JetSubCalc_bSFdn/I");
   outputTree->Branch("NJetsCSVwithSF_JetSubCalc_lSFup",&NJetsCSVwithSF_JetSubCalc_lSFup,"NJetsCSVwithSF_JetSubCalc_lSFup/I");
   outputTree->Branch("NJetsCSVwithSF_JetSubCalc_lSFdn",&NJetsCSVwithSF_JetSubCalc_lSFdn,"NJetsCSVwithSF_JetSubCalc_lSFdn/I");
   outputTree->Branch("deltaR_lepMinMlb",&deltaR_lepMinMlb,"deltaR_lepMinMlb/F");
   outputTree->Branch("deltaR_lepMinMlb_bSFup",&deltaR_lepMinMlb_bSFup,"deltaR_lepMinMlb_bSFup/F");
   outputTree->Branch("deltaR_lepMinMlb_bSFdn",&deltaR_lepMinMlb_bSFdn,"deltaR_lepMinMlb_bSFdn/F");
   outputTree->Branch("deltaR_lepMinMlb_lSFup",&deltaR_lepMinMlb_lSFup,"deltaR_lepMinMlb_lSFup/F");
   outputTree->Branch("deltaR_lepMinMlb_lSFdn",&deltaR_lepMinMlb_lSFdn,"deltaR_lepMinMlb_lSFdn/F");
   outputTree->Branch("deltaR_lepMinMlj",&deltaR_lepMinMlj,"deltaR_lepMinMlj/F");
   outputTree->Branch("minMleppBjet",&minMleppBjet,"minMleppBjet/F");
   outputTree->Branch("minMleppBjet_bSFup",&minMleppBjet_bSFup,"minMleppBjet_bSFup/F");
   outputTree->Branch("minMleppBjet_bSFdn",&minMleppBjet_bSFdn,"minMleppBjet_bSFdn/F");
   outputTree->Branch("minMleppBjet_lSFup",&minMleppBjet_lSFup,"minMleppBjet_lSFup/F");
   outputTree->Branch("minMleppBjet_lSFdn",&minMleppBjet_lSFdn,"minMleppBjet_lSFdn/F");
   outputTree->Branch("minMleppJet",&minMleppJet,"mixnMleppJet/F");
   outputTree->Branch("minDR_lepJet",&minDR_lepJet,"minDR_lepJet/F");
   outputTree->Branch("ptRel_lepJet",&ptRel_lepJet,"ptRel_lepJet/F");
   outputTree->Branch("deltaR_lepJets",&deltaR_lepJets);
   outputTree->Branch("deltaR_lepBJets",&deltaR_lepBJets);
   outputTree->Branch("deltaR_lepBJets_bSFup",&deltaR_lepBJets_bSFup);
   outputTree->Branch("deltaR_lepBJets_bSFdn",&deltaR_lepBJets_bSFdn);
   outputTree->Branch("deltaR_lepBJets_lSFup",&deltaR_lepBJets_lSFup);
   outputTree->Branch("deltaR_lepBJets_lSFdn",&deltaR_lepBJets_lSFdn);

   // AK4 gen
   outputTree->Branch("genPt_MultiLepCalc",&genPt_MultiLepCalc);
   outputTree->Branch("genEta_MultiLepCalc",&genEta_MultiLepCalc);
   outputTree->Branch("genPhi_MultiLepCalc",&genPhi_MultiLepCalc);
   outputTree->Branch("genEnergy_MultiLepCalc",&genEnergy_MultiLepCalc);
   outputTree->Branch("genStatus_MultiLepCalc",&genStatus_MultiLepCalc);
   outputTree->Branch("genID_MultiLepCalc",&genID_MultiLepCalc);
   outputTree->Branch("genJetPt_MultiLepCalc",&genJetPt_MultiLepCalc);                  
   outputTree->Branch("genJetEta_MultiLepCalc",&genJetEta_MultiLepCalc);                     
   outputTree->Branch("genJetPhi_MultiLepCalc",&genJetPhi_MultiLepCalc);                  
   outputTree->Branch("genJetEnergy_MultiLepCalc",&genJetEnergy_MultiLepCalc);  
   if(isTT){                
   outputTree->Branch("genJetPtNoClean_MultiLepCalc",&genJetPtNoClean_MultiLepCalc);                  
   outputTree->Branch("genJetEtaNoClean_MultiLepCalc",&genJetEtaNoClean_MultiLepCalc);                     
   outputTree->Branch("genJetPhiNoClean_MultiLepCalc",&genJetPhiNoClean_MultiLepCalc);                  
   outputTree->Branch("genJetEnergyNoClean_MultiLepCalc",&genJetEnergyNoClean_MultiLepCalc); 
   }                 
   outputTree->Branch("genIndex_MultiLepCalc",&genIndex_MultiLepCalc);
   outputTree->Branch("genMotherID_MultiLepCalc",&genMotherID_MultiLepCalc);
   outputTree->Branch("genMotherIndex_MultiLepCalc",&genMotherIndex_MultiLepCalc);

   // jet truth
   outputTree->Branch("HadronicVHtID_JetSubCalc",&HadronicVHtID_JetSubCalc);
   outputTree->Branch("HadronicVHtPt_JetSubCalc",&HadronicVHtPt_JetSubCalc);
   outputTree->Branch("HadronicVHtEta_JetSubCalc",&HadronicVHtEta_JetSubCalc);
   outputTree->Branch("HadronicVHtPhi_JetSubCalc",&HadronicVHtPhi_JetSubCalc);
   outputTree->Branch("HadronicVHtEnergy_JetSubCalc",&HadronicVHtEnergy_JetSubCalc);

   outputTree->Branch("isHTgt500Njetge9",&isHTgt500Njetge9,"isHTgt500Njetge9/I");
   outputTree->Branch("BJetLeadPt",&BJetLeadPt,"BJetLeadPt/F");
   outputTree->Branch("BJetLeadPt_bSFup",&BJetLeadPt_bSFup,"BJetLeadPt_bSFup/F");
   outputTree->Branch("BJetLeadPt_bSFdn",&BJetLeadPt_bSFdn,"BJetLeadPt_bSFdn/F");
   outputTree->Branch("BJetLeadPt_lSFup",&BJetLeadPt_lSFup,"BJetLeadPt_lSFup/F");
   outputTree->Branch("BJetLeadPt_lSFdn",&BJetLeadPt_lSFdn,"BJetLeadPt_lSFdn/F");
   outputTree->Branch("WJetLeadPt",&WJetLeadPt,"WJetLeadPt/F");
   outputTree->Branch("TJetLeadPt",&TJetLeadPt,"TJetLeadPt/F"); 

   outputTree->Branch("genTtbarIdCategory_TTbarMassCalc",&genTtbarIdCategory_TTbarMassCalc);    

   //B-tag calibration reweight for SF
   outputTree->Branch("CalibReaderRewgt",&CalibReaderRewgt,"CalibReaderRewgt/F");
   outputTree->Branch("CalibReaderRewgt_JESup",&CalibReaderRewgt_JESup,"CalibReaderRewgt_JESup/F");
   outputTree->Branch("CalibReaderRewgt_JESdn",&CalibReaderRewgt_JESdn,"CalibReaderRewgt_JESdn/F");
   outputTree->Branch("CalibReaderRewgt_JESAbsoluteMPFBiasup",&CalibReaderRewgt_JESAbsoluteMPFBiasup,"CalibReaderRewgt_JESAbsoluteMPFBiasup/F");
   outputTree->Branch("CalibReaderRewgt_JESAbsoluteMPFBiasdn",&CalibReaderRewgt_JESAbsoluteMPFBiasdn,"CalibReaderRewgt_JESAbsoluteMPFBiasdn/F");
   outputTree->Branch("CalibReaderRewgt_JESAbsoluteScaleup",&CalibReaderRewgt_JESAbsoluteScaleup,"CalibReaderRewgt_JESAbsoluteScaleup/F");
   outputTree->Branch("CalibReaderRewgt_JESAbsoluteScaledn",&CalibReaderRewgt_JESAbsoluteScaledn,"CalibReaderRewgt_JESAbsoluteScaledn/F");
   outputTree->Branch("CalibReaderRewgt_JESAbsoluteStatup",&CalibReaderRewgt_JESAbsoluteStatup,"CalibReaderRewgt_JESAbsoluteStatup/F");
   outputTree->Branch("CalibReaderRewgt_JESAbsoluteStatdn",&CalibReaderRewgt_JESAbsoluteStatdn,"CalibReaderRewgt_JESAbsoluteStatdn/F");
   outputTree->Branch("CalibReaderRewgt_JESFlavorQCDup",&CalibReaderRewgt_JESFlavorQCDup,"CalibReaderRewgt_JESFlavorQCDup/F");
   outputTree->Branch("CalibReaderRewgt_JESFlavorQCDdn",&CalibReaderRewgt_JESFlavorQCDdn,"CalibReaderRewgt_JESFlavorQCDdn/F");
   outputTree->Branch("CalibReaderRewgt_JESFragmentationup",&CalibReaderRewgt_JESFragmentationup,"CalibReaderRewgt_JESFragmentationup/F");
   outputTree->Branch("CalibReaderRewgt_JESFragmentationdn",&CalibReaderRewgt_JESFragmentationdn,"CalibReaderRewgt_JESFragmentationdn/F");
   outputTree->Branch("CalibReaderRewgt_JESPileUpDataMCup",&CalibReaderRewgt_JESPileUpDataMCup,"CalibReaderRewgt_JESPileUpDataMCup/F");
   outputTree->Branch("CalibReaderRewgt_JESPileUpDataMCdn",&CalibReaderRewgt_JESPileUpDataMCdn,"CalibReaderRewgt_JESPileUpDataMCdn/F");
   outputTree->Branch("CalibReaderRewgt_JESPileUpPtBBup",&CalibReaderRewgt_JESPileUpPtBBup,"CalibReaderRewgt_JESPileUpPtBBup/F");
   outputTree->Branch("CalibReaderRewgt_JESPileUpPtBBdn",&CalibReaderRewgt_JESPileUpPtBBdn,"CalibReaderRewgt_JESPileUpPtBBdn/F");
   outputTree->Branch("CalibReaderRewgt_JESPileUpPtEC1up",&CalibReaderRewgt_JESPileUpPtEC1up,"CalibReaderRewgt_JESPileUpPtEC1up/F");
   outputTree->Branch("CalibReaderRewgt_JESPileUpPtEC1dn",&CalibReaderRewgt_JESPileUpPtEC1dn,"CalibReaderRewgt_JESPileUpPtEC1dn/F");
   outputTree->Branch("CalibReaderRewgt_JESPileUpPtEC2up",&CalibReaderRewgt_JESPileUpPtEC2up,"CalibReaderRewgt_JESPileUpPtEC2up/F");
   outputTree->Branch("CalibReaderRewgt_JESPileUpPtEC2dn",&CalibReaderRewgt_JESPileUpPtEC2dn,"CalibReaderRewgt_JESPileUpPtEC2dn/F");
   outputTree->Branch("CalibReaderRewgt_JESPileUpPtHFup",&CalibReaderRewgt_JESPileUpPtHFup,"CalibReaderRewgt_JESPileUpPtHFup/F");
   outputTree->Branch("CalibReaderRewgt_JESPileUpPtHFdn",&CalibReaderRewgt_JESPileUpPtHFdn,"CalibReaderRewgt_JESPileUpPtHFdn/F");
   outputTree->Branch("CalibReaderRewgt_JESPileUpPtRefup",&CalibReaderRewgt_JESPileUpPtRefup,"CalibReaderRewgt_JESPileUpPtRefup/F");
   outputTree->Branch("CalibReaderRewgt_JESPileUpPtRefdn",&CalibReaderRewgt_JESPileUpPtRefdn,"CalibReaderRewgt_JESPileUpPtRefdn/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeBalup",&CalibReaderRewgt_JESRelativeBalup,"CalibReaderRewgt_JESRelativeBalup/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeBaldn",&CalibReaderRewgt_JESRelativeBaldn,"CalibReaderRewgt_JESRelativeBaldn/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeFSRup",&CalibReaderRewgt_JESRelativeFSRup,"CalibReaderRewgt_JESRelativeFSRup/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeFSRdn",&CalibReaderRewgt_JESRelativeFSRdn,"CalibReaderRewgt_JESRelativeFSRdn/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeJEREC1up",&CalibReaderRewgt_JESRelativeJEREC1up,"CalibReaderRewgt_JESRelativeJEREC1up/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeJEREC1dn",&CalibReaderRewgt_JESRelativeJEREC1dn,"CalibReaderRewgt_JESRelativeJEREC1dn/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeJEREC2up",&CalibReaderRewgt_JESRelativeJEREC2up,"CalibReaderRewgt_JESRelativeJEREC2up/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeJEREC2dn",&CalibReaderRewgt_JESRelativeJEREC2dn,"CalibReaderRewgt_JESRelativeJEREC2dn/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeJERHFup",&CalibReaderRewgt_JESRelativeJERHFup,"CalibReaderRewgt_JESRelativeJERHFup/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeJERHFdn",&CalibReaderRewgt_JESRelativeJERHFdn,"CalibReaderRewgt_JESRelativeJERHFdn/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativePtBBup",&CalibReaderRewgt_JESRelativePtBBup,"CalibReaderRewgt_JESRelativePtBBup/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativePtBBdn",&CalibReaderRewgt_JESRelativePtBBdn,"CalibReaderRewgt_JESRelativePtBBdn/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativePtEC1up",&CalibReaderRewgt_JESRelativePtEC1up,"CalibReaderRewgt_JESRelativePtEC1up/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativePtEC1dn",&CalibReaderRewgt_JESRelativePtEC1dn,"CalibReaderRewgt_JESRelativePtEC1dn/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativePtEC2up",&CalibReaderRewgt_JESRelativePtEC2up,"CalibReaderRewgt_JESRelativePtEC2up/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativePtEC2dn",&CalibReaderRewgt_JESRelativePtEC2dn,"CalibReaderRewgt_JESRelativePtEC2dn/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativePtHFup",&CalibReaderRewgt_JESRelativePtHFup,"CalibReaderRewgt_JESRelativePtHFup/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativePtHFdn",&CalibReaderRewgt_JESRelativePtHFdn,"CalibReaderRewgt_JESRelativePtHFdn/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeStatECup",&CalibReaderRewgt_JESRelativeStatECup,"CalibReaderRewgt_JESRelativeStatECup/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeStatECdn",&CalibReaderRewgt_JESRelativeStatECdn,"CalibReaderRewgt_JESRelativeStatECdn/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeStatFSRup",&CalibReaderRewgt_JESRelativeStatFSRup,"CalibReaderRewgt_JESRelativeStatFSRup/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeStatFSRdn",&CalibReaderRewgt_JESRelativeStatFSRdn,"CalibReaderRewgt_JESRelativeStatFSRdn/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeStatHFup",&CalibReaderRewgt_JESRelativeStatHFup,"CalibReaderRewgt_JESRelativeStatHFup/F");
   outputTree->Branch("CalibReaderRewgt_JESRelativeStatHFdn",&CalibReaderRewgt_JESRelativeStatHFdn,"CalibReaderRewgt_JESRelativeStatHFdn/F");
   outputTree->Branch("CalibReaderRewgt_JESSinglePionECALup",&CalibReaderRewgt_JESSinglePionECALup,"CalibReaderRewgt_JESSinglePionECALup/F");
   outputTree->Branch("CalibReaderRewgt_JESSinglePionECALdn",&CalibReaderRewgt_JESSinglePionECALdn,"CalibReaderRewgt_JESSinglePionECALdn/F");
   outputTree->Branch("CalibReaderRewgt_JESSinglePionHCALup",&CalibReaderRewgt_JESSinglePionHCALup,"CalibReaderRewgt_JESSinglePionHCALup/F");
   outputTree->Branch("CalibReaderRewgt_JESSinglePionHCALdn",&CalibReaderRewgt_JESSinglePionHCALdn,"CalibReaderRewgt_JESSinglePionHCALdn/F");
   outputTree->Branch("CalibReaderRewgt_JESTimePtEtaup",&CalibReaderRewgt_JESTimePtEtaup,"CalibReaderRewgt_JESTimePtEtaup/F");
   outputTree->Branch("CalibReaderRewgt_JESTimePtEtadn",&CalibReaderRewgt_JESTimePtEtadn,"CalibReaderRewgt_JESTimePtEtadn/F");

   outputTree->Branch("CalibReaderRewgt_LFup",&CalibReaderRewgt_LFup,"CalibReaderRewgt_LFup/F");
   outputTree->Branch("CalibReaderRewgt_LFdn",&CalibReaderRewgt_LFdn,"CalibReaderRewgt_LFdn/F");
   outputTree->Branch("CalibReaderRewgt_HFstat1up",&CalibReaderRewgt_HFstat1up,"CalibReaderRewgt_HFstat1up/F");
   outputTree->Branch("CalibReaderRewgt_HFstat1dn",&CalibReaderRewgt_HFstat1dn,"CalibReaderRewgt_HFstat1dn/F");
   outputTree->Branch("CalibReaderRewgt_HFstat2up",&CalibReaderRewgt_HFstat2up,"CalibReaderRewgt_HFstat2up/F");
   outputTree->Branch("CalibReaderRewgt_HFstat2dn",&CalibReaderRewgt_HFstat2dn,"CalibReaderRewgt_HFstat2dn/F");
   outputTree->Branch("CalibReaderRewgt_CFerr1up",&CalibReaderRewgt_CFerr1up,"CalibReaderRewgt_CFerr1up/F");
   outputTree->Branch("CalibReaderRewgt_CFerr1dn",&CalibReaderRewgt_CFerr1dn,"CalibReaderRewgt_CFerr1dn/F");
   outputTree->Branch("CalibReaderRewgt_CFerr2up",&CalibReaderRewgt_CFerr2up,"CalibReaderRewgt_CFerr2up/F");
   outputTree->Branch("CalibReaderRewgt_CFerr2dn",&CalibReaderRewgt_CFerr2dn,"CalibReaderRewgt_CFerr2dn/F");
   outputTree->Branch("CalibReaderRewgt_HFup",&CalibReaderRewgt_HFup,"CalibReaderRewgt_HFup/F");
   outputTree->Branch("CalibReaderRewgt_HFdn",&CalibReaderRewgt_HFdn,"CalibReaderRewgt_HFdn/F");
   outputTree->Branch("CalibReaderRewgt_LFstat1up",&CalibReaderRewgt_LFstat1up,"CalibReaderRewgt_LFstat1up/F");
   outputTree->Branch("CalibReaderRewgt_LFstat1dn",&CalibReaderRewgt_LFstat1dn,"CalibReaderRewgt_LFstat1dn/F");
   outputTree->Branch("CalibReaderRewgt_LFstat2up",&CalibReaderRewgt_LFstat2up,"CalibReaderRewgt_LFstat2up/F");
   outputTree->Branch("CalibReaderRewgt_LFstat2dn",&CalibReaderRewgt_LFstat2dn,"CalibReaderRewgt_LFstat2dn/F");
     
   //HOT tagger
   outputTree->Branch("topJet1Index_HOTTaggerCalc",&topJet1Index_HOTTaggerCalc);
   outputTree->Branch("topJet2Index_HOTTaggerCalc",&topJet2Index_HOTTaggerCalc);
   outputTree->Branch("topJet3Index_HOTTaggerCalc",&topJet3Index_HOTTaggerCalc);
   outputTree->Branch("topNAK4_HOTTaggerCalc",&topNAK4_HOTTaggerCalc,"topNAK4_HOTTaggerCalc/I");
   outputTree->Branch("topNtops_HOTTaggerCalc",&topNtops_HOTTaggerCalc,"topNtops_HOTTaggerCalc/I");
   outputTree->Branch("topDRmax_HOTTaggerCalc",&topDRmax_HOTTaggerCalc);
   outputTree->Branch("topDThetaMax_HOTTaggerCalc",&topDThetaMax_HOTTaggerCalc);
   outputTree->Branch("topDThetaMin_HOTTaggerCalc",&topDThetaMin_HOTTaggerCalc);
   outputTree->Branch("topDiscriminator_HOTTaggerCalc",&topDiscriminator_HOTTaggerCalc);
   outputTree->Branch("topEta_HOTTaggerCalc",&topEta_HOTTaggerCalc);
   outputTree->Branch("topMass_HOTTaggerCalc",&topMass_HOTTaggerCalc);
   outputTree->Branch("topNconstituents_HOTTaggerCalc",&topNconstituents_HOTTaggerCalc);
   outputTree->Branch("topPhi_HOTTaggerCalc",&topPhi_HOTTaggerCalc);
   outputTree->Branch("topPt_HOTTaggerCalc",&topPt_HOTTaggerCalc);
   outputTree->Branch("topType_HOTTaggerCalc",&topType_HOTTaggerCalc);
   outputTree->Branch("NresolvedTops1pFake",&NresolvedTops1pFake,"NresolvedTops1pFake/I");
   outputTree->Branch("NresolvedTops2pFake",&NresolvedTops2pFake,"NresolvedTops2pFake/I");
   outputTree->Branch("NresolvedTops5pFake",&NresolvedTops5pFake,"NresolvedTops5pFake/I");
   outputTree->Branch("NresolvedTops10pFake",&NresolvedTops10pFake,"NresolvedTops10pFake/I");
  
  // ----------------------------------------------------------------------------
  // Define and initialize objects / cuts / efficiencies
  // ----------------------------------------------------------------------------

   // basic cuts
   float metCut=20;
   float htCut=0;
   int   nAK8jetsCut=0;
   float lepPtCut=25;
   float elEtaCut=2.1;
   float muEtaCut=2.1;
   int   njetsCut=4;
   float jetPtCut=30;
   float jetEtaCut=2.4;
   float ak8EtaCut=2.4;
   float ak8PtCut=200;

   // counters
   int n_vectorSizeMismatch = 0;

   int n_jetstotal     = 0;
   int n_jetsnearlep   = 0;
   int n_jetspassed    = 0;
   int npass_ThreeJets = 0;
   int npass_trigger   = 0;
   int npass_mu500     = 0;
   int npass_met       = 0;
   int npass_ht        = 0;
   int npass_Njets     = 0;
   int npass_nAK8jets  = 0;
   int npass_nHjets    = 0;
   int npass_lepPt     = 0;
   int npass_ElEta     = 0;
   int npass_MuEta     = 0;
   int npass_all       = 0;
   int Nelectrons      = 0;
   int Nmuons          = 0;

   // Lorentz vectors
   TLorentzVector jet_lv;
   TLorentzVector bjet_lv;
   TLorentzVector wjet1_lv;
   TLorentzVector tjet1_lv;
   TLorentzVector lepton_lv;
   TLorentzVector ak8_lv;

   // Polynominals for WJets HT scaling
   TF1 *poly2 = new TF1("poly2","max([6],[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x)",100,5000);
   poly2->SetParameter(0,    0.998174);  
   poly2->SetParameter(1, 8.40861e-05); 
   poly2->SetParameter(2,-6.63274e-07);
   poly2->SetParameter(3, 4.09272e-10); 
   poly2->SetParameter(4,-9.50233e-14); 
   poly2->SetParameter(5, 7.59648e-18); 
   poly2->SetParameter(6,0.402806);
   
   TF1 *poly2U = new TF1("poly2U","max([6],[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x)",100,5000);
   poly2U->SetParameter(0,    0.998164);  
   poly2U->SetParameter(1, 8.51652e-05); 
   poly2U->SetParameter(2,-6.62919e-07);
   poly2U->SetParameter(3,  4.1205e-10); 
   poly2U->SetParameter(4,-9.57795e-14); 
   poly2U->SetParameter(5, 7.67371e-18); 
   poly2U->SetParameter(6,0.454456);
   
   TF1 *poly2D = new TF1("poly2D","max([6],[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x)",100,5000);
   poly2D->SetParameter(0,    0.998183);  
   poly2D->SetParameter(1, 8.30069e-05); 
   poly2D->SetParameter(2,-6.63629e-07);
   poly2D->SetParameter(3, 4.06495e-10); 
   poly2D->SetParameter(4,-9.42671e-14); 
   poly2D->SetParameter(5, 7.51924e-18); 
   poly2D->SetParameter(6,0.351156);

   // W/t TAGGING EFFICIENCIES NEED TO BE UPDATED!!!!!
   // W tagging efficiencies. Assumes each signal mass uses the same pT bins but has unique values.
   std::vector<float> ptRangeTpTp, ptRangeTTbar;
   float ptminTTbar[12] = {175,200,250,300,350,400,450,500,550,600,700,800};
   for (int i=0;i<12;++i) ptRangeTTbar.push_back(ptminTTbar[i]);

   std::vector<std::vector<float>> SignalEff;
   std::vector<std::vector<float>> SignalEffPuppi;
   
   float TTbarEff[12] = {0.468, 0.882, 0.941, 0.928, 0.915, 0.903, 0.892, 0.878, 0.868, 0.853, 0.828, 0.796};
   float STtEff[12] = {0.437, 0.876, 0.941, 0.927, 0.914, 0.900, 0.890, 0.855, 0.864, 0.842, 0.827, 0.849};
   float STtWEff[12] = {0.453, 0.885, 0.948, 0.938, 0.930, 0.927, 0.921, 0.916, 0.920, 0.910, 0.895, 0.875};
   float WVEff[12] = {0.440, 0.882, 0.945, 0.937, 0.924, 0.926, 0.909, 0.913, 0.892, 0.902, 0.881, 0.870};

   float STtEffPuppi[12] = {0.520, 0.901, 0.946, 0.933, 0.923, 0.906, 0.889, 0.869, 0.882, 0.855, 0.841, 0.875};
   float STtWEffPuppi[12] = {0.543, 0.909, 0.952, 0.945, 0.940, 0.938, 0.932, 0.929, 0.933, 0.924, 0.916, 0.900};
   float TTbarEffPuppi[12] = {0.562, 0.905, 0.944, 0.934, 0.923, 0.912, 0.899, 0.888, 0.883, 0.872, 0.841, 0.816};
   float WVEffPuppi[12] = {0.527, 0.909, 0.954, 0.946, 0.939, 0.939, 0.922, 0.929, 0.910, 0.928, 0.901, 0.914};

   // Top tagging efficiencies
   std::vector<float> ptRangeTpTpTop, ptRangeTTbarTop;
   float ptminTTbarTop[9] = {400,450,500,550,600,700,800,1000,1200};
   for (int i=0;i<9;++i) ptRangeTTbarTop.push_back(ptminTTbarTop[i]);
   float ptminTpTpTop[9] = {400,450,500,550,600,700,800,1000,1200};//X53X53
   for (int i=0;i<9;++i) ptRangeTpTpTop.push_back(ptminTpTpTop[i]);
   float TTbarEffTop[9] = {0.710731908673,0.911246812099,0.946258231276,0.950948087531,0.952308954524,0.947643078335,0.94133549608,0.929334428924,0.924661246612};
   float STEffTop[9] = {0.691643703728,0.906228313671,0.941468696261,0.948002294016,0.947021636542,0.941775836972,0.933049300061,0.936046511628,0.933862433862};

   std::vector<float> ptRangeTopTagginSF;
   float ptminTopTagginSF[5] = {300,400,480,600,1100};
   for (int i=0;i<5;++i) ptRangeTopTagginSF.push_back(ptminTopTagginSF[i]);
   float topTagginSF[5]   = {0.98331112,0.96821666,0.95967776,1.02111010,1.02111010};
   float topTagginSFup[5] = {1.02416270,0.99713147,0.99533176,1.06593850,1.11076690};
   float topTagginSFdn[5] = {0.94245958,0.93930185,0.92402375,0.97628158,0.93145317};
   
   
  // ----------------------------------------------------------------------------
  // RUN THE EVENT LOOP
  // ----------------------------------------------------------------------------

   cout << "RUN CONFIG: isMC = " << isMC << endl;
   cout << "isSig = " << isSig << ", SigMass = " << SigMass << endl;
   cout << "For W's: isTT = " << isTT << ", isSTt = " << isSTt << ", isSTtW = " << isSTtW << endl;
   cout << "Fot jets & PDF: isTOP = " << isTOP << ", isMadgraphBkg = " << isMadgraphBkg << endl;
   cout << "Pileup index: " << pileupIndex << endl;
   
   Long64_t nentries = inputTree->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = inputTree->GetEntry(jentry);   nbytes += nb;
      if (Cut(ientry) != 1) continue;
      
      //      if (ientry > 5000) continue;
      
      if(jentry % 1000 ==0) std::cout<<"Completed "<<jentry<<" out of "<<nentries<<" events"<<std::endl;

      // ----------------------------------------------------------------------------
      // Filter input file by mass or decay
      // ----------------------------------------------------------------------------

      if(isTTincMtt0to700 && ttbarMass_TTbarMassCalc >= 700.) continue;
      if(isTTincMtt0to1000 && ttbarMass_TTbarMassCalc >= 1000.) continue;
      if(isTTincMtt700to1000 && (ttbarMass_TTbarMassCalc < 700. || ttbarMass_TTbarMassCalc >= 1000.)) continue;
      if(isTTincMtt1000toInf && ttbarMass_TTbarMassCalc < 1000.) continue;

      // ----------------------------------------------------------------------------
      // Assign bool depending on gen level HT and Jet multiplicity for ttbar sample stitching
      // ----------------------------------------------------------------------------
      // HT calculated from jets with pT>30 and |eta|<2.4 > 500 GeV
      // Jet multiplicity (jet pT>30) >= 9
      isHTgt500Njetge9 = 0;
      if(isTT){
      double genHT = 0;
      int Ngenjet = 0;
      for(unsigned int ijet=0; ijet < genJetPtNoClean_MultiLepCalc->size(); ijet++){
        if(genJetPtNoClean_MultiLepCalc->at(ijet) > 30) Ngenjet+=1;
        if(genJetPtNoClean_MultiLepCalc->at(ijet) > 30 && fabs(genJetEtaNoClean_MultiLepCalc->at(ijet)) < 2.4) genHT+=genJetPtNoClean_MultiLepCalc->at(ijet);
      }
      if(genHT>500 && Ngenjet>=9) {isHTgt500Njetge9 = 1;}
      }
      
      // ----------------------------------------------------------------------------
      // Assign as electron or muon event
      // ----------------------------------------------------------------------------
      // **** If we had saveLooseLeps == True in TTsinglelep_MC_cfg.py then this needs work!
      // **** Choose an electron and muon ID for "tight" and find events with only 1 tight el or mu by that definition
      isElectron = 0;
      isMuon = 0;
      if(elPt_MultiLepCalc->size()==0 && muPt_MultiLepCalc->size()==0){std::cout << "got no leptons, something wrong" << std::endl; continue;}
      if(elPt_MultiLepCalc->size()>0 && muPt_MultiLepCalc->size()>0) std::cout << "got el and mu, something wrong" << std::endl;
      if(elPt_MultiLepCalc->size()>1 || muPt_MultiLepCalc->size()>1) std::cout << "got more than one el or mu, something wrong" << std::endl;
      if(elPt_MultiLepCalc->size()>0 && muPt_MultiLepCalc->size()==0) {isElectron = 1; isMuon = 0;}
      if(elPt_MultiLepCalc->size()==0 && muPt_MultiLepCalc->size()>0) {isElectron = 0; isMuon = 1;}
      if(isElectron==0 && isMuon==0){std::cout << "got no leptons, something wrong" << std::endl; continue;}

      if(isSM && isElectron == 1) continue;
      if(isSE && isMuon == 1) continue;

      double leppt = 0;
      double lepeta = 0;
      if(isElectron){leppt = elPt_MultiLepCalc->at(0); lepeta = elEta_MultiLepCalc->at(0);}
      if(isMuon){leppt = muPt_MultiLepCalc->at(0); lepeta = muEta_MultiLepCalc->at(0);}

      // ----------------------------------------------------------------------------
      // Pileup weight calculation
      // ----------------------------------------------------------------------------

      pileupWeight = 1.0;
      pileupWeightUp = 1.0;
      pileupWeightDown = 1.0;
	
      if(isMC){
	if(nTrueInteractions_MultiLepCalc > 99) nTrueInteractions_MultiLepCalc = 99;
        if(nTrueInteractions_MultiLepCalc > 79 && isSig) nTrueInteractions_MultiLepCalc = 79;
	if(nTrueInteractions_MultiLepCalc < 0) nTrueInteractions_MultiLepCalc = 0;
	if(pileupIndex < 0 || pileupIndex > 60){
	  std::cout << "I don't know this pileup sample, using TTToSemiLeptonic's" << std::endl;
	  pileupIndex = 26;
	}
	pileupWeight = pileupweight[pileupIndex][nTrueInteractions_MultiLepCalc];
	pileupWeightUp = pileupweightUp[pileupIndex][nTrueInteractions_MultiLepCalc];
	pileupWeightDown = pileupweightDn[pileupIndex][nTrueInteractions_MultiLepCalc];
      }

      // ----------------------------------------------------------------------------
      // Assign Lepton scale factor or efficiency weights, save trigger pass/fail
      // ----------------------------------------------------------------------------

      DataPastTrigger = 0;
      MCPastTrigger = 0;
      EGammaGsfSF = 1.0;
      lepIdSF = 1.0;
      triggerSF = 1.0;
      isoSF = 1.0;
      if(isMC){ //MC triggers check
	if(isElectron){
	  std::string string_ORa = "Ele28_eta2p1_WPTight_Gsf_HT150";
	  std::string string_ORb = "Ele35_WPTight_Gsf";
	  for(unsigned int itrig=0; itrig < vsSelMCTriggersEl_MultiLepCalc->size(); itrig++){
	    if(vsSelMCTriggersEl_MultiLepCalc->at(itrig).find(string_ORa) != std::string::npos && viSelMCTriggersEl_MultiLepCalc->at(itrig) > 0) MCPastTrigger = 1;
	    if(vsSelMCTriggersEl_MultiLepCalc->at(itrig).find(string_ORb) != std::string::npos && viSelMCTriggersEl_MultiLepCalc->at(itrig) > 0) MCPastTrigger = 1;
	  }
	  // Gsf Tracking scale factor: http://fcouderc.web.cern.ch/fcouderc/EGamma/scaleFactors/Moriond17/approval/RECO/passingRECO/egammaEffi.txt_egammaPlots.pdf
	  if (leppt < 45) {
	    if (lepeta < -2.0) EGammaGsfSF = 0.977;
	    else if (lepeta < -1.566) EGammaGsfSF = 0.982;
	    else if (lepeta < -1.442) EGammaGsfSF = 0.948;
	    else if (lepeta < -1.0) EGammaGsfSF = 0.969;
	    else if (lepeta < -0.5) EGammaGsfSF = 0.977;
	    else if (lepeta < 0.5) EGammaGsfSF = 0.970;
	    else if (lepeta < 1.0) EGammaGsfSF = 0.972;
	    else if (lepeta < 1.442) EGammaGsfSF = 0.970;
	    else if (lepeta < 1.566) EGammaGsfSF = 0.958;
	    else EGammaGsfSF = 0.980; }
	  else if (leppt < 75) {
	    if (lepeta < -2.0) EGammaGsfSF = 0.984;
	    else if (lepeta < -1.566) EGammaGsfSF = 0.982;
	    else if (lepeta < -1.442) EGammaGsfSF = 0.971;
	    else if (lepeta < -1.0) EGammaGsfSF = 0.976;
	    else if (lepeta < 0.0) EGammaGsfSF = 0.980;
	    else if (lepeta < 0.5) EGammaGsfSF = 0.978;
	    else if (lepeta < 1.0) EGammaGsfSF = 0.979;
	    else if (lepeta < 1.442) EGammaGsfSF = 0.977;
	    else if (lepeta < 1.566) EGammaGsfSF = 0.964;
	    else if (lepeta < 2.0) EGammaGsfSF = 0.983;
	    else EGammaGsfSF = 0.984; }
	  else if (leppt < 100) {
	    if (lepeta < -1.566) EGammaGsfSF = 0.997;
	    else if (lepeta < -1.442) EGammaGsfSF = 1.003;
	    else if (lepeta < -1.0) EGammaGsfSF = 0.996;
	    else if (lepeta < 1.0) EGammaGsfSF = 0.992;
	    else if (lepeta < 1.442) EGammaGsfSF = 0.996;
	    else if (lepeta < 1.566) EGammaGsfSF = 1.003;
	    else EGammaGsfSF = 0.997; }
	  else {
	    if (lepeta < -1.566) EGammaGsfSF = 0.990;
	    else if (lepeta < -1.442) EGammaGsfSF = 1.010;
	    else if (lepeta < -1.0) EGammaGsfSF = 0.985;
	    else if (lepeta < -0.5) EGammaGsfSF = 0.988;
	    else if (lepeta < 0.5) EGammaGsfSF = 0.994;
	    else if (lepeta < 1.0) EGammaGsfSF = 0.988;
	    else if (lepeta < 1.442) EGammaGsfSF = 0.985;
	    else if (lepeta < 1.566) EGammaGsfSF = 1.010;
	    else EGammaGsfSF = 0.990; }
          //Scale Factor 2: https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/2017_ElectronMVA90noiso_2D.pdf
	  if (leppt < 20) {
	    if (lepeta < -2.0) lepIdSF = 0.943;
	    else if (lepeta < -1.566) lepIdSF = 0.957;
	    else if (lepeta < -1.442) lepIdSF = 1.000;
	    else if (lepeta < -0.8) lepIdSF = 1.008;
	    else if (lepeta < 0.0) lepIdSF = 0.993;
	    else if (lepeta < 0.8) lepIdSF = 0.992;
	    else if (lepeta < 1.442) lepIdSF = 0.999;
	    else if (lepeta < 1.566) lepIdSF = 1.000;
	    else if (lepeta < 2.0) lepIdSF = 0.978;
	    else lepIdSF = 0.930; }
	  else if (leppt < 35) {
	    if (lepeta < -2.0) lepIdSF = 0.926;
	    else if (lepeta < -1.566) lepIdSF = 0.937;
	    else if (lepeta < -1.442) lepIdSF = 1.000;
	    else if (lepeta < -0.8) lepIdSF = 0.964;
	    else if (lepeta < 0.0) lepIdSF = 0.981;
	    else if (lepeta < 0.8) lepIdSF = 0.981;
	    else if (lepeta < 1.442) lepIdSF = 0.963;
	    else if (lepeta < 1.566) lepIdSF = 1.000;
	    else if (lepeta < 2.0) lepIdSF = 0.943;
	    else lepIdSF = 0.918; }
	  else if (leppt < 50) {
	    if (lepeta < -2.0) lepIdSF = 0.941;
	    else if (lepeta < -1.566) lepIdSF = 0.953;
	    else if (lepeta < -1.442) lepIdSF = 1.000;
	    else if (lepeta < -0.8) lepIdSF = 0.962;
	    else if (lepeta < 0.0) lepIdSF = 0.972;
	    else if (lepeta < 0.8) lepIdSF = 0.974;
	    else if (lepeta < 1.442) lepIdSF = 0.965;
	    else if (lepeta < 1.566) lepIdSF = 1.000;
	    else if (lepeta < 2.0) lepIdSF = 0.955;
	    else lepIdSF = 0.933; }
	  else if (leppt < 100) {
	    if (lepeta < -2.0) lepIdSF = 0.948;
	    else if (lepeta < -1.566) lepIdSF = 0.967;
	    else if (lepeta < -1.442) lepIdSF = 1.000;
	    else if (lepeta < -0.8) lepIdSF = 0.968;
	    else if (lepeta < 0.0) lepIdSF = 0.979;
	    else if (lepeta < 0.8) lepIdSF = 0.975;
	    else if (lepeta < 1.442) lepIdSF = 0.970;
	    else if (lepeta < 1.566) lepIdSF = 1.000;
	    else if (lepeta < 2.0) lepIdSF = 0.971;
	    else lepIdSF = 0.938; }
	  else if (leppt < 200) {
	    if (lepeta < -2.0) lepIdSF = 0.983;
	    else if (lepeta < -1.566) lepIdSF = 0.969;
	    else if (lepeta < -1.442) lepIdSF = 1.000;
	    else if (lepeta < -0.8) lepIdSF = 0.979;
	    else if (lepeta < 0.0) lepIdSF = 0.983;
	    else if (lepeta < 0.8) lepIdSF = 0.988;
	    else if (lepeta < 1.442) lepIdSF = 0.993;
	    else if (lepeta < 1.566) lepIdSF = 1.000;
	    else if (lepeta < 2.0) lepIdSF = 0.990;
	    else lepIdSF = 0.939; }
	  else {
	    if (lepeta < -2.0) lepIdSF = 0.922;
	    else if (lepeta < -1.566) lepIdSF = 0.985;
	    else if (lepeta < -1.442) lepIdSF = 1.000;
	    else if (lepeta < -0.8) lepIdSF = 1.007;
	    else if (lepeta < 0.0) lepIdSF = 0.993;
	    else if (lepeta < 0.8) lepIdSF = 0.959;
	    else if (lepeta < 1.442) lepIdSF = 1.013;
	    else if (lepeta < 1.566) lepIdSF = 1.000;
	    else if (lepeta < 2.0) lepIdSF = 0.949;
	    else lepIdSF = 1.057; }
	  // mini isolation scale factors:  https://indico.cern.ch/event/820185/contributions/3427809/attachments/1845276/3027274/quick_update.pdf
	  if (leppt < 30){
            if (fabs(lepeta) < 0.8) isoSF = 0.99784;
            else if (fabs(lepeta) < 1.442) isoSF = 1.00109;
            else if (fabs(lepeta) < 1.566) isoSF = 1.02922;
            else if (fabs(lepeta) < 2) isoSF = 0.99356;
            else isoSF = 0.99678;}
	  else if (leppt < 40){
            if (fabs(lepeta) < 0.8) isoSF = 0.99691;
            else if (fabs(lepeta) < 1.442) isoSF = 0.99896;
            else if (fabs(lepeta) < 1.566) isoSF = 1.01636;
            else if (fabs(lepeta) < 2) isoSF = 0.99690;
            else isoSF = 0.99483;}
	  else if (leppt < 50){
            if (fabs(lepeta) < 0.8) isoSF = 0.99899;
            else if (fabs(lepeta) < 1.442) isoSF = 0.99797;
            else if (fabs(lepeta) < 1.566) isoSF = 1.00763;
            else if (fabs(lepeta) < 2) isoSF = 0.99797;
            else isoSF = 0.99593;}
	  else if (leppt < 60){
            if (fabs(lepeta) < 0.8) isoSF = 0.99899;
            else if (fabs(lepeta) < 1.442) isoSF = 0.99899;
            else if (fabs(lepeta) < 1.566) isoSF = 1.02320;
            else if (fabs(lepeta) < 2) isoSF = 0.99899;
            else isoSF = 1.00000;}
	  else if (leppt < 100){
            if (fabs(lepeta) < 0.8) isoSF = 1.00000;
            else if (fabs(lepeta) < 1.442) isoSF = 1.00000;
            else if (fabs(lepeta) < 1.566) isoSF = 1.03027;
            else if (fabs(lepeta) < 2) isoSF = 1.00201;
            else isoSF = 1.00000;}
	  else if (leppt < 200){
            if (fabs(lepeta) < 0.8) isoSF = 1.00201;
            else if (fabs(lepeta) < 1.442) isoSF = 1.00201;
            else if (fabs(lepeta) < 1.566) isoSF = 1.01139;
            else if (fabs(lepeta) < 2) isoSF = 1.00100;
            else isoSF = 1.00100;}
	  else{
            if (fabs(lepeta) < 0.8) isoSF = 1.00100;
            else if (fabs(lepeta) < 1.442) isoSF = 1.00000;
            else if (fabs(lepeta) < 1.566) isoSF = 0.99598;
            else if (fabs(lepeta) < 2) isoSF = 1.00100;
            else isoSF = 0.99900;}
	  // Trigger Scale Factors, SF2017B_Bkg_LepPtEta_EOR.png & SF2017CDEF_Bkg_LepPtEta_EOR.png
	  float trigSFB = 1.0;
	  float trigSFCDEF = 1.0;
	  if (fabs(lepeta) < 0.8){
	    if (leppt < 50) {trigSFB = 1.0; trigSFCDEF = 1.0;}
	    else if (leppt < 55) {trigSFB = 0.800; trigSFCDEF = 1.009;}
	    else if (leppt < 60) {trigSFB = 0.797; trigSFCDEF = 1.000;}
	    else if (leppt < 70) {trigSFB = 0.796; trigSFCDEF = 1.003;}
	    else if (leppt < 100) {trigSFB = 0.795; trigSFCDEF = 1.016;}
	    else if (leppt < 200) {trigSFB = 0.780; trigSFCDEF = 1.006;}
	    else {trigSFB = 0.785; trigSFCDEF = 0.988;}
	  }else if (fabs(lepeta) < 1.442){
            if (leppt < 50) {trigSFB = 1.0; trigSFCDEF = 1.0;}
            else if (leppt < 55) {trigSFB = 0.824; trigSFCDEF = 1.007;}
            else if (leppt < 60) {trigSFB = 0.795; trigSFCDEF = 1.024;}
            else if (leppt < 70) {trigSFB = 0.727; trigSFCDEF = 1.015;}
            else if (leppt < 100) {trigSFB = 0.764; trigSFCDEF = 0.991;}
            else if (leppt < 200) {trigSFB = 0.783; trigSFCDEF = 0.999;}
            else {trigSFB = 0.756; trigSFCDEF = 0.962;}
	  }else if (fabs(lepeta) < 1.566) {trigSFB = 1.0; trigSFCDEF = 1.0;}
	  else if (fabs(lepeta) < 2.0){
            if (leppt < 50) {trigSFB = 1.0; trigSFCDEF = 1.0;}
            else if (leppt < 55) {trigSFB = 0.764; trigSFCDEF = 0.952;}
            else if (leppt < 60) {trigSFB = 0.685; trigSFCDEF = 0.984;}
            else if (leppt < 70) {trigSFB = 0.764; trigSFCDEF = 0.972;}
            else if (leppt < 100) {trigSFB = 0.780; trigSFCDEF = 0.940;}
            else if (leppt < 200) {trigSFB = 0.693; trigSFCDEF = 0.938;}
            else {trigSFB = 0.562; trigSFCDEF = 0.726;}
	  }else {
            if (leppt < 50) {trigSFB = 1.0; trigSFCDEF = 1.0;}
            else if (leppt < 55) {trigSFB = 0.713; trigSFCDEF = 1.022;}
            else if (leppt < 60) {trigSFB = 0.773; trigSFCDEF = 1.027;}
            else if (leppt < 70) {trigSFB = 0.670; trigSFCDEF = 1.031;}
            else if (leppt < 100) {trigSFB = 0.868; trigSFCDEF = 1.088;}
            else if (leppt < 200) {trigSFB = 0.828; trigSFCDEF = 1.041;}
            else {trigSFB = 0.562; trigSFCDEF = 0.814;}
	  }
	  triggerSF = (4.823*trigSFB + 36.734*trigSFCDEF)/41.557;
	}
	if(isMuon){
	  std::string string_ORb = "IsoMu24";
	  std::string string_ORc = "Mu50";
	  for(unsigned int itrig=0; itrig < vsSelMCTriggersMu_MultiLepCalc->size(); itrig++){
	    if(vsSelMCTriggersMu_MultiLepCalc->at(itrig).find(string_ORb) != std::string::npos && viSelMCTriggersMu_MultiLepCalc->at(itrig) > 0) MCPastTrigger = 1;
	    if(vsSelMCTriggersMu_MultiLepCalc->at(itrig).find(string_ORc) != std::string::npos && viSelMCTriggersMu_MultiLepCalc->at(itrig) > 0) MCPastTrigger = 1;
	  }	  // MiniIsoTight/Tight
	  // https://wiwong.web.cern.ch/wiwong/Muon_Eff_Plots/2017_Efficiency20_miniTight_Tight_abseta/20_miniIsoTight_SF.json
	  if(leppt < 30){
            if(fabs(lepeta) < 0.9) isoSF= 1.027698;
            else if(fabs(lepeta) <  1.2) isoSF= 1.026264;
            else if(fabs(lepeta) <  2.1) isoSF= 1.018457;
            else if(fabs(lepeta) <  2.4) isoSF= 1.008373;
          }
          else if(leppt < 40){
            if(fabs(lepeta) < 0.9) isoSF= 1.013940;
            else if(fabs(lepeta) <  1.2) isoSF= 1.013229;
            else if(fabs(lepeta) <  2.1) isoSF= 1.009491;
            else if(fabs(lepeta) <  2.4) isoSF= 1.004073;
          }
          else if(leppt < 50){
            if(fabs(lepeta) < 0.9) isoSF= 1.004970;
            else if(fabs(lepeta) <  1.2) isoSF= 1.004221;
            else if(fabs(lepeta) <  2.1) isoSF= 1.003037;
            else if(fabs(lepeta) <  2.4) isoSF= 1.001450;
          }
          else if(leppt < 60){
            if(fabs(lepeta) < 0.9) isoSF= 1.002595;
            else if(fabs(lepeta) <  1.2) isoSF= 1.002115;
            else if(fabs(lepeta) <  2.1) isoSF= 1.001620;
            else if(fabs(lepeta) <  2.4) isoSF= 1.001054;
          }
          else if(leppt < 120){
            if(fabs(lepeta) < 0.9) isoSF= 1.000577;
            else if(fabs(lepeta) <  1.2) isoSF= 1.000272;
            else if(fabs(lepeta) <  2.1) isoSF= 1.001289;
            else if(fabs(lepeta) <  2.4) isoSF= 1.001343;
          }
	  else{
	    if(fabs(lepeta) < 0.9) isoSF= 0.999416;
            else if(fabs(lepeta) <  1.2) isoSF= 0.998657;
            else if(fabs(lepeta) <  2.1) isoSF= 0.998686;
            else if(fabs(lepeta) <  2.4) isoSF= 0.996060;
	  }
	  
	  // Cut based tight id
	  // Scale Factor 3: https://twiki.cern.ch/twiki/pub/CMS/MuonReferenceEffs2017/RunBCDEF_SF_ID.json
	  if (fabs(lepeta) < 0.90) {
	    if (leppt < 25.0) lepIdSF = 0.9910777627756951;
	    else if (leppt < 30.0) lepIdSF = 0.987410468262084;
	    else if (leppt < 40.0) lepIdSF = 0.9907753279135898;
	    else if (leppt < 50.0) lepIdSF = 0.9892483588952047;
	    else if (leppt < 60.0) lepIdSF = 0.9855545160334763;
	    else lepIdSF = 0.9898057377093389; 
	  }
	  else if (fabs(lepeta) < 1.20) {
	    if (leppt < 25.0) lepIdSF = 0.9927389275515244;
	    else if (leppt < 30.0) lepIdSF = 0.985063939762512;
	    else if (leppt < 40.0) lepIdSF = 0.9865359464182247;
	    else if (leppt < 50.0) lepIdSF = 0.984913093101493;
	    else if (leppt < 60.0) lepIdSF = 0.9839056384760008;
	    else lepIdSF = 0.984060403143468; 
	  }  
	  else if (fabs(lepeta) < 2.10) {
	    if (leppt < 25.0) lepIdSF = 0.9924252719877384;
	    else if (leppt < 30.0) lepIdSF = 0.9890884461284933;
	    else if (leppt < 40.0) lepIdSF = 0.9946469069883841;
	    else if (leppt < 50.0) lepIdSF = 0.9926528825155183;
	    else if (leppt < 60.0) lepIdSF = 0.9906364222943529;
	    else lepIdSF = 0.9920464322143979; 
	  }
	  else {
	    if (leppt < 25.0) lepIdSF = 0.9758095839531763;
	    else if (leppt < 30.0) lepIdSF = 0.9745153594179884;
	    else if (leppt < 40.0) lepIdSF = 0.9787410500158746;
	    else if (leppt < 50.0) lepIdSF = 0.978189122919501;
	    else if (leppt < 60.0) lepIdSF = 0.9673568416097894;
	    else lepIdSF = 0.9766311856731202; 
	  }

	  float triggerSFB = 1.0;
          float triggerSFCDEF = 1.0;
	  if (fabs(lepeta) < 0.90){
	    if (leppt < 50.0){
	      triggerSFB = 1.0;
	      triggerSFCDEF = 1.027;
	    }
	    else if (leppt < 55.0){
	      triggerSFB = 0.872;
	      triggerSFCDEF = 1.014;
	    }
	    else if (leppt < 60.0){
	      triggerSFB = 0.978;
	      triggerSFCDEF = 1.026;
	    }
	    else if (leppt < 70.0){
	      triggerSFB = 1.003;
	      triggerSFCDEF = 1.021;
	    }
	    else if (leppt < 100){
	      triggerSFB = 0.990;
	      triggerSFCDEF = 1.013;
	    }
	    else if (leppt < 200){
	      triggerSFB = 0.979;
	      triggerSFCDEF = 1.014;
	    }
	    else{
	      triggerSFB = 0.949;
	      triggerSFCDEF = 1.006;
	    }
	  }
	  else if (fabs(lepeta) < 1.20){
	    if (leppt < 50.0){
	      triggerSFB = 1.0;
	      triggerSFCDEF = 1.020;
	    }
	    else if (leppt < 55.0){
	      triggerSFB = 0.897;
	      triggerSFCDEF = 1.015;
	    }
	    else if (leppt < 60.0){
	      triggerSFB = 1.024;
	      triggerSFCDEF = 0.995;
	    }
	    else if (leppt < 70.0){
	      triggerSFB = 0.958;
	      triggerSFCDEF = 0.995;
	    }
	    else if (leppt < 100){
	      triggerSFB = 0.990;
	      triggerSFCDEF = 1.004;
	    }
	    else if (leppt < 200){
	      triggerSFB = 0.951;
	      triggerSFCDEF = 1.000;
	    }
	    else{
	      triggerSFB = 0.884;
	      triggerSFCDEF = 0.994;
	    }
	  }
	  else if (fabs(lepeta) < 2.10){
	    if (leppt < 50.0){
	      triggerSFB = 1.0;
	      triggerSFCDEF = 1.052;
	    }
	    else if (leppt < 55.0){
	      triggerSFB = 0.800;
	      triggerSFCDEF = 1.064;
	    }
	    else if (leppt < 60.0){
	      triggerSFB = 0.936;
	      triggerSFCDEF = 1.052;
	    }
	    else if (leppt < 70.0){
	      triggerSFB = 0.993;
	      triggerSFCDEF = 1.038;
	    }
	    else if (leppt < 100){
	      triggerSFB = 0.984;
	      triggerSFCDEF = 1.041;
	    }
	    else if (leppt < 200){
	      triggerSFB = 0.963;
	      triggerSFCDEF = 1.024;
	    }
	    else{
	      triggerSFB = 0.991;
	      triggerSFCDEF = 1.018;
	    }
	  }
	  else{
	    if (leppt < 50.0){
	      triggerSFB = 1.0;
	      triggerSFCDEF = 1.109;
	    }
	    else  if (leppt < 55.0){
	      triggerSFB = 1.0;
	      triggerSFCDEF = 1.061;
	    }
	    else if (leppt < 60.0){
	      triggerSFB = 0.751;
	      triggerSFCDEF = 1.086;
	    }
	    else if (leppt < 70.0){
	      triggerSFB = 0.804;
	      triggerSFCDEF = 1.113;
	    }
	    else if (leppt < 100){
	      triggerSFB = 0.915;
	      triggerSFCDEF = 1.105;
	    }
	    else if (leppt < 200){
	      triggerSFB = 1.032;
	      triggerSFCDEF = 1.146;
	    }
	    else{
	      triggerSFB = 0.835;
	      triggerSFCDEF = 1.136;
	    }
	  }
	  triggerSF = (4.823*triggerSFB+36.734*triggerSFCDEF)/41.557;
	 
	}
	DataPastTrigger = 1;
        // Trigger SF Muon 
      }

      else{ //Data triggers check
	if(isElectron){
	  std::string string_ORa = "Ele28_eta2p1_WPTight_Gsf_HT150";
	  std::string string_ORb = "Ele35_WPTight_Gsf";
	  for(unsigned int itrig=0; itrig < vsSelTriggersEl_MultiLepCalc->size(); itrig++){
	    if(vsSelTriggersEl_MultiLepCalc->at(itrig).find(string_ORa) != std::string::npos && viSelTriggersEl_MultiLepCalc->at(itrig) > 0) DataPastTrigger = 1;
	    if(vsSelTriggersEl_MultiLepCalc->at(itrig).find(string_ORb) != std::string::npos && viSelTriggersEl_MultiLepCalc->at(itrig) > 0) DataPastTrigger = 1;
	  }
	}



	if(isMuon){
	  std::string string_ORa = "IsoMu24";
	  std::string string_ORb = "Mu50";
	  for(unsigned int itrig=0; itrig < vsSelTriggersMu_MultiLepCalc->size(); itrig++){
	    if(vsSelTriggersMu_MultiLepCalc->at(itrig).find(string_ORa) != std::string::npos && viSelTriggersMu_MultiLepCalc->at(itrig) > 0) DataPastTrigger = 1;
	    if(vsSelTriggersMu_MultiLepCalc->at(itrig).find(string_ORb) != std::string::npos && viSelTriggersMu_MultiLepCalc->at(itrig) > 0) DataPastTrigger = 1;
	  }
	}	MCPastTrigger = 1;
      }
      
      if(isMC && MCPastTrigger) npass_trigger+=1;
      if(!isMC && DataPastTrigger) npass_trigger+=1;

      // ----------------------------------------------------------------------------
      // Generator-level HT correction
      // ----------------------------------------------------------------------------      

      HTSF_Pol = 1;
      HTSF_PolUp = 1;
      HTSF_PolDn = 1;

      if(isMadgraphBkg){
	// Piece-wise splice with a flat line. Uncertainty from upper/lower error bar fits
	HTSF_Pol = poly2->Eval(HTfromHEPUEP_MultiLepCalc);
	HTSF_PolUp = poly2U->Eval(HTfromHEPUEP_MultiLepCalc);
	HTSF_PolDn = poly2D->Eval(HTfromHEPUEP_MultiLepCalc);
      }

      // ----------------------------------------------------------------------------
      // Lepton 4-vectors, calculate MT and electron trigger presel value
      // ----------------------------------------------------------------------------

      // Set lepton 4-vectors
      double lepM;
      double lepphi;
      if (isMuon){ 
	lepM = 0.105658367;
	lepphi = muPhi_MultiLepCalc->at(0);
	lepton_lv.SetPtEtaPhiM(muPt_MultiLepCalc->at(0),muEta_MultiLepCalc->at(0),muPhi_MultiLepCalc->at(0),lepM);
      }
      else{
	lepM = 0.00051099891;
	lepphi = elPhi_MultiLepCalc->at(0);
	lepton_lv.SetPtEtaPhiM(elPt_MultiLepCalc->at(0),elEta_MultiLepCalc->at(0),elPhi_MultiLepCalc->at(0),lepM);
      }      

      MT_lepMet = sqrt(2*leppt*corr_met_MultiLepCalc*(1 - cos(lepphi - corr_met_phi_MultiLepCalc)));

      // ----------------------------------------------------------------------------
      // Loop over AK4 jets for calculations and pt ordering pair
      // ----------------------------------------------------------------------------

      NJets_JetSubCalc = 0;
      AK4HT = 0;
      vector<pair<double,int>> jetptindpair;

      CalibReaderRewgt = 1.0;
      // only in JECup, JECdn samples
      CalibReaderRewgt_JESup = 1.0;
      CalibReaderRewgt_JESdn = 1.0;      
      CalibReaderRewgt_JESAbsoluteMPFBiasup = 1.0;      
      CalibReaderRewgt_JESAbsoluteMPFBiasdn = 1.0;      
      CalibReaderRewgt_JESAbsoluteScaleup = 1.0;      
      CalibReaderRewgt_JESAbsoluteScaledn = 1.0;      
      CalibReaderRewgt_JESAbsoluteStatup = 1.0;      
      CalibReaderRewgt_JESAbsoluteStatdn = 1.0;      
      CalibReaderRewgt_JESFlavorQCDup = 1.0;      
      CalibReaderRewgt_JESFlavorQCDdn = 1.0;      
      CalibReaderRewgt_JESFragmentationup = 1.0;      
      CalibReaderRewgt_JESFragmentationdn = 1.0;      
      CalibReaderRewgt_JESPileUpDataMCup = 1.0;      
      CalibReaderRewgt_JESPileUpDataMCdn = 1.0;      
      CalibReaderRewgt_JESPileUpPtBBup = 1.0;      
      CalibReaderRewgt_JESPileUpPtBBdn = 1.0;      
      CalibReaderRewgt_JESPileUpPtEC1up = 1.0;      
      CalibReaderRewgt_JESPileUpPtEC1dn = 1.0;      
      CalibReaderRewgt_JESPileUpPtEC2up = 1.0;      
      CalibReaderRewgt_JESPileUpPtEC2dn = 1.0;      
      CalibReaderRewgt_JESPileUpPtHFup = 1.0;
      CalibReaderRewgt_JESPileUpPtHFdn = 1.0;
      CalibReaderRewgt_JESPileUpPtRefup = 1.0;
      CalibReaderRewgt_JESPileUpPtRefdn = 1.0;
      CalibReaderRewgt_JESRelativeBalup = 1.0;
      CalibReaderRewgt_JESRelativeBaldn = 1.0;
      CalibReaderRewgt_JESRelativeFSRup = 1.0;
      CalibReaderRewgt_JESRelativeFSRdn = 1.0;
      CalibReaderRewgt_JESRelativeJEREC1up = 1.0;
      CalibReaderRewgt_JESRelativeJEREC1dn = 1.0;
      CalibReaderRewgt_JESRelativeJEREC2up = 1.0;
      CalibReaderRewgt_JESRelativeJEREC2dn = 1.0;
      CalibReaderRewgt_JESRelativeJERHFup = 1.0;
      CalibReaderRewgt_JESRelativeJERHFdn = 1.0;
      CalibReaderRewgt_JESRelativePtBBup = 1.0;
      CalibReaderRewgt_JESRelativePtBBdn = 1.0;
      CalibReaderRewgt_JESRelativePtEC1up = 1.0;
      CalibReaderRewgt_JESRelativePtEC1dn = 1.0;
      CalibReaderRewgt_JESRelativePtEC2up = 1.0;
      CalibReaderRewgt_JESRelativePtEC2dn = 1.0;
      CalibReaderRewgt_JESRelativePtHFup = 1.0;
      CalibReaderRewgt_JESRelativePtHFdn = 1.0;
      CalibReaderRewgt_JESRelativeStatECup = 1.0;
      CalibReaderRewgt_JESRelativeStatECdn = 1.0;
      CalibReaderRewgt_JESRelativeStatFSRup = 1.0;
      CalibReaderRewgt_JESRelativeStatFSRdn = 1.0;
      CalibReaderRewgt_JESRelativeStatHFup = 1.0;
      CalibReaderRewgt_JESRelativeStatHFdn = 1.0;
      CalibReaderRewgt_JESSinglePionECALup = 1.0;
      CalibReaderRewgt_JESSinglePionECALdn = 1.0;
      CalibReaderRewgt_JESSinglePionHCALup = 1.0;
      CalibReaderRewgt_JESSinglePionHCALdn = 1.0;
      CalibReaderRewgt_JESTimePtEtaup = 1.0;
      CalibReaderRewgt_JESTimePtEtadn = 1.0;            
      
      // b-flavor only
      CalibReaderRewgt_LFup = 1.0;
      CalibReaderRewgt_LFdn = 1.0;
      CalibReaderRewgt_HFstat1up = 1.0;
      CalibReaderRewgt_HFstat1dn = 1.0;
      CalibReaderRewgt_HFstat2up = 1.0;
      CalibReaderRewgt_HFstat2dn = 1.0;
      // c-flavor only
      CalibReaderRewgt_CFerr1up = 1.0;
      CalibReaderRewgt_CFerr1dn = 1.0;
      CalibReaderRewgt_CFerr2up = 1.0;
      CalibReaderRewgt_CFerr2dn = 1.0;
      // udsg-flavor only
      CalibReaderRewgt_HFup = 1.0;
      CalibReaderRewgt_HFdn = 1.0;
      CalibReaderRewgt_LFstat1up = 1.0;
      CalibReaderRewgt_LFstat1dn = 1.0;
      CalibReaderRewgt_LFstat2up = 1.0;
      CalibReaderRewgt_LFstat2dn = 1.0;
      for(unsigned int ijet=0; ijet < theJetPt_JetSubCalc->size(); ijet++){
	// ----------------------------------------------------------------------------
	// Basic cuts
	// ----------------------------------------------------------------------------

	if(theJetPt_JetSubCalc->at(ijet) < jetPtCut || fabs(theJetEta_JetSubCalc->at(ijet)) > jetEtaCut) continue;

	// ----------------------------------------------------------------------------
	// Counts and pt ordering pair
	// ----------------------------------------------------------------------------

	jetptindpair.push_back(std::make_pair(theJetPt_JetSubCalc->at(ijet),ijet));
	NJets_JetSubCalc+=1;
	AK4HT+=theJetPt_JetSubCalc->at(ijet);
        if(isMC){
          // ----------------------------------------------------------------------------
          // Discriminator reweighting b-tagging Scale factor
          // ----------------------------------------------------------------------------

          int flav = theJetHFlav_JetSubCalc->at(ijet);
          double ak4jetpt = theJetPt_JetSubCalc->at(ijet);
          double jeteta = theJetEta_JetSubCalc->at(ijet);
          double jetcsv = theJetDeepCSVb_JetSubCalc->at(ijet)+theJetDeepCSVbb_JetSubCalc->at(ijet);

          BTagEntryForLJMet::JetFlavor jf;
          if(abs(flav) == 5) jf = BTagEntryForLJMet::FLAV_B;
          else if(abs(flav) == 4) jf = BTagEntryForLJMet::FLAV_C;
          else jf = BTagEntryForLJMet::FLAV_UDSG;
          if(isNOMshift){
            CalibReaderRewgt *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
            if(abs(flav) == 5){
              CalibReaderRewgt_LFup *= reader.eval_auto_bounds("up_lf",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFdn *= reader.eval_auto_bounds("down_lf",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFstat1up *= reader.eval_auto_bounds("up_hfstats1",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFstat1dn *= reader.eval_auto_bounds("down_hfstats1",jf,jeteta,ak4jetpt,jetcsv);
              CalibReaderRewgt_HFstat2up *= reader.eval_auto_bounds("up_hfstats2",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFstat2dn *= reader.eval_auto_bounds("down_hfstats2",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_CFerr1up *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_CFerr1dn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_CFerr2up *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv);
              CalibReaderRewgt_CFerr2dn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFup *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFdn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFstat1up *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFstat1dn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFstat2up *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFstat2dn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
            }else if(abs(flav) == 4){
              CalibReaderRewgt_CFerr1up *= reader.eval_auto_bounds("up_cferr1",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_CFerr1dn *= reader.eval_auto_bounds("down_cferr1",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_CFerr2up *= reader.eval_auto_bounds("up_cferr2",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_CFerr2dn *= reader.eval_auto_bounds("down_cferr2",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFup *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFdn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFstat1up *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFstat1dn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFstat2up *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFstat2dn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFup *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFdn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFstat1up *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFstat1dn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFstat2up *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFstat2dn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
            }else{
              CalibReaderRewgt_HFup *= reader.eval_auto_bounds("up_hf",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFdn *= reader.eval_auto_bounds("down_hf",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFstat1up *= reader.eval_auto_bounds("up_lfstats1",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFstat1dn *= reader.eval_auto_bounds("down_lfstats1",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFstat2up *= reader.eval_auto_bounds("up_lfstats2",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFstat2dn *= reader.eval_auto_bounds("down_lfstats2",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_CFerr1up *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_CFerr1dn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_CFerr2up *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_CFerr2dn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFup *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_LFdn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFstat1up *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFstat1dn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFstat2up *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
              CalibReaderRewgt_HFstat2dn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 
            }
          }else if(isJECshift){
            if(jf!=1){            
                CalibReaderRewgt_JESup                  *= reader.eval_auto_bounds("up_jes",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESdn                  *= reader.eval_auto_bounds("down_jes",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESAbsoluteMPFBiasup   *= reader.eval_auto_bounds("up_jesAbsoluteMPFBias",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESAbsoluteMPFBiasdn   *= reader.eval_auto_bounds("down_jesAbsoluteMPFBias",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESAbsoluteScaleup     *= reader.eval_auto_bounds("up_jesAbsoluteScale",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESAbsoluteScaledn     *= reader.eval_auto_bounds("down_jesAbsoluteScale",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESAbsoluteStatup      *= reader.eval_auto_bounds("up_jesAbsoluteStat",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESAbsoluteStatdn      *= reader.eval_auto_bounds("down_jesAbsoluteStat",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESFlavorQCDup         *= reader.eval_auto_bounds("up_jesFlavorQCD",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESFlavorQCDdn         *= reader.eval_auto_bounds("down_jesFlavorQCD",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESFragmentationup     *= reader.eval_auto_bounds("up_jesFragmentation",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESFragmentationdn     *= reader.eval_auto_bounds("down_jesFragmentation",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESPileUpDataMCup      *= reader.eval_auto_bounds("up_jesPileUpDataMC",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESPileUpDataMCdn      *= reader.eval_auto_bounds("down_jesPileUpDataMC",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESPileUpPtBBup        *= reader.eval_auto_bounds("up_jesPileUpPtBB",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESPileUpPtBBdn        *= reader.eval_auto_bounds("down_jesPileUpPtBB",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESPileUpPtEC1up       *= reader.eval_auto_bounds("up_jesPileUpPtEC1",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESPileUpPtEC1dn       *= reader.eval_auto_bounds("down_jesPileUpPtEC1",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESPileUpPtEC2up       *= reader.eval_auto_bounds("up_jesPileUpPtEC2",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESPileUpPtEC2dn       *= reader.eval_auto_bounds("down_jesPileUpPtEC2",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESPileUpPtHFup        *= reader.eval_auto_bounds("up_jesPileUpPtHF",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESPileUpPtHFdn        *= reader.eval_auto_bounds("down_jesPileUpPtHF",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESPileUpPtRefup       *= reader.eval_auto_bounds("up_jesPileUpPtRef",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESPileUpPtRefdn       *= reader.eval_auto_bounds("down_jesPileUpPtRef",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeBalup       *= reader.eval_auto_bounds("up_jesRelativeBal",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeBaldn       *= reader.eval_auto_bounds("down_jesRelativeBal",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeFSRup       *= reader.eval_auto_bounds("up_jesRelativeFSR",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeFSRdn       *= reader.eval_auto_bounds("down_jesRelativeFSR",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeJEREC1up    *= reader.eval_auto_bounds("up_jesRelativeJEREC1",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeJEREC1dn    *= reader.eval_auto_bounds("down_jesRelativeJEREC1",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeJEREC2up    *= reader.eval_auto_bounds("up_jesRelativeJEREC2",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeJEREC2dn    *= reader.eval_auto_bounds("down_jesRelativeJEREC2",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeJERHFup     *= reader.eval_auto_bounds("up_jesRelativeJERHF",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeJERHFdn     *= reader.eval_auto_bounds("down_jesRelativeJERHF",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativePtBBup      *= reader.eval_auto_bounds("up_jesRelativePtBB",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativePtBBdn      *= reader.eval_auto_bounds("down_jesRelativePtBB",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativePtEC1up     *= reader.eval_auto_bounds("up_jesRelativePtEC1",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativePtEC1dn     *= reader.eval_auto_bounds("down_jesRelativePtEC1",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativePtEC2up     *= reader.eval_auto_bounds("up_jesRelativePtEC2",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativePtEC2dn     *= reader.eval_auto_bounds("down_jesRelativePtEC2",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativePtHFup      *= reader.eval_auto_bounds("up_jesRelativePtHF",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativePtHFdn      *= reader.eval_auto_bounds("down_jesRelativePtHF",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeStatECup    *= reader.eval_auto_bounds("up_jesRelativeStatEC",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeStatECdn    *= reader.eval_auto_bounds("down_jesRelativeStatEC",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeStatFSRup   *= reader.eval_auto_bounds("up_jesRelativeStatFSR",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeStatFSRdn   *= reader.eval_auto_bounds("down_jesRelativeStatFSR",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeStatHFup    *= reader.eval_auto_bounds("up_jesRelativeStatHF",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESRelativeStatHFdn    *= reader.eval_auto_bounds("down_jesRelativeStatHF",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESSinglePionECALup    *= reader.eval_auto_bounds("up_jesSinglePionECAL",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESSinglePionECALdn    *= reader.eval_auto_bounds("down_jesSinglePionECAL",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESSinglePionHCALup    *= reader.eval_auto_bounds("up_jesSinglePionHCAL",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESSinglePionHCALdn    *= reader.eval_auto_bounds("down_jesSinglePionHCAL",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESTimePtEtaup         *= reader.eval_auto_bounds("up_jesTimePtEta",jf,jeteta,ak4jetpt,jetcsv); 
                CalibReaderRewgt_JESTimePtEtadn         *= reader.eval_auto_bounds("down_jesTimePtEta",jf,jeteta,ak4jetpt,jetcsv); 
            }
            else if(jf==1){ /// In case of c-jet, give it central
                CalibReaderRewgt_JESup *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv);
                CalibReaderRewgt_JESdn *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv);
            }
            std::cout<<"jf : "<<jf<<"  jeteta : "<<jeteta<<"  ak4jetpt : "<<ak4jetpt<<"  jetcsv : "<<jetcsv<<std::endl;
            std::cout<<"up_jes : "<<CalibReaderRewgt_JESup<<std::endl;
            std::cout<<"down_jes : "<<CalibReaderRewgt_JESdn<<std::endl;	    
            std::cout<<"##############################################################"<<std::endl;
          }else if(isJERshift){
                CalibReaderRewgt *= reader.eval_auto_bounds("central",jf,jeteta,ak4jetpt,jetcsv); 	        
      
          }       
      }
      }
 
      // ----------------------------------------------------------------------------
      // Apply kinematic cuts
      // ----------------------------------------------------------------------------
	                
      int isPastHTCut = 0;
      if(AK4HT >= htCut){npass_ht+=1;isPastHTCut=1;}
      
      int isPastMETcut = 0;
      if(corr_met_MultiLepCalc > metCut){npass_met+=1;isPastMETcut=1;}

      int isPastLepPtCut = 0;
      if(leppt > lepPtCut){npass_lepPt+=1;isPastLepPtCut=1;}
      
      int isPastElEtaCut = 0;
      int isPastMuEtaCut = 0;
      if(isElectron && fabs(lepeta) < elEtaCut){npass_ElEta+=1;isPastElEtaCut=1;}
      if(isMuon && fabs(lepeta) < muEtaCut){npass_MuEta+=1;isPastMuEtaCut=1;}
      if(isElectron){Nelectrons+=1;}
      if(isMuon){Nmuons+=1;}
      
      int isPastNjetsCut = 0;
      if(NJets_JetSubCalc >= njetsCut){npass_Njets+=1;isPastNjetsCut=1;}

      // ----------------------------------------------------------------------------
      // Skip failing events
      // ----------------------------------------------------------------------------
            
      if(!(isPastMETcut && isPastHTCut && isPastNjetsCut && isPastLepPtCut && (isPastElEtaCut || isPastMuEtaCut))) continue;
      npass_all+=1;
      
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////// ONLY ON SELECTED EVENTS ////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////      

      AK4HTpMETpLepPt = 0; //ST
      AK4HTpMETpLepPt = AK4HT + corr_met_MultiLepCalc + leppt;
      
      // ----------------------------------------------------------------------------
      // Combine lepton variables into one set
      // ----------------------------------------------------------------------------

      if(isElectron){
	leptonPt_MultiLepCalc = leppt;
	leptonEta_MultiLepCalc = lepeta;
	leptonPhi_MultiLepCalc = elPhi_MultiLepCalc->at(0);
	leptonEnergy_MultiLepCalc = elEnergy_MultiLepCalc->at(0);
	leptonMiniIso_MultiLepCalc = elMiniIso_MultiLepCalc->at(0);
	leptonMVAValue_MultiLepCalc = -99.9; //elMVAValue_MultiLepCalc->at(0);
      }
      if(isMuon){
	leptonPt_MultiLepCalc = leppt;
	leptonEta_MultiLepCalc = lepeta;
	leptonPhi_MultiLepCalc = muPhi_MultiLepCalc->at(0);
	leptonEnergy_MultiLepCalc = muEnergy_MultiLepCalc->at(0);
	leptonMiniIso_MultiLepCalc = muMiniIso_MultiLepCalc->at(0);
	leptonMVAValue_MultiLepCalc = -99.9;
      }

      // ----------------------------------------------------------------------------
      // Apply pt ordering to AK4 vectors
      // ----------------------------------------------------------------------------

      std::sort(jetptindpair.begin(), jetptindpair.end(), comparepair);
      theJetPt_JetSubCalc_PtOrdered.clear();
      theJetEta_JetSubCalc_PtOrdered.clear();
      theJetPhi_JetSubCalc_PtOrdered.clear();
      theJetEnergy_JetSubCalc_PtOrdered.clear();
      theJetDeepCSVb_JetSubCalc_PtOrdered.clear();
      theJetDeepCSVbb_JetSubCalc_PtOrdered.clear();
      theJetDeepCSVc_JetSubCalc_PtOrdered.clear();
      theJetDeepCSVudsg_JetSubCalc_PtOrdered.clear();
      theJetHFlav_JetSubCalc_PtOrdered.clear();
      theJetPFlav_JetSubCalc_PtOrdered.clear();
      theJetBTag_JetSubCalc_PtOrdered.clear();
      theJetBTag_bSFup_JetSubCalc_PtOrdered.clear();
      theJetBTag_bSFdn_JetSubCalc_PtOrdered.clear();
      theJetBTag_lSFup_JetSubCalc_PtOrdered.clear();
      theJetBTag_lSFdn_JetSubCalc_PtOrdered.clear();
      for(unsigned int ijet=0; ijet < jetptindpair.size(); ijet++){
      	theJetPt_JetSubCalc_PtOrdered.push_back(theJetPt_JetSubCalc->at(jetptindpair[ijet].second));
      	theJetEta_JetSubCalc_PtOrdered.push_back(theJetEta_JetSubCalc->at(jetptindpair[ijet].second));
      	theJetPhi_JetSubCalc_PtOrdered.push_back(theJetPhi_JetSubCalc->at(jetptindpair[ijet].second));
      	theJetEnergy_JetSubCalc_PtOrdered.push_back(theJetEnergy_JetSubCalc->at(jetptindpair[ijet].second));
      	theJetDeepCSVb_JetSubCalc_PtOrdered.push_back(theJetDeepCSVb_JetSubCalc->at(jetptindpair[ijet].second));
		theJetDeepCSVbb_JetSubCalc_PtOrdered.push_back(theJetDeepCSVbb_JetSubCalc->at(jetptindpair[ijet].second));
		theJetDeepCSVc_JetSubCalc_PtOrdered.push_back(theJetDeepCSVc_JetSubCalc->at(jetptindpair[ijet].second));
		theJetDeepCSVudsg_JetSubCalc_PtOrdered.push_back(theJetDeepCSVudsg_JetSubCalc->at(jetptindpair[ijet].second));
      	theJetHFlav_JetSubCalc_PtOrdered.push_back(theJetHFlav_JetSubCalc->at(jetptindpair[ijet].second));
      	theJetPFlav_JetSubCalc_PtOrdered.push_back(theJetPFlav_JetSubCalc->at(jetptindpair[ijet].second));
		theJetBTag_JetSubCalc_PtOrdered.push_back(theJetBTag_JetSubCalc->at(jetptindpair[ijet].second));
	theJetBTag_bSFup_JetSubCalc_PtOrdered.push_back(theJetBTag_bSFup_JetSubCalc->at(jetptindpair[ijet].second));
	theJetBTag_bSFdn_JetSubCalc_PtOrdered.push_back(theJetBTag_bSFdn_JetSubCalc->at(jetptindpair[ijet].second));
	theJetBTag_lSFup_JetSubCalc_PtOrdered.push_back(theJetBTag_lSFup_JetSubCalc->at(jetptindpair[ijet].second));
	theJetBTag_lSFdn_JetSubCalc_PtOrdered.push_back(theJetBTag_lSFdn_JetSubCalc->at(jetptindpair[ijet].second));
      }

      // ----------------------------------------------------------------------------
      // AK4 Jet - lepton associations
      // ----------------------------------------------------------------------------
         
      BJetLeadPt = -99;
      BJetLeadPt_bSFup = -99;
      BJetLeadPt_bSFdn = -99;
      BJetLeadPt_lSFup = -99;
      BJetLeadPt_lSFdn = -99;
      minMleppBjet = 1e8;
      minMleppBjet_bSFup = 1e8;
      minMleppBjet_bSFdn = 1e8;
      minMleppBjet_lSFup = 1e8;
      minMleppBjet_lSFdn = 1e8;
      minMleppJet = 1e8;      
      minDPhi_MetJet = 1e8;
      minDR_lepJet = 1e8;
      ptRel_lepJet = -99;
      deltaR_lepJets.clear();
      deltaR_lepBJets.clear();
      deltaR_lepBJets_bSFup.clear();
      deltaR_lepBJets_bSFdn.clear();
      deltaR_lepBJets_lSFup.clear();
      deltaR_lepBJets_lSFdn.clear();
      deltaR_lepMinMlb = 1e8;
      deltaR_lepMinMlb_bSFup = 1e8;
      deltaR_lepMinMlb_bSFdn = 1e8;
      deltaR_lepMinMlb_lSFup = 1e8;
      deltaR_lepMinMlb_lSFdn = 1e8;
      deltaR_lepMinMlj = 1e8;
      NJetsCSV_JetSubCalc = 0;
      NJetsCSVwithSF_JetSubCalc = 0;
      NJetsCSVwithSF_JetSubCalc_bSFup = 0;
      NJetsCSVwithSF_JetSubCalc_bSFdn = 0;
      NJetsCSVwithSF_JetSubCalc_lSFup = 0;
      NJetsCSVwithSF_JetSubCalc_lSFdn = 0;
      TLorentzVector nu;
      nu.SetPtEtaPhiE(corr_met_MultiLepCalc,0,corr_met_phi_MultiLepCalc,corr_met_MultiLepCalc);

      for(unsigned int ijet=0; ijet < theJetPt_JetSubCalc_PtOrdered.size(); ijet++){
        jet_lv.SetPtEtaPhiE(theJetPt_JetSubCalc_PtOrdered.at(ijet),theJetEta_JetSubCalc_PtOrdered.at(ijet),theJetPhi_JetSubCalc_PtOrdered.at(ijet),theJetEnergy_JetSubCalc_PtOrdered.at(ijet));
	if((lepton_lv + jet_lv).M() < minMleppJet) {
	  minMleppJet = fabs((lepton_lv + jet_lv).M());
	  deltaR_lepMinMlj = jet_lv.DeltaR(lepton_lv);
	}

	deltaR_lepJets.push_back(lepton_lv.DeltaR(jet_lv));

   	if(theJetDeepCSVb_JetSubCalc_PtOrdered.at(ijet) + theJetDeepCSVbb_JetSubCalc_PtOrdered.at(ijet) > 0.4941){
          NJetsCSV_JetSubCalc += 1;
        }

	if(theJetBTag_JetSubCalc_PtOrdered.at(ijet) == 1){
	  NJetsCSVwithSF_JetSubCalc += 1;
          if(theJetPt_JetSubCalc_PtOrdered.at(ijet) > BJetLeadPt) BJetLeadPt = theJetPt_JetSubCalc_PtOrdered.at(ijet);
          deltaR_lepBJets.push_back(lepton_lv.DeltaR(jet_lv));
	  
          if((lepton_lv + jet_lv).M() < minMleppBjet) {
            minMleppBjet = fabs( (lepton_lv + jet_lv).M() );
	    deltaR_lepMinMlb = jet_lv.DeltaR(lepton_lv);
          }
	}
	if(theJetBTag_bSFup_JetSubCalc_PtOrdered.at(ijet) == 1){
	  NJetsCSVwithSF_JetSubCalc_bSFup += 1;
          if(theJetPt_JetSubCalc_PtOrdered.at(ijet) > BJetLeadPt_bSFup) BJetLeadPt_bSFup = theJetPt_JetSubCalc_PtOrdered.at(ijet);
          deltaR_lepBJets_bSFup.push_back(lepton_lv.DeltaR(jet_lv));
	  
          if((lepton_lv + jet_lv).M() < minMleppBjet_bSFup) {
            minMleppBjet_bSFup = fabs( (lepton_lv + jet_lv).M() );
	    deltaR_lepMinMlb_bSFup = jet_lv.DeltaR(lepton_lv);
          }
	}
	if(theJetBTag_bSFdn_JetSubCalc_PtOrdered.at(ijet) == 1){
	  NJetsCSVwithSF_JetSubCalc_bSFdn += 1;
          if(theJetPt_JetSubCalc_PtOrdered.at(ijet) > BJetLeadPt_bSFdn) BJetLeadPt_bSFdn = theJetPt_JetSubCalc_PtOrdered.at(ijet);
          deltaR_lepBJets_bSFdn.push_back(lepton_lv.DeltaR(jet_lv));
	  
          if((lepton_lv + jet_lv).M() < minMleppBjet_bSFdn) {
            minMleppBjet_bSFdn = fabs( (lepton_lv + jet_lv).M() );
	    deltaR_lepMinMlb_bSFdn = jet_lv.DeltaR(lepton_lv);
          }
	}
	if(theJetBTag_lSFup_JetSubCalc_PtOrdered.at(ijet) == 1){
	  NJetsCSVwithSF_JetSubCalc_lSFup += 1;
          if(theJetPt_JetSubCalc_PtOrdered.at(ijet) > BJetLeadPt_lSFup) BJetLeadPt_lSFup = theJetPt_JetSubCalc_PtOrdered.at(ijet);
          deltaR_lepBJets_lSFup.push_back(lepton_lv.DeltaR(jet_lv));
	  
          if((lepton_lv + jet_lv).M() < minMleppBjet_lSFup) {
            minMleppBjet_lSFup = fabs( (lepton_lv + jet_lv).M() );
	    deltaR_lepMinMlb_lSFup = jet_lv.DeltaR(lepton_lv);
          }
	}
	if(theJetBTag_lSFdn_JetSubCalc_PtOrdered.at(ijet) == 1){
	  NJetsCSVwithSF_JetSubCalc_lSFdn += 1;
          if(theJetPt_JetSubCalc_PtOrdered.at(ijet) > BJetLeadPt_lSFdn) BJetLeadPt_lSFdn = theJetPt_JetSubCalc_PtOrdered.at(ijet);
          deltaR_lepBJets_lSFdn.push_back(lepton_lv.DeltaR(jet_lv));
	  
          if((lepton_lv + jet_lv).M() < minMleppBjet_lSFdn) {
            minMleppBjet_lSFdn = fabs( (lepton_lv + jet_lv).M() );
	    deltaR_lepMinMlb_lSFdn = jet_lv.DeltaR(lepton_lv);
          }
	}
	
	if(jet_lv.DeltaPhi(nu) < minDPhi_MetJet) minDPhi_MetJet = jet_lv.DeltaPhi(nu);	  

 	if(deltaR_lepJets[ijet] < minDR_lepJet) {
	  minDR_lepJet = deltaR_lepJets[ijet];
	  ptRel_lepJet = lepton_lv.P()*(jet_lv.Vect().Cross(lepton_lv.Vect()).Mag()/jet_lv.P()/lepton_lv.P());
	}
      }

      // ----------------------------------------------------------------------------
      // 13TeV Top pT reweighting -- TTbarMassCalc top vectors are the wrong tops....
      // ----------------------------------------------------------------------------
      
      genTopPt = -999;
      genAntiTopPt = -999;
      topPtWeight13TeV = 1.0;
      
      if(isTT){
	int top1index = -1;
	int top2index = -1;
	for(unsigned int itop=0; itop < allTopsStatus_TTbarMassCalc->size(); itop++){
	  if(allTopsStatus_TTbarMassCalc->at(itop) != 62) continue;
	  if(top1index < 0){
	    top1index = itop;
	    continue;
	  }
	  if(top2index < 0){
	    top2index = itop;
	    continue;
	  }
	}

	float SFtop1 = TMath::Exp(0.0615-0.0005*allTopsPt_TTbarMassCalc->at(top1index));
	float SFtop2 = TMath::Exp(0.0615-0.0005*allTopsPt_TTbarMassCalc->at(top2index));
	topPtWeight13TeV = TMath::Sqrt(SFtop1*SFtop2);

	if(allTopsID_TTbarMassCalc->at(top1index) == 6){
	  genTopPt = allTopsPt_TTbarMassCalc->at(top1index);
	  genAntiTopPt = allTopsPt_TTbarMassCalc->at(top2index);
	}
	else{
	  genAntiTopPt = allTopsPt_TTbarMassCalc->at(top1index);
	  genTopPt = allTopsPt_TTbarMassCalc->at(top2index);
	}
      }

      // ----------------------------------------------------------------------------
      // HOT tagger variables -- SCALE FACTORS TO BE ADDED!!!
      // ----------------------------------------------------------------------------
     NresolvedTops1pFake = 0;
     NresolvedTops2pFake = 0;
     NresolvedTops5pFake = 0;
     NresolvedTops10pFake = 0;
     for(unsigned int itop=0; itop < topDiscriminator_HOTTaggerCalc->size(); itop++){
        if(topDiscriminator_HOTTaggerCalc->at(itop) > 0.75) NresolvedTops10pFake+=1;
        if(topDiscriminator_HOTTaggerCalc->at(itop) > 0.85) NresolvedTops5pFake+=1;
        if(topDiscriminator_HOTTaggerCalc->at(itop) > 0.92) NresolvedTops2pFake+=1;
        if(topDiscriminator_HOTTaggerCalc->at(itop) > 0.95) NresolvedTops1pFake+=1;
        }

      // ----------------------------------------------------------------------------
      // PDF and Matrix Element energy scale weights
      // ----------------------------------------------------------------------------

      std::vector<double> renorm;
      std::vector<double> pdf;
      renormWeights.clear();
      pdfWeights.clear();
      pdfNewWeights.clear();
      pdfNewNominalWeight = 1.0;
      if(isSig){
	pdfNewNominalWeight = NewPDFweights_MultiLepCalc->at(0);
	// SEEMS TO APPLY TO ALL B2G MG+PYTHIA SIGNALS. NNLO 4-FLAVOR PDF
	for(unsigned int i = 0; i < LHEweightids_MultiLepCalc->size(); i++){
	  if(i > 0 && i < 101) pdfNewWeights.push_back(NewPDFweights_MultiLepCalc->at(i));
	  if(LHEweightids_MultiLepCalc->at(i) > 1 && LHEweightids_MultiLepCalc->at(i) < 10){
	    if(LHEweightids_MultiLepCalc->at(i) == 6 || LHEweightids_MultiLepCalc->at(i) == 8) continue;
	    renormWeights.push_back(LHEweights_MultiLepCalc->at(i));
	  }
	  if(LHEweightids_MultiLepCalc->at(i) > 474 && LHEweightids_MultiLepCalc->at(i) < 575){
	    pdfWeights.push_back(LHEweights_MultiLepCalc->at(i));
	  }
	}
      }
      else if(isMadgraphBkg){
	// SEEMS TO APPLY TO OTHER MG+PYTHIA BACKGROUNDS. LEADING ORDER 5-FLAVOR PDF
	for(unsigned int i = 0; i < LHEweightids_MultiLepCalc->size(); i++){
	  if(LHEweightids_MultiLepCalc->at(i) > 1 && LHEweightids_MultiLepCalc->at(i) < 10){
	    if(LHEweightids_MultiLepCalc->at(i) == 6 || LHEweightids_MultiLepCalc->at(i) == 8) continue;
	    renorm.push_back(LHEweights_MultiLepCalc->at(i));
	    renormWeights.push_back(LHEweights_MultiLepCalc->at(i));
	  }
	  if(LHEweightids_MultiLepCalc->at(i) > 10 && LHEweightids_MultiLepCalc->at(i) < 111){
	    pdf.push_back(LHEweights_MultiLepCalc->at(i));
	    pdfWeights.push_back(LHEweights_MultiLepCalc->at(i));
	  }
	  if(LHEweightids_MultiLepCalc->at(i) == 111 || LHEweightids_MultiLepCalc->at(i) == 112){
	    alphaSWeights.push_back(LHEweights_MultiLepCalc->at(i));
	  }
	}
      }
      else{
	// SEEMS TO APPLY TO ALL POWHEG AND MC@NLO BACKGROUNDS. NLO PDFs
	for(unsigned int i = 0; i < LHEweightids_MultiLepCalc->size(); i++){
	  if(LHEweightids_MultiLepCalc->at(i) > 1001 && LHEweightids_MultiLepCalc->at(i) < 1010){
	    if(LHEweightids_MultiLepCalc->at(i) == 1006 || LHEweightids_MultiLepCalc->at(i) == 1008) continue;
	    renorm.push_back(LHEweights_MultiLepCalc->at(i));
	    renormWeights.push_back(LHEweights_MultiLepCalc->at(i));
	  }
	  if(LHEweightids_MultiLepCalc->at(i) > 2000 && LHEweightids_MultiLepCalc->at(i) < 2101){
	    pdf.push_back(LHEweights_MultiLepCalc->at(i));
	    pdfWeights.push_back(LHEweights_MultiLepCalc->at(i));
	  }
	  if(LHEweightids_MultiLepCalc->at(i) == 2101 || LHEweightids_MultiLepCalc->at(i) == 2102){
	    alphaSWeights.push_back(LHEweights_MultiLepCalc->at(i));
	  }
	}
      }

      if(renormWeights.size() == 0){
	if(isVV){
	  renormWeights.push_back(0.85);
	  renormWeights.push_back(1.15);
	  renormWeights.push_back(0.85);
	  renormWeights.push_back(0.85);
	  renormWeights.push_back(1.15);
	  renormWeights.push_back(1.15);
	}else{
	  for(int irn = 0; irn < 6; irn++){
	    renormWeights.push_back(1.0);
	  }
	}
      }
      if(pdfWeights.size() == 0){
	for(int ipdf = 0; ipdf < 100; ipdf++){
	  pdfWeights.push_back(1.0);
	}
      }

      // ----------------------------------------------------------------------------
      // DONE!! Write the tree
      // ----------------------------------------------------------------------------
      
      outputTree->Fill();
   }
   std::cout<<"Nelectrons      = "<<Nelectrons<<" / "<<nentries<<std::endl;
   std::cout<<"Npassed_ElEta   = "<<npass_ElEta<<" / "<<nentries<<std::endl;
   std::cout<<"Nmuons          = "<<Nmuons<<" / "<<nentries<<std::endl;
   std::cout<<"Npassed_MuEta   = "<<npass_MuEta<<" / "<<nentries<<std::endl;
   std::cout<<"Npassed_nAK8Jets= "<<npass_nAK8jets<<" / "<<nentries<<std::endl;
   std::cout<<"Npassed_Trigger = "<<npass_trigger<<" / "<<nentries<<std::endl;
   std::cout<<"Npassed_MET     = "<<npass_met<<" / "<<nentries<<std::endl;
   std::cout<<"Npassed_lepPt   = "<<npass_lepPt<<" / "<<nentries<<std::endl;
   std::cout<<"Npassed_HT      = "<<npass_ht<<" / "<<nentries<<std::endl;
   std::cout<<"npass_Njets     = "<<npass_Njets<<" / "<<nentries<<std::endl;
   std::cout<<"Npassed_ALL     = "<<npass_all<<" / "<<nentries<<std::endl;

   outputTree->Write();
   delete outputTree;
   delete poly2;
   delete poly2U;
   delete poly2D;

}
