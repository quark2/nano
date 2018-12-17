#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "nano/analysis/interface/topEventSelectionSL.h"
#include "nano/analysis/interface/hadAnalyser.h"
#include "nano/analysis/interface/HadTruthEvents.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"


using namespace std;


class testForJECAnalyser : public topEventSelectionSL {
private:
  
public:
  Float_t Jet_pt_jer_nom[35];
  
  Float_t Jet_mass_jer_up[35];
  Float_t Jet_mass_jer_dn[35];
  Float_t Jet_mass_jes_up[35];
  Float_t Jet_mass_jes_dn[35];
  
  Float_t Jet_pt_jer_up[35];
  Float_t Jet_pt_jer_dn[35];
  Float_t Jet_pt_jes_up[35];
  Float_t Jet_pt_jes_dn[35];
  
  TBranch *b_Jet_pt_jer_nom;
  
  TBranch *b_Jet_mass_jer_up;
  TBranch *b_Jet_mass_jer_dn;
  TBranch *b_Jet_mass_jes_up;
  TBranch *b_Jet_mass_jes_dn;
  
  TBranch *b_Jet_pt_jer_up;
  TBranch *b_Jet_pt_jer_dn;
  TBranch *b_Jet_pt_jes_up;
  TBranch *b_Jet_pt_jes_dn;
  
  std::vector<Float_t> b_jetM, b_jetPt, b_jetEta, b_jetPhi;
  std::vector<Float_t> b_jetPtNom;
  std::vector<Float_t> b_jetM_Unc_pp, b_jetPt_Unc_pp;
  std::vector<Float_t> b_jetM_Unc_an, b_jetPt_Unc_an;
  
  std::vector<Int_t> b_GenMatched;
  
  Long64_t m_nIdxEntry;

  void setOutput(std::string outputName);

  testForJECAnalyser(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t sle = false, Bool_t slm = false, UInt_t unFlag = 0) :
    topEventSelectionSL(tree, had, hadTruth, isMC, sle, slm, unFlag)
  {
    //fChain->SetBranchAddress("Jet_pt_nom", Jet_pt_jer_nom, &b_Jet_pt_jer_nom);
    fChain->SetBranchAddress("JetPreSmeared_pt", Jet_pt_jer_nom, &b_Jet_pt_jer_nom);
    
    fChain->SetBranchAddress("Jet_mass_jerUp", Jet_mass_jer_up, &b_Jet_mass_jer_up);
    fChain->SetBranchAddress("Jet_mass_jerDown", Jet_mass_jer_dn, &b_Jet_mass_jer_dn);
    fChain->SetBranchAddress("Jet_mass_jesTotalUp", Jet_mass_jes_up, &b_Jet_mass_jes_up);
    fChain->SetBranchAddress("Jet_mass_jesTotalDown", Jet_mass_jes_dn, &b_Jet_mass_jes_dn);
    
    fChain->SetBranchAddress("Jet_pt_jerUp", Jet_pt_jer_up, &b_Jet_pt_jer_up);
    fChain->SetBranchAddress("Jet_pt_jerDown", Jet_pt_jer_dn, &b_Jet_pt_jer_dn);
    fChain->SetBranchAddress("Jet_pt_jesTotalUp", Jet_pt_jes_up, &b_Jet_pt_jes_up);
    fChain->SetBranchAddress("Jet_pt_jesTotalDown", Jet_pt_jes_dn, &b_Jet_pt_jes_dn);
    
    //cut_JetID = -1;
    //cut_JetPt = 0;
    //cut_JetEta = 100000000000;
    //cut_JetConeSizeOverlap = 0.0;
    
    cut_ElectronPt = 35;
    cut_ElectronEta = 2.1;
    cut_ElectronIDType = Electron_cutBased;
    cut_ElectronIDCut = 4;
    cut_ElectronSCEtaLower = 1.4442;
    cut_ElectronSCEtaUpper = 10000000000;
    cut_ElectronRelIso03All = 0.0588;
    
    cut_MuonIDType = Muon_tightId;
    cut_MuonPt = 26;
    cut_MuonEta = 2.4;
    cut_MuonRelIso04All = 0.06;
    
    cut_VetoElectronPt = 15;
    cut_VetoElectronEta = 2.5;
    cut_VetoElectronIDType = Electron_cutBased;
    cut_VetoElectronIDCut = 1;
    cut_VetoElectronSCEtaLower = 1.4442;
    cut_VetoElectronSCEtaUpper = 10000000000;
    cut_VetoElectronRelIso03All = 10000000000;
    
    cut_VetoMuonIDType = NULL;
    cut_VetoMuonPt = 10;
    cut_VetoMuonEta = 2.4;
    cut_VetoMuonRelIso04All = 0.2;
    cut_VetoMuonRelIso04All = 1048576;
    
    cut_GenJetPt = 30;
    cut_GenJetEta = 2.4;
    cut_GenJetConeSizeOverlap = 0.4;
    
    cut_JetID = 1;
    cut_JetPt = 40;
    cut_JetEta = 4.7;
    cut_JetConeSizeOverlap = 0.4;
    cut_JetConeSizeOverlap = 0.0;
    
    cut_BJetID = 1;
    cut_BJetPt = 40;
    cut_BJetEta = 2.4;
    cut_BJetConeSizeOverlap = 0.4;
    cut_BJetConeSizeOverlap = 0.0;
    cut_BJetTypeBTag = Jet_btagCMVA;
    cut_BJetBTagCut = 0.9432;
  }
  void resetBranch();
  
  virtual void Loop();
  virtual bool additionalConditionForJet(UInt_t nIdx, Float_t &fJetPt, Float_t &fJetEta, Float_t &fJetPhi, Float_t &fJetMass);
  
  ~testForJECAnalyser() {};
};


void testForJECAnalyser::resetBranch() {
  b_jetM.clear();
  b_jetPt.clear();
  b_jetEta.clear();
  b_jetPhi.clear();
  
  b_jetPtNom.clear();
  
  b_jetM_Unc_pp.clear();
  b_jetPt_Unc_pp.clear();
  b_jetM_Unc_an.clear();
  b_jetPt_Unc_an.clear();
  
  b_GenMatched.clear();
}


bool testForJECAnalyser::additionalConditionForJet(UInt_t nIdx, Float_t &fJetPt, Float_t &fJetEta, Float_t &fJetPhi, Float_t &fJetMass) 
{
  Float_t fJetMassUncPP, fJetPtUncPP;
  Int_t nIdxGen;
  
  b_jetM.push_back(Jet_mass[ nIdx ]);
  b_jetPt.push_back(Jet_pt[ nIdx ]);
  b_jetEta.push_back(fJetEta);
  b_jetPhi.push_back(fJetPhi);
  
  b_jetPtNom.push_back(Jet_pt_jer_nom[ nIdx ]);
  
  if ( ( m_unFlag & OptFlag_JER_Up ) != 0 ) {
    fJetMassUncPP = Jet_mass_jer_up[ nIdx ];
    fJetPtUncPP = Jet_pt_jer_up[ nIdx ];
  } else if ( ( m_unFlag & OptFlag_JER_Dn ) != 0 ) {
    fJetMassUncPP = Jet_mass_jer_dn[ nIdx ];
    fJetPtUncPP = Jet_pt_jer_dn[ nIdx ];
  } else if ( ( m_unFlag & OptFlag_JES_Up ) != 0 ) {
    fJetMassUncPP = Jet_mass_jes_up[ nIdx ];
    fJetPtUncPP = Jet_pt_jes_up[ nIdx ];
  } else if ( ( m_unFlag & OptFlag_JES_Dn ) != 0 ) {
    fJetMassUncPP = Jet_mass_jes_dn[ nIdx ];
    fJetPtUncPP = Jet_pt_jes_dn[ nIdx ];
  } else {
    fJetMassUncPP = Jet_mass[ nIdx ];
    //fJetPtUncPP = Jet_pt[ nIdx ];
    fJetPtUncPP = Jet_pt_jer_nom[ nIdx ];
  }
  
  b_jetM_Unc_pp.push_back(fJetMassUncPP);
  b_jetPt_Unc_pp.push_back(fJetPtUncPP);
  b_jetM_Unc_an.push_back(fJetMass);
  b_jetPt_Unc_an.push_back(fJetPt);
  
  //nIdxGen = Jet_genJetIdx[ nIdx ];
  nIdxGen = GetMatchGenJet(nIdx);
  b_GenMatched.push_back(( nIdxGen >= 0 ? 1 : 0 ));
  
  if ( nIdxGen >= 0 && abs(( fJetPt - fJetPtUncPP ) / fJetPtUncPP) > 0.000001 ) {
    //printf("%lf - %lf; ", ( fJetPt - fJetPtUncPP ) / fJetPtUncPP, ( fJetMass - fJetMassUncPP ) / fJetMassUncPP);
    printf("%i, %i, %lf, %lf, %lf, %lf\n", (int)m_nIdxEntry, (int)nIdx, Jet_pt[ nIdx ], fJetPt, fJetPtUncPP, ( fJetPt - fJetPtUncPP ) / fJetPtUncPP);
  }
  
  return true;
}


void testForJECAnalyser::setOutput(std::string outputName) {
  m_output = TFile::Open(outputName.c_str(), "recreate");
  
  m_tree = new TTree("event", "event");
  
  m_tree->Branch("jetMass", &b_jetM);
  m_tree->Branch("jetPt", &b_jetPt);
  m_tree->Branch("jetEta", &b_jetEta);
  m_tree->Branch("jetPhi", &b_jetPhi);
  
  m_tree->Branch("jetPtNom", &b_jetPtNom);
  
  m_tree->Branch("jetMass_Unc_pp", &b_jetM_Unc_pp);
  m_tree->Branch("jetPt_Unc_pp", &b_jetPt_Unc_pp);
  m_tree->Branch("jetMass_Unc_an", &b_jetM_Unc_an);
  m_tree->Branch("jetPt_Unc_an", &b_jetPt_Unc_an);
  
  m_tree->Branch("GenMatched", &b_GenMatched);
}


void testForJECAnalyser::Loop() {
  Long64_t nentries, nbytes, nb;
  
  if (fChain == 0) return;
  
  nentries = fChain->GetEntries();
  nbytes = 0;

  // Events loop
  for (Long64_t iev=0; iev < nentries; iev++) {
    //GetEntry(iev);
    if ( LoadTree(iev) < 0 ) break;
    nb = fChain->GetEntry(iev);
    nbytes += nb;
    
    Reset();
    resetBranch();
    
    m_nIdxEntry = iev;
    
    jetSelection();
    //EventSelection();
    
    m_tree->Fill();
    
    //if ( iev % 1000 == 0 ) printf("Progress: %i / %i\n", (int)iev, (int)nentries);
  }
}


int main(int argc, char *argv[]) {
  // ORG: /xrootd/store/group/nanoAOD/run2_2016v5/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180610_143635/0000/nanoAOD_197.root
  //std::string strSrc = "org_Skim.root";
  std::string strSrc = "/cms/ldap_home/quark2930/Work/top/singletop/nanoAOD_949_Oct18_01/src/nano/nanoAOD/prod/nanoAOD.root";
  std::string strDst = "res.root";
  
  auto inFile = TFile::Open(strSrc.c_str());
  auto inTree = (TTree*) inFile->Get("Events");
  
  int nFlag = 0;
  
  if ( strcmp(argv[ 1 ], "JERUp") == 0 ) {
    printf("Test on JERUp\n");
    nFlag = topObjectSelection::OptFlag_JER_Up;
  } else if ( strcmp(argv[ 1 ], "JERDn") == 0 ) {
    printf("Test on JERDn\n");
    nFlag = topObjectSelection::OptFlag_JER_Dn;
  } else if ( strcmp(argv[ 1 ], "JESUp") == 0 ) {
    printf("Test on JESUp\n");
    nFlag = topObjectSelection::OptFlag_JES_Up;
  } else if ( strcmp(argv[ 1 ], "JESDn") == 0 ) {
    printf("Test on JESDn\n");
    nFlag = topObjectSelection::OptFlag_JES_Dn;
  }
  
  testForJECAnalyser ana(inTree, 0, 0, true, false, false, nFlag);
  ana.setOutput(strDst);
  ana.Loop();
  
  return 0;
}


