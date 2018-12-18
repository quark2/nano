#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "DataFormats/Math/interface/deltaR.h"

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
  Float_t Jet_pt_jec[35];
  
  Int_t nJetPreSmeared;
  Float_t Jet_pt_jer_nom[35];
  Float_t Jet_eta_jer_nom[35];
  Float_t Jet_phi_jer_nom[35];
  Float_t Jet_mass_jer_nom[35];
  
  Float_t Jet_mass_jer_up[35];
  Float_t Jet_mass_jer_dn[35];
  Float_t Jet_mass_jes_up[35];
  Float_t Jet_mass_jes_dn[35];
  
  Float_t Jet_pt_jer_up[35];
  Float_t Jet_pt_jer_dn[35];
  Float_t Jet_pt_jes_up[35];
  Float_t Jet_pt_jes_dn[35];
  
  TBranch *b_Jet_pt_jec_nom;
  
  TBranch *b_nJetPreSmeared;
  TBranch *b_Jet_pt_jer_nom;
  TBranch *b_Jet_eta_jer_nom;
  TBranch *b_Jet_phi_jer_nom;
  TBranch *b_Jet_mass_jer_nom;
  
  TBranch *b_Jet_mass_jer_up;
  TBranch *b_Jet_mass_jer_dn;
  TBranch *b_Jet_mass_jes_up;
  TBranch *b_Jet_mass_jes_dn;
  
  TBranch *b_Jet_pt_jer_up;
  TBranch *b_Jet_pt_jer_dn;
  TBranch *b_Jet_pt_jes_up;
  TBranch *b_Jet_pt_jes_dn;
  
  std::vector<Float_t> b_jetM, b_jetPt, b_jetEta, b_jetPhi;
  std::vector<Float_t> b_jetPtJEC;
  std::vector<Float_t> b_jetM_Unc_pp, b_jetPt_Unc_pp;
  std::vector<Float_t> b_jetM_Unc_an, b_jetPt_Unc_an;
  
  std::vector<Int_t> b_GenMatched;
  
  Long64_t m_nIdxEntry;

  void setOutput(std::string outputName);

  testForJECAnalyser(TTree *tree=0, TTree *had=0, TTree *hadTruth=0, Bool_t isMC = false, Bool_t sle = false, Bool_t slm = false, UInt_t unFlag = 0) :
    topEventSelectionSL(tree, had, hadTruth, isMC, sle, slm, unFlag)
  {
    // Loading the jets to which JEC is applied in PAT
    fChain->SetBranchAddress("JetPreJEC_pt", Jet_pt_jec, &b_Jet_pt_jec_nom);
    
    // Loading the smeared jets in PAT
    fChain->SetBranchAddress("nJetPreSmeared",     &nJetPreSmeared,  &b_nJetPreSmeared);
    fChain->SetBranchAddress("JetPreSmeared_pt",   Jet_pt_jer_nom,   &b_Jet_pt_jer_nom);
    fChain->SetBranchAddress("JetPreSmeared_eta",  Jet_eta_jer_nom,  &b_Jet_eta_jer_nom);
    fChain->SetBranchAddress("JetPreSmeared_phi",  Jet_phi_jer_nom,  &b_Jet_phi_jer_nom);
    fChain->SetBranchAddress("JetPreSmeared_mass", Jet_mass_jer_nom, &b_Jet_mass_jer_nom);
    
    // For uncertainty study; not related to nominal smearing
    fChain->SetBranchAddress("Jet_mass_jerUp", Jet_mass_jer_up, &b_Jet_mass_jer_up);
    fChain->SetBranchAddress("Jet_mass_jerDown", Jet_mass_jer_dn, &b_Jet_mass_jer_dn);
    fChain->SetBranchAddress("Jet_mass_jesTotalUp", Jet_mass_jes_up, &b_Jet_mass_jes_up);
    fChain->SetBranchAddress("Jet_mass_jesTotalDown", Jet_mass_jes_dn, &b_Jet_mass_jes_dn);
    
    // For uncertainty study; not related to nominal smearing
    fChain->SetBranchAddress("Jet_pt_jerUp", Jet_pt_jer_up, &b_Jet_pt_jer_up);
    fChain->SetBranchAddress("Jet_pt_jerDown", Jet_pt_jer_dn, &b_Jet_pt_jer_dn);
    fChain->SetBranchAddress("Jet_pt_jesTotalUp", Jet_pt_jes_up, &b_Jet_pt_jes_up);
    fChain->SetBranchAddress("Jet_pt_jesTotalDown", Jet_pt_jes_dn, &b_Jet_pt_jes_dn);
    
    // Easing cuts on jets
    cut_JetID = -1;
    cut_JetPt = 0;
    cut_JetEta = 104857600;
    cut_JetConeSizeOverlap = 0.0;
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
  
  b_jetPtJEC.clear();
  
  b_jetM_Unc_pp.clear();
  b_jetPt_Unc_pp.clear();
  b_jetM_Unc_an.clear();
  b_jetPt_Unc_an.clear();
  
  b_GenMatched.clear();
}


bool testForJECAnalyser::additionalConditionForJet(UInt_t nIdx, Float_t &fJetPt, Float_t &fJetEta, Float_t &fJetPhi, Float_t &fJetMass) 
{
  Float_t fJetMassUncPP, fJetPtUncPP;
  
  b_jetM.push_back(Jet_mass[ nIdx ]);
  b_jetPt.push_back(Jet_pt[ nIdx ]);
  b_jetEta.push_back(fJetEta);
  b_jetPhi.push_back(fJetPhi);
  
  b_jetPtJEC.push_back(Jet_pt_jec[ nIdx ]);
  
  // Seeking the matching smeared jet in PAT
  // Because it is sorted by pT after smearing, 
  // we don't know what is smeared one by scaling method and what is by stochastic method directly.
  // This deltaR method is the only way...
  Float_t fDRMin = 10485760;
  Int_t nIdxMatchedPre = -1;
  
  for ( Int_t i = 0 ; i < (Int_t)nJetPreSmeared ; i++ ) {
    Float_t fDR = deltaR(fJetEta, fJetPhi, Jet_eta_jer_nom[ i ], Jet_phi_jer_nom[ i ]);
    
    if ( fDRMin > fDR ) {
      fDRMin = fDR;
      nIdxMatchedPre = i;
    }
  }
  
  if ( nIdxMatchedPre < 0 ) {
    printf("Something is wrong on this jet...! - %i in %i-th event\n", nIdx, (int)m_nIdxEntry);
  }
  
  // For jet uncertainty study...
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
    // BUUUT THIS IS FOR SMEARER STUDY
    fJetMassUncPP = Jet_mass_jer_nom[ nIdxMatchedPre ];
    fJetPtUncPP = Jet_pt_jer_nom[ nIdxMatchedPre ];
  }
  
  b_jetM_Unc_pp.push_back(fJetMassUncPP);
  b_jetPt_Unc_pp.push_back(fJetPtUncPP);
  b_jetM_Unc_an.push_back(fJetMass);
  b_jetPt_Unc_an.push_back(fJetPt);
  
  // Checking whether it has matched gen jet (so that scaling method is applied) or not
  JME::JetParameters jetPars = {{JME::Binning::JetPt, Jet_pt[ nIdx ]},
                                {JME::Binning::JetEta, Jet_eta[ nIdx ]},
                                {JME::Binning::Rho, fixedGridRhoFastjetAll}};
  const double jetRes = jetResObj.getResolution(jetPars); // Note: this is relative resolution.
  b_GenMatched.push_back(( GetMatchGenJet(nIdx, jetRes * Jet_pt[ nIdx ]) >= 0 ? 1 : 0 ));
  
  return true;
}


void testForJECAnalyser::setOutput(std::string outputName) {
  m_output = TFile::Open(outputName.c_str(), "recreate");
  
  m_tree = new TTree("event", "event");
  
  m_tree->Branch("jetMass", &b_jetM);
  m_tree->Branch("jetPt", &b_jetPt);
  m_tree->Branch("jetEta", &b_jetEta);
  m_tree->Branch("jetPhi", &b_jetPhi);
  
  m_tree->Branch("jetPtJEC", &b_jetPtJEC);
  
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
  std::string strSrc = "nanoAOD.root";
  std::string strDst = "res.root";
  
  auto inFile = TFile::Open(strSrc.c_str());
  auto inTree = (TTree*) inFile->Get("Events");
  
  int nFlag = 0;
  
  // For jet uncertainty study
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


