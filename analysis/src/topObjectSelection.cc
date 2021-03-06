#include "nano/analysis/interface/topObjectSelection.h"

using std::vector;

topObjectSelection::topObjectSelection(TTree *tree, TTree *had, TTree *hadTruth, Bool_t isMC, Bool_t _isDilep, Bool_t _isSemiLep) :
  nanoBase(tree, had, hadTruth, isMC),
  isDilep(_isDilep),
  isSemiLep(_isSemiLep)
{}

vector<TParticle> topObjectSelection::elecSelection() {
  vector<TParticle> elecs; 
  for (UInt_t i = 0; i < nElectron; ++i){
    if (Electron_pt[i] < 20) continue;
    if (isSemiLep) { if (Electron_pt[i] < 35) continue; }
    if (std::abs(Electron_eta[i]) > 2.1) continue;
    if (Electron_cutBased[i] < 4) continue; 
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( std::abs(el_scEta) > 1.4442 &&  std::abs(el_scEta) < 1.566 ) continue;
    //if ( std::abs(el_scEta) >= 1.566 ) continue; // For AN-2017/083; it must be turned off when QCD-studying
    if ( Electron_pfRelIso03_all[ i ] > 0.0588 ) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);

    auto elec = TParticle();
    elec.SetPdgCode(11*Electron_charge[i]*-1);
    elec.SetMomentum(mom);
    elec.SetWeight(el_scEta);
    elec.SetFirstMother(i);
    elecs.push_back(elec);
    b_isolep = Electron_pfRelIso03_all[ i ];
  }
  return elecs;
}

vector<TParticle> topObjectSelection::muonSelection() {
  vector<TParticle> muons; 
  for (UInt_t i = 0; i < nMuon; ++i){
    if (!Muon_tightId[i]) continue;
    if (Muon_pt[i] < 20) continue;
    if (isSemiLep) { if (Muon_pt[i] < 26) continue; }
    if (std::abs(Muon_eta[i]) > 2.4) continue;
    //if (isSemiLep) { if (std::abs(Muon_eta[i]) > 2.1) continue; }
    if (Muon_pfRelIso04_all[i] > 0.06) continue;
    if (!Muon_globalMu[i]) continue;
    if (!Muon_isPFcand[i]) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);
    muon.SetFirstMother(i);
    muons.push_back(muon);
    b_isolep = Muon_pfRelIso04_all[ i ];
  }
  return muons;
}

/// TODO: Implement veto electron selection (currently empty on TTbarXSecSynchronization)
vector<TParticle> topObjectSelection::vetoElecSelection() {
  vector<TParticle> elecs; 
  for (UInt_t i = 0; i < nElectron; ++i){
    if (Electron_pt[i] < 15) continue;
    if (std::abs(Electron_eta[i]) > 2.5) continue;
    if (Electron_cutBased[i] < 1) continue; 
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( std::abs(el_scEta) > 1.4442 &&  std::abs(el_scEta) < 1.566 ) continue;
    if ( std::abs(el_scEta) >= 1.566 ) continue; // For AN-2017/083; it must be turned off when QCD-studying
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);
    auto elec = TParticle();
    elec.SetPdgCode(11*Electron_charge[i]*-1);
    elec.SetMomentum(mom);
    elec.SetWeight(el_scEta);
    elec.SetFirstMother(i);
    elecs.push_back(elec);
  }
  return elecs;
}

vector<TParticle> topObjectSelection::vetoMuonSelection() {
  vector<TParticle> muons; 
  for (UInt_t i = 0; i < nMuon; ++i){
    //if (!Muon_looseId[i]) continue;

    // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Loose_Muon
    if (!Muon_isPFcand[i]) continue;
    if (!(Muon_globalMu[i] || Muon_trackerMu[i])) continue;

    
    if (Muon_pt[i] < 10) continue;
    if (std::abs(Muon_eta[i]) > 2.4) continue;
    if (Muon_pfRelIso04_all[i] > 0.2) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);
    muon.SetFirstMother(i);
    muons.push_back(muon);
  }
  return muons;
}

vector<TParticle> topObjectSelection::genJetSelection() {
  vector<TParticle> jets;
  for (UInt_t i = 0; i < nJet; ++i){
    if (GenJet_pt[i] < 30) continue;
    if (std::abs(GenJet_eta[i]) > 2.4) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(GenJet_pt[i], GenJet_eta[i], GenJet_phi[i], GenJet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps){
        if (mom.TLorentzVector::DeltaR(lep) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    auto jet = TParticle();
    jet.SetMomentum(mom);
    jet.SetFirstMother(i);
    jets.push_back(jet);
  }
  return jets;
}

vector<TParticle> topObjectSelection::jetSelection(std::vector<Float_t> *csvVal) {
  vector<TParticle> jets;
  //float Jet_SF_CSV[19] = {1.0,};
  for (UInt_t i = 0; i < nJet; ++i){
    // For AN-2017/083
    if ( std::abs(Jet_eta[i]) > 4.7 ) continue;
    if ( !( 2.7 <= std::abs(Jet_eta[i]) && std::abs(Jet_eta[i]) < 3.0 ) ) {
      if (Jet_pt[i] < 40) continue;
    } else {
      if (Jet_pt[i] < 50) continue;
    }
    //if (Jet_pt[i] < 40) continue;
    //if (std::abs(Jet_eta[i]) > 4.7) continue;
    if (Jet_jetId[i] < 1) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps){
        if (mom.TLorentzVector::DeltaR(lep) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    auto jet = TParticle();
    jet.SetMomentum(mom);
    jet.SetFirstMother(i);
    jets.push_back(jet);
    //b_btagCSVV2 = Jet_btagCSVV2[i];
    //BTagEntry::JetFlavor JF = BTagEntry::FLAV_UDSG;
    BTagEntry::JetFlavor JF;
    if (abs(Jet_hadronFlavour[i]) == 5) JF = BTagEntry::FLAV_B;
    //else if (abs(Jet_hadronFlavour[i]) == 4) JF = BTagEntry::FLAV_C;
    auto bjetSF = m_btagSF.eval_auto_bounds("central", JF , Jet_eta[i], Jet_pt[i], Jet_btagCSVV2[i]);
    b_btagweight *= bjetSF;
    if ( csvVal != NULL ) csvVal->push_back(Jet_btagCSVV2[ i ]);
  }
  
  return jets;
}

vector<TParticle> topObjectSelection::bjetSelection() {
  vector<TParticle> bjets;
  for (UInt_t i = 0; i < nJet; ++i ) {
    // For AN-2017/083
    /*if ( !( 2.7 <= std::abs(Jet_eta[i]) && std::abs(Jet_eta[i]) < 3.0 ) ) {
      if (Jet_pt[i] < 40) continue;
    } else {
      if (Jet_pt[i] < 50) continue;
    }*/
    // According to p. 9, AN-2017/083, one can find the following eta cut
    if (std::abs(Jet_eta[i]) > 2.4) continue;
    if (Jet_pt[i] < 40) continue;
    if (Jet_jetId[i] < 1) continue;
    //if (Jet_btagCSVV2[i] < 0.8484) continue;
    if (Jet_btagCSVV2[i] < 0.9535) {
      if ( b_maxBDiscr_nonb < Jet_btagCSVV2[ i ] ) {
        b_maxBDiscr_nonb = Jet_btagCSVV2[ i ];
      }
      
      continue;
    }
    // For AN-2017/083
    // See p. 6, AN-2017/056, or 
    // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Boosted_event_topologies
    //if (Jet_btagCMVA[i] < 0.9432) {continue;}
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps) {
        if (mom.TLorentzVector::DeltaR(lep) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    auto bjet = TParticle();
    bjet.SetMomentum(mom);
    bjet.SetFirstMother(i);
    bjets.push_back(bjet);
  }
  return bjets;
}
