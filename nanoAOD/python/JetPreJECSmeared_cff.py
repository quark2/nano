import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *


jetJECTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
  #src = cms.InputTag("linkedObjects","jets"),
  src = cms.InputTag("slimmedJetsNewJEC"),
  cut = cms.string(""), #we should not filter on cross linked collections
  name = cms.string("JetPreJEC"),
  doc  = cms.string("slimmedJets, i.e. ak4 PFJets CHS with JECs ALREADY applied"),
  singleton = cms.bool(False), # the number of entries is variable
  extension = cms.bool(False), # this is the main table for the jets
  variables = cms.PSet(P4Vars), 
)

jetSmearTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
  #src = cms.InputTag("linkedObjects","jets"),
  src = cms.InputTag("slimmedJetsSmeared"),
  cut = cms.string(""), #we should not filter on cross linked collections
  name = cms.string("JetPreSmeared"),
  doc  = cms.string("slimmedJets, i.e. ak4 PFJets CHS with JECs ALREADY applied and also SMEARED"),
  singleton = cms.bool(False), # the number of entries is variable
  extension = cms.bool(False), # this is the main table for the jets
  variables = cms.PSet(P4Vars), 
)

jetJECTable.variables.pt.precision=20
jetSmearTable.variables.pt.precision=20

jetPreJECSmearTables = cms.Sequence(jetJECTable+jetSmearTable)


