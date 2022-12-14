import FWCore.ParameterSet.Config as cms

myelectrons = cms.EDAnalyzer('Electron_Analyzer',
        electrons = cms.InputTag("slimmedElectrons"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        testBool = cms.bool(False)
)
