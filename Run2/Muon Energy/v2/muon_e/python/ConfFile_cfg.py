# Import CMS-specific Python classes and functions.
import FWCore.ParameterSet.Config as cms

# Create a CMSSW process object. "Demo" is the name of the process.
process = cms.Process("Demo")

# Control how the message logging is handled during the job execution.
# Here we load the MessageLogger module and change just one parameter.
process.load("FWCore.MessageService.MessageLogger_cfi")
# Make the Framework to report every 5 events instead of each event.
# https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_X/FWCore/MessageService/python/MessageLogger_cfi.py#L29
process.MessageLogger.cerr.FwkReport.reportEvery = 5

# Set the number of events to be processed in the CMSSW job.
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

# The first module (also an object by itself) to be attached to our process object.
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://eospublic.cern.ch//eos/opendata/cms/Run2015D/SingleElectron/MINIAOD/08Jun2016-v1/10000/001A703B-B52E-E611-BA13-0025905A60B6.root'
    )
)

#NOTE: There are muons in a "SingleElectron" dataset. This is, of course, expected.

process.demo = cms.EDAnalyzer('muon_e',
       muons = cms.InputTag("slimmedMuons")
       #muons = cms.InputTag("slimmedElectrons") # produces no output because there is no slimmedElectrons tag for muon objects.
)

process.hltHighLevel = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
    HLTPaths = cms.vstring('HLT_Ele27_WPLoose_Gsf_v*'),               # provide list of HLT paths (or patterns) you want
    eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
    andOr = cms.bool(True),              # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(True)        # throw exception on unknown path names
)


process.p = cms.Path(process.hltHighLevel+process.demo) #applies the filter. #the execution of the trigger path will stop if the hltHighLevel filter module throws a False result.
#process.p = cms.Path(process.demo+process.hltHighLevel) #does not apply the filter.
#process.p = cms.Path(process.demo) #does not apply the filter.
