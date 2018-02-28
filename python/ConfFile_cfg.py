import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
# process.load("Configuration.StandardSequences.Services_cff")
process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string('flatTuple.root')
                                   )



process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(True),
                     SkipEvent = cms.untracked.vstring('ProductNotFound'),
                     allowUnscheduled = cms.untracked.bool(True)
                     )
                     
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2017', '')


# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/CMSSW_9_2_1/src/btag_ntracks/TprimeBToTH_M-3000_Width-30p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8.root'
    )
)

process.demo = cms.EDAnalyzer('HitAnalyzer')
# save only events passing the full path SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ), 
# save PAT Layer 1 output; you need a '*' to 
# unpack the list of commands 'patEventContent' outputCommands = cms.untracked.vstring('drop *', 'keep *_patMuonAnalyzer_*_*') )

process.p = cms.Path(process.demo)
