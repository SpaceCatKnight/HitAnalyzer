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
# process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v2', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_v10', '')

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
    # fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIISpring16reHLT80/ZprimeToTTJet_M-4000_TuneCUETP8M1_13TeV-amcatnlo-pythia8/AODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/00E593DA-7139-E611-86E4-0CC47A4D99A4.root')
    fileNames = cms.untracked.vstring(
        # 'file:/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/CMSSW_9_2_1/src/btag_ntracks/TprimeBToTH_M-3000_Width-30p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8.root'
        # 'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/ZprimeToTT_M-2500_W-25_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/GEN-SIM-RECO/Asympt25nsReco_MCRUN2_74_V9-v1/10000/64461CF9-9176-E511-95F8-A0000420FE80.root'
        'file:TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_GEN-SIM-RECODEBUG_92X_upgrade2017_realistic_v10_ext1-v2.root'
    )
)
process.demo = cms.EDAnalyzer('HitAnalyzer',
    Verbosity = cms.untracked.bool(False),
    phase1 = cms.untracked.bool(True),
    src = cms.InputTag("siPixelClusters"),
)

process.p = cms.Path(process.demo)
