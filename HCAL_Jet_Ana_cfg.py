import FWCore.ParameterSet.Config as cms
import sys
from Configuration.Eras.Era_Run3_2023_cff import Run3_2023

process = cms.Process('Test',Run3_2023)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# geometry
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_postBPix_v2', '')

process.options = cms.untracked.PSet(
  numberOfThreads=cms.untracked.uint32(1)
)

f = open(sys.argv[1], "r")
my_list = f.readlines()
f.close()

if len(sys.argv) == 3:
    nFile = int(sys.argv[2])
    my_list = my_list[0:nFile]

OutputFile = sys.argv[1].split("/")[-1]
OutputFile = OutputFile.split(".")[0]
OutputFile = OutputFile + "_TTree.root"

# source
process.source = cms.Source(
    "PoolSource",
    fileNames  = cms.untracked.vstring(
	#"file:results_temp/reco_data_RAW2DIGI_RECO.root" 
	my_list
	),
    #firstEvent = cms.untracked.uint32(30)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.TFileService = cms.Service("TFileService", fileName = cms.string(OutputFile) )

#TightCut requires each event has two CaloJet match to GenJet
#Set to False if using single particle gun
process.myAna = cms.EDAnalyzer("HCAL_Jet_Ana", TightCut = cms.untracked.bool(False))

process.path = cms.Path(process.myAna)
