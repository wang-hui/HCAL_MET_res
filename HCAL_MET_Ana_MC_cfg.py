import FWCore.ParameterSet.Config as cms
import sys
from Configuration.StandardSequences.Eras import eras
process = cms.Process("Test", eras.Run2_2018)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# geometry
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# global tag
globalTag = '106X_upgrade2018_realistic_v11_L1v1'
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = globalTag

process.options = cms.untracked.PSet(
  numberOfThreads=cms.untracked.uint32(1)
)

f = open(sys.argv[2], "r")
my_list = f.readlines()
f.close()

if len(sys.argv) == 4:
    nFile = int(sys.argv[3])
    my_list = my_list[0:nFile]

OutputFile = sys.argv[2].split("/")[-1]
OutputFile = OutputFile.split(".")[0]
OutputFile = OutputFile + "_plots.root"
#OutputFile = "test_plots.root"

# source
process.source = cms.Source(
    "PoolSource",
    fileNames  = cms.untracked.vstring(
	#"file:results_temp/reco_data_RAW2DIGI_RECO.root" 
	my_list
	),
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.TFileService = cms.Service("TFileService", fileName = cms.string(OutputFile) )
process.myAna = cms.EDAnalyzer("HCAL_MET_Ana", PrintChannel = cms.untracked.bool(False), IsMC = cms.untracked.bool(True), RunMod = cms.untracked.string("Zmumu"), IsHighPtZ = cms.untracked.bool(True), LowPtMatching = cms.untracked.bool(False))

process.path = cms.Path(process.myAna)
