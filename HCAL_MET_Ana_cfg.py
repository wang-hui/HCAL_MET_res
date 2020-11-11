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
globalTag = '102X_dataRun2_v12'
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = globalTag

f = open(sys.argv[2], "r")
my_list = f.readlines()
f.close()

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
process.myAna = cms.EDAnalyzer("HCAL_MET_Ana", print_channel = cms.untracked.bool(False))

process.path = cms.Path(process.myAna)
