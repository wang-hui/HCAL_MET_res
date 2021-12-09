import sys
import ROOT as rt
import DataFormats.FWLite as fw

f = open(sys.argv[1], "r")
InputFileList = f.readlines()
f.close()

if len(sys.argv) == 3:
    nFile = int(sys.argv[2])
    InputFileList = InputFileList[0:nFile]
for i, InputFile in enumerate(InputFileList):
    InputFileList[i] = InputFile.strip()

OutputFileName = sys.argv[1].split("/")[-1]
OutputFileName = OutputFileName.split(".")[0]
OutputFileName = OutputFileName + "_quick_test.root"

MaxEvents = -1
OutputFile = rt.TFile(OutputFileName, "RECREATE")

GenPionEta_h = rt.TH1F("GenPionEta_h", "GenPionEta_h", 70, -3.5, 3.5)
GenJetEta_h = rt.TH1F("GenJetEta_h", "GenJetEta_h", 70, -3.5, 3.5)
CaloJetEta_h = rt.TH1F("CaloJetEta_h", "CaloJetEta_h", 70, -3.5, 3.5)

GenPionPt_h = rt.TH1F("GenPionPt_h", "GenPionPt_h", 100, 0, 1000)
GenJetPt_h = rt.TH1F("GenJetPt_h", "GenJetPt_h", 100, 0, 1000)
CaloJetPt_h = rt.TH1F("CaloJetPt_h", "CaloJetPt_h", 100, 0, 1000)

GenPartHandle = fw.Handle("vector<reco::GenParticle>")
GenJetHandle = fw.Handle("vector<reco::GenJet>")
CaloJetHandle = fw.Handle("vector<reco::CaloJet>")

for InputFile in InputFileList:
    print InputFile
    Events = fw.Events(InputFile, maxEvents=MaxEvents)
    for Event in Events:
        Event.getByLabel("genParticles", GenPartHandle)
        GenParts = GenPartHandle.product()

        Event.getByLabel("ak4GenJetsNoNu", GenJetHandle)
        GenJets = GenJetHandle.product()

        Event.getByLabel("ak4CaloJets", CaloJetHandle)
        CaloJets = CaloJetHandle.product()

        for GenPart in GenParts:
            #print GenPart.pdgId()
            if abs(GenPart.pdgId()) == 211:
                GenPionEta_h.Fill(GenPart.eta())
                GenPionPt_h.Fill(GenPart.pt())

        for GenJet in GenJets:
            GenJetEta_h.Fill(GenJet.eta())
            GenJetPt_h.Fill(GenJet.pt())

        for CaloJet in CaloJets:
            CaloJetEta_h.Fill(CaloJet.eta())
            CaloJetPt_h.Fill(CaloJet.pt())

OutputFile.cd()
OutputFile.Write()
OutputFile.Close()
