import ROOT as rt
import DataFormats.FWLite as fw

FileList = "FileList/UL_DoublePion_E-50_RECO_pCaloHits_noPU_DLPHIN.list"
MaxFiles = -1
MaxEvents = -1

OutputFile = rt.TFile("HCAL_quick_test.root","RECREATE")

GenPionEta_h = rt.TH1F("GenPionEta_h", "GenPionEta_h", 70, -3.5, 3.5)
GenJetEta_h = rt.TH1F("GenJetEta_h", "GenJetEta_h", 70, -3.5, 3.5)
CaloJetEta_h = rt.TH1F("CaloJetEta_h", "CaloJetEta_h", 70, -3.5, 3.5)

GenPartHandle = fw.Handle("vector<reco::GenParticle>")
GenJetHandle = fw.Handle("vector<reco::GenJet>")
CaloJetHandle = fw.Handle("vector<reco::CaloJet>")

f = open(FileList, "r")
#InputFiles = f.readlines()
InputFiles = f.read().splitlines()
f.close()
MaxFilesTemp = len(InputFiles)
if MaxFiles > 0:
    MaxFilesTemp = MaxFiles

for InputFile in InputFiles[:MaxFilesTemp]:
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

        for GenJet in GenJets:
            GenJetEta_h.Fill(GenJet.eta())

        for CaloJet in CaloJets:
            CaloJetEta_h.Fill(CaloJet.eta())

OutputFile.cd()
OutputFile.Write()
OutputFile.Close()
