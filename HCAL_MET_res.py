#!/usr/bin/env python

#case convention
#lowercase with underline for function names
#canmel case for variable/object names

import sys
import math
import array
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

#===================== define input/output here ===========================
OutputFile = ""      #if set to "", then OutputFile = <FileListName> + _plots.root
MaxEvents = 0        #set to 0 to run all events in each file
MaxFiles = -1        #set to -1 to run all files in file list
#MaxFiles = 10

#====================== define custome cuts here===========================
LepPtCut = 10
EtaCut = 2.4
MuIdCut = "mediumId" #options: "looseId", "mediumId", "tightId"
MiniIsoCut = 0.2

ResXEdges = array.array('d',[2.5 * (x*x+x) for x in range(10)])
ResXBins = len(ResXEdges) - 1

class ExampleAnalysis(Module):
    def __init__(self):
        self.writeHistFile=True

    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)

        #==================== 1D hist ============================
        self.ZmumuCand_mass_h = ROOT.TH1F('ZmumuCand_mass_h', 'ZmumuCand_mass_h', 100, 0, 200)
        self.addObject(self.ZmumuCand_mass_h)
        self.MET_h = ROOT.TH1F('MET_h', 'MET_h', 100, 0, 200)
        self.addObject(self.MET_h)
        self.MET_with_Z_h = ROOT.TH1F('MET_with_Z_h', 'MET_with_Z_h', 100, 0, 200)
        self.addObject(self.MET_with_Z_h)
        self.METphi_with_Z_h = ROOT.TH1F('METphi_with_Z_h', 'METphi_with_Z_h', 100, -3.2, 3.2)
        self.addObject(self.METphi_with_Z_h)
        self.METPara_h = ROOT.TH1F('METPara_h', 'METPara_h', 100, -200, 200)
        self.addObject(self.METPara_h)
        self.METVert_h = ROOT.TH1F('METVert_h', 'METVert_h', 100, -200, 200)
        self.addObject(self.METVert_h)

        self.CaloMET_h = ROOT.TH1F('CaloMET_h', 'CaloMET_h', 100, 0, 200)
        self.addObject(self.CaloMET_h)
        self.CaloMET_with_Z_h = ROOT.TH1F('CaloMET_with_Z_h', 'CaloMET_with_Z_h', 100, 0, 200)
        self.addObject(self.CaloMET_with_Z_h)
        self.CaloMETphi_with_Z_h = ROOT.TH1F('CaloMETphi_with_Z_h', 'CaloMETphi_with_Z_h', 100, -3.2, 3.2)
        self.addObject(self.CaloMETphi_with_Z_h)
        self.CaloMETPara_h = ROOT.TH1F('CaloMETPara_h', 'CaloMETPara_h', 100, -200, 200)
        self.addObject(self.CaloMETPara_h)
        self.CaloMETVert_h = ROOT.TH1F('CaloMETVert_h', 'CaloMETVert_h', 100, -200, 200)
        self.addObject(self.CaloMETVert_h)
        #===================== 2D hist =============================
        self.UPara_ratio_vs_Zpt_h = ROOT.TH2F('UPara_ratio_vs_Zpt_h', 'UPara_ratio_vs_Zpt_h', 20, 0, 200, 200, -10, 10)
        self.addObject(self.UPara_ratio_vs_Zpt_h)
        self.UPara_vs_Zpt_h = ROOT.TH2F('UPara_vs_Zpt_h', 'UPara_vs_Zpt_h', ResXBins, ResXEdges, 200, -300, 200)
        self.addObject(self.UPara_vs_Zpt_h)
        self.UVert_vs_Zpt_h = ROOT.TH2F('UVert_vs_Zpt_h', 'UVert_vs_Zpt_h', ResXBins, ResXEdges, 200, -150, 150)
        self.addObject(self.UVert_vs_Zpt_h)

        self.UCaloPara_ratio_vs_Zpt_h = ROOT.TH2F('UCaloPara_ratio_vs_Zpt_h', 'UCaloPara_ratio_vs_Zpt_h', 20, 0, 200, 200, -10, 10)
        self.addObject(self.UCaloPara_ratio_vs_Zpt_h)
        self.UCaloPara_vs_Zpt_h = ROOT.TH2F('UCaloPara_vs_Zpt_h', 'UCaloPara_vs_Zpt_h', ResXBins, ResXEdges, 200, -300, 200)
        self.addObject(self.UCaloPara_vs_Zpt_h)
        self.UCaloVert_vs_Zpt_h = ROOT.TH2F('UCaloVert_vs_Zpt_h', 'UCaloVert_vs_Zpt_h', ResXBins, ResXEdges, 200, -150, 150)
        self.addObject(self.UCaloVert_vs_Zpt_h)
        #===================== baseline test =======================
        self.BaseLineTest_h=ROOT.TH1F('BaseLineTest_h', '0: all. 1: pass filter', 5, 0, 5)
        self.addObject(self.BaseLineTest_h)

    def sel_muons(self, Muons):
        SelMuons = []
        for Muon in Muons:
            if Muon.pt > LepPtCut and abs(Muon.eta) < EtaCut and getattr(Muon, MuIdCut) and Muon.miniPFRelIso_all < MiniIsoCut:
                SelMuons.append(Muon)
        return SelMuons

    def sel_zmumu_cand(self, SelMuons):
        Mu1 = SelMuons[0]
        Mu2 = SelMuons[1]
        ZmumuCand = None
        if Mu1.pt > 20 and Mu1.charge + Mu2.charge == 0:        
            ZmumuCand = Mu1.p4() + Mu2.p4()
        return ZmumuCand

    def pass_event_filter(self, Flags):
        PassEventFilter = (
            Flags.goodVertices
            and Flags.HBHENoiseFilter
            and Flags.HBHENoiseIsoFilter
            and Flags.EcalDeadCellTriggerPrimitiveFilter
            and Flags.BadPFMuonFilter
            and Flags.ecalBadCalibFilter
            and Flags.globalSuperTightHalo2016Filter
            and Flags.eeBadScFilter)
        return PassEventFilter

    def analyze(self, event):
        #================ read needed collections/objects here ====================
        Muons = Collection(event, "Muon")
        Flags = Object(event, "Flag")
        MET = Object(event, "MET")
        CaloMET = Object(event, "CaloMET")

        #================ analyze each event ======================================
        self.BaseLineTest_h.Fill(0)
        self.MET_h.Fill(MET.pt)
        self.CaloMET_h.Fill(CaloMET.pt)

        #print event.event, CaloMET.pt

        if self.pass_event_filter(Flags):
            self.BaseLineTest_h.Fill(1)
        
            SelMuons = self.sel_muons(Muons)
            if len(SelMuons) == 2:
                ZmumuCand = self.sel_zmumu_cand(SelMuons)
                if ZmumuCand is not None:
                    ZmumuCandMass = ZmumuCand.M()
                    self.ZmumuCand_mass_h.Fill(ZmumuCandMass)
                    if ZmumuCandMass > 81 and ZmumuCandMass < 101:
                        ZmumuCandPt = ZmumuCand.Pt()
                        self.MET_with_Z_h.Fill(MET.pt)
                        self.METphi_with_Z_h.Fill(MET.phi)
                        dPhiZMET = MET.phi - ZmumuCand.Phi()
                        METPara = MET.pt * math.cos(dPhiZMET)
                        METVert = MET.pt * math.sin(dPhiZMET)
                        UPara = - METPara - ZmumuCandPt
                        UVert = - METVert

                        self.METPara_h.Fill(METPara)
                        self.METVert_h.Fill(METVert)
                        self.UPara_ratio_vs_Zpt_h.Fill(ZmumuCandPt, - UPara / ZmumuCandPt)
                        self.UPara_vs_Zpt_h.Fill(ZmumuCandPt, UPara)
                        self.UVert_vs_Zpt_h.Fill(ZmumuCandPt, UVert)

                        self.CaloMET_with_Z_h.Fill(CaloMET.pt)
                        self.CaloMETphi_with_Z_h.Fill(CaloMET.phi)
                        dPhiZCaloMET = CaloMET.phi - ZmumuCand.Phi()
                        CaloMETPara = CaloMET.pt * math.cos(dPhiZCaloMET)
                        CaloMETVert = CaloMET.pt * math.sin(dPhiZCaloMET)
                        UCaloPara = - CaloMETPara - ZmumuCandPt
                        UCaloVert = - CaloMETVert

                        self.CaloMETPara_h.Fill(CaloMETPara)
                        self.CaloMETVert_h.Fill(CaloMETVert)
                        self.UCaloPara_ratio_vs_Zpt_h.Fill(ZmumuCandPt, - UCaloPara / ZmumuCandPt)
                        self.UCaloPara_vs_Zpt_h.Fill(ZmumuCandPt, UCaloPara)
                        self.UCaloVert_vs_Zpt_h.Fill(ZmumuCandPt, UCaloVert)
        #return true to move on to the next module. return false to go to the next event
        return True

def read_file_list(FileList, MaxFiles):
    f=open(FileList, "r")
    InputFiles = f.readlines()
    f.close()

    nInputFiles = len(InputFiles)
    if MaxFiles > nInputFiles:
        print "MaxFiles", MaxFiles, "> nInputFiles", nInputFiles
        quit()

    nOutputFiles = MaxFiles
    if MaxFiles == -1: nOutputFiles = nInputFiles

    for i in range(nInputFiles): InputFiles[i] = InputFiles[i].strip()
    print InputFiles[0:nOutputFiles]
    return InputFiles[0:nOutputFiles]

FileList = sys.argv[1]
if OutputFile == "":
    OutputFile = FileList.split("/")[-1]
    OutputFile = OutputFile.split(".")[0]
    OutputFile = OutputFile + "_plots.root"
preselection=""
p=PostProcessor(".",read_file_list(FileList, MaxFiles),cut=preselection,branchsel=None,modules=[ExampleAnalysis()],noOut=True,histFileName=OutputFile,histDirName="plots",maxEntries=MaxEvents)
p.run()
