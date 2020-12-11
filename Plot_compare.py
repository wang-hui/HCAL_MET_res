#lowercase with underline for function/class names
#canmel case for variables

import ROOT as rt

BaseName = "default recHit"
BaseFileList = ["2018_DYJetsToMuMu_M-50_RECO_origin_recHit_plots.root"]
BaseHistList = ["myCaloMETBE_h", "myCaloMETBE_HB_h", "myCaloMETBE_HE_h"]

Comp1Name = "DLPHIN no SF"
Comp1FileList = BaseFileList
Comp1FileList = ["2018_DYJetsToMuMu_M-50_RECO_DLPHIN_energy_plots.root"]
Comp1HistList = BaseHistList
#Comp1HistList = ["myCaloMETBE_HB_h"]

Comp2Name = "DLPHIN no SF"
Comp2FileList = BaseFileList
Comp2FileList = ["DoubleMuon_Run2018A_Run_315512_RECO_DLPHIN_no_SF_plots.root"]
Comp2HistList = BaseHistList
#Comp2HistList = ["myCaloMETBE_HE_h"]

ShapeComp = True

YTitle = "Events"
XTitle = "MET [GeV]"

FileDir = ""
HistDir = "myAna/"

class MyStruct:
    def __init__(self, Name, FileList, HistList, Color, StructList):
        self.Name = Name
        self.FileList = FileList
        self.HistList = HistList
        self.Color = Color
        StructList.append(self)

StructList = []
Base = MyStruct(BaseName, BaseFileList, BaseHistList, rt.kBlack, StructList)
Comp1 = MyStruct(Comp1Name, Comp1FileList, Comp1HistList, rt.kRed, StructList)
#Comp2 = MyStruct(Comp2Name, Comp2FileList, Comp2HistList, rt.kBlue, StructList)

rt.TH1.AddDirectory(rt.kFALSE)
#rt.TH1.__init__._creates = False

for i in range(len(BaseFileList)):
    for j in range(len(BaseHistList)):
        MyCanvas = rt.TCanvas("MyCanvas", "MyCanvas", 600, 600)
        rt.gStyle.SetOptStat(rt.kFALSE)

        MyLeg = rt.TLegend(0.5,0.65,0.9,0.9)

        PadUp = rt.TPad("PadUp", "PadUp", 0, 0.3, 1, 1.0)
        PadUp.SetBottomMargin(0.01)
        PadUp.Draw()

        PadDown = rt.TPad("PadDown", "PadDown", 0, 0, 1, 0.3)
        PadDown.SetTopMargin(0.03)
        PadDown.SetBottomMargin(0.3)
        PadDown.SetGrid()
        PadDown.Draw()

        OutName = ""
        BaseHist = None
        MaxY = 0

        for k in range(len(StructList)):
            iFileName = StructList[k].FileList[i]
            iHistName = StructList[k].HistList[j]
            OutName = OutName + iHistName + "_"
            iFile = rt.TFile.Open(FileDir + iFileName)
            iHist = iFile.Get(HistDir + iHistName)
            iHist.SetLineColor(StructList[k].Color)
            iHist.Sumw2()
            if ShapeComp: iHist.Scale(1/iHist.GetEntries())
            print k, iHist
            MyLeg.AddEntry(iHist, StructList[k].Name, "l")

            MyCanvas.cd()
            PadUp.cd()

            if k == 0:
                iHist.GetXaxis().SetRangeUser(0, 50)
                iHist.GetYaxis().SetTitle(YTitle)
                #iHist.SetTitle("")
                iHist.Draw("hist")
                BaseHist = iHist.Clone()

                MyCanvas.cd()
                PadDown.cd()
                BaseFrame = iHist.Clone()
                BaseFrame.Reset()
                BaseFrame.SetTitle("")

                BaseFrame.GetYaxis().SetTitle("Ratio")
                BaseFrame.GetYaxis().SetTitleOffset(0.4)
                BaseFrame.GetYaxis().SetTitleSize(0.1)
                BaseFrame.GetYaxis().SetLabelSize(0.08)
                BaseFrame.GetYaxis().SetRangeUser(0, 2)

                BaseFrame.GetXaxis().SetTitle(XTitle)
                BaseFrame.GetXaxis().SetTitleOffset(0.8)
                BaseFrame.GetXaxis().SetTitleSize(0.1)
                BaseFrame.GetXaxis().SetLabelSize(0.08)
                BaseFrame.Draw()

                MyLine = rt.TLine(0, 1.0, 50, 1.0)
                MyLine.Draw()
            else:
                iHist.Draw("histsame")

                MyCanvas.cd()
                PadDown.cd()
                RatioHist = iHist.Clone()
                rt.SetOwnership(RatioHist, False)
                RatioHist.Divide(BaseHist)
                RatioHist.Draw("same")

        MyCanvas.cd()
        PadUp.cd()
        MyLeg.Draw()
        MyCanvas.SaveAs("plots_temp/" + OutName + ".png")

