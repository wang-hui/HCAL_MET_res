#lowercase with underline for function/class names
#canmel case for variables

import ROOT as rt

BaseName = "CaloMET"
BaseFileList = ["UL_DYJetsToEE_M-50_RECO_PU_DLPHIN_plots.root"]
BaseHistList = ["CaloMET_h"]
#BaseHistList = ["LeadingCaloJet_nConstituents_h"]

Comp1Name = "CaloMETBE"
Comp1FileList = BaseFileList
#Comp1FileList = ["UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_DLPHIN_1dHB_2dHE_plots.root"]
Comp1HistList = BaseHistList
Comp1HistList = ["CaloMETBE_h"]

Comp2Name = "HBHE recHits MET"
Comp2FileList = BaseFileList
#Comp2FileList = ["UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_DLPHIN_1dHB_2dHE_truncate_plots.root"]
Comp2HistList = BaseHistList
Comp2HistList = ["myCaloMETBE_h"]

ShapeComp = True

YTitle = "Events"
XTitle = "MET [GeV]"
#XTitle = "nConstituents"

XMin = 0
XMax = 100
YScale = 2.0

FileDir = "results/"
HistDir = "myAna/"

if ShapeComp: YTitle = "A.U."

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
Comp2 = MyStruct(Comp2Name, Comp2FileList, Comp2HistList, rt.kBlue, StructList)

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

        OutName = BaseFileList[i].replace(".root", "_")
        BaseHist = None

        for k in range(len(StructList)):
            iFileName = StructList[k].FileList[i]
            iHistName = StructList[k].HistList[j]
            OutName = OutName + iHistName.replace("_h", "_")
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
                XMinTemp = iHist.GetXaxis().GetXmin()
                XMaxTemp = iHist.GetXaxis().GetXmax()
                if XMax > 0:
                    XMinTemp = XMin
                    XMaxTemp = XMax
                iHist.GetXaxis().SetRangeUser(XMinTemp, XMaxTemp)
                iHist.GetYaxis().SetTitle(YTitle)
                iHist.SetMaximum(iHist.GetMaximum() * YScale)
                iHist.SetTitle("")
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

                MyLine = rt.TLine(XMinTemp, 1.0, XMaxTemp, 1.0)
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

