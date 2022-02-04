#lowercase with underline for function/class names
#canmel case for variables

import math
import ROOT as rt
#rt.TH1.SetDefaultSumw2()
rt.TH1.AddDirectory(False)
rt.gStyle.SetOptStat(False)
rt.gROOT.SetBatch(True)

import CMS_lumi
#import tdrstyle
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "(13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPeriod = 0
iPos = 11

PlotsFolder = "plots/"

SubdetRangeList = ["0 < |#eta_{track}| < 1.2", "1.5 < |#eta_{track}| < 2.3", "0 < |#eta_{track}| < 2.3"]
SubdetMaxList = [300, 800, 1000]
SubdetHistDict = {
"HB" : ["raw1etaB11One", "raw2etaB11One", "raw0etaB11One"],
"HE" : ["raw1etaB13One", "raw2etaB13One", "raw0etaB13One"],
#"HB+HE" : ["raw1etaB14One", "raw2etaB14One", "raw0etaB14One"],
}
HistColorList = [rt.kBlue-4, rt.kBlack, rt.kRed]
HistLegList = ["M0", "M2", "MAHI"]
HistMarkerList = [20, 21, 22]

File = rt.TFile.Open("IsoTrackPlots.root")

for iSubdet, (Subdet, HistNameList) in enumerate(SubdetHistDict.items()):
    MyCanvas = rt.TCanvas("MyCanvas", "MyCanvas", 600, 600)
    MyCanvas.cd()
    MyCanvas.SetTopMargin(0.08)
    MyCanvas.SetRightMargin(0.1)
    MyCanvas.SetLeftMargin(0.15)

    MyLeg = rt.TLegend(0.7,0.7,0.89,0.91)
    MyLeg.SetBorderSize(0)

    HistStrList = []

    for iHist, HistName in enumerate(HistNameList):
        #print iHist, HistName
        Hist = File.Get(HistName)

        Hist.SetMarkerStyle(HistMarkerList[iHist])
        Hist.SetMarkerColor(HistColorList[iHist])
        Hist.SetLineColor(HistColorList[iHist])

        MyFit = Hist.GetListOfFunctions().FindObject("gaus")
        MyFit.SetLineColor(HistColorList[iHist])
        MyFit.SetLineWidth(2)
        #MyFit.SetLineStyle(2) # 2 = "- - -"
        C = MyFit.GetParameter(0)
        Mu = MyFit.GetParameter(1)
        Sigma = MyFit.GetParameter(2)
        print HistLegList[iHist], Mu, "+-", Sigma
        Str = "#bf{" + HistLegList[iHist] + " #sigma/#mu = %.3f}" % (Sigma/Mu)
        HistStrList.append(Str)

        MyLeg.AddEntry(Hist, HistLegList[iHist], "pe")
        if iHist == 0:
            Hist.SetTitle("")
            Hist.GetYaxis().SetTitleOffset(1.3)
            Hist.SetMaximum(SubdetMaxList[iSubdet])
            Hist.Draw("pex0")
        else:
            Hist.Draw("pex0same")

    MyLeg.Draw("same")
    MyLatex = rt.TLatex()
    MyLatex.SetTextSize(0.04)
    MyLatex.SetNDC()
    MyLatex.DrawLatex(0.55,0.65,"#bf{20 < p_{track} < 30 GeV}")
    MyLatex.DrawLatex(0.55,0.6,"#bf{" + SubdetRangeList[iSubdet]+ "}")
    for iStr, Str in enumerate(HistStrList):
        MyLatex.SetTextColor(HistColorList[iStr])
        MyLatex.DrawLatex(0.6,0.52-0.04*iStr,Str)
    #MyLatex.SetTextColor(rt.kBlack)
    #MyLatex.DrawLatex(0.15,0.91,"CMS #bf{Preliminary}")
    #MyLatex.SetTextAlign(31)
    #MyLatex.DrawLatex(0.9,0.91,"#bf{13TeV}")
    CMS_lumi.CMS_lumi(MyCanvas, iPeriod, iPos)

    MyCanvas.SaveAs(PlotsFolder + "IsoTrackRes_20to30GeV_" + Subdet + ".png")
    MyCanvas.SaveAs(PlotsFolder + "IsoTrackRes_20to30GeV_" + Subdet + ".pdf")
        
