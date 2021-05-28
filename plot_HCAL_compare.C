int plot_HCAL_compare()
{
    bool plot_CaloJet_vs_GenJet = true;
    bool plot_CaloJet_vs_GenJet_pull = false;
    bool plot_MET_response = false;
    bool plot_MET_resolution = false;

/*
    std::vector<TString> hist_list =
    {
        //"LeadingCaloJet_vs_LeadingGenJet", "LeadingCaloJet_vs_LeadingGenJet", "LeadingCaloJet_vs_LeadingGenJet",
        //"LeadingCaloJet_vs_LeadingGenJet_HB", "LeadingCaloJet_vs_LeadingGenJet_HB", "LeadingCaloJet_vs_LeadingGenJet_HB",
        "LeadingCaloJet_vs_LeadingGenJet_HE", "LeadingCaloJet_vs_LeadingGenJet_HE", "LeadingCaloJet_vs_LeadingGenJet_HE",
        //"CaloJet_vs_GenJet", "CaloJet_vs_GenJet", "CaloJet_vs_GenJet",
        //"CaloJet_vs_GenJet_etaL", "CaloJet_vs_GenJet_etaL", "CaloJet_vs_GenJet_etaL",
        //"CaloJet_vs_GenJet_etaM", "CaloJet_vs_GenJet_etaM", "CaloJet_vs_GenJet_etaM",
        //"CaloJet_vs_GenJet_etaH", "CaloJet_vs_GenJet_etaH", "CaloJet_vs_GenJet_etaH",
    };
*/

    std::vector<TString> hist_list(4,
                                    //"LeadingCaloJet_vs_LeadingGenJet"
                                    //"LeadingCaloJet_vs_LeadingGenJet_HB"
                                    //"LeadingCaloJet_vs_LeadingGenJet_HE"
                                    //"DiJet_CaloJet_vs_GenJet"
                                    "DiJet_CaloJet_vs_GenJet_HB"
                                    //"DiJet_CaloJet_vs_GenJet_HE"
                                    );

    std::vector<TString> file_list = {"UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_origin_recHit_plots", "UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_M0_energy_plots"};
    file_list = {"UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_origin_recHit_plots", "UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_M0_energy_plots", "UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_DLPHIN_1dHB_2dHE_plots", "UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_DLPHIN_1dHB_2dHE_scaled_plots"};

    std::vector<TString> leg_list = {"MAHI", "M0"};
    leg_list = {"MAHI", "M0", "DLPHIN", "DLPHIN scaled"};

    std::vector<int> color_list = {kRed, kBlue};
    color_list = {kBlack, kRed, kBlue, kViolet};

    std::vector<TProfile*> SD_list, px_list;

    bool overlay_px = false;
    bool do_fit = true;
    bool scale_SD = true;

    TString hist_folder = "";
    hist_folder = "myAna";
    int rebin_x = 5;
    float px_scale = 0.5;
    float px_shift = 0;

    float SD_min = 0;
    float SD_max = 120;

    TString x_title = "gen pt [GeV]";
    TString y_title = "reco pt [GeV]";
    TString SD_title = "#sigma_{reco pt} x gen pt / reco pt";

    if (plot_CaloJet_vs_GenJet_pull)
    {
        hist_list =
        {
            "CaloJet_vs_GenJet_pull", "CaloJet_vs_GenJet_pull", //"CaloJet_vs_GenJet_pull",
            //"CaloJet_vs_GenJet_etaL_pull", "CaloJet_vs_GenJet_etaL_pull", "CaloJet_vs_GenJet_etaL_pull",
            //"CaloJet_vs_GenJet_etaM_pull", "CaloJet_vs_GenJet_etaM_pull", "CaloJet_vs_GenJet_etaM_pull",
            //"CaloJet_vs_GenJet_etaH_pull", "CaloJet_vs_GenJet_etaH_pull", "CaloJet_vs_GenJet_etaH_pull",
        };
        y_title = "pull";
        SD_title = "pull";
        overlay_px = true;
        do_fit = false; 
    }

    if (plot_MET_response)
    {
        file_list = {"UL_EGamma_Run2018B_Run_317392_RECO_origin_recHit_plots", "UL_EGamma_Run2018B_Run_317392_RECO_DLPHIN_1dHB_2dHE_plots", "UL_EGamma_Run2018B_Run_317392_RECO_DLPHIN_1dHB_2dHE_respCorr_plots"};
        leg_list = {"MAHI", "DLPHIN no respCorr", "DLPHIN with respCorr"};
        hist_list =
        {
            "UPara_ratio_vs_Zpt", "UPara_ratio_vs_Zpt", "UPara_ratio_vs_Zpt",
        };
        px_scale = 1;

        x_title = "Z pt [GeV]";
        y_title = "- Upara / Z pt";
        overlay_px = true;
        do_fit = false; 
    }

    if (plot_MET_resolution)
    {
        file_list = {"UL_EGamma_Run2018B_Run_317392_RECO_origin_recHit_plots", "UL_EGamma_Run2018B_Run_317392_RECO_DLPHIN_1dHB_2dHE_plots", "UL_EGamma_Run2018B_Run_317392_RECO_DLPHIN_1dHB_2dHE_respCorr_plots"};
        leg_list = {"MAHI", "DLPHIN no respCorr", "DLPHIN with respCorr"};
        hist_list =
        {
            "UPara_vs_Zpt", "UPara_vs_Zpt", "UPara_vs_Zpt",
            //"UVert_vs_Zpt", "UVert_vs_Zpt", "UVert_vs_Zpt",
        };
        px_scale = 1;

        x_title = "Z pt [GeV]";
        y_title = "Upara [GeV]";
        SD_title = "#sigma_{Upara}";
        //y_title = "Uvert [GeV]";
        //SD_title = "#sigma_{Uvert}";
        do_fit = false;
        scale_SD = false;
        SD_min = 0;
        SD_max = 100;
        rebin_x = 1;
        px_shift = 5;
    }

    for(int i = 0; i < hist_list.size(); i++)
    {
        TString file_name = file_list.at(i);
        TString hist_name = hist_list.at(i);

        TFile *f1 = new TFile("results/" + file_name + ".root");
        TH2F *h1 = (TH2F*)f1->Get(hist_folder + "/" + hist_name + "_h");

        TCanvas* mycanvas = new TCanvas("mycanvas", "mycanvas", 600, 600);
        gStyle->SetOptStat(kFALSE);

        auto xmin = h1->GetXaxis()->GetXmin();
        auto xmax = h1->GetXaxis()->GetXmax();
        auto ymin = h1->GetYaxis()->GetXmin();
        auto ymax = h1->GetYaxis()->GetXmax();

        h1->Draw("colz");
        h1->SetTitle(hist_name);
        //h1->SetTitle("");
        h1->RebinX(rebin_x);
        h1->GetXaxis()->SetTitle(x_title);
        h1->GetXaxis()->SetRangeUser(xmin, xmax);
        h1->GetYaxis()->SetTitle(y_title);
        h1->GetYaxis()->SetRangeUser(ymin, ymax);
        gPad->SetLogz();

        mycanvas->SetLeftMargin(0.15);
        mycanvas->SetRightMargin(0.15);
        mycanvas->SaveAs("plots_temp/" + file_name + "_" + hist_name + ".png");
        //mycanvas->SaveAs("plots_temp/" + file_name + "_"  + hist_name + ".pdf");

        TProfile *px = h1->ProfileX();
        if(!overlay_px) px->BuildOptions(0, 0, "s");
        px->SetTitle(hist_name);
        //px->SetTitle("");
        //px->GetYaxis()->SetNdivisions(512);
        px->GetXaxis()->SetRangeUser(xmin + px_shift, xmax * px_scale);
        //px->GetYaxis()->SetRangeUser(0.5, 2);
        px->SetLineColor(kRed);
        //px->SetLineWidth(2);
        //px->SetMarkerStyle(8);
        px->GetXaxis()->SetTitle(x_title);
        px->GetYaxis()->SetTitle(y_title);
        px->Draw();

        px_list.push_back(px);

        gPad->SetGrid();

        if(do_fit)
        {
            px->Fit("pol1", "w");
            TF1 *f = (TF1*)px->GetListOfFunctions()->FindObject("pol1");
            if(f)
            {
                f->SetLineColor(kBlue);
                f->SetLineWidth(1);
                f->SetLineStyle(2); // 2 = "- - -"
            }

            TLine *l=new TLine(xmin + px_shift, ymin + px_shift, xmax * px_scale, ymax * px_scale);
            l->SetLineColor(kBlack);
            if(!hist_name.Contains("err"))
            {l->Draw("same");}
        }
        mycanvas->SaveAs("plots_temp/" + file_name + "_"  + hist_name + "_profile.png");
        //mycanvas->SaveAs("plots_temp/" + file_name + "_"  + hist_name + "_profile.pdf");

        TProfile* SD_h = (TProfile*)px->Clone(hist_name + "_SD");
        SD_h->Reset();
        for(int i = 1; i<= px->GetNbinsX(); i++)
        {
            auto x_center = px->GetBinCenter(i);
            auto error = px->GetBinError(i);
            auto y_center = px->GetBinContent(i);
            auto SD = error;
            if (scale_SD && y_center > 0) SD = error * x_center / y_center;
            SD_h->SetBinContent(i, SD);
            SD_h->SetBinEntries(i, 1);
        }
        //SD_h->SetTitle("");
        SD_h->GetXaxis()->SetTitle(x_title);
        SD_h->GetYaxis()->SetTitle(SD_title);
        SD_h->GetXaxis()->SetRangeUser(xmin + px_shift, xmax * px_scale);
        SD_h->GetYaxis()->SetRangeUser(SD_min, SD_max);
        SD_h->Draw("hist");
        SD_list.push_back(SD_h);
        mycanvas->SaveAs("plots_temp/" + file_name + "_"  + hist_name + "_SD.png");
        //mycanvas->SaveAs("plots_temp/" + file_name + "_"  + hist_name + "_SD.pdf");
    }

    TCanvas* mycanvas = new TCanvas("mycanvas", "mycanvas", 600, 600);
    gStyle->SetOptStat(kFALSE);
    mycanvas->SetLeftMargin(0.15);
    //mycanvas->SetRightMargin(0.1);

    TLegend* leg = new TLegend(0.65,0.7,0.9,0.9);

    for(int i = 0; i < SD_list.size(); i++)
    {
        //std::cout << SD_list.at(i) << std::endl;
        auto my_hist = SD_list.at(i);
        if(overlay_px) my_hist = px_list.at(i);
        my_hist->SetLineColor(color_list.at(i));
        if(i==0) my_hist->Draw("hist");
        else my_hist->Draw("histsame");
        leg->AddEntry(my_hist,leg_list.at(i),"l");
    }
    leg->Draw("same");
    gPad->SetGrid();

    TString postfix = "_SD";
    if(overlay_px) postfix = "_px";

    mycanvas->SaveAs("plots_temp/" + hist_list.at(0) + postfix + ".png");
    //mycanvas->SaveAs("plots_temp/" + hist_list.at(0) + postfix + ".pdf");

    return 0;
}
