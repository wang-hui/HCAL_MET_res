int plot_HCAL_compare_paper()
{
    bool plot_MET_response = true;
    bool plot_MET_resolution = false;

    std::vector<TString> hist_list, file_list, leg_list;
    std::vector<int> color_list = {kRed, kBlue, kGreen+1, kYellow+1};
    std::vector<TProfile*> SD_list, px_list;
    std::vector<float> response_list;

    bool do_SD = false;
    bool scale_SD = false;
    bool print_test_value = true;
    float test_value = 50;

    TString hist_folder = "";
    hist_folder = "myAna";
    int rebin_x = 5;
    float px_min = 0;
    float px_max = 0;

    float y_min = 0;
    float y_max = 120;

    TString x_title = "gen pt [GeV]";
    TString y_title = "reco pt [GeV]";
    TString SD_title = "#sigma_{reco pt} x gen pt / reco pt";

    if (plot_MET_response)
    {
        file_list =
        {
            //"DYJetsToMuMu_M-50_Zpt-150toInf_RECO_mahi_energy_plots",
            //"DYJetsToMuMu_M-50_Zpt-150toInf_RECO_M0_energy_plots",
            "DoubleMuon_Run2018B_Run_317392_RECO_mahi_energy_plots",
            "DoubleMuon_Run2018B_Run_317392_RECO_M0_energy_plots",
            //"DoubleMuon_Run2018B_Run_317392_RECO_DLPHIN_respCorr_ZeroOut_plots",
        };
        leg_list = {"MAHI", "M0"};
        hist_list.assign(2,
                "CaloMETBE_UPara_ratio_vs_Zpt"
                //"CaloMETBE_UPara_ratio_vs_Zpt_PUL"
                //"CaloMETBE_UPara_ratio_vs_Zpt_PUH"
                );
        x_title = "Z pt [GeV]";
        //y_title = "- Upara / Z pt";
        y_title = "response";
        do_SD = false;
        y_min = 0.0;
        y_max = 1.0;
        rebin_x = 1;
        px_min = 0;
        px_max = -40;
    }

    if (plot_MET_resolution)
    {
        file_list =
        {
            //"DYJetsToMuMu_M-50_Zpt-150toInf_RECO_mahi_energy_plots",
            //"DYJetsToMuMu_M-50_Zpt-150toInf_RECO_M0_energy_plots",
            "DoubleMuon_Run2018B_Run_317392_RECO_mahi_energy_plots",
            "DoubleMuon_Run2018B_Run_317392_RECO_M0_energy_plots",
            //"DoubleMuon_Run2018B_Run_317392_RECO_DLPHIN_respCorr_ZeroOut_plots",
        };
        leg_list = {"MAHI", "M0"};
        hist_list.assign(2,
                //"CaloMETBE_UPara_vs_Zpt"
                //"CaloMETBE_UPara_vs_Zpt_PUL"
                //"CaloMETBE_UPara_vs_Zpt_PUH"
                "CaloMETBE_UVert_vs_Zpt"
                //"CaloMETBE_UVert_vs_Zpt_PUL"
                //"CaloMETBE_UVert_vs_Zpt_PUH"
                );
        //response_list = {0.736737, 0.764297};         //all PU
        //response_list = {0.732903, 0.759424};         //low PU
        //response_list = {0.739548, 0.770508};         //high PU
        response_list = {0.482794, 0.516279};           //data

        x_title = "Z pt [GeV]";
        //y_title = "Upara [GeV]";
        //SD_title = "#sigma_{Upara} / response";
        //y_max = 80;
        y_title = "Uvert [GeV]";
        SD_title = "#sigma_{Uvert} / response";
        y_max = 60;
        do_SD = true;
        scale_SD = true;
        y_min = 20;
        rebin_x = 1;
        px_min = 0;
        px_max = -90;
    }

    for(int i_hist = 0; i_hist < hist_list.size(); i_hist++)
    {
        TString file_name = file_list.at(i_hist);
        TString hist_name = hist_list.at(i_hist);

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
        h1->SetTitle("");
        h1->RebinX(rebin_x);
        h1->GetXaxis()->SetTitle(x_title);
        h1->GetXaxis()->SetRangeUser(xmin + px_min, xmax + px_max);
        h1->GetYaxis()->SetTitle(y_title);
        h1->GetYaxis()->SetRangeUser(ymin, ymax);
        gPad->SetLogz();

        mycanvas->SetLeftMargin(0.15);
        mycanvas->SetRightMargin(0.15);
        //mycanvas->SaveAs("plots_temp/" + file_name + "_" + hist_name + ".png");
        //mycanvas->SaveAs("plots_temp/" + file_name + "_"  + hist_name + ".pdf");

        TProfile *px = h1->ProfileX();
        if(do_SD) px->BuildOptions(0, 0, "s");

        px->SetLineColor(kRed);
        px->GetYaxis()->SetTitle(y_title);
        //px->SetLineWidth(2);
        //px->SetMarkerStyle(8);
        px->Draw("same");
        if(print_test_value)
        {
            std::cout << file_name << std::endl;
            std::cout << "profile at x = " << test_value << " is " << px->GetBinContent(px->GetXaxis()->FindBin(test_value)) << std::endl;
        }
        px_list.push_back(px);

        gPad->SetGrid();

        mycanvas->SaveAs("plots_temp/" + file_name + "_"  + hist_name + "_profile.png");
        //mycanvas->SaveAs("plots_temp/" + file_name + "_"  + hist_name + "_profile.pdf");

        if(do_SD)
        {
            auto response = response_list.at(i_hist);
            TProfile* SD_h = (TProfile*)px->Clone(hist_name + "_SD");
            SD_h->Reset();
            for(int i_bin = 1; i_bin<= px->GetNbinsX(); i_bin++)
            {
                auto x_center = px->GetBinCenter(i_bin);
                auto error = px->GetBinError(i_bin);
                auto y_center = px->GetBinContent(i_bin);
                auto SD = error;
                if (scale_SD) SD = error / response;
                SD_h->SetBinContent(i_bin, SD);
                SD_h->SetBinEntries(i_bin, 1);
            }
            //SD_h->SetTitle("");
            SD_h->GetXaxis()->SetTitle(x_title);
            SD_h->GetYaxis()->SetTitle(SD_title);
            SD_h->GetXaxis()->SetRangeUser(xmin + px_min, xmax + px_max);
            SD_h->Draw("hist");
            if(print_test_value)
            {
                std::cout << "SD at x = " << test_value << " is " << SD_h->GetBinContent(SD_h->GetXaxis()->FindBin(test_value)) << std::endl;
            }

            SD_list.push_back(SD_h);
            mycanvas->SaveAs("plots_temp/" + file_name + "_"  + hist_name + "_SD.png");
            //mycanvas->SaveAs("plots_temp/" + file_name + "_"  + hist_name + "_SD.pdf");
        }
    }

    TCanvas* mycanvas = new TCanvas("mycanvas", "mycanvas", 600, 600);
    gStyle->SetOptStat(kFALSE);
    mycanvas->SetLeftMargin(0.15);
    //mycanvas->SetRightMargin(0.1);

    TLegend* leg = new TLegend(0.65,0.7,0.9,0.9);

    for(int i_px = 0; i_px < px_list.size(); i_px++)
    {
        //std::cout << SD_list.at(i_px) << std::endl;
        auto my_hist = px_list.at(i_px);
        if(do_SD) my_hist = SD_list.at(i_px);
        my_hist->SetLineColor(color_list.at(i_px));
        my_hist->GetYaxis()->SetRangeUser(y_min, y_max);
        if(do_SD)
        {
            if(i_px==0) my_hist->Draw("hist");
            else my_hist->Draw("histsame");
        }
        else
        {
            if(i_px==0) my_hist->Draw("");
            else my_hist->Draw("same");
        }
        leg->AddEntry(my_hist,leg_list.at(i_px),"l");
    }
    leg->Draw("same");
    gPad->SetGrid();

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    //latex.DrawLatex(0.15,0.91,"CMS #bf{Simulation 2018}");
    latex.DrawLatex(0.15,0.91,"CMS #bf{Preliminary 2018}");
    latex.SetTextAlign(31);  //align at right bottom
    latex.DrawLatex(0.9,0.91,"#bf{13TeV}");

    TString postfix = "_px";
    if(do_SD) postfix = "_SD";

    mycanvas->SaveAs("plots_temp/" + hist_list.at(0) + postfix + ".png");
    mycanvas->SaveAs("plots_temp/" + hist_list.at(0) + postfix + ".pdf");

    return 0;
}
