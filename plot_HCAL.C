int plot_HCAL()
{
    TString file_name = "UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_origin_recHit_plots";
    //file_name = "UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_M0_energy_plots";
    //file_name = "UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_DLPHIN_1dHB_2dHE_plots";
    //file_name = "UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_DLPHIN_1dHB_2dHE_truncate_plots";
    //file_name = "UL_RSGravitonToQuarkQuark_kMpl01_M_2000_RECO_DLPHIN_1dHB_2dHE_scaled_plots";
    //file_name = "UL_1TeV_pion_gun_RECO_noPU_mahi_energy_plots";

    bool plot_MET = false;
    bool plot_METphi = false;
    bool plot_MET_resolution = false;
    bool plot_MET_response = false;
    bool plot_MET_vs_PU = false;
    bool plot_METphi_vs_PU = false;
    bool plot_CaloJet_vs_GenJet = false;
    bool plot_CaloJet_vs_GenJet_pull = false;
    bool plot_CaloTowerET_vs_eta = false;
    bool plot_leading_jet_ratio = false;
    bool plot_dijet_mass = true;

    bool plot_2D = false;
    bool plot_1D = false;
    bool plot_ratio = false;

    bool plot_underflow_overflow_bins = false;
    bool plot_log = false;
    bool profile_SD = false;
    bool set_SD_range = false;

    TFile *f1 = new TFile("results/" + file_name + ".root");

    std::vector<TString> hist_list =
    {
        "ZmumuCand_mass", "MET", 
    };

    TString hist_dir = "myAna/";

    TString x_title = "";
    TString y_title = "";

    float xmin = 999;
    float xmax = 999;
    float ymin = 999;
    float ymax = 999;
    float SDmin, SDmax;

    if(plot_MET)
    {
        hist_list = 
        {
            "MET_with_Z", "CaloMET_with_Z", "METPara", "METVert", "CaloMETPara", "CaloMETVert",
        };

        plot_1D = true;
        plot_log = true;

        x_title = "MET";
    }

    if(plot_METphi)
    {
        hist_list = 
        {
            "METphi_with_Z", "CaloMETphi_with_Z",
        };

        plot_1D = true;
        plot_log = false;

        x_title = "METphi";
    }

    if(plot_MET_resolution)
    {
        hist_list = 
        {
            "UPara_vs_Zpt", "UVert_vs_Zpt",
        };

        plot_2D = true;
        plot_log = true;

        x_title = "Z pt";
        profile_SD = true;
        set_SD_range = true;
        SDmin = 0;
        SDmax = 50;
    }

    if(plot_MET_response)
    {
        hist_list = 
        {
            "UPara_ratio_vs_Zpt",
        };

        plot_2D = true;
        plot_log = false;

        ymin = 0;
        ymax = 1.2;

        x_title = "Z pt";
        y_title = "- Upara / Zpt";
    }

    if(plot_leading_jet_ratio)
    {
        hist_dir = "myAna/";

        hist_list = 
        {
            //"LeadingCaloJet_ratio"
            "DiJet_ratio", "DiJet_ratio_HB", "DiJet_ratio_HE"
        };

        plot_ratio = true;
        //plot_log = true;
        plot_underflow_overflow_bins = true;

        //ymin = 0;
        //ymax = 60;

        x_title = "CaloJet / GenJet";
        //y_title = "";
    }

    if(plot_dijet_mass)
    {
        hist_dir = "myAna/";

        hist_list = 
        {
            "DiJet_CaloJet_mass", "DiJet_CaloJet_mass_BB", "DiJet_CaloJet_mass_EE", "DiJet_CaloJet_mass_BE"
        };

        plot_ratio = true;
        //plot_log = true;
        //plot_underflow_overflow_bins = true;

        //ymin = 0;
        //ymax = 60;

        x_title = "di-jet mass [GeV]";
        //y_title = "";
    }

    if(plot_MET_vs_PU)
    {
        hist_dir = "myAna/";

        hist_list = 
        {
            "CaloMETBE_vs_PU", "myCaloMETBE_Muon_vs_PU"
        };

        plot_2D = true;
        plot_log = true;

        ymin = 0;
        ymax = 60;

        x_title = "PU";
        y_title = "MET";
    }

    if(plot_METphi_vs_PU)
    {
        hist_dir = "myAna/";

        hist_list = 
        {
            "CaloMETBE_phi_vs_PU", "myCaloMETBE_Muon_phi_vs_PU"
        };

        plot_2D = true;
        plot_log = false;

        //ymin = 0;
        //ymax = 1.2;

        x_title = "PU";
        y_title = "MET phi";
    }

    if(plot_CaloJet_vs_GenJet)
    {
        hist_dir = "myAna/";

        hist_list = 
        {
            //"CaloJet_vs_GenJet", "CaloJet_vs_GenJet_etaL", "CaloJet_vs_GenJet_etaM", "CaloJet_vs_GenJet_etaH"
            "DiJet_CaloJet_vs_GenJet_HB", "DiJet_CaloJet_vs_GenJet_HE"
        };

        plot_2D = true;
        plot_log = true;
        profile_SD = true;

        //ymin = 0;
        //ymax = 1.2;

        x_title = "gen jet pt";
        y_title = "reco jet pt";
        set_SD_range = true;
        SDmin = 0;
        SDmax = 1;
    }

    if(plot_CaloJet_vs_GenJet_pull)
    {
        hist_dir = "myAna/";

        hist_list = 
        {
            "CaloJet_vs_GenJet_pull", "CaloJet_vs_GenJet_etaL_pull", "CaloJet_vs_GenJet_etaM_pull", "CaloJet_vs_GenJet_etaH_pull",
        };

        plot_2D = true;
        plot_log = true;
        //profile_SD = true;

        //ymin = 0;
        //ymax = 1.2;

        x_title = "gen jet pt";
        y_title = "pull";
        //set_SD_range = true;
        //SDmin = 0;
        //SDmax = 1;
    }

    if(plot_CaloTowerET_vs_eta)
    {
        hist_dir = "myAna/";

        hist_list = 
        {
            "CaloTowerET_vs_eta", "CaloTowerEMET_vs_eta", "CaloTowerHadET_vs_eta",
        };

        plot_2D = true;
        plot_log = true;
        //profile_SD = true;

        ymin = 0;
        ymax = 0.1;

        x_title = "eta";
        y_title = "ET / PU";
    }

    for(int i = 0; i < hist_list.size(); i++)
    {
        TString hist_name = hist_list.at(i);

        if(plot_2D)
        {
            TString h1_name = hist_dir + hist_name + "_h";

            TH2F *h1 = (TH2F*)f1->Get(h1_name);

            TCanvas* mycanvas = new TCanvas("mycanvas", "mycanvas", 600, 600);
            gStyle->SetOptStat(kFALSE);

            if(xmin == 999) xmin = h1->GetXaxis()->GetXmin();
            if(xmax == 999) xmax = h1->GetXaxis()->GetXmax();
            if(ymin == 999) ymin = h1->GetYaxis()->GetXmin();
            if(ymax == 999) ymax = h1->GetYaxis()->GetXmax();

            h1->Draw("colz");
            //h1->Draw("surf3");
            h1->SetTitle(hist_name);
            h1->GetXaxis()->SetTitle(x_title);
            h1->GetYaxis()->SetTitle(y_title);
            //h1->RebinX(10);
            if(plot_log) gPad->SetLogz();

            mycanvas->SetLeftMargin(0.15);
            mycanvas->SetRightMargin(0.15);
            mycanvas->SaveAs("plots_temp/" + hist_name + ".png");

            //TProfile *px = h1->ProfileX("px", 1, -1, "os");
            TProfile *px = h1->ProfileX();
            if(profile_SD)px->BuildOptions(0, 0, "s");
            px->SetTitle(hist_name);
            //px->GetYaxis()->SetNdivisions(512);
            px->GetXaxis()->SetRangeUser(xmin, xmax);
            px->GetYaxis()->SetRangeUser(ymin, ymax);
            px->GetYaxis()->SetTitle(y_title);
            px->SetLineColor(kRed);
            //px->SetLineWidth(2);
            //px->SetMarkerStyle(8);
            px->Draw();
            gPad->SetGrid();
            /*px->Fit("pol1", "w");
            TF1 *f = (TF1*)px->GetListOfFunctions()->FindObject("pol1");
            if(f)
            {
                f->SetLineColor(kBlue);
                f->SetLineWidth(1);
                f->SetLineStyle(2); // 2 = "- - -"
            }*/

            //TLine *l=new TLine(xmin, ymin, xmax * 0.8, ymax * 0.8);
            //l->SetLineColor(kBlack);
            //l->Draw("same");

            mycanvas->SaveAs("plots_temp/" + hist_name + "_profile.png");

            TProfile* SD_h = (TProfile*)px->Clone(hist_name + "_SD");
            SD_h->Reset();
            for(int i = 1; i <= px->GetNbinsX(); i++)
            {
                float error = px->GetBinError(i);
                float center = px->GetBinContent(i);
                float rel_error = 0;
                if (center > 0) rel_error = error/center;
                SD_h->SetBinContent(i, error);
                SD_h->SetBinEntries(i, 1);
                //std::cout << i << ", " << SD_h->GetBinContent(i) << ", " << SD_h->GetBinError(i) << std::endl;
            }
            SD_h->GetXaxis()->SetRangeUser(xmin, xmax);
            if(set_SD_range)SD_h->GetYaxis()->SetRangeUser(SDmin, SDmax);
            SD_h->GetXaxis()->SetTitle(x_title);
            SD_h->GetYaxis()->SetTitle(y_title);
            SD_h->Draw("hist");
            mycanvas->SaveAs("plots_temp/" + hist_name + "_SD.png");
        }

        if(plot_1D)
        {
            TString h1_name = hist_dir + hist_name + "_h";
            TH1F *h1 = (TH1F*)f1->Get(h1_name);

            TCanvas* mycanvas = new TCanvas("mycanvas", "mycanvas", 600, 600);
            //gStyle->SetOptStat(kFALSE);

            h1->Draw();
            //h1->Fit("gaus","","",0.7,2);
            h1->GetXaxis()->SetTitle(x_title);
            //h1->SetTitle(hist_name);
            h1->GetYaxis()->SetTitle("Events");
            if (plot_log) gPad->SetLogy();

            //TLine *l=new TLine(1, h1->GetMinimum(),1, h1->GetMaximum());
            //l->SetLineColor(kBlack);
            //l->Draw("same");

            mycanvas->SetLeftMargin(0.15);
            mycanvas->SetRightMargin(0.1);
            mycanvas->SaveAs("plots_temp/" + hist_name + ".png");
        }

        if(plot_ratio)
        {
            TString h1_name = hist_dir + hist_name + "_h";
            TH1F *h1 = (TH1F*)f1->Get(h1_name);

            TCanvas* mycanvas = new TCanvas("mycanvas", "mycanvas", 600, 600);
            //gStyle->SetOptStat(kFALSE);

            h1->Draw("e");
            auto h1_mean = h1->GetMean();
            auto h1_std = h1->GetStdDev();
            auto max_bin = h1->GetMaximumBin();
            if(max_bin != 1)
            {h1_mean = h1->GetXaxis()->GetBinCenter(max_bin);}
            //std::cout << h1_mean << ", " << h1_std << std::endl;
            h1->Fit("gaus", "", "", h1_mean - 1.5 * h1_std, h1_mean + 1.5 * h1_std);
            auto f1 = h1->GetFunction("gaus");
            auto c = f1->GetParameter(0);
            auto mu = f1->GetParameter(1);
            auto sigma = f1->GetParameter(2);
            auto c_err = f1->GetParError(0);
            auto mu_err = f1->GetParError(1);
            auto sigma_err = f1->GetParError(2);
            //std::cout << c << ", " << mu << ", " << sigma << std::endl;

            std::stringstream s1;
            s1 << "#bf{#mu = " << std::setprecision(3) << mu << " #pm " << mu_err << "}";
            TString TS1 = s1.str();
            std::stringstream s2;
            s2 << "#bf{#sigma = " << std::setprecision(3) << sigma << " #pm " << sigma_err << "}";
            TString TS2 = s2.str();
            std::stringstream s3;
            s3 << "#bf{#sigma / #mu = " << std::setprecision(3) << sigma / mu << "}";
            TString TS3 = s3.str();

            h1->GetXaxis()->SetTitle(x_title);
            if(plot_underflow_overflow_bins) h1->GetXaxis()->SetRange(0, h1->GetNbinsX() + 1);
            else h1->GetXaxis()->SetRangeUser(0.5,1.5);
            h1->SetTitle(hist_name);
            h1->GetYaxis()->SetTitle(y_title);
            //gPad->SetLogz();

            TLine *l = new TLine(1, h1->GetMinimum(),1, h1->GetMaximum());
            l->SetLineColor(kBlack);
            l->Draw("same");

            TLatex latex;
            latex.SetTextSize(0.04);
            latex.SetNDC();
            latex.DrawLatex(0.6,0.7,TS1);
            latex.DrawLatex(0.6,0.65,TS2);
            latex.DrawLatex(0.6,0.6,TS3);

            mycanvas->SetLeftMargin(0.15);
            mycanvas->SetRightMargin(0.1);
            mycanvas->SaveAs("plots_temp/" + file_name + "_" + hist_name + ".png");
        }
    }
    return 0;
}
