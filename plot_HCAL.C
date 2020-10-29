int plot_HCAL()
{
    TFile *f1 = new TFile("DoubleMuon_Run2018A_NANOAOD_plots.root");
    //TFile *f1 = new TFile("DoubleMuon_Run2018A_Run_315512_NANO_DLPHIN_plots.root");

    bool plot_MET = false;
    bool plot_METphi = false;
    bool plot_resolution = true;
    bool plot_response = false;

    bool plot_2D = false;
    bool plot_1D = false;
    bool plot_log = false;
    bool profile_SD = false;
    bool set_SD_range = false;

    std::vector<TString> hist_list =
    {
        //"ZmumuCand_mass", "MET", 
    };

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

    if(plot_resolution)
    {
        hist_list = 
        {
            "UPara_vs_Zpt", "UVert_vs_Zpt", "UCaloPara_vs_Zpt", "UCaloVert_vs_Zpt",
        };

        plot_2D = true;
        plot_log = false;

        x_title = "Z pt";
        profile_SD = true;
        set_SD_range = true;
        SDmin = 0;
        SDmax = 50;
    }

    if(plot_response)
    {
        hist_list = 
        {
            "UPara_ratio_vs_Zpt", "UCaloPara_ratio_vs_Zpt"
        };

        plot_2D = true;
        plot_log = false;

        ymin = 0;
        ymax = 1.2;

        x_title = "Z pt";
        y_title = "- Upara / Zpt";
    }

    for(int i = 0; i < hist_list.size(); i++)
    {
        TString hist_name = hist_list.at(i);

        if(plot_2D)
        {
            TString h1_name = "plots/" + hist_name + "_h";

            TH2F *h1 = (TH2F*)f1->Get(h1_name);

            TCanvas* mycanvas = new TCanvas("mycanvas", "mycanvas", 600, 600);
            gStyle->SetOptStat(kFALSE);

            if(xmin == 999) xmin = h1->GetXaxis()->GetXmin();
            if(xmax == 999) xmax = h1->GetXaxis()->GetXmax();
            if(ymin == 999) ymin = h1->GetYaxis()->GetXmin();
            if(ymax == 999) ymax = h1->GetYaxis()->GetXmax();

            h1->Draw("colz");
            //h1->SetTitle(h1_name);
            h1->GetXaxis()->SetTitle(x_title);
            h1->GetYaxis()->SetTitle(y_title);
            //h1->RebinX(10);
            gPad->SetLogz();

            mycanvas->SetLeftMargin(0.15);
            mycanvas->SetRightMargin(0.15);
            mycanvas->SaveAs("plots_temp/" + hist_name + ".png");

            //TProfile *px = h1->ProfileX("px", 1, -1, "os");
            TProfile *px = h1->ProfileX();
            if(profile_SD)px->BuildOptions(0, 0, "s");
            //px->SetTitle(h1_name);
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

            TProfile* SD_h = (TProfile*)px->Clone(hist_name + "_SD_h");
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
            TString h1_name = "plots/" + hist_name + "_h";
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
    }
    return 0;
}
