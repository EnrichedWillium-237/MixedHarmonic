# include "TFile.h"
# include "TGraphErrors.h"
# include "TCanvas.h"
# include "TF1.h"
# include "TH1.h"
# include "TH2.h"
# include "TLegend.h"
# include "TMath.h"
# include "TPaveText.h"
# include "TStyle.h"
# include <fstream>
# include <iostream>

Bool_t close_plots = kFALSE;
Bool_t gridlines = kFALSE;

# include "style.h"

static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static const int netabins = 12;
static const double etabins[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,
                     0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
static const int nadbins = 6;
static const double adbins[] = {0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
static const int ncentbins = 8;
static const int centBins[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
static const int nanals = 16;
string AnalNames[] = {
  // 0         1             2             3             4             5             6             7
  "v1SP",   "v1SP_mid",   "v1SP_102",   "v1SP_106",   "v1SP_110",   "v1SP_114",   "v1SP_118",   "v1SP_122",
  // 8         9            10            11            12            13            14            15
  "v1SPmc", "v1SPmc_mid", "v1SPmc_102", "v1SPmc_106", "v1SPmc_110", "v1SPmc_114", "v1SPmc_118", "v1SPmc_122"
};

using namespace std;

TH1D * v1Stat_eta[nanals][ncentbins];
TH1D * v1Syst_eta[nanals][ncentbins];

TGraphErrors * grv1Stat_eta[nanals][ncentbins];
TGraphErrors * grv1Syst_eta[nanals][ncentbins];

TFile * tfin;

void intV1_eta()
{

    TH1::SetDefaultSumw2();

    tfin = new TFile("../outputs/final_outputs/v1Int_ErrCalc.root");

    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            v1Stat_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v1Stat_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v1Syst_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v1Syst_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
        }
    }

    //-- turn into TGraphErrors
    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            double neta = v1Stat_eta[i][cbin]->GetNbinsX();
            double x[netabins];
            double y[netabins];
            double xerr[netabins];
            double yerrStat[netabins];
            double xerrSyst[netabins];
            double yerrSyst[netabins];
            for (int ebin = 0; ebin<netabins; ebin++) {
                x[ebin] = v1Stat_eta[i][cbin]->GetBinCenter(ebin+1);
                y[ebin] = v1Stat_eta[i][cbin]->GetBinContent(ebin+1);
                xerr[ebin] = 0;
                yerrStat[ebin] = v1Stat_eta[i][cbin]->GetBinError(ebin+1);
                xerrSyst[ebin] = v1Stat_eta[i][cbin]->GetBinCenter(ebin+1)-v1Stat_eta[i][cbin]->GetBinLowEdge(ebin+1);
                yerrSyst[ebin] = v1Syst_eta[i][cbin]->GetBinError(ebin+1);
            }
            grv1Stat_eta[i][cbin] = new TGraphErrors(neta, x, y, xerr, yerrStat);
            grv1Syst_eta[i][cbin] = new TGraphErrors(neta, x, y, xerrSyst, yerrSyst);
        }
    }

    //-- plotting options

    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            grv1Stat_eta[i][cbin]->SetMarkerColor(kBlue);
            grv1Stat_eta[i][cbin]->SetLineColor(kBlue);
            grv1Stat_eta[i][cbin]->SetMarkerStyle(21);
            grv1Stat_eta[i][cbin]->SetMarkerSize(1.1);

            grv1Syst_eta[i][cbin]->SetMarkerColor(kBlue);
            grv1Syst_eta[i][cbin]->SetLineColor(kBlue);
            grv1Syst_eta[i][cbin]->SetMarkerStyle(21);
            grv1Syst_eta[i][cbin]->SetMarkerSize(1.1);
            grv1Syst_eta[i][cbin]->SetFillColor(kRed-7);
        }
    }

    //-- make plots

    if (!fopen("plots","r")) system("mkdir plots");
    if (!fopen("plots/intv1","r")) system("mkdir plots/intv1");
    if (!fopen("plots/intv1/intv1_eta","r")) system("mkdir plots/intv1/intv1_eta");

    int anal; // choice of analysis


    // v1odd(eta) using the HF event planes
    anal = 7;
    if (!fopen(Form("plots/intv1/intv1_eta/int%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/intv1/intv1_eta/int%s",AnalNames[anal].data()));

    TCanvas * cv1oddIntHF_eta = new TCanvas("cv1oddIntHF_eta","cv1oddIntHF_eta",1100,620);
    cv1oddIntHF_eta->Divide(4,2,0,0);
    TH1D * hv1oddIntHF_eta_tmp = new TH1D("hv1oddIntHF_eta_tmp", "", 100, -2.4, 2.4);
    hv1oddIntHF_eta_tmp->SetTitle("");
    hv1oddIntHF_eta_tmp->SetStats(0);
    hv1oddIntHF_eta_tmp->SetXTitle("#eta");
    hv1oddIntHF_eta_tmp->SetYTitle("v_{1}^{odd}");
    hv1oddIntHF_eta_tmp->GetYaxis()->SetRangeUser(-0.06, 0.06);
    hv1oddIntHF_eta_tmp->GetXaxis()->CenterTitle();
    hv1oddIntHF_eta_tmp->GetYaxis()->CenterTitle();
    hv1oddIntHF_eta_tmp->GetYaxis()->SetNdivisions(509);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1oddIntHF_eta = (TPad *) cv1oddIntHF_eta->cd(cbin+1);
        if (gridlines) padv1oddIntHF_eta->SetGrid();
        if (cbin == 3 || cbin == 7) padv1oddIntHF_eta->SetRightMargin(0.02);
        if (cbin<=3) padv1oddIntHF_eta->SetTopMargin(0.08);
        TH1D * hv1oddIntHF_eta = (TH1D *) hv1oddIntHF_eta_tmp->Clone(Form("hv1oddIntHF_eta_%c",cbin));
        if (cbin == 0) {
            hv1oddIntHF_eta->GetYaxis()->SetTitleSize(0.07);
            hv1oddIntHF_eta->GetYaxis()->SetTitleOffset(1.33);
            hv1oddIntHF_eta->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1oddIntHF_eta->GetXaxis()->SetTitleSize(0.06);
            hv1oddIntHF_eta->GetXaxis()->SetTitleOffset(1.12);
            hv1oddIntHF_eta->GetXaxis()->SetLabelSize(0.06);
            hv1oddIntHF_eta->GetXaxis()->SetLabelOffset(0.018);
            hv1oddIntHF_eta->GetYaxis()->SetTitleSize(0.06);
            hv1oddIntHF_eta->GetYaxis()->SetTitleOffset(1.50);
            hv1oddIntHF_eta->GetYaxis()->SetLabelSize(0.05);
            hv1oddIntHF_eta->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1oddIntHF_eta->GetXaxis()->SetTitleSize(0.07);
            hv1oddIntHF_eta->GetXaxis()->SetTitleOffset(1.00);
            hv1oddIntHF_eta->GetXaxis()->SetLabelSize(0.07);
            hv1oddIntHF_eta->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1oddIntHF_eta->Draw();
        grv1Stat_eta[anal][cbin]->Draw("same p");

        TPaveText * txv1oddIntHF_eta_0;
        if (cbin == 0) txv1oddIntHF_eta_0 = new TPaveText(0.25, 0.04, 0.46, 0.17,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1oddIntHF_eta_0 = new TPaveText(0.06, 0.04, 0.26, 0.17,"NDC");
        else if (cbin == 4) txv1oddIntHF_eta_0 = new TPaveText(0.25, 0.21, 0.46, 0.30,"NDC");
        else txv1oddIntHF_eta_0 = new TPaveText(0.06, 0.21, 0.26, 0.30,"NDC");
        SetTPaveTxt(txv1oddIntHF_eta_0, 18);
        txv1oddIntHF_eta_0->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1oddIntHF_eta_0->Draw();
    }
    cv1oddIntHF_eta->cd(1);
    TPaveText * txv1oddIntHF_eta_1 = new TPaveText(0.22, 0.65, 0.81, 0.86,"NDC");
    SetTPaveTxt(txv1oddIntHF_eta_1, 18);
    txv1oddIntHF_eta_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1oddIntHF_eta_1->AddText("0.3 < p_{T} < 3.0 (GeV/c)");
    txv1oddIntHF_eta_1->Draw();

    TPaveText * txv1oddIntHF_eta_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1oddIntHF_eta_2, 18);
    txv1oddIntHF_eta_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1oddIntHF_eta_2->Draw();

    cv1oddIntHF_eta->Print(Form("plots/intv1/intv1_eta/int%s/v1oddIntHF_eta_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1oddIntHF_eta->Close();



    // v1odd(eta) using the HF event planes with systematics due to antisymmetry
    TCanvas * cv1oddIntHF_eta_syst = new TCanvas("cv1oddIntHF_eta_syst","cv1oddIntHF_eta_syst",1100,620);
    cv1oddIntHF_eta_syst->Divide(4,2,0,0);
    TH1D * hv1oddIntHF_eta_syst_tmp = new TH1D("hv1oddIntHF_eta_syst_tmp", "", 100, -2.4, 2.4);
    hv1oddIntHF_eta_syst_tmp->SetTitle("");
    hv1oddIntHF_eta_syst_tmp->SetStats(0);
    hv1oddIntHF_eta_syst_tmp->SetXTitle("#eta");
    hv1oddIntHF_eta_syst_tmp->SetYTitle("v_{1}^{odd}");
    hv1oddIntHF_eta_syst_tmp->GetYaxis()->SetRangeUser(-0.06, 0.06);
    hv1oddIntHF_eta_syst_tmp->GetXaxis()->CenterTitle();
    hv1oddIntHF_eta_syst_tmp->GetYaxis()->CenterTitle();
    hv1oddIntHF_eta_syst_tmp->GetYaxis()->SetNdivisions(509);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1oddIntHF_eta_syst = (TPad *) cv1oddIntHF_eta_syst->cd(cbin+1);
        if (gridlines) padv1oddIntHF_eta_syst->SetGrid();
        if (cbin == 3 || cbin == 7) padv1oddIntHF_eta_syst->SetRightMargin(0.02);
        if (cbin<=3) padv1oddIntHF_eta_syst->SetTopMargin(0.08);
        TH1D * hv1oddIntHF_eta_syst = (TH1D *) hv1oddIntHF_eta_syst_tmp->Clone(Form("hv1oddIntHF_eta_syst_%c",cbin));
        if (cbin == 0) {
            hv1oddIntHF_eta_syst->GetYaxis()->SetTitleSize(0.07);
            hv1oddIntHF_eta_syst->GetYaxis()->SetTitleOffset(1.33);
            hv1oddIntHF_eta_syst->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1oddIntHF_eta_syst->GetXaxis()->SetTitleSize(0.06);
            hv1oddIntHF_eta_syst->GetXaxis()->SetTitleOffset(1.12);
            hv1oddIntHF_eta_syst->GetXaxis()->SetLabelSize(0.06);
            hv1oddIntHF_eta_syst->GetXaxis()->SetLabelOffset(0.018);
            hv1oddIntHF_eta_syst->GetYaxis()->SetTitleSize(0.06);
            hv1oddIntHF_eta_syst->GetYaxis()->SetTitleOffset(1.50);
            hv1oddIntHF_eta_syst->GetYaxis()->SetLabelSize(0.05);
            hv1oddIntHF_eta_syst->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1oddIntHF_eta_syst->GetXaxis()->SetTitleSize(0.07);
            hv1oddIntHF_eta_syst->GetXaxis()->SetTitleOffset(1.00);
            hv1oddIntHF_eta_syst->GetXaxis()->SetLabelSize(0.07);
            hv1oddIntHF_eta_syst->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1oddIntHF_eta_syst->Draw();
        grv1Syst_eta[anal][cbin]->Draw("same []2");
        grv1Stat_eta[anal][cbin]->Draw("same p");

        TPaveText * txv1oddIntHF_eta_syst_0;
        if (cbin == 0) txv1oddIntHF_eta_syst_0 = new TPaveText(0.25, 0.04, 0.46, 0.17,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1oddIntHF_eta_syst_0 = new TPaveText(0.06, 0.04, 0.26, 0.17,"NDC");
        else if (cbin == 4) txv1oddIntHF_eta_syst_0 = new TPaveText(0.25, 0.21, 0.46, 0.30,"NDC");
        else txv1oddIntHF_eta_syst_0 = new TPaveText(0.06, 0.21, 0.26, 0.30,"NDC");
        SetTPaveTxt(txv1oddIntHF_eta_syst_0, 18);
        txv1oddIntHF_eta_syst_0->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1oddIntHF_eta_syst_0->Draw();
    }
    cv1oddIntHF_eta_syst->cd(1);
    TPaveText * txv1oddIntHF_eta_syst_1 = new TPaveText(0.22, 0.65, 0.81, 0.86,"NDC");
    SetTPaveTxt(txv1oddIntHF_eta_syst_1, 18);
    txv1oddIntHF_eta_syst_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1oddIntHF_eta_syst_1->AddText("0.3 < p_{T} < 3.0 (GeV/c)");
    txv1oddIntHF_eta_syst_1->Draw();

    TPaveText * txv1oddIntHF_eta_syst_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1oddIntHF_eta_syst_2, 18);
    txv1oddIntHF_eta_syst_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1oddIntHF_eta_syst_2->Draw();

    cv1oddIntHF_eta_syst->Print(Form("plots/intv1/intv1_eta/int%s/v1oddIntHF_eta_syst_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1oddIntHF_eta_syst->Close();


}
