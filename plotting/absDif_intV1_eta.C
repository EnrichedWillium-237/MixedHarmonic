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

TH1D * v1p_eta[nanals][ncentbins];
TH1D * v1m_eta[nanals][ncentbins];
TH1D * v1odd_eta[nanals][ncentbins];
TH1D * v1even_eta[nanals][ncentbins];

TH1D * v112p_eta[nanals][ncentbins];
TH1D * v112m_eta[nanals][ncentbins];
TH1D * v112odd_eta[nanals][ncentbins];
TH1D * v112even_eta[nanals][ncentbins];

TH1D * v123p_eta[nanals][ncentbins];
TH1D * v123m_eta[nanals][ncentbins];
TH1D * v123odd_eta[nanals][ncentbins];
TH1D * v123even_eta[nanals][ncentbins];

TH1D * absDiffv1odd_eta[nanals][ncentbins];
TH1D * absDiffv1even_eta[nanals][ncentbins];

TH1D * runParms[nanals];

void absDif_intV1_eta()
{

    TH1::SetDefaultSumw2();

    TFile * tfin = new TFile("../outputs/final_outputs/v1Int.root");

    //-- retrieve histograms from final output file
    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            v1p_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v1p_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v1m_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v1m_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v1odd_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v1odd_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v1even_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v1even_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));

            v112p_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v112p_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v112m_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v112m_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v112odd_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v112odd_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v112even_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v112even_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));

            v123p_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v123p_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v123m_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v123m_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v123odd_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v123odd_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v123even_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v123even_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));

            absDiffv1odd_eta[i][cbin] = new TH1D(Form("absDiffv1odd_eta_%s_%d",AnalNames[i].data(),cbin), "", nadbins, adbins);
            absDiffv1even_eta[i][cbin] = new TH1D(Form("absDiffv1even_eta_%s_%d",AnalNames[i].data(),cbin), "", nadbins, adbins);
        }
    }

    // calculate absolute differences
    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            for (int j = 0; j<6; j++) {
                double v1odd_pos = v1odd_eta[i][cbin]->GetBinContent(12-j);
                double v1odd_neg = v1odd_eta[i][cbin]->GetBinContent(j+1);
                double v1odd_pos_err = v1odd_eta[i][cbin]->GetBinError(12-j);
                double v1odd_neg_err = v1odd_eta[i][cbin]->GetBinError(j+1);

                double adodd = fabs(v1odd_pos) - fabs(v1odd_neg);
                double adodd_err = sqrt( pow(v1odd_pos_err,2) + pow(v1odd_neg_err,2) );
                absDiffv1odd_eta[i][cbin]->SetBinContent(j+1, adodd);
                absDiffv1odd_eta[i][cbin]->SetBinError(j+1, adodd_err);

                double v1even_pos = v1even_eta[i][cbin]->GetBinContent(12-j);
                double v1even_neg = v1even_eta[i][cbin]->GetBinContent(j+1);
                double v1even_pos_err = v1even_eta[i][cbin]->GetBinError(12-j);
                double v1even_neg_err = v1even_eta[i][cbin]->GetBinError(j+1);

                double adeven = fabs(v1even_pos) - fabs(v1even_neg);
                double adeven_err = sqrt( pow(v1even_pos_err,2) + pow(v1even_neg_err,2) );
                absDiffv1even_eta[i][cbin]->SetBinContent(j+1, adeven);
                absDiffv1even_eta[i][cbin]->SetBinError(j+1, adeven_err);
            }
        }
    }

    //-- plotting options

    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            v1p_eta[i][cbin]->SetMarkerColor(kRed);
            v1p_eta[i][cbin]->SetLineColor(kRed);
            v1p_eta[i][cbin]->SetMarkerStyle(21);
            v1p_eta[i][cbin]->SetMarkerSize(1.1);

            v1m_eta[i][cbin]->SetMarkerColor(kBlue);
            v1m_eta[i][cbin]->SetLineColor(kBlue);
            v1m_eta[i][cbin]->SetMarkerStyle(21);
            v1m_eta[i][cbin]->SetMarkerSize(1.1);

            v1odd_eta[i][cbin]->SetMarkerColor(kBlue);
            v1odd_eta[i][cbin]->SetLineColor(kBlue);
            v1odd_eta[i][cbin]->SetMarkerStyle(21);
            v1odd_eta[i][cbin]->SetMarkerSize(1.1);

            v1even_eta[i][cbin]->SetMarkerColor(kBlue);
            v1even_eta[i][cbin]->SetLineColor(kBlue);
            v1even_eta[i][cbin]->SetMarkerStyle(21);
            v1even_eta[i][cbin]->SetMarkerSize(1.1);

            absDiffv1odd_eta[i][cbin]->SetMarkerColor(kBlack);
            absDiffv1odd_eta[i][cbin]->SetLineColor(kBlack);
            absDiffv1odd_eta[i][cbin]->SetMarkerStyle(21);
            absDiffv1odd_eta[i][cbin]->SetMarkerSize(1.1);

            absDiffv1even_eta[i][cbin]->SetMarkerColor(kBlack);
            absDiffv1even_eta[i][cbin]->SetLineColor(kBlack);
            absDiffv1even_eta[i][cbin]->SetMarkerStyle(21);
            absDiffv1even_eta[i][cbin]->SetMarkerSize(1.1);
        }
    }
    TF1 * fitleg = new TF1("fitleg","pol0",0,1); // dummy TF1 for legends
    fitleg->SetLineColor(kBlue);


    //-- make plots
    if (!fopen("plots","r")) system("mkdir plots");
    if (!fopen("plots/intv1","r")) system("mkdir plots/intv1");
    if (!fopen("plots/intv1/intv1_eta","r")) system("mkdir plots/intv1/intv1_eta");

    int anal; // choice of analysis


    // integrated v1(eta) using HF+/- for each centrality bin
    anal = 7;
    if (!fopen(Form("plots/intv1/intv1_eta/int%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/intv1/intv1_eta/int%s",AnalNames[anal].data()));

    TCanvas * cv1HFpm_eta = new TCanvas("cv1HFpm_eta","cv1HFpm_eta",1100,620);
    TH1D * hv1HFpm_eta_tmp = new TH1D("hv1HFpm_eta", "", 40, -2.4, 2.4);
    hv1HFpm_eta_tmp->SetTitle("");
    hv1HFpm_eta_tmp->SetStats(0);
    hv1HFpm_eta_tmp->SetXTitle("#eta");
    hv1HFpm_eta_tmp->SetYTitle("v_{1}");
    hv1HFpm_eta_tmp->GetYaxis()->SetRangeUser(-0.19, 0.19);
    hv1HFpm_eta_tmp->GetYaxis()->CenterTitle();
    hv1HFpm_eta_tmp->SetNdivisions(509);
    cv1HFpm_eta->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1HFpm_eta = (TPad *) cv1HFpm_eta->cd(cbin+1);
        if (gridlines) padv1HFpm_eta->SetGrid();
        if (cbin == 3 || cbin == 7) padv1HFpm_eta->SetRightMargin(0.02);
        if (cbin <= 3) padv1HFpm_eta->SetTopMargin(0.08);
        TH1D * hv1HFpm_eta = (TH1D *) hv1HFpm_eta_tmp->Clone(Form("hv1HFpm_eta_%c",cbin));
        if (cbin == 0) {
            hv1HFpm_eta->GetYaxis()->SetTitleSize(0.07);
            hv1HFpm_eta->GetYaxis()->SetTitleOffset(1.33);
            hv1HFpm_eta->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1HFpm_eta->GetXaxis()->SetTitleSize(0.06);
            hv1HFpm_eta->GetXaxis()->SetTitleOffset(1.12);
            hv1HFpm_eta->GetXaxis()->SetLabelSize(0.06);
            hv1HFpm_eta->GetXaxis()->SetLabelOffset(0.018);
            hv1HFpm_eta->GetYaxis()->SetTitleSize(0.06);
            hv1HFpm_eta->GetYaxis()->SetTitleOffset(1.50);
            hv1HFpm_eta->GetYaxis()->SetLabelSize(0.05);
            hv1HFpm_eta->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1HFpm_eta->GetXaxis()->SetTitleSize(0.07);
            hv1HFpm_eta->GetXaxis()->SetTitleOffset(1.00);
            hv1HFpm_eta->GetXaxis()->SetLabelSize(0.07);
            hv1HFpm_eta->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1HFpm_eta->Draw();
        v1p_eta[anal][cbin]->Draw("same");
        v1m_eta[anal][cbin]->Draw("same");

        TPaveText * txv1HFpm_eta;
        if (cbin == 0) txv1HFpm_eta = new TPaveText(0.25, 0.04, 0.46, 0.17,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1HFpm_eta = new TPaveText(0.06, 0.04, 0.26, 0.17,"NDC");
        else if (cbin == 4) txv1HFpm_eta = new TPaveText(0.25, 0.21, 0.46, 0.30,"NDC");
        else txv1HFpm_eta = new TPaveText(0.06, 0.21, 0.26, 0.30,"NDC");
        SetTPaveTxt(txv1HFpm_eta, 18);
        txv1HFpm_eta->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1HFpm_eta->Draw();
    }
    cv1HFpm_eta->cd(1);
    TPaveText * txv1HFpm_eta_1 = new TPaveText(0.22, 0.65, 0.81, 0.86,"NDC");
    SetTPaveTxt(txv1HFpm_eta_1, 18);
    txv1HFpm_eta_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1HFpm_eta_1->AddText("0.3 < p_{T} < 3.0 (GeV/c)");
    txv1HFpm_eta_1->Draw();

    TPaveText * txv1HFpm_eta_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1HFpm_eta_2, 18);
    txv1HFpm_eta_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1HFpm_eta_2->Draw();

    cv1HFpm_eta->cd(2);
    TLegend * legv1HFpm_eta = new TLegend(0.06, 0.67, 0.25, 0.86);
    SetLegend(legv1HFpm_eta, 18);
    legv1HFpm_eta->AddEntry(v1p_eta[anal][0]," v_{1}(HF+)","p");
    legv1HFpm_eta->AddEntry(v1m_eta[anal][0]," v_{1}(HF-)","p");
    legv1HFpm_eta->Draw();

    cv1HFpm_eta->Print(Form("plots/intv1/intv1_eta/int%s/v1_pm_eta_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFpm_eta->Close();



    // integrated v1^even(eta) using HF for each centrality bin

    TCanvas * cv1HFeven_eta = new TCanvas("cv1HFeven_eta","cv1HFeven_eta",1100,620);
    TH1D * hv1HFeven_eta_tmp = new TH1D("hv1HFeven_eta", "", 40, -2.4, 2.4);
    hv1HFeven_eta_tmp->SetTitle("");
    hv1HFeven_eta_tmp->SetStats(0);
    hv1HFeven_eta_tmp->SetXTitle("#eta");
    hv1HFeven_eta_tmp->SetYTitle("v_{1}^{even}");
    hv1HFeven_eta_tmp->GetYaxis()->SetRangeUser(-0.15, 0.0);
    hv1HFeven_eta_tmp->GetYaxis()->CenterTitle();
    hv1HFeven_eta_tmp->SetNdivisions(509);
    cv1HFeven_eta->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1HFeven_eta = (TPad *) cv1HFeven_eta->cd(cbin+1);
        if (gridlines) padv1HFeven_eta->SetGrid();
        if (cbin == 3 || cbin == 7) padv1HFeven_eta->SetRightMargin(0.02);
        if (cbin <= 3) padv1HFeven_eta->SetTopMargin(0.08);
        TH1D * hv1HFeven_eta = (TH1D *) hv1HFeven_eta_tmp->Clone(Form("hv1HFeven_eta_%c",cbin));
        if (cbin == 0) {
            hv1HFeven_eta->GetYaxis()->SetTitleSize(0.07);
            hv1HFeven_eta->GetYaxis()->SetTitleOffset(1.33);
            hv1HFeven_eta->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1HFeven_eta->GetXaxis()->SetTitleSize(0.06);
            hv1HFeven_eta->GetXaxis()->SetTitleOffset(1.12);
            hv1HFeven_eta->GetXaxis()->SetLabelSize(0.06);
            hv1HFeven_eta->GetXaxis()->SetLabelOffset(0.018);
            hv1HFeven_eta->GetYaxis()->SetTitleSize(0.06);
            hv1HFeven_eta->GetYaxis()->SetTitleOffset(1.50);
            hv1HFeven_eta->GetYaxis()->SetLabelSize(0.05);
            hv1HFeven_eta->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1HFeven_eta->GetXaxis()->SetTitleSize(0.07);
            hv1HFeven_eta->GetXaxis()->SetTitleOffset(1.00);
            hv1HFeven_eta->GetXaxis()->SetLabelSize(0.07);
            hv1HFeven_eta->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1HFeven_eta->Draw();
        v1even_eta[anal][cbin]->Draw("same");

        TPaveText * txv1HFeven_eta;
        if (cbin == 0) txv1HFeven_eta = new TPaveText(0.25, 0.04, 0.46, 0.17,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1HFeven_eta = new TPaveText(0.06, 0.04, 0.26, 0.17,"NDC");
        else if (cbin == 4) txv1HFeven_eta = new TPaveText(0.25, 0.21, 0.46, 0.30,"NDC");
        else txv1HFeven_eta = new TPaveText(0.06, 0.21, 0.26, 0.30,"NDC");
        SetTPaveTxt(txv1HFeven_eta, 18);
        txv1HFeven_eta->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1HFeven_eta->Draw();
    }
    cv1HFeven_eta->cd(1);
    TPaveText * txv1HFeven_eta_1 = new TPaveText(0.23, 0.16, 0.82, 0.37,"NDC");
    SetTPaveTxt(txv1HFeven_eta_1, 18);
    txv1HFeven_eta_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1HFeven_eta_1->AddText("0.3 < p_{T} < 3.0 (GeV/c)");
    txv1HFeven_eta_1->Draw();

    TPaveText * txv1HFeven_eta_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1HFeven_eta_2, 18);
    txv1HFeven_eta_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1HFeven_eta_2->Draw();

    cv1HFeven_eta->Print(Form("plots/intv1/intv1_eta/int%s/v1_HFeven_eta_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFeven_eta->Close();



    // integrated v1(eta) using Track+/- for each centrality bin
    anal = 15;
    if (!fopen(Form("plots/intv1/intv1_eta/int%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/intv1/intv1_eta/int%s",AnalNames[anal].data()));

    TCanvas * cv1Trkpm_eta = new TCanvas("cv1Trkpm_eta","cv1Trkpm_eta",1100,620);
    TH1D * hv1Trkpm_eta_tmp = new TH1D("hv1Trkpm_eta", "", 40, -2.4, 2.4);
    hv1Trkpm_eta_tmp->SetTitle("");
    hv1Trkpm_eta_tmp->SetStats(0);
    hv1Trkpm_eta_tmp->SetXTitle("#eta");
    hv1Trkpm_eta_tmp->SetYTitle("v_{1}");
    hv1Trkpm_eta_tmp->GetYaxis()->SetRangeUser(-0.19, 0.19);
    hv1Trkpm_eta_tmp->GetYaxis()->CenterTitle();
    hv1Trkpm_eta_tmp->SetNdivisions(509);
    cv1Trkpm_eta->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1Trkpm_eta = (TPad *) cv1Trkpm_eta->cd(cbin+1);
        if (gridlines) padv1Trkpm_eta->SetGrid();
        if (cbin == 3 || cbin == 7) padv1Trkpm_eta->SetRightMargin(0.02);
        if (cbin == 3 || cbin == 7) padv1Trkpm_eta->SetRightMargin(0.02);
        TH1D * hv1Trkpm_eta = (TH1D *) hv1Trkpm_eta_tmp->Clone(Form("hv1Trkpm_eta_%c",cbin));
        if (cbin == 0) {
            hv1Trkpm_eta->GetYaxis()->SetTitleSize(0.07);
            hv1Trkpm_eta->GetYaxis()->SetTitleOffset(1.33);
            hv1Trkpm_eta->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1Trkpm_eta->GetXaxis()->SetTitleSize(0.06);
            hv1Trkpm_eta->GetXaxis()->SetTitleOffset(1.12);
            hv1Trkpm_eta->GetXaxis()->SetLabelSize(0.06);
            hv1Trkpm_eta->GetXaxis()->SetLabelOffset(0.018);
            hv1Trkpm_eta->GetYaxis()->SetTitleSize(0.06);
            hv1Trkpm_eta->GetYaxis()->SetTitleOffset(1.50);
            hv1Trkpm_eta->GetYaxis()->SetLabelSize(0.05);
            hv1Trkpm_eta->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1Trkpm_eta->GetXaxis()->SetTitleSize(0.07);
            hv1Trkpm_eta->GetXaxis()->SetTitleOffset(1.00);
            hv1Trkpm_eta->GetXaxis()->SetLabelSize(0.07);
            hv1Trkpm_eta->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1Trkpm_eta->Draw();
        v1p_eta[anal][cbin]->Draw("same");
        v1m_eta[anal][cbin]->Draw("same");

        TPaveText * txv1Trkpm_eta;
        if (cbin == 0) txv1Trkpm_eta = new TPaveText(0.25, 0.04, 0.46, 0.17,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1Trkpm_eta = new TPaveText(0.06, 0.04, 0.26, 0.17,"NDC");
        else if (cbin == 4) txv1Trkpm_eta = new TPaveText(0.25, 0.21, 0.46, 0.30,"NDC");
        else txv1Trkpm_eta = new TPaveText(0.06, 0.21, 0.26, 0.30,"NDC");
        SetTPaveTxt(txv1Trkpm_eta, 18);
        txv1Trkpm_eta->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1Trkpm_eta->Draw();
    }
    cv1Trkpm_eta->cd(1);
    TPaveText * txv1Trkpm_eta_1 = new TPaveText(0.22, 0.65, 0.81, 0.86,"NDC");
    SetTPaveTxt(txv1Trkpm_eta_1, 18);
    txv1Trkpm_eta_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1Trkpm_eta_1->AddText("0.3 < p_{T} < 3.0 (GeV/c)");
    txv1Trkpm_eta_1->Draw();

    TPaveText * txv1Trkpm_eta_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1Trkpm_eta_2, 18);
    txv1Trkpm_eta_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1Trkpm_eta_2->Draw();

    cv1Trkpm_eta->cd(2);
    TLegend * legv1Trkpm_eta = new TLegend(0.06, 0.67, 0.25, 0.86);
    SetLegend(legv1Trkpm_eta, 18);
    legv1Trkpm_eta->AddEntry(v1p_eta[anal][0]," Tracker+","p");
    legv1Trkpm_eta->AddEntry(v1m_eta[anal][0]," Tracker-","p");
    legv1Trkpm_eta->Draw();

    cv1Trkpm_eta->Print(Form("plots/intv1/intv1_eta/int%s/v1_pm_eta_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1Trkpm_eta->Close();



    // integrated v1^even(eta) using the tracker (actually called v1odd in output file)

    TCanvas * cv1Trkeven_eta = new TCanvas("cv1Trkeven_eta","cv1Trkeven_eta",1100,620);
    TH1D * hv1Trkeven_eta_tmp = new TH1D("hv1Trkeven_eta", "", 40, -2.4, 2.4);
    hv1Trkeven_eta_tmp->SetTitle("");
    hv1Trkeven_eta_tmp->SetStats(0);
    hv1Trkeven_eta_tmp->SetXTitle("#eta");
    hv1Trkeven_eta_tmp->SetYTitle("v_{1}^{even}");
    hv1Trkeven_eta_tmp->GetYaxis()->SetRangeUser(-0.06, 0.06);
    hv1Trkeven_eta_tmp->GetYaxis()->CenterTitle();
    hv1Trkeven_eta_tmp->SetNdivisions(509);
    cv1Trkeven_eta->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1Trkeven_eta = (TPad *) cv1Trkeven_eta->cd(cbin+1);
        if (gridlines) padv1Trkeven_eta->SetGrid();
        if (cbin == 3 || cbin == 7) padv1Trkeven_eta->SetRightMargin(0.02);
        if (cbin <= 3) padv1Trkeven_eta->SetTopMargin(0.08);
        TH1D * hv1Trkeven_eta = (TH1D *) hv1Trkeven_eta_tmp->Clone(Form("hv1Trkeven_eta_%c",cbin));
        if (cbin == 0) {
            hv1Trkeven_eta->GetYaxis()->SetTitleSize(0.07);
            hv1Trkeven_eta->GetYaxis()->SetTitleOffset(1.33);
            hv1Trkeven_eta->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1Trkeven_eta->GetXaxis()->SetTitleSize(0.06);
            hv1Trkeven_eta->GetXaxis()->SetTitleOffset(1.12);
            hv1Trkeven_eta->GetXaxis()->SetLabelSize(0.06);
            hv1Trkeven_eta->GetXaxis()->SetLabelOffset(0.018);
            hv1Trkeven_eta->GetYaxis()->SetTitleSize(0.06);
            hv1Trkeven_eta->GetYaxis()->SetTitleOffset(1.50);
            hv1Trkeven_eta->GetYaxis()->SetLabelSize(0.05);
            hv1Trkeven_eta->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1Trkeven_eta->GetXaxis()->SetTitleSize(0.07);
            hv1Trkeven_eta->GetXaxis()->SetTitleOffset(1.00);
            hv1Trkeven_eta->GetXaxis()->SetLabelSize(0.07);
            hv1Trkeven_eta->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1Trkeven_eta->Draw();
        v1odd_eta[anal][cbin]->Draw("same");

        TPaveText * txv1Trkeven_eta;
        if (cbin == 0) txv1Trkeven_eta = new TPaveText(0.25, 0.04, 0.46, 0.17,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1Trkeven_eta = new TPaveText(0.06, 0.04, 0.26, 0.17,"NDC");
        else if (cbin == 4) txv1Trkeven_eta = new TPaveText(0.25, 0.21, 0.46, 0.30,"NDC");
        else txv1Trkeven_eta = new TPaveText(0.06, 0.21, 0.26, 0.30,"NDC");
        SetTPaveTxt(txv1Trkeven_eta, 18);
        txv1Trkeven_eta->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1Trkeven_eta->Draw();
    }
    cv1Trkeven_eta->cd(1);
    TPaveText * txv1Trkeven_eta_1 = new TPaveText(0.22, 0.65, 0.81, 0.86,"NDC");
    SetTPaveTxt(txv1Trkeven_eta_1, 18);
    txv1Trkeven_eta_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1Trkeven_eta_1->AddText("0.3 < p_{T} < 3.0 (GeV/c)");
    txv1Trkeven_eta_1->Draw();

    TPaveText * txv1Trkeven_eta_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1Trkeven_eta_2, 18);
    txv1Trkeven_eta_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1Trkeven_eta_2->Draw();

    cv1Trkeven_eta->Print(Form("plots/intv1/intv1_eta/int%s/v1_odd_eta_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1Trkeven_eta->Close();



    // |v1odd(+eta)| - |v1odd(-eta)| odd for the HF
    anal = 7;
    if (!fopen(Form("plots/intv1/intv1_eta/int%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/intv1/intv1_eta/int%s",AnalNames[anal].data()));

    TCanvas * cv1HFoddAbsDiff_eta = new TCanvas("cv1HFoddAbsDiff_eta","cv1HFoddAbsDiff_eta",1100,620);
    TH1D * hv1HFoddAbsDiff_eta_tmp = new TH1D("hv1HFoddAbsDiff_eta_tmp", "", 100, 0, 2.5);
    hv1HFoddAbsDiff_eta_tmp->SetTitle("");
    hv1HFoddAbsDiff_eta_tmp->SetStats(0);
    hv1HFoddAbsDiff_eta_tmp->SetXTitle("|#eta|");
    hv1HFoddAbsDiff_eta_tmp->SetYTitle("|v_{1}(+#eta)| - |v_{1}(-#eta)|");
    hv1HFoddAbsDiff_eta_tmp->GetXaxis()->SetRangeUser(0, 2.4);
    hv1HFoddAbsDiff_eta_tmp->GetYaxis()->SetRangeUser(-.003, 0.003);
    hv1HFoddAbsDiff_eta_tmp->SetNdivisions(509);
    hv1HFoddAbsDiff_eta_tmp->GetXaxis()->CenterTitle();
    hv1HFoddAbsDiff_eta_tmp->GetYaxis()->CenterTitle();
    cv1HFoddAbsDiff_eta->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1HFoddAbsDiff_eta = (TPad *) cv1HFoddAbsDiff_eta->cd(cbin+1);
        if (gridlines) padv1HFoddAbsDiff_eta->SetGrid();
        if (cbin == 3 || cbin == 7) padv1HFoddAbsDiff_eta->SetRightMargin(0.02);
        if (cbin <= 3) padv1HFoddAbsDiff_eta->SetTopMargin(0.08);
        TH1D * hv1HFoddAbsDiff_eta = (TH1D *) hv1HFoddAbsDiff_eta_tmp->Clone(Form("hv1HFoddAbsDiff_eta_%c",cbin));
        if (cbin == 0) {
            hv1HFoddAbsDiff_eta->GetYaxis()->SetTitleSize(0.07);
            hv1HFoddAbsDiff_eta->GetYaxis()->SetTitleOffset(1.33);
            hv1HFoddAbsDiff_eta->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1HFoddAbsDiff_eta->GetXaxis()->SetTitleSize(0.06);
            hv1HFoddAbsDiff_eta->GetXaxis()->SetTitleOffset(1.12);
            hv1HFoddAbsDiff_eta->GetXaxis()->SetLabelSize(0.06);
            hv1HFoddAbsDiff_eta->GetXaxis()->SetLabelOffset(0.018);
            hv1HFoddAbsDiff_eta->GetYaxis()->SetTitleSize(0.06);
            hv1HFoddAbsDiff_eta->GetYaxis()->SetTitleOffset(1.50);
            hv1HFoddAbsDiff_eta->GetYaxis()->SetLabelSize(0.05);
            hv1HFoddAbsDiff_eta->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1HFoddAbsDiff_eta->GetXaxis()->SetTitleSize(0.07);
            hv1HFoddAbsDiff_eta->GetXaxis()->SetTitleOffset(1.00);
            hv1HFoddAbsDiff_eta->GetXaxis()->SetLabelSize(0.07);
            hv1HFoddAbsDiff_eta->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1HFoddAbsDiff_eta->Draw();
        absDiffv1odd_eta[anal][cbin]->Draw("same");

        TF1 * fit1 = new TF1("fit1", "pol0", 0, 2.4);
        fit1->SetLineColor(kBlue);
        absDiffv1odd_eta[anal][cbin]->Fit(fit1,"QR");
        double par0 = fit1->GetParameter(0);
        double par0E = fit1->GetParError(0);
        double par0Chi2 = fit1->GetChisquare();

        TPaveText * txtxv1HFoddAbsDiff_eta_fit;
        if (cbin == 0) txtxv1HFoddAbsDiff_eta_fit = new TPaveText(0.24, 0.07, 0.74, 0.27,"NDC");
        else if (cbin >= 1 && cbin <= 3) txtxv1HFoddAbsDiff_eta_fit = new TPaveText(0.08, 0.07, 0.56, 0.27,"NDC");
        else if (cbin == 4) txtxv1HFoddAbsDiff_eta_fit = new TPaveText(0.24, 0.21, 0.74, 0.39,"NDC");
        else txtxv1HFoddAbsDiff_eta_fit = new TPaveText(0.08, 0.21, 0.56, 0.39,"NDC");
        SetTPaveTxt(txtxv1HFoddAbsDiff_eta_fit, 16);
        txtxv1HFoddAbsDiff_eta_fit->AddText(Form("mean: %0.4f #pm %0.4f",par0,par0E));
        txtxv1HFoddAbsDiff_eta_fit->AddText(Form("#chi^{2}: %0.4f",par0Chi2));
        txtxv1HFoddAbsDiff_eta_fit->Draw();

        TPaveText * txv1HFoddAbsDiff_eta;
        if (cbin == 0) txv1HFoddAbsDiff_eta = new TPaveText(0.75, 0.78, 0.93, 0.87,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1HFoddAbsDiff_eta = new TPaveText(0.68, 0.78, 0.86, 0.87,"NDC");
        else if (cbin == 4) txv1HFoddAbsDiff_eta = new TPaveText(0.75, 0.86, 0.93, 0.95,"NDC");
        else txv1HFoddAbsDiff_eta = new TPaveText(0.68, 0.86, 0.86, 0.95,"NDC");
        SetTPaveTxt(txv1HFoddAbsDiff_eta, 18);
        txv1HFoddAbsDiff_eta->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1HFoddAbsDiff_eta->Draw();
    }
    cv1HFoddAbsDiff_eta->cd(1);
    TPaveText * txv1HFoddAbsDiff_eta_1 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1HFoddAbsDiff_eta_1, 18);
    txv1HFoddAbsDiff_eta_1->AddText("#bf{CMS} #it{Preliminary}");
    txv1HFoddAbsDiff_eta_1->Draw();

    cv1HFoddAbsDiff_eta->Print(Form("plots/intv1/intv1_eta/int%s/AbsDiff_v1HFodd_eta_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFoddAbsDiff_eta->Close();



    // |v1odd(+eta)| - |v1odd(-eta)| even for the HF
    anal = 7;
    if (!fopen(Form("plots/intv1/intv1_eta/int%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/intv1/intv1_eta/int%s",AnalNames[anal].data()));

    TCanvas * cv1HFevenAbsDiff_eta = new TCanvas("cv1HFevenAbsDiff_eta","cv1HFevenAbsDiff_eta",1100,620);
    TH1D * hv1HFevenAbsDiff_eta_tmp = new TH1D("hv1HFevenAbsDiff_eta_tmp", "", 100, 0, 2.5);
    hv1HFevenAbsDiff_eta_tmp->SetTitle("");
    hv1HFevenAbsDiff_eta_tmp->SetStats(0);
    hv1HFevenAbsDiff_eta_tmp->SetXTitle("|#eta|");
    hv1HFevenAbsDiff_eta_tmp->SetYTitle("|v_{1}(+#eta)| - |v_{1}(-#eta)|");
    hv1HFevenAbsDiff_eta_tmp->GetXaxis()->SetRangeUser(0, 2.4);
    hv1HFevenAbsDiff_eta_tmp->GetYaxis()->SetRangeUser(-0.003, 0.003);
    hv1HFevenAbsDiff_eta_tmp->SetNdivisions(509);
    hv1HFevenAbsDiff_eta_tmp->GetXaxis()->CenterTitle();
    hv1HFevenAbsDiff_eta_tmp->GetYaxis()->CenterTitle();
    cv1HFevenAbsDiff_eta->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1HFevenAbsDiff_eta = (TPad *) cv1HFevenAbsDiff_eta->cd(cbin+1);
        if (gridlines) padv1HFevenAbsDiff_eta->SetGrid();
        if (cbin == 3 || cbin == 7) padv1HFevenAbsDiff_eta->SetRightMargin(0.02);
        if (cbin <= 3) padv1HFevenAbsDiff_eta->SetTopMargin(0.08);
        TH1D * hv1HFevenAbsDiff_eta = (TH1D *) hv1HFevenAbsDiff_eta_tmp->Clone(Form("hv1HFevenAbsDiff_eta_%c",cbin));
        if (cbin == 0) {
            hv1HFevenAbsDiff_eta->GetYaxis()->SetTitleSize(0.07);
            hv1HFevenAbsDiff_eta->GetYaxis()->SetTitleOffset(1.33);
            hv1HFevenAbsDiff_eta->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1HFevenAbsDiff_eta->GetXaxis()->SetTitleSize(0.06);
            hv1HFevenAbsDiff_eta->GetXaxis()->SetTitleOffset(1.12);
            hv1HFevenAbsDiff_eta->GetXaxis()->SetLabelSize(0.06);
            hv1HFevenAbsDiff_eta->GetXaxis()->SetLabelOffset(0.018);
            hv1HFevenAbsDiff_eta->GetYaxis()->SetTitleSize(0.06);
            hv1HFevenAbsDiff_eta->GetYaxis()->SetTitleOffset(1.50);
            hv1HFevenAbsDiff_eta->GetYaxis()->SetLabelSize(0.05);
            hv1HFevenAbsDiff_eta->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1HFevenAbsDiff_eta->GetXaxis()->SetTitleSize(0.07);
            hv1HFevenAbsDiff_eta->GetXaxis()->SetTitleOffset(1.00);
            hv1HFevenAbsDiff_eta->GetXaxis()->SetLabelSize(0.07);
            hv1HFevenAbsDiff_eta->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1HFevenAbsDiff_eta->Draw();
        absDiffv1even_eta[anal][cbin]->Draw("same");

        TF1 * fit1 = new TF1("fit1", "pol0", 0, 2.4);
        fit1->SetLineColor(kBlue);
        absDiffv1even_eta[anal][cbin]->Fit(fit1,"QR");
        double par0 = fit1->GetParameter(0);
        double par0E = fit1->GetParError(0);
        double par0Chi2 = fit1->GetChisquare();

        TPaveText * txtxv1HFoddAbsDiff_eta_fit;
        if (cbin == 0) txtxv1HFoddAbsDiff_eta_fit = new TPaveText(0.24, 0.07, 0.74, 0.27,"NDC");
        else if (cbin >= 1 && cbin <= 3) txtxv1HFoddAbsDiff_eta_fit = new TPaveText(0.08, 0.07, 0.56, 0.27,"NDC");
        else if (cbin == 4) txtxv1HFoddAbsDiff_eta_fit = new TPaveText(0.24, 0.21, 0.74, 0.39,"NDC");
        else txtxv1HFoddAbsDiff_eta_fit = new TPaveText(0.08, 0.21, 0.56, 0.39,"NDC");
        SetTPaveTxt(txtxv1HFoddAbsDiff_eta_fit, 16);
        txtxv1HFoddAbsDiff_eta_fit->AddText(Form("mean: %0.4f #pm %0.4f",par0,par0E));
        txtxv1HFoddAbsDiff_eta_fit->AddText(Form("#chi^{2}: %0.4f",par0Chi2));
        txtxv1HFoddAbsDiff_eta_fit->Draw();

        TPaveText * txv1HFoddAbsDiff_eta;
        if (cbin == 0) txv1HFoddAbsDiff_eta = new TPaveText(0.75, 0.78, 0.93, 0.87,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1HFoddAbsDiff_eta = new TPaveText(0.68, 0.78, 0.86, 0.87,"NDC");
        else if (cbin == 4) txv1HFoddAbsDiff_eta = new TPaveText(0.75, 0.86, 0.93, 0.95,"NDC");
        else txv1HFoddAbsDiff_eta = new TPaveText(0.68, 0.86, 0.86, 0.95,"NDC");
        SetTPaveTxt(txv1HFoddAbsDiff_eta, 18);
        txv1HFoddAbsDiff_eta->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1HFoddAbsDiff_eta->Draw();
    }
    cv1HFevenAbsDiff_eta->cd(1);
    TPaveText * txv1HFevenAbsDiff_eta_1 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1HFevenAbsDiff_eta_1, 18);
    txv1HFevenAbsDiff_eta_1->AddText("#bf{CMS} #it{Preliminary}");
    txv1HFevenAbsDiff_eta_1->Draw();

    cv1HFevenAbsDiff_eta->Print(Form("plots/intv1/intv1_eta/int%s/AbsDiff_v1HFeven_eta_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFevenAbsDiff_eta->Close();



    //-- plot both v1odd(eta) and it's absolute differences
    // 0 - 40% centrality
    anal = 7;
    TCanvas * cv1HFodd_and_absdiff_eta_0to40 = new TCanvas("cv1HFodd_and_absdiff_eta_0to40","cv1HFodd_and_absdiff_eta_0to40",1100,700);
    TH1D * hv1HFodd_and_absdiff_eta_0to40_tmp = new TH1D("hv1HFodd_and_absdiff_eta_0to40", "", 40, -2.4, 2.4);
    hv1HFodd_and_absdiff_eta_0to40_tmp->SetTitle("");
    hv1HFodd_and_absdiff_eta_0to40_tmp->SetStats(0);
    hv1HFodd_and_absdiff_eta_0to40_tmp->SetXTitle("#eta");
    hv1HFodd_and_absdiff_eta_0to40_tmp->SetYTitle("v_{1}^{odd}");
    hv1HFodd_and_absdiff_eta_0to40_tmp->GetYaxis()->SetRangeUser(-0.06, 0.06);
    hv1HFodd_and_absdiff_eta_0to40_tmp->GetYaxis()->CenterTitle();
    hv1HFodd_and_absdiff_eta_0to40_tmp->SetNdivisions(509);
    cv1HFodd_and_absdiff_eta_0to40->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1HFodd_and_absdiff_eta_0to40 = (TPad *) cv1HFodd_and_absdiff_eta_0to40->cd(cbin+1);
        if (gridlines) padv1HFodd_and_absdiff_eta_0to40->SetGrid();
        if (cbin == 3 || cbin == 7) padv1HFodd_and_absdiff_eta_0to40->SetRightMargin(0.02);
        if (cbin <= 3) {
            padv1HFodd_and_absdiff_eta_0to40->SetTopMargin(0.08);
            padv1HFodd_and_absdiff_eta_0to40->SetBottomMargin(0.14);
        }
        if (cbin >= 4) padv1HFodd_and_absdiff_eta_0to40->SetTopMargin(0.15);
        TH1D * hv1HFodd_and_absdiff_eta_0to40 = (TH1D *) hv1HFodd_and_absdiff_eta_0to40_tmp->Clone(Form("hv1HFodd_and_absdiff_eta_0to40_%c",cbin));
        if (cbin == 0) {
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetTitleSize(0.06);
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetTitleOffset(1.12);
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetLabelSize(0.06);
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetLabelOffset(0.018);
            hv1HFodd_and_absdiff_eta_0to40->GetYaxis()->SetTitleSize(0.07);
            hv1HFodd_and_absdiff_eta_0to40->GetYaxis()->SetTitleOffset(1.33);
            hv1HFodd_and_absdiff_eta_0to40->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin >= 1 && cbin <= 3) {
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetTitleSize(0.07);
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetTitleOffset(1.00);
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetLabelSize(0.07);
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetLabelOffset(0.008);
        }
        if (cbin >= 4) {
            hv1HFodd_and_absdiff_eta_0to40->SetYTitle("");
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetRangeUser(0, 2.4);
            hv1HFodd_and_absdiff_eta_0to40->GetYaxis()->SetRangeUser(-0.003, 0.003);
        }
        if (cbin == 4) {
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetTitleSize(0.06);
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetTitleOffset(1.12);
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetLabelSize(0.06);
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetLabelOffset(0.018);
            hv1HFodd_and_absdiff_eta_0to40->GetYaxis()->SetTitleSize(0.06);
            hv1HFodd_and_absdiff_eta_0to40->GetYaxis()->SetTitleOffset(1.50);
            hv1HFodd_and_absdiff_eta_0to40->GetYaxis()->SetLabelSize(0.05);
            hv1HFodd_and_absdiff_eta_0to40->GetYaxis()->SetLabelOffset(0.010);
            hv1HFodd_and_absdiff_eta_0to40->SetXTitle("|#eta|");
        }
        if (cbin >=5) {
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetTitleSize(0.07);
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetTitleOffset(1.00);
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetLabelSize(0.07);
            hv1HFodd_and_absdiff_eta_0to40->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1HFodd_and_absdiff_eta_0to40->Draw();

        TPaveText * txv1HFodd_and_absdiff_eta_0to40;
        if (cbin <= 3) {
            v1odd_eta[anal][cbin]->Draw("same");
            if (cbin == 0) txv1HFodd_and_absdiff_eta_0to40 = new TPaveText(0.24, 0.17, 0.46, 0.24,"NDC");
            else if (cbin >= 1 && cbin <= 3) txv1HFodd_and_absdiff_eta_0to40 = new TPaveText(0.05, 0.17, 0.25, 0.24,"NDC");
            SetTPaveTxt(txv1HFodd_and_absdiff_eta_0to40, 18);
            txv1HFodd_and_absdiff_eta_0to40->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
            txv1HFodd_and_absdiff_eta_0to40->Draw();
        } else {
            absDiffv1odd_eta[anal][cbin-4]->Draw("same");

            TF1 * fit1 = new TF1("fit1", "pol0", 0, 2.4);
            fit1->SetLineColor(kBlue);
            absDiffv1odd_eta[anal][cbin-4]->Fit(fit1,"QR");
            double par0 = fit1->GetParameter(0);
            double par0E = fit1->GetParError(0);
            double par0Chi2 = fit1->GetChisquare();

            TPaveText * tx1HFodd_and_absdiff_eta_fit;
            if (cbin == 4) tx1HFodd_and_absdiff_eta_fit = new TPaveText(0.22, 0.21, 0.72, 0.35,"NDC");
            else tx1HFodd_and_absdiff_eta_fit = new TPaveText(0.06, 0.21, 0.54, 0.35,"NDC");
            SetTPaveTxt(tx1HFodd_and_absdiff_eta_fit, 16);
            tx1HFodd_and_absdiff_eta_fit->AddText(Form("mean: %0.4f #pm %0.4f",par0,par0E));
            tx1HFodd_and_absdiff_eta_fit->AddText(Form("#chi^{2}: %0.4f",par0Chi2));
            tx1HFodd_and_absdiff_eta_fit->Draw();
        }

    }
    cv1HFodd_and_absdiff_eta_0to40->cd(1);
    TPaveText * txv1HFodd_and_absdiff_eta_0to40_1 = new TPaveText(0.22, 0.68, 0.81, 0.87,"NDC");
    SetTPaveTxt(txv1HFodd_and_absdiff_eta_0to40_1, 18);
    txv1HFodd_and_absdiff_eta_0to40_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1HFodd_and_absdiff_eta_0to40_1->AddText("0.3 < p_{T} < 3.0 (GeV/c)");
    txv1HFodd_and_absdiff_eta_0to40_1->Draw();

    TPaveText * txv1HFoddAbsDiff_eta_0to40_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1HFoddAbsDiff_eta_0to40_2, 18);
    txv1HFoddAbsDiff_eta_0to40_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1HFoddAbsDiff_eta_0to40_2->Draw();

    cv1HFodd_and_absdiff_eta_0to40->cd(5);
    TLegend * legv1HFoddAbsDiff_eta_0to40 = new TLegend(0.24, 0.67, 0.51, 0.84);
    SetLegend(legv1HFoddAbsDiff_eta_0to40, 18);
    legv1HFoddAbsDiff_eta_0to40->SetHeader("|v_{1}(+#eta)| - |v_{1}(-#eta)|");
    legv1HFoddAbsDiff_eta_0to40->AddEntry(fitleg," pol0 fit","l");
    legv1HFoddAbsDiff_eta_0to40->Draw();

    cv1HFodd_and_absdiff_eta_0to40->Print(Form("plots/intv1/intv1_eta/int%s/v1odd_with_absdiff_eta_0to40cent_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFodd_and_absdiff_eta_0to40->Close();


    // 40 - 80% centrality
    TCanvas * cv1HFodd_and_absdiff_eta_40to80 = new TCanvas("cv1HFodd_and_absdiff_eta_40to80","cv1HFodd_and_absdiff_eta_40to80",1100,700);
    TH1D * hv1HFodd_and_absdiff_eta_40to80_tmp = new TH1D("hv1HFodd_and_absdiff_eta_40to80", "", 40, -2.4, 2.4);
    hv1HFodd_and_absdiff_eta_40to80_tmp->SetTitle("");
    hv1HFodd_and_absdiff_eta_40to80_tmp->SetStats(0);
    hv1HFodd_and_absdiff_eta_40to80_tmp->SetXTitle("#eta");
    hv1HFodd_and_absdiff_eta_40to80_tmp->SetYTitle("v_{1}^{odd}");
    hv1HFodd_and_absdiff_eta_40to80_tmp->GetYaxis()->SetRangeUser(-0.06, 0.06);
    hv1HFodd_and_absdiff_eta_40to80_tmp->GetYaxis()->CenterTitle();
    hv1HFodd_and_absdiff_eta_40to80_tmp->SetNdivisions(509);
    cv1HFodd_and_absdiff_eta_40to80->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1HFodd_and_absdiff_eta_40to80 = (TPad *) cv1HFodd_and_absdiff_eta_40to80->cd(cbin+1);
        if (gridlines) padv1HFodd_and_absdiff_eta_40to80->SetGrid();
        if (cbin == 3 || cbin == 7) padv1HFodd_and_absdiff_eta_40to80->SetRightMargin(0.02);
        if (cbin <= 3) {
            padv1HFodd_and_absdiff_eta_40to80->SetTopMargin(0.08);
            padv1HFodd_and_absdiff_eta_40to80->SetBottomMargin(0.14);
        }
        if (cbin >= 4) padv1HFodd_and_absdiff_eta_40to80->SetTopMargin(0.15);
        TH1D * hv1HFodd_and_absdiff_eta_40to80 = (TH1D *) hv1HFodd_and_absdiff_eta_40to80_tmp->Clone(Form("hv1HFodd_and_absdiff_eta_40to80_%c",cbin));
        if (cbin == 0) {
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetTitleSize(0.06);
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetTitleOffset(1.12);
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetLabelSize(0.06);
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetLabelOffset(0.018);
            hv1HFodd_and_absdiff_eta_40to80->GetYaxis()->SetTitleSize(0.07);
            hv1HFodd_and_absdiff_eta_40to80->GetYaxis()->SetTitleOffset(1.33);
            hv1HFodd_and_absdiff_eta_40to80->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin >= 1 && cbin <= 3) {
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetTitleSize(0.07);
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetTitleOffset(1.00);
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetLabelSize(0.07);
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetLabelOffset(0.008);
        }
        if (cbin >= 4) {
            hv1HFodd_and_absdiff_eta_40to80->SetYTitle("");
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetRangeUser(0, 2.4);
            hv1HFodd_and_absdiff_eta_40to80->GetYaxis()->SetRangeUser(-0.004, 0.004);
        }
        if (cbin == 4) {
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetTitleSize(0.06);
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetTitleOffset(1.12);
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetLabelSize(0.06);
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetLabelOffset(0.018);
            hv1HFodd_and_absdiff_eta_40to80->GetYaxis()->SetTitleSize(0.06);
            hv1HFodd_and_absdiff_eta_40to80->GetYaxis()->SetTitleOffset(1.50);
            hv1HFodd_and_absdiff_eta_40to80->GetYaxis()->SetLabelSize(0.05);
            hv1HFodd_and_absdiff_eta_40to80->GetYaxis()->SetLabelOffset(0.010);
            hv1HFodd_and_absdiff_eta_40to80->SetXTitle("|#eta|");
        }
        if (cbin >=5) {
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetTitleSize(0.07);
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetTitleOffset(1.00);
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetLabelSize(0.07);
            hv1HFodd_and_absdiff_eta_40to80->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1HFodd_and_absdiff_eta_40to80->Draw();

        TPaveText * txv1HFodd_and_absdiff_eta_40to80;
        if (cbin <= 3) {
            v1odd_eta[anal][cbin+4]->Draw("same");
            if (cbin == 0) txv1HFodd_and_absdiff_eta_40to80 = new TPaveText(0.24, 0.17, 0.46, 0.24,"NDC");
            else if (cbin >= 1 && cbin <= 3) txv1HFodd_and_absdiff_eta_40to80 = new TPaveText(0.05, 0.17, 0.25, 0.24,"NDC");
            SetTPaveTxt(txv1HFodd_and_absdiff_eta_40to80, 18);
            txv1HFodd_and_absdiff_eta_40to80->AddText(Form("%d-%d%%",centBins[cbin+4],centBins[cbin+5]));
            txv1HFodd_and_absdiff_eta_40to80->Draw();
        } else {
            absDiffv1odd_eta[anal][cbin]->Draw("same");

            TF1 * fit1 = new TF1("fit1", "pol0", 0, 2.4);
            fit1->SetLineColor(kBlue);
            absDiffv1odd_eta[anal][cbin]->Fit(fit1,"QR");
            double par0 = fit1->GetParameter(0);
            double par0E = fit1->GetParError(0);
            double par0Chi2 = fit1->GetChisquare();

            TPaveText * tx1HFodd_and_absdiff_eta_fit;
            if (cbin == 4) tx1HFodd_and_absdiff_eta_fit = new TPaveText(0.22, 0.21, 0.72, 0.35,"NDC");
            else tx1HFodd_and_absdiff_eta_fit = new TPaveText(0.06, 0.21, 0.54, 0.35,"NDC");
            SetTPaveTxt(tx1HFodd_and_absdiff_eta_fit, 16);
            tx1HFodd_and_absdiff_eta_fit->AddText(Form("mean: %0.4f #pm %0.4f",par0,par0E));
            tx1HFodd_and_absdiff_eta_fit->AddText(Form("#chi^{2}: %0.4f",par0Chi2));
            tx1HFodd_and_absdiff_eta_fit->Draw();
        }

    }
    cv1HFodd_and_absdiff_eta_40to80->cd(1);
    TPaveText * txv1HFodd_and_absdiff_eta_40to80_1 = new TPaveText(0.22, 0.68, 0.81, 0.87,"NDC");
    SetTPaveTxt(txv1HFodd_and_absdiff_eta_40to80_1, 18);
    txv1HFodd_and_absdiff_eta_40to80_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1HFodd_and_absdiff_eta_40to80_1->AddText("0.3 < p_{T} < 3.0 (GeV/c)");
    txv1HFodd_and_absdiff_eta_40to80_1->Draw();

    TPaveText * txv1HFoddAbsDiff_eta_40to80_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1HFoddAbsDiff_eta_40to80_2, 18);
    txv1HFoddAbsDiff_eta_40to80_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1HFoddAbsDiff_eta_40to80_2->Draw();

    cv1HFodd_and_absdiff_eta_40to80->cd(5);
    TLegend * legv1HFoddAbsDiff_eta_40to80 = new TLegend(0.24, 0.67, 0.51, 0.84);
    SetLegend(legv1HFoddAbsDiff_eta_40to80, 18);
    legv1HFoddAbsDiff_eta_40to80->SetHeader("|v_{1}(+#eta)| - |v_{1}(-#eta)|");
    legv1HFoddAbsDiff_eta_40to80->AddEntry(fitleg," pol0 fit","l");
    legv1HFoddAbsDiff_eta_40to80->Draw();

    cv1HFodd_and_absdiff_eta_40to80->Print(Form("plots/intv1/intv1_eta/int%s/v1odd_with_absdiff_eta_40to80cent_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFodd_and_absdiff_eta_40to80->Close();


}
