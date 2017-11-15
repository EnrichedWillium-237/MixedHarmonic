# include "TFile.h"
# include "TGraphErrors.h"
# include "TCanvas.h"
# include "TH1D.h"
# include "TH2D.h"
# include "TLegend.h"
# include "TMath.h"
# include "TPaveText.h"
# include "TStyle.h"
# include <fstream>
# include <iostream>

Bool_t close_plots = kTRUE;
Bool_t gridlines = kFALSE;

# include "style.h"

static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static const int netabins = 12;
static const double etabins[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,
                     0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
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

TH1D * v1p_eta[nanals][ncentbins][nptbins];
TH1D * v1m_eta[nanals][ncentbins][nptbins];
TH1D * v1odd_eta[nanals][ncentbins][nptbins];
TH1D * v1even_eta[nanals][ncentbins][nptbins];

TH1D * v112p_eta[nanals][ncentbins][nptbins];
TH1D * v112m_eta[nanals][ncentbins][nptbins];
TH1D * v112odd_eta[nanals][ncentbins][nptbins];
TH1D * v112even_eta[nanals][ncentbins][nptbins];

TH1D * v123p_eta[nanals][ncentbins][nptbins];
TH1D * v123m_eta[nanals][ncentbins][nptbins];
TH1D * v123odd_eta[nanals][ncentbins][nptbins];
TH1D * v123even_eta[nanals][ncentbins][nptbins];

TH1D * runParms[nanals];

TFile * tfin;

void diffV1_eta()
{

    TH1::SetDefaultSumw2();

    tfin = new TFile("../outputs/final_outputs/v1Diff.root","read");

    //-- retrieve histograms from final output file
    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            for (int pbin = 0; pbin<nptbins; pbin++) {
                string tag0 = Form("cent_%d-%d/pt_%0.1f-%0.1f",centBins[cbin],centBins[cbin+1],ptbins[pbin],ptbins[pbin+1]);
                string tag1 = Form("%s_c%d_p%d",AnalNames[i].data(),cbin,pbin);

                v1p_eta[i][cbin][pbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%s/v1p_eta_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v1m_eta[i][cbin][pbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%s/v1m_eta_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v1odd_eta[i][cbin][pbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%s/v1odd_eta_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v1even_eta[i][cbin][pbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%s/v1even_eta_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                /*
                v112p_eta[i][cbin][pbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%s/v112p_eta_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v112m_eta[i][cbin][pbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%s/v112m_eta_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v112odd_eta[i][cbin][pbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%s/v112odd_eta_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v112even_eta[i][cbin][pbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%s/v112even_eta_%s",AnalNames[i].data(),tag0.data(),tag1.data()));

                v123p_eta[i][cbin][pbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%s/v123p_eta_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v123m_eta[i][cbin][pbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%s/v123m_eta_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v123odd_eta[i][cbin][pbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%s/v123odd_eta_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v123even_eta[i][cbin][pbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%s/v123even_eta_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                */
            }
        }
    }

    //-- plotting options

    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            for (int pbin = 0; pbin<nptbins; pbin++) {
                v1p_eta[i][cbin][pbin]->SetMarkerColor(kRed);
                v1p_eta[i][cbin][pbin]->SetLineColor(kRed);
                v1p_eta[i][cbin][pbin]->SetMarkerStyle(21);
                v1p_eta[i][cbin][pbin]->SetMarkerSize(1.1);

                v1m_eta[i][cbin][pbin]->SetMarkerColor(kBlue);
                v1m_eta[i][cbin][pbin]->SetLineColor(kBlue);
                v1m_eta[i][cbin][pbin]->SetMarkerStyle(21);
                v1m_eta[i][cbin][pbin]->SetMarkerSize(1.1);

                v1odd_eta[i][cbin][pbin]->SetMarkerColor(kRed);
                v1odd_eta[i][cbin][pbin]->SetLineColor(kRed);
                v1odd_eta[i][cbin][pbin]->SetMarkerStyle(21);
                v1odd_eta[i][cbin][pbin]->SetMarkerSize(1.1);

                v1even_eta[i][cbin][pbin]->SetMarkerColor(kBlue);
                v1even_eta[i][cbin][pbin]->SetLineColor(kBlue);
                v1even_eta[i][cbin][pbin]->SetMarkerStyle(21);
                v1even_eta[i][cbin][pbin]->SetMarkerSize(1.1);
            }
        }
    }


    //-- make plots
    if (!fopen("plots","r")) system("mkdir plots");
    if (!fopen("plots/diffv1","r")) system("mkdir plots/diffv1");
    if (!fopen("plots/diffv1/diffv1_eta","r")) system("mkdir plots/diffv1/diffv1_eta");

    int centMarkerStyle[] = {21, 24, 20, 25, 33, 27, 34, 28};
    float centMarkerSize[] = {1.1, 1.2, 1.2, 1.1, 1.7, 1.7, 1.5};

    int anal; // choice of analysis


    // differential v1(eta) using HF+/- for each centrality bin
    anal = 7;
    if (!fopen(Form("plots/diffv1/diffv1_eta/diff%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/diffv1/diffv1_eta/diff%s",AnalNames[anal].data()));

    TCanvas * cv1HFpm_eta[ncentbins];
    TH1D * hv1HFpm_eta_tmp = new TH1D("hv1HFpm_eta_tmp", "", 40, -2.5, 2.5);
    hv1HFpm_eta_tmp->SetTitle("");
    hv1HFpm_eta_tmp->SetStats(0);
    hv1HFpm_eta_tmp->SetXTitle("#eta");
    hv1HFpm_eta_tmp->SetYTitle("v_{1}");
    hv1HFpm_eta_tmp->GetYaxis()->SetRangeUser(-0.3, 0.3);
    hv1HFpm_eta_tmp->SetNdivisions(509);
    for (int cbin = 0; cbin<ncentbins; cbin++) {

        cv1HFpm_eta[cbin] = new TCanvas(Form("cv1HFpm_eta_cent%d-%d",centBins[cbin],centBins[cbin+1]),"cv1HFpm_eta",1100,850);
        cv1HFpm_eta[cbin]->Divide(4,3,0,0);

        for (int pbin = 0; pbin<netabins; pbin++) {
            TPad * padv1HFpm_eta = (TPad *) cv1HFpm_eta[cbin]->cd(pbin+1);
            if (gridlines) padv1HFpm_eta->SetGrid();
            if (pbin == 3 || pbin == 7 || pbin == 11) padv1HFpm_eta->SetRightMargin(0.02);
            TH1D * hv1HFpm_eta = (TH1D *) hv1HFpm_eta_tmp->Clone(Form("hv1HFpm_eta_%c_%d",cbin,pbin));
            // hv1HFpm_eta->GetYaxis()->SetRangeUser(-0.14-0.1*cbin,0.14+0.1*cbin);
            if (pbin == 0 || pbin == 4) {
                hv1HFpm_eta->GetYaxis()->CenterTitle();
                hv1HFpm_eta->GetYaxis()->SetTitleSize(0.07);
                hv1HFpm_eta->GetYaxis()->SetTitleOffset(1.34);
                hv1HFpm_eta->GetYaxis()->SetLabelSize(0.06);
            } else if (pbin == 8) {
                hv1HFpm_eta->GetXaxis()->SetTitleSize(0.06);
                hv1HFpm_eta->GetXaxis()->SetTitleOffset(1.14);
                hv1HFpm_eta->GetYaxis()->CenterTitle();
                hv1HFpm_eta->GetYaxis()->SetTitleSize(0.06);
                hv1HFpm_eta->GetYaxis()->SetTitleOffset(1.48);
                hv1HFpm_eta->GetYaxis()->SetLabelSize(0.05);
            } else if (pbin>=9) {
                hv1HFpm_eta->GetXaxis()->SetTitleSize(0.07);
                hv1HFpm_eta->GetXaxis()->SetTitleOffset(1.00);
                hv1HFpm_eta->GetXaxis()->SetLabelSize(0.06);
                hv1HFpm_eta->GetXaxis()->SetLabelOffset(0.005);
            }
            hv1HFpm_eta->Draw();
            v1p_eta[anal][cbin][pbin]->Draw("same");
            v1m_eta[anal][cbin][pbin]->Draw("same");

            TPaveText * txv1HFpm_eta;
            if (pbin == 0 || pbin == 4) txv1HFpm_eta = new TPaveText(0.23, 0.05, 0.62, 0.15,"NDC");
            else if (pbin == 8) txv1HFpm_eta = new TPaveText(0.23, 0.19, 0.62, 0.29,"NDC");
            else if (pbin>=9) txv1HFpm_eta = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
            else txv1HFpm_eta = new TPaveText(0.05, 0.05, 0.50, 0.15,"NDC");
            SetTPaveTxt(txv1HFpm_eta, 18);
            txv1HFpm_eta->AddText(Form("%0.2f<p_{T}<%0.2f (GeV/c)",ptbins[pbin],ptbins[pbin+1]));
            txv1HFpm_eta->Draw();
        }
        cv1HFpm_eta[cbin]->cd(1);
        TPaveText * txv1HFpm_eta_cent = new TPaveText(0.24, 0.84, 0.38, 0.95,"NDC");
        SetTPaveTxt(txv1HFpm_eta_cent, 18);
        txv1HFpm_eta_cent->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1HFpm_eta_cent->Draw();

        cv1HFpm_eta[cbin]->cd(2);
        TLegend * legv1HFpm_eta = new TLegend(0.07, 0.74, 0.28, 0.94);
        SetLegend(legv1HFpm_eta, 18);
        legv1HFpm_eta->AddEntry(v1p_eta[anal][cbin][0],"HF+","p");
        legv1HFpm_eta->AddEntry(v1m_eta[anal][cbin][0],"HF-","p");
        legv1HFpm_eta->Draw();

        cv1HFpm_eta[cbin]->Print(Form("plots/diffv1/diffv1_eta/diff%s/v1_pm_eta_%s_cent%d-%d.png",AnalNames[anal].data(),AnalNames[anal].data(),centBins[cbin],centBins[cbin+1]),"png");
        if (close_plots) cv1HFpm_eta[cbin]->Close();
    }



    // differential v1(eta) using HF+/- with all centralities on the same plot

    TCanvas * cv1HFpm_etaCent;
    TH1D * hv1HFpm_etaCent_tmp = new TH1D("hv1HFpm_etaCent_tmp", "", 40, -2.5, 2.5);
    hv1HFpm_etaCent_tmp->SetTitle("");
    hv1HFpm_etaCent_tmp->SetStats(0);
    hv1HFpm_etaCent_tmp->SetXTitle("#eta");
    hv1HFpm_etaCent_tmp->SetYTitle("v_{1}");
    hv1HFpm_etaCent_tmp->GetYaxis()->SetRangeUser(-0.26,0.26);
    hv1HFpm_etaCent_tmp->SetNdivisions(509);

    cv1HFpm_etaCent = new TCanvas("cv1HFpm_etaCent","cv1HFpm_etaCent",1100,850);
    cv1HFpm_etaCent->Divide(4,3,0,0);

    TLegend * legv1HFpm_etaCent_0 = new TLegend(0.23, 0.65, 0.45, 0.94);
    SetLegend(legv1HFpm_etaCent_0, 16);

    for (int pbin = 0; pbin<netabins; pbin++) {
        TPad * padv1HFpm_etaCent = (TPad *) cv1HFpm_etaCent->cd(pbin+1);
        if (gridlines) padv1HFpm_etaCent->SetGrid();
        if (pbin == 3 || pbin == 7 || pbin == 11) padv1HFpm_etaCent->SetRightMargin(0.02);
        TH1D * hv1HFpm_etaCent = (TH1D *) hv1HFpm_etaCent_tmp->Clone(Form("hv1HFpm_etaCent_%d",pbin));
        if (pbin == 0 || pbin == 4) {
            hv1HFpm_etaCent->GetYaxis()->CenterTitle();
            hv1HFpm_etaCent->GetYaxis()->SetTitleSize(0.07);
            hv1HFpm_etaCent->GetYaxis()->SetTitleOffset(1.34);
            hv1HFpm_etaCent->GetYaxis()->SetLabelSize(0.06);
        } else if (pbin == 8) {
            hv1HFpm_etaCent->GetXaxis()->SetTitleSize(0.06);
            hv1HFpm_etaCent->GetXaxis()->SetTitleOffset(1.14);
            hv1HFpm_etaCent->GetYaxis()->CenterTitle();
            hv1HFpm_etaCent->GetYaxis()->SetTitleSize(0.06);
            hv1HFpm_etaCent->GetYaxis()->SetTitleOffset(1.48);
            hv1HFpm_etaCent->GetYaxis()->SetLabelSize(0.05);
        } else if (pbin>=9) {
            hv1HFpm_etaCent->GetXaxis()->SetTitleSize(0.07);
            hv1HFpm_etaCent->GetXaxis()->SetTitleOffset(1.00);
            hv1HFpm_etaCent->GetXaxis()->SetLabelSize(0.06);
            hv1HFpm_etaCent->GetXaxis()->SetLabelOffset(0.005);
        }
        hv1HFpm_etaCent->Draw();
        TH1D * v1p_eta_tmp[ncentbins];
        TH1D * v1m_eta_tmp[ncentbins];
        for (int cbin = 0; cbin<5; cbin++) {
            v1p_eta_tmp[cbin] = (TH1D *) v1p_eta[anal][cbin][pbin]->Clone(Form("v1p_eta_tmp%d",cbin));
            v1m_eta_tmp[cbin] = (TH1D *) v1m_eta[anal][cbin][pbin]->Clone(Form("v1m_eta_tmp%d",cbin));
            v1p_eta_tmp[cbin]->SetMarkerStyle(centMarkerStyle[cbin]);
            v1p_eta_tmp[cbin]->SetMarkerSize(centMarkerSize[cbin]);
            v1m_eta_tmp[cbin]->SetMarkerStyle(centMarkerStyle[cbin]);
            v1m_eta_tmp[cbin]->SetMarkerSize(centMarkerSize[cbin]);
            v1p_eta_tmp[cbin]->Draw("same");
            v1m_eta_tmp[cbin]->Draw("same");
            if (pbin == 0) {
                legv1HFpm_etaCent_0->AddEntry(v1m_eta_tmp[cbin],Form("%d-%d%%",centBins[cbin],centBins[cbin+1]),"p");
            }
        }

        TPaveText * txv1HFpm_etaCent;
        if (pbin == 0 || pbin == 4) txv1HFpm_etaCent = new TPaveText(0.23, 0.05, 0.62, 0.15,"NDC");
        else if (pbin == 8) txv1HFpm_etaCent = new TPaveText(0.23, 0.19, 0.62, 0.29,"NDC");
        else if (pbin>=9) txv1HFpm_etaCent = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
        else txv1HFpm_etaCent = new TPaveText(0.05, 0.05, 0.50, 0.15,"NDC");
        SetTPaveTxt(txv1HFpm_etaCent, 18);
        txv1HFpm_etaCent->AddText(Form("%0.2f<p_{T}<%0.2f (GeV/c)",ptbins[pbin],ptbins[pbin+1]));
        txv1HFpm_etaCent->Draw();
    }

    cv1HFpm_etaCent->cd(1);
    legv1HFpm_etaCent_0->Draw();

    cv1HFpm_etaCent->cd(2);
    TLegend * legv1HFpm_etaCent_1 = new TLegend(0.07, 0.74, 0.28, 0.94);
    SetLegend(legv1HFpm_etaCent_1, 18);
    legv1HFpm_etaCent_1->AddEntry(v1p_eta[anal][0][0],"HF+","p");
    legv1HFpm_etaCent_1->AddEntry(v1m_eta[anal][0][0],"HF-","p");
    legv1HFpm_etaCent_1->Draw();

    cv1HFpm_etaCent->Print(Form("plots/diffv1/diffv1_eta/diff%s/v1_pm_eta_CentScan_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFpm_etaCent->Close();



    // differential v1(eta) using HF odd/even for each centrality bin
    anal = 7;

    TCanvas * cv1HFoddeven_eta[ncentbins];
    TH1D * hv1HFoddeven_eta_tmp = new TH1D("hv1HFoddeven_eta_tmp", "", 40, -2.5, 2.5);
    hv1HFoddeven_eta_tmp->SetTitle("");
    hv1HFoddeven_eta_tmp->SetStats(0);
    hv1HFoddeven_eta_tmp->SetXTitle("#eta");
    hv1HFoddeven_eta_tmp->SetYTitle("v_{1}");
    hv1HFoddeven_eta_tmp->GetYaxis()->SetRangeUser(-0.21, 0.21);
    hv1HFoddeven_eta_tmp->SetNdivisions(509);
    for (int cbin = 0; cbin<ncentbins; cbin++) {

        cv1HFoddeven_eta[cbin] = new TCanvas(Form("cv1HFoddeven_eta_cent%d-%d",centBins[cbin],centBins[cbin+1]),"cv1HFoddeven_eta",1100,850);
        cv1HFoddeven_eta[cbin]->Divide(4,3,0,0);

        for (int pbin = 0; pbin<netabins; pbin++) {
            TPad * padv1HFoddeven_eta = (TPad *) cv1HFoddeven_eta[cbin]->cd(pbin+1);
            if (gridlines) padv1HFoddeven_eta->SetGrid();
            if (pbin == 3 || pbin == 7 || pbin == 11) padv1HFoddeven_eta->SetRightMargin(0.02);
            TH1D * hv1HFoddeven_eta = (TH1D *) hv1HFoddeven_eta_tmp->Clone(Form("hv1HFoddeven_eta_%c_%d",cbin,pbin));
            // hv1HFoddeven_eta->GetYaxis()->SetRangeUser(-0.14-0.06*cbin,0.14+0.06*cbin);
            if (pbin == 0 || pbin == 4) {
                hv1HFoddeven_eta->GetYaxis()->CenterTitle();
                hv1HFoddeven_eta->GetYaxis()->SetTitleSize(0.07);
                hv1HFoddeven_eta->GetYaxis()->SetTitleOffset(1.34);
                hv1HFoddeven_eta->GetYaxis()->SetLabelSize(0.06);
            } else if (pbin == 8) {
                hv1HFoddeven_eta->GetXaxis()->SetTitleSize(0.06);
                hv1HFoddeven_eta->GetXaxis()->SetTitleOffset(1.14);
                hv1HFoddeven_eta->GetYaxis()->CenterTitle();
                hv1HFoddeven_eta->GetYaxis()->SetTitleSize(0.06);
                hv1HFoddeven_eta->GetYaxis()->SetTitleOffset(1.48);
                hv1HFoddeven_eta->GetYaxis()->SetLabelSize(0.05);
            } else if (pbin>=9) {
                hv1HFoddeven_eta->GetXaxis()->SetTitleSize(0.07);
                hv1HFoddeven_eta->GetXaxis()->SetTitleOffset(1.00);
                hv1HFoddeven_eta->GetXaxis()->SetLabelSize(0.06);
                hv1HFoddeven_eta->GetXaxis()->SetLabelOffset(0.005);
            }
            hv1HFoddeven_eta->Draw();
            v1odd_eta[anal][cbin][pbin]->Draw("same");
            v1even_eta[anal][cbin][pbin]->Draw("same");

            TPaveText * txv1HFoddeven_eta;
            if (pbin == 0 || pbin == 4) txv1HFoddeven_eta = new TPaveText(0.23, 0.05, 0.62, 0.15,"NDC");
            else if (pbin == 8) txv1HFoddeven_eta = new TPaveText(0.23, 0.19, 0.62, 0.29,"NDC");
            else if (pbin>=9) txv1HFoddeven_eta = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
            else txv1HFoddeven_eta = new TPaveText(0.05, 0.05, 0.50, 0.15,"NDC");
            SetTPaveTxt(txv1HFoddeven_eta, 18);
            txv1HFoddeven_eta->AddText(Form("%0.2f<p_{T}<%0.2f (GeV/c)",ptbins[pbin],ptbins[pbin+1]));
            txv1HFoddeven_eta->Draw();
        }
        cv1HFoddeven_eta[cbin]->cd(1);
        TPaveText * txv1HFoddeven_eta_cent = new TPaveText(0.24, 0.84, 0.38, 0.95,"NDC");
        SetTPaveTxt(txv1HFoddeven_eta_cent, 18);
        txv1HFoddeven_eta_cent->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1HFoddeven_eta_cent->Draw();

        cv1HFoddeven_eta[cbin]->cd(2);
        TLegend * legv1HFoddeven_eta = new TLegend(0.07, 0.74, 0.28, 0.94);
        SetLegend(legv1HFoddeven_eta, 18);
        legv1HFoddeven_eta->AddEntry(v1odd_eta[anal][cbin][0],"v_{1}^{odd}","p");
        legv1HFoddeven_eta->AddEntry(v1even_eta[anal][cbin][0],"v_{1}^{even}","p");
        legv1HFoddeven_eta->Draw();

        cv1HFoddeven_eta[cbin]->Print(Form("plots/diffv1/diffv1_eta/diff%s/v1_oddeven_eta_%s_cent%d-%d.png",AnalNames[anal].data(),AnalNames[anal].data(),centBins[cbin],centBins[cbin+1]),"png");
        if (close_plots) cv1HFoddeven_eta[cbin]->Close();
    }



    // differential v1(eta) using HF odd/even with all centralities on the same plot

    TCanvas * cv1HFoddeven_etaCent;
    TH1D * hv1HFoddeven_etaCent_tmp = new TH1D("hv1HFoddeven_etaCent_tmp", "", 40, -2.5, 2.5);
    hv1HFoddeven_etaCent_tmp->SetTitle("");
    hv1HFoddeven_etaCent_tmp->SetStats(0);
    hv1HFoddeven_etaCent_tmp->SetXTitle("#eta");
    hv1HFoddeven_etaCent_tmp->SetYTitle("v_{1}");
    hv1HFoddeven_etaCent_tmp->GetYaxis()->SetRangeUser(-0.16, 0.16);
    hv1HFoddeven_etaCent_tmp->SetNdivisions(509);

    cv1HFoddeven_etaCent = new TCanvas("cv1HFoddeven_etaCent","cv1HFoddeven_etaCent",1100,850);
    cv1HFoddeven_etaCent->Divide(4,3,0,0);

    TLegend * legv1HFoddeven_etaCent_0 = new TLegend(0.23, 0.65, 0.45, 0.94);
    SetLegend(legv1HFoddeven_etaCent_0, 16);

    for (int pbin = 0; pbin<netabins; pbin++) {
        TPad * padv1HFoddeven_etaCent = (TPad *) cv1HFoddeven_etaCent->cd(pbin+1);
        if (gridlines) padv1HFoddeven_etaCent->SetGrid();
        if (pbin == 3 || pbin == 7 || pbin == 11) padv1HFoddeven_etaCent->SetRightMargin(0.02);
        TH1D * hv1HFoddeven_etaCent = (TH1D *) hv1HFoddeven_etaCent_tmp->Clone(Form("hv1HFoddeven_etaCent_%d",pbin));
        if (pbin == 0 || pbin == 4) {
            hv1HFoddeven_etaCent->GetYaxis()->CenterTitle();
            hv1HFoddeven_etaCent->GetYaxis()->SetTitleSize(0.07);
            hv1HFoddeven_etaCent->GetYaxis()->SetTitleOffset(1.34);
            hv1HFoddeven_etaCent->GetYaxis()->SetLabelSize(0.06);
        } else if (pbin == 8) {
            hv1HFoddeven_etaCent->GetXaxis()->SetTitleSize(0.06);
            hv1HFoddeven_etaCent->GetXaxis()->SetTitleOffset(1.14);
            hv1HFoddeven_etaCent->GetYaxis()->CenterTitle();
            hv1HFoddeven_etaCent->GetYaxis()->SetTitleSize(0.06);
            hv1HFoddeven_etaCent->GetYaxis()->SetTitleOffset(1.48);
            hv1HFoddeven_etaCent->GetYaxis()->SetLabelSize(0.05);
        } else if (pbin>=9) {
            hv1HFoddeven_etaCent->GetXaxis()->SetTitleSize(0.07);
            hv1HFoddeven_etaCent->GetXaxis()->SetTitleOffset(1.00);
            hv1HFoddeven_etaCent->GetXaxis()->SetLabelSize(0.06);
            hv1HFoddeven_etaCent->GetXaxis()->SetLabelOffset(0.005);
        }
        hv1HFoddeven_etaCent->Draw();
        TH1D * v1odd_eta_tmp[ncentbins];
        TH1D * v1even_eta_tmp[ncentbins];
        for (int cbin = 0; cbin<5; cbin++) {
            v1odd_eta_tmp[cbin] = (TH1D *) v1odd_eta[anal][cbin][pbin]->Clone(Form("v1odd_eta_tmp%d",cbin));
            v1even_eta_tmp[cbin] = (TH1D *) v1even_eta[anal][cbin][pbin]->Clone(Form("v1even_eta_tmp%d",cbin));
            v1odd_eta_tmp[cbin]->SetMarkerStyle(centMarkerStyle[cbin]);
            v1odd_eta_tmp[cbin]->SetMarkerSize(centMarkerSize[cbin]);
            v1even_eta_tmp[cbin]->SetMarkerStyle(centMarkerStyle[cbin]);
            v1even_eta_tmp[cbin]->SetMarkerSize(centMarkerSize[cbin]);
            v1odd_eta_tmp[cbin]->Draw("same");
            v1even_eta_tmp[cbin]->Draw("same");
            if (pbin == 0) {
                legv1HFoddeven_etaCent_0->AddEntry(v1even_eta_tmp[cbin],Form("%d-%d%%",centBins[cbin],centBins[cbin+1]),"p");
            }
        }

        TPaveText * txv1HFoddeven_etaCent;
        if (pbin == 0 || pbin == 4) txv1HFoddeven_etaCent = new TPaveText(0.23, 0.05, 0.62, 0.15,"NDC");
        else if (pbin == 8) txv1HFoddeven_etaCent = new TPaveText(0.23, 0.19, 0.62, 0.29,"NDC");
        else if (pbin>=9) txv1HFoddeven_etaCent = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
        else txv1HFoddeven_etaCent = new TPaveText(0.05, 0.05, 0.50, 0.15,"NDC");
        SetTPaveTxt(txv1HFoddeven_etaCent, 18);
        txv1HFoddeven_etaCent->AddText(Form("%0.2f<p_{T}<%0.2f (GeV/c)",ptbins[pbin],ptbins[pbin+1]));
        txv1HFoddeven_etaCent->Draw();
    }

    cv1HFoddeven_etaCent->cd(1);
    legv1HFoddeven_etaCent_0->Draw();

    cv1HFoddeven_etaCent->cd(2);
    TLegend * legv1HFoddeven_etaCent_1 = new TLegend(0.07, 0.74, 0.28, 0.94);
    SetLegend(legv1HFoddeven_etaCent_1, 18);
    legv1HFoddeven_etaCent_1->AddEntry(v1odd_eta[anal][0][0],"v_{1}^{odd}","p");
    legv1HFoddeven_etaCent_1->AddEntry(v1even_eta[anal][0][0],"v_{1}^{even}","p");
    legv1HFoddeven_etaCent_1->Draw();

    cv1HFoddeven_etaCent->Print(Form("plots/diffv1/diffv1_eta/diff%s/v1_oddeven_eta_CentScan_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFoddeven_etaCent->Close();



    // differential v1(eta) using tracker+/- for each centrality bin
    anal = 15;
    if (!fopen(Form("plots/diffv1/diffv1_eta/diff%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/diffv1/diffv1_eta/diff%s",AnalNames[anal].data()));

    TCanvas * cv1Trackpm_eta[ncentbins];
    TH1D * hv1Trackpm_eta_tmp = new TH1D("hv1Trackpm_eta_tmp", "", 40, -2.5, 2.5);
    hv1Trackpm_eta_tmp->SetTitle("");
    hv1Trackpm_eta_tmp->SetStats(0);
    hv1Trackpm_eta_tmp->SetXTitle("#eta");
    hv1Trackpm_eta_tmp->SetYTitle("v_{1}");
    hv1Trackpm_eta_tmp->GetYaxis()->SetRangeUser(-0.1, 0.2);
    hv1Trackpm_eta_tmp->SetNdivisions(509);
    for (int cbin = 0; cbin<ncentbins; cbin++) {

        cv1Trackpm_eta[cbin] = new TCanvas(Form("cv1Trackpm_eta_cent%d-%d",centBins[cbin],centBins[cbin+1]),"cv1Trackpm_eta",1000,800);
        cv1Trackpm_eta[cbin]->Divide(4,3,0,0);

        for (int pbin = 0; pbin<netabins; pbin++) {
            TPad * padv1Trackpm_eta = (TPad *) cv1Trackpm_eta[cbin]->cd(pbin+1);
            if (gridlines) padv1Trackpm_eta->SetGrid();
            if (pbin == 3 || pbin == 7 || pbin == 11) padv1Trackpm_eta->SetRightMargin(0.02);
            TH1D * hv1Trackpm_eta = (TH1D *) hv1Trackpm_eta_tmp->Clone(Form("hv1Trackpm_eta_%c_%d",cbin,pbin));
            // hv1Trackpm_eta->GetYaxis()->SetRangeUser(-0.05-0.03, 0.15+0.03*cbin);
            if (pbin == 0 || pbin == 4) {
                hv1Trackpm_eta->GetYaxis()->CenterTitle();
                hv1Trackpm_eta->GetYaxis()->SetTitleSize(0.07);
                hv1Trackpm_eta->GetYaxis()->SetTitleOffset(1.34);
                hv1Trackpm_eta->GetYaxis()->SetLabelSize(0.06);
            } else if (pbin == 8) {
                hv1Trackpm_eta->GetXaxis()->SetTitleSize(0.06);
                hv1Trackpm_eta->GetXaxis()->SetTitleOffset(1.14);
                hv1Trackpm_eta->GetYaxis()->CenterTitle();
                hv1Trackpm_eta->GetYaxis()->SetTitleSize(0.06);
                hv1Trackpm_eta->GetYaxis()->SetTitleOffset(1.48);
                hv1Trackpm_eta->GetYaxis()->SetLabelSize(0.05);
            } else if (pbin>=9) {
                hv1Trackpm_eta->GetXaxis()->SetTitleSize(0.07);
                hv1Trackpm_eta->GetXaxis()->SetTitleOffset(1.00);
                hv1Trackpm_eta->GetXaxis()->SetLabelSize(0.06);
                hv1Trackpm_eta->GetXaxis()->SetLabelOffset(0.005);
            }
            hv1Trackpm_eta->Draw();
            v1p_eta[anal][cbin][pbin]->Draw("same");
            v1m_eta[anal][cbin][pbin]->Draw("same");

            TPaveText * txv1Trackpm_eta;
            if (pbin == 0 || pbin == 4) txv1Trackpm_eta = new TPaveText(0.23, 0.05, 0.62, 0.15,"NDC");
            else if (pbin == 8) txv1Trackpm_eta = new TPaveText(0.23, 0.19, 0.62, 0.29,"NDC");
            else if (pbin>=9) txv1Trackpm_eta = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
            else txv1Trackpm_eta = new TPaveText(0.05, 0.05, 0.50, 0.15,"NDC");
            SetTPaveTxt(txv1Trackpm_eta, 18);
            txv1Trackpm_eta->AddText(Form("%0.2f<p_{T}<%0.2f (GeV/c)",ptbins[pbin],ptbins[pbin+1]));
            txv1Trackpm_eta->Draw();
        }
        cv1Trackpm_eta[cbin]->cd(1);
        TPaveText * txv1Trackpm_eta_cent = new TPaveText(0.24, 0.84, 0.38, 0.95,"NDC");
        SetTPaveTxt(txv1Trackpm_eta_cent, 18);
        txv1Trackpm_eta_cent->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1Trackpm_eta_cent->Draw();

        cv1Trackpm_eta[cbin]->cd(2);
        TLegend * legv1Trackpm_eta = new TLegend(0.07, 0.74, 0.28, 0.94);
        SetLegend(legv1Trackpm_eta, 18);
        legv1Trackpm_eta->AddEntry(v1p_eta[anal][cbin][0],"Track+","p");
        legv1Trackpm_eta->AddEntry(v1m_eta[anal][cbin][0],"Track-","p");
        legv1Trackpm_eta->Draw();

        cv1Trackpm_eta[cbin]->Print(Form("plots/diffv1/diffv1_eta/diff%s/v1_pm_eta_%s_cent%d-%d.png",AnalNames[anal].data(),AnalNames[anal].data(),centBins[cbin],centBins[cbin+1]),"png");
        if (close_plots) cv1Trackpm_eta[cbin]->Close();
    }


}
