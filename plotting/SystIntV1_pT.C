# include "TFile.h"
# include "TGraphErrors.h"
# include "TCanvas.h"
# include "TF1.h"
# include "TH1D.h"
# include "TH2D.h"
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

TH1D * v1p_pt[nanals][ncentbins];
TH1D * v1m_pt[nanals][ncentbins];
TH1D * v1odd_pt[nanals][ncentbins];
TH1D * v1even_pt[nanals][ncentbins];

TH1D * v112p_pt[nanals][ncentbins];
TH1D * v112m_pt[nanals][ncentbins];
TH1D * v112odd_pt[nanals][ncentbins];
TH1D * v112even_pt[nanals][ncentbins];

TH1D * v123p_pt[nanals][ncentbins];
TH1D * v123m_pt[nanals][ncentbins];
TH1D * v123odd_pt[nanals][ncentbins];
TH1D * v123even_pt[nanals][ncentbins];

TH1D * absDiffv1_pt[nanals][ncentbins];

TH1D * runParms[nanals];

void SystIntV1_pT()
{

    TH1::SetDefaultSumw2();

    TFile * tfin = new TFile("../outputs/final_outputs/v1Int.root");

    //-- retrieve histograms from final output file
    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            v1p_pt[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%d-%d/v1p_pt_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v1m_pt[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%d-%d/v1m_pt_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v1odd_pt[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%d-%d/v1odd_pt_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v1even_pt[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%d-%d/v1even_pt_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));

            v112p_pt[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%d-%d/v112p_pt_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v112m_pt[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%d-%d/v112m_pt_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v112odd_pt[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%d-%d/v112odd_pt_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v112even_pt[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%d-%d/v112even_pt_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));

            v123p_pt[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%d-%d/v123p_pt_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v123m_pt[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%d-%d/v123m_pt_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v123odd_pt[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%d-%d/v123odd_pt_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v123even_pt[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%d-%d/v123even_pt_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));

            absDiffv1_pt[i][cbin] = new TH1D(Form("absDiffv1_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);
        }
    }

    // calculate ratios
    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            for (int j = 0; j<nptbins; j++) {
                double v1_pos = v1p_pt[i][cbin]->GetBinContent(j+1);
                double v1_neg = v1m_pt[i][cbin]->GetBinContent(j+1);
                double v1_pos_err = v1p_pt[i][cbin]->GetBinError(j+1);
                double v1_neg_err = v1m_pt[i][cbin]->GetBinError(j+1);

                double adv1 = fabs(v1_pos) - fabs(v1_neg);
                double adv1_err = sqrt( pow(v1_pos_err,2) + pow(v1_neg_err,2) );
                absDiffv1_pt[i][cbin]->SetBinContent(j+1, adv1);
                absDiffv1_pt[i][cbin]->SetBinError(j+1, adv1_err);
            }
        }
    }

    //-- plotting options

    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            absDiffv1_pt[i][cbin]->SetMarkerColor(kBlack);
            absDiffv1_pt[i][cbin]->SetLineColor(kBlack);
            absDiffv1_pt[i][cbin]->SetMarkerStyle(21);
            absDiffv1_pt[i][cbin]->SetMarkerSize(1.1);
        }
    }


    //-- make plots
    if (!fopen("plots","r")) system("mkdir plots");
    if (!fopen("plots/intv1","r")) system("mkdir plots/intv1");
    if (!fopen("plots/intv1/intv1_pt","r")) system("mkdir plots/intv1/intv1_pt");

    int anal; // choice of analysis

    TLine * lnetaRatio = new TLine(0.0, 1.0, 12, 1.0);
    lnetaRatio->SetLineWidth(1);


    // v1odd(+eta)/v1odd(-eta) for the HF
    anal = 7;
    if (!fopen(Form("plots/intv1/intv1_pt/int%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/intv1/intv1_pt/int%s",AnalNames[anal].data()));

    TCanvas * cv1HFoddAbsDiff_pt = new TCanvas("cv1HFoddAbsDiff_pt","cv1HFoddAbsDiff_pt",1100,620);
    TH1D * hv1HFoddAbsDiff_pt_tmp = new TH1D("hv1HFoddAbsDiff_pt_tmp", "", 100, 0, 12);
    hv1HFoddAbsDiff_pt_tmp->SetTitle("");
    hv1HFoddAbsDiff_pt_tmp->SetStats(0);
    hv1HFoddAbsDiff_pt_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1HFoddAbsDiff_pt_tmp->SetYTitle("|v_{1}(+#eta)| - |v_{1}(-#eta)|");
    hv1HFoddAbsDiff_pt_tmp->GetXaxis()->SetRangeUser(0, 12);
    hv1HFoddAbsDiff_pt_tmp->GetYaxis()->SetRangeUser(-0.05, 0.05);
    hv1HFoddAbsDiff_pt_tmp->SetNdivisions(509);
    hv1HFoddAbsDiff_pt_tmp->GetXaxis()->CenterTitle();
    hv1HFoddAbsDiff_pt_tmp->GetYaxis()->CenterTitle();
    cv1HFoddAbsDiff_pt->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1HFoddAbsDiff_pt = (TPad *) cv1HFoddAbsDiff_pt->cd(cbin+1);
        if (gridlines) padv1HFoddAbsDiff_pt->SetGrid();
        if (cbin == 3 || cbin == 7) padv1HFoddAbsDiff_pt->SetRightMargin(0.02);
        if (cbin <= 3) padv1HFoddAbsDiff_pt->SetTopMargin(0.08);
        TH1D * hv1HFoddAbsDiff_pt = (TH1D *) hv1HFoddAbsDiff_pt_tmp->Clone(Form("hv1HFoddAbsDiff_pt_%c",cbin));
        if (cbin == 0) {
            hv1HFoddAbsDiff_pt->GetYaxis()->SetTitleSize(0.07);
            hv1HFoddAbsDiff_pt->GetYaxis()->SetTitleOffset(1.33);
            hv1HFoddAbsDiff_pt->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1HFoddAbsDiff_pt->GetXaxis()->SetTitleSize(0.06);
            hv1HFoddAbsDiff_pt->GetXaxis()->SetTitleOffset(1.12);
            hv1HFoddAbsDiff_pt->GetXaxis()->SetLabelSize(0.06);
            hv1HFoddAbsDiff_pt->GetXaxis()->SetLabelOffset(0.018);
            hv1HFoddAbsDiff_pt->GetYaxis()->SetTitleSize(0.06);
            hv1HFoddAbsDiff_pt->GetYaxis()->SetTitleOffset(1.50);
            hv1HFoddAbsDiff_pt->GetYaxis()->SetLabelSize(0.05);
            hv1HFoddAbsDiff_pt->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >= 5) {
            hv1HFoddAbsDiff_pt->GetXaxis()->SetTitleSize(0.07);
            hv1HFoddAbsDiff_pt->GetXaxis()->SetTitleOffset(1.00);
            hv1HFoddAbsDiff_pt->GetXaxis()->SetLabelSize(0.07);
            hv1HFoddAbsDiff_pt->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1HFoddAbsDiff_pt->Draw();
        lnetaRatio->Draw();
        absDiffv1_pt[anal][cbin]->Draw("same");

        TF1 * fit1 = new TF1("fit1", "pol0", 0, 12);
        fit1->SetLineColor(kBlue);
        absDiffv1_pt[anal][cbin]->Fit(fit1,"QR");
        double par0 = fit1->GetParameter(0);
        double par0E = fit1->GetParError(0);
        double par0Chi2 = fit1->GetChisquare();

        TPaveText * txtxv1HFoddAbsDiff_pt_fit;
        if (cbin == 0) txtxv1HFoddAbsDiff_pt_fit = new TPaveText(0.24, 0.07, 0.74, 0.27,"NDC");
        else if (cbin >= 1 && cbin <= 3) txtxv1HFoddAbsDiff_pt_fit = new TPaveText(0.08, 0.07, 0.56, 0.27,"NDC");
        else if (cbin == 4) txtxv1HFoddAbsDiff_pt_fit = new TPaveText(0.24, 0.21, 0.74, 0.39,"NDC");
        else txtxv1HFoddAbsDiff_pt_fit = new TPaveText(0.08, 0.21, 0.56, 0.39,"NDC");
        SetTPaveTxt(txtxv1HFoddAbsDiff_pt_fit, 16);
        txtxv1HFoddAbsDiff_pt_fit->AddText(Form("mean: %0.4f #pm %0.4f",par0,par0E));
        txtxv1HFoddAbsDiff_pt_fit->AddText(Form("#chi^{2}: %0.4f",par0Chi2));
        txtxv1HFoddAbsDiff_pt_fit->Draw();

        TPaveText * txv1HFoddAbsDiff_pt;
        if (cbin == 0) txv1HFoddAbsDiff_pt = new TPaveText(0.75, 0.78, 0.93, 0.87,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1HFoddAbsDiff_pt = new TPaveText(0.68, 0.78, 0.86, 0.87,"NDC");
        else if (cbin == 4) txv1HFoddAbsDiff_pt = new TPaveText(0.75, 0.86, 0.93, 0.95,"NDC");
        else txv1HFoddAbsDiff_pt = new TPaveText(0.68, 0.86, 0.86, 0.95,"NDC");
        SetTPaveTxt(txv1HFoddAbsDiff_pt, 18);
        txv1HFoddAbsDiff_pt->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1HFoddAbsDiff_pt->Draw();
    }
    cv1HFoddAbsDiff_pt->cd(1);
    TPaveText * txv1HFoddAbsDiff_pt_1 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1HFoddAbsDiff_pt_1, 18);
    txv1HFoddAbsDiff_pt_1->AddText("#bf{CMS} #it{Preliminary}");
    txv1HFoddAbsDiff_pt_1->Draw();

    cv1HFoddAbsDiff_pt->Print(Form("plots/intv1/intv1_pt/int%s/v1HFoddAbsDiff_pt_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFoddAbsDiff_pt->Close();



    // v1odd(+eta)/v1odd(-eta) for the Tracker mom-cons
    anal = 15;
    if (!fopen(Form("plots/intv1/intv1_pt/int%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/intv1/intv1_pt/int%s",AnalNames[anal].data()));

    TCanvas * cv1TrkevenRatio_pt = new TCanvas("cv1TrkevenRatio_pt","cv1TrkevenRatio_pt",1100,620);
    TH2D * hv1TrkevenRatio_pt_tmp = new TH2D("hv1TrkevenRatio_pt_tmp", "", 100, 0, 12, 100, -0.5, 2.5);
    hv1TrkevenRatio_pt_tmp->SetTitle("");
    hv1TrkevenRatio_pt_tmp->SetStats(0);
    hv1TrkevenRatio_pt_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1TrkevenRatio_pt_tmp->SetYTitle("|v_{1}(+#eta)| - |v_{1}(-#eta)|");
    hv1TrkevenRatio_pt_tmp->GetXaxis()->SetRangeUser(0, 12);
    hv1TrkevenRatio_pt_tmp->GetYaxis()->SetRangeUser(-0.1, 0.1);
    hv1TrkevenRatio_pt_tmp->SetNdivisions(509);
    hv1TrkevenRatio_pt_tmp->GetXaxis()->CenterTitle();
    hv1TrkevenRatio_pt_tmp->GetYaxis()->CenterTitle();
    cv1TrkevenRatio_pt->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1TrkevenRatio_pt = (TPad *) cv1TrkevenRatio_pt->cd(cbin+1);
        if (gridlines) padv1TrkevenRatio_pt->SetGrid();
        if (cbin == 3 || cbin == 7) padv1TrkevenRatio_pt->SetRightMargin(0.02);
        if (cbin <= 3) padv1TrkevenRatio_pt->SetTopMargin(0.08);
        TH1D * hv1TrkevenRatio_pt = (TH1D *) hv1TrkevenRatio_pt_tmp->Clone(Form("hv1TrkevenRatio_pt_%c",cbin));
        if (cbin == 0) {
            hv1TrkevenRatio_pt->GetYaxis()->SetTitleSize(0.07);
            hv1TrkevenRatio_pt->GetYaxis()->SetTitleOffset(1.33);
            hv1TrkevenRatio_pt->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1TrkevenRatio_pt->GetXaxis()->SetTitleSize(0.06);
            hv1TrkevenRatio_pt->GetXaxis()->SetTitleOffset(1.12);
            hv1TrkevenRatio_pt->GetXaxis()->SetLabelSize(0.06);
            hv1TrkevenRatio_pt->GetXaxis()->SetLabelOffset(0.018);
            hv1TrkevenRatio_pt->GetYaxis()->SetTitleSize(0.06);
            hv1TrkevenRatio_pt->GetYaxis()->SetTitleOffset(1.50);
            hv1TrkevenRatio_pt->GetYaxis()->SetLabelSize(0.05);
            hv1TrkevenRatio_pt->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1TrkevenRatio_pt->GetXaxis()->SetTitleSize(0.07);
            hv1TrkevenRatio_pt->GetXaxis()->SetTitleOffset(1.00);
            hv1TrkevenRatio_pt->GetXaxis()->SetLabelSize(0.07);
            hv1TrkevenRatio_pt->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1TrkevenRatio_pt->Draw();
        lnetaRatio->Draw();
        absDiffv1_pt[anal][cbin]->Draw("same");

        TF1 * fit1 = new TF1("fit1", "pol0", 0, 12);
        fit1->SetLineColor(kBlue);
        absDiffv1_pt[anal][cbin]->Fit(fit1,"QR");
        double par0 = fit1->GetParameter(0);
        double par0E = fit1->GetParError(0);
        double par0Chi2 = fit1->GetChisquare();

        TPaveText * txtxv1TrkevenRatio_pt_fit;
        if (cbin == 0) txtxv1TrkevenRatio_pt_fit = new TPaveText(0.24, 0.07, 0.74, 0.27,"NDC");
        else if (cbin >= 1 && cbin <= 3) txtxv1TrkevenRatio_pt_fit = new TPaveText(0.08, 0.07, 0.56, 0.27,"NDC");
        else if (cbin == 4) txtxv1TrkevenRatio_pt_fit = new TPaveText(0.24, 0.21, 0.74, 0.39,"NDC");
        else txtxv1TrkevenRatio_pt_fit = new TPaveText(0.08, 0.21, 0.56, 0.39,"NDC");
        SetTPaveTxt(txtxv1TrkevenRatio_pt_fit, 16);
        txtxv1TrkevenRatio_pt_fit->AddText(Form("mean: %0.4f #pm %0.4f",par0,par0E));
        txtxv1TrkevenRatio_pt_fit->AddText(Form("#chi^{2}: %0.4f",par0Chi2));
        txtxv1TrkevenRatio_pt_fit->Draw();

        TPaveText * txv1TrkevenRatio_pt;
        if (cbin == 0) txv1TrkevenRatio_pt = new TPaveText(0.75, 0.78, 0.93, 0.87,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1TrkevenRatio_pt = new TPaveText(0.68, 0.78, 0.86, 0.87,"NDC");
        else if (cbin == 4) txv1TrkevenRatio_pt = new TPaveText(0.75, 0.86, 0.93, 0.95,"NDC");
        else txv1TrkevenRatio_pt = new TPaveText(0.68, 0.86, 0.86, 0.95,"NDC");
        SetTPaveTxt(txv1TrkevenRatio_pt, 18);
        txv1TrkevenRatio_pt->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1TrkevenRatio_pt->Draw();
    }
    cv1TrkevenRatio_pt->cd(1);
    TPaveText * txv1TrkevenRatio_pt_1 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1TrkevenRatio_pt_1, 18);
    txv1TrkevenRatio_pt_1->AddText("#bf{CMS} #it{Preliminary}");
    txv1TrkevenRatio_pt_1->Draw();

    cv1TrkevenRatio_pt->Print(Form("plots/intv1/intv1_pt/int%s/v1TrkevenRatio_pt_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1TrkevenRatio_pt->Close();

}
