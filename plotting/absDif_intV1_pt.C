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

void absDif_intV1_pt()
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

    // calculate absolute differences
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
            v1p_pt[i][cbin]->SetMarkerColor(kRed);
            v1p_pt[i][cbin]->SetLineColor(kRed);
            v1p_pt[i][cbin]->SetMarkerStyle(21);
            v1p_pt[i][cbin]->SetMarkerSize(1.1);

            v1m_pt[i][cbin]->SetMarkerColor(kBlue);
            v1m_pt[i][cbin]->SetLineColor(kBlue);
            v1m_pt[i][cbin]->SetMarkerStyle(21);
            v1m_pt[i][cbin]->SetMarkerSize(1.1);

            v1odd_pt[i][cbin]->SetMarkerColor(kBlue);
            v1odd_pt[i][cbin]->SetLineColor(kBlue);
            v1odd_pt[i][cbin]->SetMarkerStyle(21);
            v1odd_pt[i][cbin]->SetMarkerSize(1.1);

            v1even_pt[i][cbin]->SetMarkerColor(kBlue);
            v1even_pt[i][cbin]->SetLineColor(kBlue);
            v1even_pt[i][cbin]->SetMarkerStyle(21);
            v1even_pt[i][cbin]->SetMarkerSize(1.1);

            absDiffv1_pt[i][cbin]->SetMarkerColor(kBlack);
            absDiffv1_pt[i][cbin]->SetLineColor(kBlack);
            absDiffv1_pt[i][cbin]->SetMarkerStyle(21);
            absDiffv1_pt[i][cbin]->SetMarkerSize(1.1);
        }
    }
    TF1 * fitleg = new TF1("fitleg","pol0",0,1); // dummy TF1 for legends
    fitleg->SetLineColor(kBlue);


    //-- make plots
    if (!fopen("plots","r")) system("mkdir plots");
    if (!fopen("plots/intv1","r")) system("mkdir plots/intv1");
    if (!fopen("plots/intv1/intv1_pt","r")) system("mkdir plots/intv1/intv1_pt");

    int anal; // choice of analysis


    // integrated v1(pT) using HF+/- for each centrality bin
    anal = 7;
    if (!fopen(Form("plots/intv1/intv1_pt/int%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/intv1/intv1_pt/int%s",AnalNames[anal].data()));

    TCanvas * cv1HFpm_pt = new TCanvas("cv1HFpm_pt","cv1HFpm_pt",1100,620);
    TH1D * hv1HFpm_pt_tmp = new TH1D("hv1HFpm_pt", "", 40, 0, 12);
    hv1HFpm_pt_tmp->SetTitle("");
    hv1HFpm_pt_tmp->SetStats(0);
    hv1HFpm_pt_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1HFpm_pt_tmp->SetYTitle("v_{1}");
    hv1HFpm_pt_tmp->GetYaxis()->SetRangeUser(-0.3, 0.3);
    hv1HFpm_pt_tmp->GetXaxis()->CenterTitle();
    hv1HFpm_pt_tmp->SetNdivisions(509);
    cv1HFpm_pt->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1HFpm_pt = (TPad *) cv1HFpm_pt->cd(cbin+1);
        if (gridlines) padv1HFpm_pt->SetGrid();
        if (cbin == 3 || cbin == 7) padv1HFpm_pt->SetRightMargin(0.02);
        if (cbin <= 3) padv1HFpm_pt->SetTopMargin(0.08);
        TH1D * hv1HFpm_pt = (TH1D *) hv1HFpm_pt_tmp->Clone(Form("hv1HFpm_pt_%c",cbin));
        if (cbin == 0) {
            hv1HFpm_pt->GetYaxis()->SetTitleSize(0.07);
            hv1HFpm_pt->GetYaxis()->SetTitleOffset(1.33);
            hv1HFpm_pt->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1HFpm_pt->GetXaxis()->SetTitleSize(0.06);
            hv1HFpm_pt->GetXaxis()->SetTitleOffset(1.12);
            hv1HFpm_pt->GetXaxis()->SetLabelSize(0.06);
            hv1HFpm_pt->GetXaxis()->SetLabelOffset(0.018);
            hv1HFpm_pt->GetYaxis()->SetTitleSize(0.06);
            hv1HFpm_pt->GetYaxis()->SetTitleOffset(1.50);
            hv1HFpm_pt->GetYaxis()->SetLabelSize(0.05);
            hv1HFpm_pt->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1HFpm_pt->GetXaxis()->SetTitleSize(0.07);
            hv1HFpm_pt->GetXaxis()->SetTitleOffset(1.00);
            hv1HFpm_pt->GetXaxis()->SetLabelSize(0.07);
            hv1HFpm_pt->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1HFpm_pt->Draw();
        v1p_pt[anal][cbin]->Draw("same");
        v1m_pt[anal][cbin]->Draw("same");

        TPaveText * txv1HFpm_pt;
        if (cbin == 0) txv1HFpm_pt = new TPaveText(0.25, 0.04, 0.46, 0.17,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1HFpm_pt = new TPaveText(0.06, 0.04, 0.26, 0.17,"NDC");
        else if (cbin == 4) txv1HFpm_pt = new TPaveText(0.25, 0.21, 0.46, 0.30,"NDC");
        else txv1HFpm_pt = new TPaveText(0.06, 0.21, 0.26, 0.30,"NDC");
        SetTPaveTxt(txv1HFpm_pt, 18);
        txv1HFpm_pt->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1HFpm_pt->Draw();
    }
    cv1HFpm_pt->cd(1);
    TPaveText * txv1HFpm_pt_1 = new TPaveText(0.22, 0.65, 0.81, 0.86,"NDC");
    SetTPaveTxt(txv1HFpm_pt_1, 18);
    txv1HFpm_pt_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1HFpm_pt_1->AddText("#pm#eta < #pm2.4");
    txv1HFpm_pt_1->Draw();

    TPaveText * txv1HFpm_pt_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1HFpm_pt_2, 18);
    txv1HFpm_pt_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1HFpm_pt_2->Draw();

    cv1HFpm_pt->cd(2);
    TLegend * legv1HFpm_pt = new TLegend(0.06, 0.67, 0.25, 0.86);
    SetLegend(legv1HFpm_pt, 18);
    legv1HFpm_pt->AddEntry(v1p_pt[anal][0]," HF+","p");
    legv1HFpm_pt->AddEntry(v1m_pt[anal][0]," HF-","p");
    legv1HFpm_pt->Draw();

    cv1HFpm_pt->Print(Form("plots/intv1/intv1_pt/int%s/v1_pm_pt_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFpm_pt->Close();



    // integrated v1^odd(pT) using HF for each centrality bin

    TCanvas * cv1HFodd_pt = new TCanvas("cv1HFodd_pt","cv1HFodd_pt",1100,620);
    TH1D * hv1HFodd_pt_tmp = new TH1D("hv1HFodd_pt", "", 40, 0, 12);
    hv1HFodd_pt_tmp->SetTitle("");
    hv1HFodd_pt_tmp->SetStats(0);
    hv1HFodd_pt_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1HFodd_pt_tmp->SetYTitle("v_{1}^{odd}");
    hv1HFodd_pt_tmp->GetYaxis()->SetRangeUser(-0.04, 0.04);
    hv1HFodd_pt_tmp->GetYaxis()->CenterTitle();
    hv1HFodd_pt_tmp->SetNdivisions(509);
    cv1HFodd_pt->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1HFodd_pt = (TPad *) cv1HFodd_pt->cd(cbin+1);
        if (gridlines) padv1HFodd_pt->SetGrid();
        if (cbin == 3 || cbin == 7) padv1HFodd_pt->SetRightMargin(0.02);
        if (cbin <= 3) padv1HFodd_pt->SetTopMargin(0.08);
        TH1D * hv1HFodd_pt = (TH1D *) hv1HFodd_pt_tmp->Clone(Form("hv1HFodd_pt_%c",cbin));
        if (cbin == 0) {
            hv1HFodd_pt->GetYaxis()->SetTitleSize(0.07);
            hv1HFodd_pt->GetYaxis()->SetTitleOffset(1.33);
            hv1HFodd_pt->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1HFodd_pt->GetXaxis()->SetTitleSize(0.06);
            hv1HFodd_pt->GetXaxis()->SetTitleOffset(1.12);
            hv1HFodd_pt->GetXaxis()->SetLabelSize(0.06);
            hv1HFodd_pt->GetXaxis()->SetLabelOffset(0.018);
            hv1HFodd_pt->GetYaxis()->SetTitleSize(0.06);
            hv1HFodd_pt->GetYaxis()->SetTitleOffset(1.50);
            hv1HFodd_pt->GetYaxis()->SetLabelSize(0.05);
            hv1HFodd_pt->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1HFodd_pt->GetXaxis()->SetTitleSize(0.07);
            hv1HFodd_pt->GetXaxis()->SetTitleOffset(1.00);
            hv1HFodd_pt->GetXaxis()->SetLabelSize(0.07);
            hv1HFodd_pt->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1HFodd_pt->Draw();
        v1odd_pt[anal][cbin]->Draw("same");

        TPaveText * txv1HFodd_pt;
        if (cbin == 0) txv1HFodd_pt = new TPaveText(0.25, 0.04, 0.46, 0.17,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1HFodd_pt = new TPaveText(0.06, 0.04, 0.26, 0.17,"NDC");
        else if (cbin == 4) txv1HFodd_pt = new TPaveText(0.25, 0.21, 0.46, 0.30,"NDC");
        else txv1HFodd_pt = new TPaveText(0.06, 0.21, 0.26, 0.30,"NDC");
        SetTPaveTxt(txv1HFodd_pt, 18);
        txv1HFodd_pt->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1HFodd_pt->Draw();
    }
    cv1HFodd_pt->cd(1);
    TPaveText * txv1HFodd_pt_1 = new TPaveText(0.22, 0.65, 0.81, 0.86,"NDC");
    SetTPaveTxt(txv1HFodd_pt_1, 18);
    txv1HFodd_pt_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1HFodd_pt_1->AddText("|#eta| < 2.4");
    txv1HFodd_pt_1->Draw();

    TPaveText * txv1HFodd_eta_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1HFodd_eta_2, 18);
    txv1HFodd_eta_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1HFodd_eta_2->Draw();

    cv1HFodd_pt->Print(Form("plots/intv1/intv1_pt/int%s/v1_odd_pt_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFodd_pt->Close();



    // integrated v1^even(pT) using HF for each centrality bin

    TCanvas * cv1HFeven_pt = new TCanvas("cv1HFeven_pt","cv1HFeven_pt",1100,620);
    TH1D * hv1HFeven_pt_tmp = new TH1D("hv1HFeven_pt", "", 40, 0, 12);
    hv1HFeven_pt_tmp->SetTitle("");
    hv1HFeven_pt_tmp->SetStats(0);
    hv1HFeven_pt_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1HFeven_pt_tmp->SetYTitle("v_{1}^{even}");
    hv1HFeven_pt_tmp->GetYaxis()->SetRangeUser(-0.25, 0.0);
    hv1HFeven_pt_tmp->GetYaxis()->CenterTitle();
    hv1HFeven_pt_tmp->SetNdivisions(509);
    cv1HFeven_pt->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1HFeven_pt = (TPad *) cv1HFeven_pt->cd(cbin+1);
        if (gridlines) padv1HFeven_pt->SetGrid();
        if (cbin == 3 || cbin == 7) padv1HFeven_pt->SetRightMargin(0.02);
        if (cbin <= 3) padv1HFeven_pt->SetTopMargin(0.08);
        TH1D * hv1HFeven_pt = (TH1D *) hv1HFeven_pt_tmp->Clone(Form("hv1HFeven_pt_%c",cbin));
        if (cbin == 0) {
            hv1HFeven_pt->GetYaxis()->SetTitleSize(0.07);
            hv1HFeven_pt->GetYaxis()->SetTitleOffset(1.33);
            hv1HFeven_pt->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1HFeven_pt->GetXaxis()->SetTitleSize(0.06);
            hv1HFeven_pt->GetXaxis()->SetTitleOffset(1.12);
            hv1HFeven_pt->GetXaxis()->SetLabelSize(0.06);
            hv1HFeven_pt->GetXaxis()->SetLabelOffset(0.018);
            hv1HFeven_pt->GetYaxis()->SetTitleSize(0.06);
            hv1HFeven_pt->GetYaxis()->SetTitleOffset(1.50);
            hv1HFeven_pt->GetYaxis()->SetLabelSize(0.05);
            hv1HFeven_pt->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1HFeven_pt->GetXaxis()->SetTitleSize(0.07);
            hv1HFeven_pt->GetXaxis()->SetTitleOffset(1.00);
            hv1HFeven_pt->GetXaxis()->SetLabelSize(0.07);
            hv1HFeven_pt->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1HFeven_pt->Draw();
        v1even_pt[anal][cbin]->Draw("same");

        TPaveText * txv1HFeven_pt;
        if (cbin == 0) txv1HFeven_pt = new TPaveText(0.25, 0.04, 0.46, 0.17,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1HFeven_pt = new TPaveText(0.06, 0.04, 0.26, 0.17,"NDC");
        else if (cbin == 4) txv1HFeven_pt = new TPaveText(0.25, 0.21, 0.46, 0.30,"NDC");
        else txv1HFeven_pt = new TPaveText(0.06, 0.21, 0.26, 0.30,"NDC");
        SetTPaveTxt(txv1HFeven_pt, 18);
        txv1HFeven_pt->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1HFeven_pt->Draw();
    }
    cv1HFeven_pt->cd(1);
    TPaveText * txv1HFeven_pt_1 = new TPaveText(0.24, 0.17, 0.69, 0.33,"NDC");
    SetTPaveTxt(txv1HFeven_pt_1, 18);
    txv1HFeven_pt_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1HFeven_pt_1->AddText("|#eta| < 2.4");
    txv1HFeven_pt_1->Draw();

    TPaveText * txv1HFeven_pt_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1HFeven_pt_2, 18);
    txv1HFeven_pt_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1HFeven_pt_2->Draw();

    cv1HFeven_pt->Print(Form("plots/intv1/intv1_pt/int%s/v1_even_pt_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFeven_pt->Close();



    // integrated v1(pT) using Track+/- for each centrality bin
    anal = 15;
    if (!fopen(Form("plots/intv1/intv1_pt/int%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/intv1/intv1_pt/int%s",AnalNames[anal].data()));

    TCanvas * cv1Trkpm_pt = new TCanvas("cv1Trkpm_pt","cv1Trkpm_pt",1100,620);
    TH1D * hv1Trkpm_pt_tmp = new TH1D("hv1Trkpm_pt", "", 40, 0, 12);
    hv1Trkpm_pt_tmp->SetTitle("");
    hv1Trkpm_pt_tmp->SetStats(0);
    hv1Trkpm_pt_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1Trkpm_pt_tmp->SetYTitle("v_{1}");
    hv1Trkpm_pt_tmp->GetYaxis()->SetRangeUser(-0.06, 0.26);
    hv1Trkpm_pt_tmp->GetYaxis()->CenterTitle();
    hv1Trkpm_pt_tmp->SetNdivisions(509);
    cv1Trkpm_pt->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1Trkpm_pt = (TPad *) cv1Trkpm_pt->cd(cbin+1);
        if (gridlines) padv1Trkpm_pt->SetGrid();
        if (cbin == 3 || cbin == 7) padv1Trkpm_pt->SetRightMargin(0.02);
        if (cbin <= 3) padv1Trkpm_pt->SetTopMargin(0.08);
        TH1D * hv1Trkpm_pt = (TH1D *) hv1Trkpm_pt_tmp->Clone(Form("hv1Trkpm_pt_%c",cbin));
        if (cbin == 0) {
            hv1Trkpm_pt->GetYaxis()->SetTitleSize(0.07);
            hv1Trkpm_pt->GetYaxis()->SetTitleOffset(1.33);
            hv1Trkpm_pt->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1Trkpm_pt->GetXaxis()->SetTitleSize(0.06);
            hv1Trkpm_pt->GetXaxis()->SetTitleOffset(1.12);
            hv1Trkpm_pt->GetXaxis()->SetLabelSize(0.06);
            hv1Trkpm_pt->GetXaxis()->SetLabelOffset(0.018);
            hv1Trkpm_pt->GetYaxis()->SetTitleSize(0.06);
            hv1Trkpm_pt->GetYaxis()->SetTitleOffset(1.50);
            hv1Trkpm_pt->GetYaxis()->SetLabelSize(0.05);
            hv1Trkpm_pt->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1Trkpm_pt->GetXaxis()->SetTitleSize(0.07);
            hv1Trkpm_pt->GetXaxis()->SetTitleOffset(1.00);
            hv1Trkpm_pt->GetXaxis()->SetLabelSize(0.07);
            hv1Trkpm_pt->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1Trkpm_pt->Draw();
        v1p_pt[anal][cbin]->Draw("same");
        v1m_pt[anal][cbin]->Draw("same");

        TPaveText * txv1Trkpm_pt;
        if (cbin == 0) txv1Trkpm_pt = new TPaveText(0.75, 0.05, 0.93, 0.14,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1Trkpm_pt = new TPaveText(0.69, 0.05, 0.87, 0.14,"NDC");
        else if (cbin == 4) txv1Trkpm_pt = new TPaveText(0.75, 0.20, 0.93, 0.29,"NDC");
        else txv1Trkpm_pt = new TPaveText(0.69, 0.20, 0.87, 0.29,"NDC");
        SetTPaveTxt(txv1Trkpm_pt, 18);
        txv1Trkpm_pt->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1Trkpm_pt->Draw();
    }
    cv1Trkpm_pt->cd(1);
    TPaveText * txv1Trkpm_pt_1 = new TPaveText(0.22, 0.65, 0.81, 0.86,"NDC");
    SetTPaveTxt(txv1Trkpm_pt_1, 18);
    txv1Trkpm_pt_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1Trkpm_pt_1->AddText("|#eta| < 2.4");
    txv1Trkpm_pt_1->Draw();

    TPaveText * txv1Trkpm_pt_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1Trkpm_pt_2, 18);
    txv1Trkpm_pt_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1Trkpm_pt_2->Draw();

    cv1Trkpm_pt->cd(2);
    TLegend * legv1Trkpm_pt = new TLegend(0.06, 0.67, 0.25, 0.86);
    SetLegend(legv1Trkpm_pt, 18);
    legv1Trkpm_pt->AddEntry(v1p_pt[anal][0],"Tracker+","p");
    legv1Trkpm_pt->AddEntry(v1m_pt[anal][0],"Tracker-","p");
    legv1Trkpm_pt->Draw();

    cv1Trkpm_pt->Print(Form("plots/intv1/intv1_pt/int%s/v1_pm_pt_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1Trkpm_pt->Close();



    // integrated v1^even(pT) using the tracker (called v1odd in output file)

    TCanvas * cv1Trkeven_pt = new TCanvas("cv1Trkeven_pt","cv1Trkeven_pt",1100,620);
    TH1D * hv1Trkeven_pt_tmp = new TH1D("hv1Trkeven_pt", "", 40, 0, 12);
    hv1Trkeven_pt_tmp->SetTitle("");
    hv1Trkeven_pt_tmp->SetStats(0);
    hv1Trkeven_pt_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1Trkeven_pt_tmp->SetYTitle("v_{1}^{even}");
    hv1Trkeven_pt_tmp->GetYaxis()->SetRangeUser(-0.06, 0.26);
    hv1Trkeven_pt_tmp->GetYaxis()->CenterTitle();
    hv1Trkeven_pt_tmp->SetNdivisions(509);
    cv1Trkeven_pt->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1Trkeven_pt = (TPad *) cv1Trkeven_pt->cd(cbin+1);
        if (gridlines) padv1Trkeven_pt->SetGrid();
        if (cbin == 3 || cbin == 7) padv1Trkeven_pt->SetRightMargin(0.02);
        if (cbin <= 3) padv1Trkeven_pt->SetTopMargin(0.08);
        TH1D * hv1Trkeven_pt = (TH1D *) hv1Trkeven_pt_tmp->Clone(Form("hv1Trkeven_pt_%c",cbin));
        if (cbin == 0) {
            hv1Trkeven_pt->GetYaxis()->SetTitleSize(0.07);
            hv1Trkeven_pt->GetYaxis()->SetTitleOffset(1.33);
            hv1Trkeven_pt->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin == 4) {
            hv1Trkeven_pt->GetXaxis()->SetTitleSize(0.06);
            hv1Trkeven_pt->GetXaxis()->SetTitleOffset(1.12);
            hv1Trkeven_pt->GetXaxis()->SetLabelSize(0.06);
            hv1Trkeven_pt->GetXaxis()->SetLabelOffset(0.018);
            hv1Trkeven_pt->GetYaxis()->SetTitleSize(0.06);
            hv1Trkeven_pt->GetYaxis()->SetTitleOffset(1.50);
            hv1Trkeven_pt->GetYaxis()->SetLabelSize(0.05);
            hv1Trkeven_pt->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1Trkeven_pt->GetXaxis()->SetTitleSize(0.07);
            hv1Trkeven_pt->GetXaxis()->SetTitleOffset(1.00);
            hv1Trkeven_pt->GetXaxis()->SetLabelSize(0.07);
            hv1Trkeven_pt->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1Trkeven_pt->Draw();
        v1odd_pt[anal][cbin]->Draw("same");

        TPaveText * txv1Trkeven_pt;
        if (cbin == 0) txv1Trkeven_pt = new TPaveText(0.75, 0.05, 0.93, 0.14,"NDC");
        else if (cbin >= 1 && cbin <= 3) txv1Trkeven_pt = new TPaveText(0.69, 0.05, 0.87, 0.14,"NDC");
        else if (cbin == 4) txv1Trkeven_pt = new TPaveText(0.75, 0.20, 0.93, 0.29,"NDC");
        else txv1Trkeven_pt = new TPaveText(0.69, 0.20, 0.87, 0.29,"NDC");
        SetTPaveTxt(txv1Trkeven_pt, 18);
        txv1Trkeven_pt->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1Trkeven_pt->Draw();
    }
    cv1Trkeven_pt->cd(1);
    TPaveText * txv1Trkeven_pt_1 = new TPaveText(0.22, 0.65, 0.81, 0.86,"NDC");
    SetTPaveTxt(txv1Trkeven_pt_1, 18);
    txv1Trkeven_pt_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1Trkeven_pt_1->AddText("|#eta| < 2.4");
    txv1Trkeven_pt_1->Draw();

    TPaveText * txv1Trkeven_pt_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1Trkeven_pt_2, 18);
    txv1Trkeven_pt_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1Trkeven_pt_2->Draw();

    cv1Trkeven_pt->Print(Form("plots/intv1/intv1_pt/int%s/v1_even_pt_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1Trkeven_pt->Close();



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
    TH1D * hv1TrkevenRatio_pt_tmp = new TH1D("hv1TrkevenRatio_pt_tmp", "", 100, 0, 12);
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



    //-- plot both v1odd(pT) and it's absolute differences
    // 0 - 40% centrality
    anal = 7;
    TCanvas * cv1HFodd_and_absdiff_pt_0to40 = new TCanvas("cv1HFodd_and_absdiff_pt_0to40","cv1HFodd_and_absdiff_pt_0to40",1100,700);
    TH1D * hv1HFodd_and_absdiff_pt_0to40_tmp = new TH1D("hv1HFodd_and_absdiff_pt_0to40", "", 40, 0, 12);
    hv1HFodd_and_absdiff_pt_0to40_tmp->SetTitle("");
    hv1HFodd_and_absdiff_pt_0to40_tmp->SetStats(0);
    hv1HFodd_and_absdiff_pt_0to40_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1HFodd_and_absdiff_pt_0to40_tmp->SetYTitle("v_{1}^{odd}");
    hv1HFodd_and_absdiff_pt_0to40_tmp->GetYaxis()->SetRangeUser(-0.02, 0.02);
    hv1HFodd_and_absdiff_pt_0to40_tmp->GetYaxis()->CenterTitle();
    hv1HFodd_and_absdiff_pt_0to40_tmp->SetNdivisions(509);
    cv1HFodd_and_absdiff_pt_0to40->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1HFodd_and_absdiff_pt_0to40 = (TPad *) cv1HFodd_and_absdiff_pt_0to40->cd(cbin+1);
        if (gridlines) padv1HFodd_and_absdiff_pt_0to40->SetGrid();
        if (cbin == 3 || cbin == 7) padv1HFodd_and_absdiff_pt_0to40->SetRightMargin(0.02);
        if (cbin <= 3) {
            padv1HFodd_and_absdiff_pt_0to40->SetTopMargin(0.08);
            padv1HFodd_and_absdiff_pt_0to40->SetBottomMargin(0.12);
        }
        if (cbin >= 4) padv1HFodd_and_absdiff_pt_0to40->SetTopMargin(0.14);
        TH1D * hv1HFodd_and_absdiff_pt_0to40 = (TH1D *) hv1HFodd_and_absdiff_pt_0to40_tmp->Clone(Form("hv1HFodd_and_absdiff_pt_0to40_%c",cbin));
        if (cbin == 0) {
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetTitleSize(0.06);
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetTitleOffset(1.12);
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetLabelSize(0.06);
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetLabelOffset(0.018);
            hv1HFodd_and_absdiff_pt_0to40->GetYaxis()->SetTitleSize(0.07);
            hv1HFodd_and_absdiff_pt_0to40->GetYaxis()->SetTitleOffset(1.33);
            hv1HFodd_and_absdiff_pt_0to40->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin >= 1 && cbin <= 3) {
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetTitleSize(0.07);
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetTitleOffset(1.00);
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetLabelSize(0.07);
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetLabelOffset(0.008);
        }
        if (cbin >= 4) {
            hv1HFodd_and_absdiff_pt_0to40->SetYTitle("");
            hv1HFodd_and_absdiff_pt_0to40->GetYaxis()->SetRangeUser(-0.03, 0.03);
        }
        if (cbin == 4) {
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetTitleSize(0.06);
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetTitleOffset(1.12);
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetLabelSize(0.06);
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetLabelOffset(0.018);
            hv1HFodd_and_absdiff_pt_0to40->GetYaxis()->SetTitleSize(0.06);
            hv1HFodd_and_absdiff_pt_0to40->GetYaxis()->SetTitleOffset(1.50);
            hv1HFodd_and_absdiff_pt_0to40->GetYaxis()->SetLabelSize(0.05);
            hv1HFodd_and_absdiff_pt_0to40->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetTitleSize(0.07);
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetTitleOffset(1.00);
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetLabelSize(0.07);
            hv1HFodd_and_absdiff_pt_0to40->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1HFodd_and_absdiff_pt_0to40->Draw();

        TPaveText * txv1HFodd_and_absdiff_pt_0to40;
        if (cbin <= 3) {
            // v1p_pt[anal][cbin]->Draw("same");
            // v1m_pt[anal][cbin]->Draw("same");
            v1odd_pt[anal][cbin]->Draw("same");
            if (cbin == 0) txv1HFodd_and_absdiff_pt_0to40 = new TPaveText(0.24, 0.17, 0.46, 0.24,"NDC");
            else if (cbin >= 1 && cbin <= 3) txv1HFodd_and_absdiff_pt_0to40 = new TPaveText(0.05, 0.17, 0.25, 0.24,"NDC");
            SetTPaveTxt(txv1HFodd_and_absdiff_pt_0to40, 18);
            txv1HFodd_and_absdiff_pt_0to40->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
            txv1HFodd_and_absdiff_pt_0to40->Draw();
        } else {
            absDiffv1_pt[anal][cbin-4]->Draw("same");

            TF1 * fit1 = new TF1("fit1", "pol0", 0, 12);
            fit1->SetLineColor(kBlue);
            absDiffv1_pt[anal][cbin-4]->Fit(fit1,"QR");
            double par0 = fit1->GetParameter(0);
            double par0E = fit1->GetParError(0);
            double par0Chi2 = fit1->GetChisquare();

            TPaveText * tx1HFodd_and_absdiff_pt_fit;
            if (cbin == 4) tx1HFodd_and_absdiff_pt_fit = new TPaveText(0.22, 0.21, 0.72, 0.35,"NDC");
            else tx1HFodd_and_absdiff_pt_fit = new TPaveText(0.06, 0.21, 0.54, 0.35,"NDC");
            SetTPaveTxt(tx1HFodd_and_absdiff_pt_fit, 16);
            tx1HFodd_and_absdiff_pt_fit->AddText(Form("mean: %0.4f #pm %0.4f",par0,par0E));
            tx1HFodd_and_absdiff_pt_fit->AddText(Form("#chi^{2}: %0.4f",par0Chi2));
            tx1HFodd_and_absdiff_pt_fit->Draw();
        }

    }
    cv1HFodd_and_absdiff_pt_0to40->cd(1);
    TPaveText * txv1HFodd_and_absdiff_pt_0to40_1 = new TPaveText(0.22, 0.68, 0.81, 0.87,"NDC");
    SetTPaveTxt(txv1HFodd_and_absdiff_pt_0to40_1, 18);
    txv1HFodd_and_absdiff_pt_0to40_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1HFodd_and_absdiff_pt_0to40_1->AddText("|#eta| < 2.4 (GeV/c)");
    txv1HFodd_and_absdiff_pt_0to40_1->Draw();

    TPaveText * txv1HFoddAbsDiff_pt_0to40_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1HFoddAbsDiff_pt_0to40_2, 18);
    txv1HFoddAbsDiff_pt_0to40_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1HFoddAbsDiff_pt_0to40_2->Draw();

    cv1HFodd_and_absdiff_pt_0to40->cd(5);
    TLegend * legv1HFoddAbsDiff_pt_0to40 = new TLegend(0.24, 0.67, 0.51, 0.84);
    SetLegend(legv1HFoddAbsDiff_pt_0to40, 18);
    legv1HFoddAbsDiff_pt_0to40->SetHeader("|v_{1}(+#eta)| - |v_{1}(-#eta)|");
    legv1HFoddAbsDiff_pt_0to40->AddEntry(fitleg," pol0 fit","l");
    legv1HFoddAbsDiff_pt_0to40->Draw();

    cv1HFodd_and_absdiff_pt_0to40->Print(Form("plots/intv1/intv1_pt/int%s/v1odd_with_absdiff_pt_0to40cent_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFodd_and_absdiff_pt_0to40->Close();



    // 40 - 80% centrality
    TCanvas * cv1HFodd_and_absdiff_pt_40to80 = new TCanvas("cv1HFodd_and_absdiff_pt_40to80","cv1HFodd_and_absdiff_pt_40to80",1100,700);
    TH1D * hv1HFodd_and_absdiff_pt_40to80_tmp = new TH1D("hv1HFodd_and_absdiff_pt_40to80", "", 40, 0, 12);
    hv1HFodd_and_absdiff_pt_40to80_tmp->SetTitle("");
    hv1HFodd_and_absdiff_pt_40to80_tmp->SetStats(0);
    hv1HFodd_and_absdiff_pt_40to80_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1HFodd_and_absdiff_pt_40to80_tmp->SetYTitle("v_{1}^{odd}");
    hv1HFodd_and_absdiff_pt_40to80_tmp->GetYaxis()->SetRangeUser(-0.04, 0.04);
    hv1HFodd_and_absdiff_pt_40to80_tmp->GetYaxis()->CenterTitle();
    hv1HFodd_and_absdiff_pt_40to80_tmp->SetNdivisions(509);
    cv1HFodd_and_absdiff_pt_40to80->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1HFodd_and_absdiff_pt_40to80 = (TPad *) cv1HFodd_and_absdiff_pt_40to80->cd(cbin+1);
        if (gridlines) padv1HFodd_and_absdiff_pt_40to80->SetGrid();
        if (cbin == 3 || cbin == 7) padv1HFodd_and_absdiff_pt_40to80->SetRightMargin(0.02);
        if (cbin <= 3) {
            padv1HFodd_and_absdiff_pt_40to80->SetTopMargin(0.08);
            padv1HFodd_and_absdiff_pt_40to80->SetBottomMargin(0.12);
        }
        if (cbin >= 4) padv1HFodd_and_absdiff_pt_40to80->SetTopMargin(0.14);
        TH1D * hv1HFodd_and_absdiff_pt_40to80 = (TH1D *) hv1HFodd_and_absdiff_pt_40to80_tmp->Clone(Form("hv1HFodd_and_absdiff_pt_40to80_%c",cbin));
        if (cbin == 0) {
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetTitleSize(0.06);
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetTitleOffset(1.12);
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetLabelSize(0.06);
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetLabelOffset(0.018);
            hv1HFodd_and_absdiff_pt_40to80->GetYaxis()->SetTitleSize(0.07);
            hv1HFodd_and_absdiff_pt_40to80->GetYaxis()->SetTitleOffset(1.33);
            hv1HFodd_and_absdiff_pt_40to80->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin >= 1 && cbin <= 3) {
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetTitleSize(0.07);
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetTitleOffset(1.00);
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetLabelSize(0.07);
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetLabelOffset(0.008);
        }
        if (cbin >= 4) {
            hv1HFodd_and_absdiff_pt_40to80->SetYTitle("");
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetRangeUser(0, 12);
            hv1HFodd_and_absdiff_pt_40to80->GetYaxis()->SetRangeUser(-0.06, 0.06);
        }
        if (cbin == 4) {
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetTitleSize(0.06);
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetTitleOffset(1.12);
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetLabelSize(0.06);
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetLabelOffset(0.018);
            hv1HFodd_and_absdiff_pt_40to80->GetYaxis()->SetTitleSize(0.06);
            hv1HFodd_and_absdiff_pt_40to80->GetYaxis()->SetTitleOffset(1.50);
            hv1HFodd_and_absdiff_pt_40to80->GetYaxis()->SetLabelSize(0.05);
            hv1HFodd_and_absdiff_pt_40to80->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->CenterTitle();
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetTitleSize(0.07);
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetTitleOffset(1.00);
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetLabelSize(0.07);
            hv1HFodd_and_absdiff_pt_40to80->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1HFodd_and_absdiff_pt_40to80->Draw();

        TPaveText * txv1HFodd_and_absdiff_pt_40to80;
        if (cbin <= 3) {
            // v1p_pt[anal][cbin+4]->Draw("same");
            // v1m_pt[anal][cbin+4]->Draw("same");
            v1odd_pt[anal][cbin+4]->Draw("same");
            if (cbin == 0) txv1HFodd_and_absdiff_pt_40to80 = new TPaveText(0.24, 0.17, 0.46, 0.24,"NDC");
            else if (cbin >= 1 && cbin <= 3) txv1HFodd_and_absdiff_pt_40to80 = new TPaveText(0.05, 0.17, 0.25, 0.24,"NDC");
            SetTPaveTxt(txv1HFodd_and_absdiff_pt_40to80, 18);
            txv1HFodd_and_absdiff_pt_40to80->AddText(Form("%d-%d%%",centBins[cbin+4],centBins[cbin+5]));
            txv1HFodd_and_absdiff_pt_40to80->Draw();
        } else {
            absDiffv1_pt[anal][cbin]->Draw("same");

            TF1 * fit1 = new TF1("fit1", "pol0", 0, 12);
            fit1->SetLineColor(kBlue);
            absDiffv1_pt[anal][cbin]->Fit(fit1,"QR");
            double par0 = fit1->GetParameter(0);
            double par0E = fit1->GetParError(0);
            double par0Chi2 = fit1->GetChisquare();

            TPaveText * tx1HFodd_and_absdiff_pt_fit;
            if (cbin == 4) tx1HFodd_and_absdiff_pt_fit = new TPaveText(0.22, 0.21, 0.72, 0.35,"NDC");
            else tx1HFodd_and_absdiff_pt_fit = new TPaveText(0.06, 0.21, 0.54, 0.35,"NDC");
            SetTPaveTxt(tx1HFodd_and_absdiff_pt_fit, 16);
            tx1HFodd_and_absdiff_pt_fit->AddText(Form("mean: %0.4f #pm %0.4f",par0,par0E));
            tx1HFodd_and_absdiff_pt_fit->AddText(Form("#chi^{2}: %0.4f",par0Chi2));
            tx1HFodd_and_absdiff_pt_fit->Draw();
        }

    }
    cv1HFodd_and_absdiff_pt_40to80->cd(1);
    TPaveText * txv1HFodd_and_absdiff_pt_40to80_1 = new TPaveText(0.22, 0.68, 0.81, 0.87,"NDC");
    SetTPaveTxt(txv1HFodd_and_absdiff_pt_40to80_1, 18);
    txv1HFodd_and_absdiff_pt_40to80_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    txv1HFodd_and_absdiff_pt_40to80_1->AddText("|#eta| < 2.4");
    txv1HFodd_and_absdiff_pt_40to80_1->Draw();

    TPaveText * txv1HFoddAbsDiff_pt_40to80_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1HFoddAbsDiff_pt_40to80_2, 18);
    txv1HFoddAbsDiff_pt_40to80_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1HFoddAbsDiff_pt_40to80_2->Draw();

    cv1HFodd_and_absdiff_pt_40to80->cd(5);
    TLegend * legv1HFoddAbsDiff_pt_40to80 = new TLegend(0.24, 0.67, 0.51, 0.84);
    SetLegend(legv1HFoddAbsDiff_pt_40to80, 18);
    legv1HFoddAbsDiff_pt_40to80->SetHeader("|v_{1}(+#eta)| - |v_{1}(-#eta)|");
    legv1HFoddAbsDiff_pt_40to80->AddEntry(fitleg," pol0 fit","l");
    legv1HFoddAbsDiff_pt_40to80->Draw();

    cv1HFodd_and_absdiff_pt_40to80->Print(Form("plots/intv1/intv1_pt/int%s/v1odd_with_absdiff_pt_40to80cent_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFodd_and_absdiff_pt_40to80->Close();



    //-- plot both v1odd(pT) and it's absolute differences
    // 0 - 40% centrality
    anal = 15;
    TCanvas * cv1Trkeven_and_absdiff_pt_0to40 = new TCanvas("cv1Trkeven_and_absdiff_pt_0to40","cv1Trkeven_and_absdiff_pt_0to40",1100,700);
    TH1D * hv1Trkeven_and_absdiff_pt_0to40_tmp = new TH1D("hv1Trkeven_and_absdiff_pt_0to40", "", 40, 0, 12);
    hv1Trkeven_and_absdiff_pt_0to40_tmp->SetTitle("");
    hv1Trkeven_and_absdiff_pt_0to40_tmp->SetStats(0);
    hv1Trkeven_and_absdiff_pt_0to40_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1Trkeven_and_absdiff_pt_0to40_tmp->SetYTitle("v_{1}^{even}");
    hv1Trkeven_and_absdiff_pt_0to40_tmp->GetYaxis()->SetRangeUser(-0.05, 0.26);
    hv1Trkeven_and_absdiff_pt_0to40_tmp->GetYaxis()->CenterTitle();
    hv1Trkeven_and_absdiff_pt_0to40_tmp->SetNdivisions(509);
    cv1Trkeven_and_absdiff_pt_0to40->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1Trkeven_and_absdiff_pt_0to40 = (TPad *) cv1Trkeven_and_absdiff_pt_0to40->cd(cbin+1);
        if (gridlines) padv1Trkeven_and_absdiff_pt_0to40->SetGrid();
        if (cbin == 3 || cbin == 7) padv1Trkeven_and_absdiff_pt_0to40->SetRightMargin(0.02);
        if (cbin <= 3) {
            padv1Trkeven_and_absdiff_pt_0to40->SetTopMargin(0.08);
            padv1Trkeven_and_absdiff_pt_0to40->SetBottomMargin(0.14);
        }
        if (cbin >= 4) padv1Trkeven_and_absdiff_pt_0to40->SetTopMargin(0.15);
        TH1D * hv1Trkeven_and_absdiff_pt_0to40 = (TH1D *) hv1Trkeven_and_absdiff_pt_0to40_tmp->Clone(Form("hv1Trkeven_and_absdiff_pt_0to40_%c",cbin));
        if (cbin == 0) {
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->CenterTitle();
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetTitleSize(0.06);
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetTitleOffset(1.12);
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetLabelSize(0.06);
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetLabelOffset(0.018);
            hv1Trkeven_and_absdiff_pt_0to40->GetYaxis()->SetTitleSize(0.07);
            hv1Trkeven_and_absdiff_pt_0to40->GetYaxis()->SetTitleOffset(1.33);
            hv1Trkeven_and_absdiff_pt_0to40->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin >= 1 && cbin <= 3) {
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->CenterTitle();
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetTitleSize(0.07);
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetTitleOffset(1.00);
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetLabelSize(0.07);
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetLabelOffset(0.008);
        }
        if (cbin >= 4) {
            hv1Trkeven_and_absdiff_pt_0to40->SetYTitle("");
            hv1Trkeven_and_absdiff_pt_0to40->GetYaxis()->SetRangeUser(-0.15, 0.15);
        }
        if (cbin == 4) {
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->CenterTitle();
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetTitleSize(0.06);
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetTitleOffset(1.12);
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetLabelSize(0.06);
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetLabelOffset(0.018);
            hv1Trkeven_and_absdiff_pt_0to40->GetYaxis()->SetTitleSize(0.06);
            hv1Trkeven_and_absdiff_pt_0to40->GetYaxis()->SetTitleOffset(1.50);
            hv1Trkeven_and_absdiff_pt_0to40->GetYaxis()->SetLabelSize(0.05);
            hv1Trkeven_and_absdiff_pt_0to40->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->CenterTitle();
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetTitleSize(0.07);
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetTitleOffset(1.00);
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetLabelSize(0.07);
            hv1Trkeven_and_absdiff_pt_0to40->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1Trkeven_and_absdiff_pt_0to40->Draw();

        TPaveText * txv1Trkeven_and_absdiff_pt_0to40;
        if (cbin <= 3) {
            v1p_pt[anal][cbin]->Draw("same");
            v1m_pt[anal][cbin]->Draw("same");
            // v1odd_pt[anal][cbin]->Draw("same");
            if (cbin == 0) txv1Trkeven_and_absdiff_pt_0to40 = new TPaveText(0.75, 0.80, 0.97, 0.87,"NDC");
            else if (cbin >= 1 && cbin <= 3) txv1Trkeven_and_absdiff_pt_0to40 = new TPaveText(0.72, 0.80, 0.92, 0.87,"NDC");
            SetTPaveTxt(txv1Trkeven_and_absdiff_pt_0to40, 18);
            txv1Trkeven_and_absdiff_pt_0to40->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
            txv1Trkeven_and_absdiff_pt_0to40->Draw();
        } else {
            absDiffv1_pt[anal][cbin-4]->Draw("same");

            TF1 * fit1 = new TF1("fit1", "pol0", 0, 12);
            fit1->SetLineColor(kBlue);
            absDiffv1_pt[anal][cbin-4]->Fit(fit1,"QR");
            double par0 = fit1->GetParameter(0);
            double par0E = fit1->GetParError(0);
            double par0Chi2 = fit1->GetChisquare();

            TPaveText * tx1Trkeven_and_absdiff_pt_fit;
            if (cbin == 4) tx1Trkeven_and_absdiff_pt_fit = new TPaveText(0.22, 0.21, 0.72, 0.35,"NDC");
            else tx1Trkeven_and_absdiff_pt_fit = new TPaveText(0.06, 0.21, 0.54, 0.35,"NDC");
            SetTPaveTxt(tx1Trkeven_and_absdiff_pt_fit, 16);
            tx1Trkeven_and_absdiff_pt_fit->AddText(Form("mean: %0.4f #pm %0.4f",par0,par0E));
            tx1Trkeven_and_absdiff_pt_fit->AddText(Form("#chi^{2}: %0.4f",par0Chi2));
            tx1Trkeven_and_absdiff_pt_fit->Draw();
        }

    }
    cv1Trkeven_and_absdiff_pt_0to40->cd(1);
    // TPaveText * txv1Trkeven_and_absdiff_pt_0to40_1 = new TPaveText(0.22, 0.68, 0.81, 0.87,"NDC");
    // SetTPaveTxt(txv1Trkeven_and_absdiff_pt_0to40_1, 18);
    // txv1Trkeven_and_absdiff_pt_0to40_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    // txv1Trkeven_and_absdiff_pt_0to40_1->AddText("|#eta| < 2.4 (GeV/c)");
    // txv1Trkeven_and_absdiff_pt_0to40_1->Draw();

    TLegend * legv1Trkeven_and_absdiff_pt_0to40 = new TLegend(0.22, 0.71, 0.62, 0.87);
    SetLegend(legv1Trkeven_and_absdiff_pt_0to40, 18);
    legv1Trkeven_and_absdiff_pt_0to40->AddEntry(v1p_pt[anal][0]," 0 < #eta < 2.4","p");
    legv1Trkeven_and_absdiff_pt_0to40->AddEntry(v1m_pt[anal][0],"-2.4 < #eta < 0","p");
    legv1Trkeven_and_absdiff_pt_0to40->Draw();

    TPaveText * txv1TrkevenAbsDiff_pt_0to40_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1TrkevenAbsDiff_pt_0to40_2, 18);
    txv1TrkevenAbsDiff_pt_0to40_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1TrkevenAbsDiff_pt_0to40_2->Draw();

    cv1Trkeven_and_absdiff_pt_0to40->cd(5);
    TLegend * legv1TrkevenAbsDiff_pt_0to40 = new TLegend(0.24, 0.67, 0.51, 0.84);
    SetLegend(legv1TrkevenAbsDiff_pt_0to40, 18);
    legv1TrkevenAbsDiff_pt_0to40->SetHeader("|v_{1}(+#eta)| - |v_{1}(-#eta)|");
    legv1TrkevenAbsDiff_pt_0to40->AddEntry(fitleg," pol0 fit","l");
    legv1TrkevenAbsDiff_pt_0to40->Draw();

    cv1Trkeven_and_absdiff_pt_0to40->Print(Form("plots/intv1/intv1_pt/int%s/v1even_with_absdiff_pt_0to40cent_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1Trkeven_and_absdiff_pt_0to40->Close();



    // 40 - 80% centrality
    TCanvas * cv1Trkeven_and_absdiff_pt_40to80 = new TCanvas("cv1Trkeven_and_absdiff_pt_40to80","cv1Trkeven_and_absdiff_pt_40to80",1100,700);
    TH1D * hv1Trkeven_and_absdiff_pt_40to80_tmp = new TH1D("hv1Trkeven_and_absdiff_pt_40to80", "", 40, 0, 12);
    hv1Trkeven_and_absdiff_pt_40to80_tmp->SetTitle("");
    hv1Trkeven_and_absdiff_pt_40to80_tmp->SetStats(0);
    hv1Trkeven_and_absdiff_pt_40to80_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1Trkeven_and_absdiff_pt_40to80_tmp->SetYTitle("v_{1}^{even}");
    hv1Trkeven_and_absdiff_pt_40to80_tmp->GetYaxis()->SetRangeUser(-0.06, 0.26);
    hv1Trkeven_and_absdiff_pt_40to80_tmp->GetYaxis()->CenterTitle();
    hv1Trkeven_and_absdiff_pt_40to80_tmp->SetNdivisions(509);
    cv1Trkeven_and_absdiff_pt_40to80->Divide(4,2,0,0);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TPad * padv1Trkeven_and_absdiff_pt_40to80 = (TPad *) cv1Trkeven_and_absdiff_pt_40to80->cd(cbin+1);
        if (gridlines) padv1Trkeven_and_absdiff_pt_40to80->SetGrid();
        if (cbin == 3 || cbin == 7) padv1Trkeven_and_absdiff_pt_40to80->SetRightMargin(0.02);
        if (cbin <= 3) {
            padv1Trkeven_and_absdiff_pt_40to80->SetTopMargin(0.08);
            padv1Trkeven_and_absdiff_pt_40to80->SetBottomMargin(0.14);
        }
        if (cbin >= 4) padv1Trkeven_and_absdiff_pt_40to80->SetTopMargin(0.15);
        TH1D * hv1Trkeven_and_absdiff_pt_40to80 = (TH1D *) hv1Trkeven_and_absdiff_pt_40to80_tmp->Clone(Form("hv1Trkeven_and_absdiff_pt_40to80_%c",cbin));
        if (cbin == 0) {
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->CenterTitle();
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetTitleSize(0.06);
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetTitleOffset(1.12);
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetLabelSize(0.06);
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetLabelOffset(0.018);
            hv1Trkeven_and_absdiff_pt_40to80->GetYaxis()->SetTitleSize(0.07);
            hv1Trkeven_and_absdiff_pt_40to80->GetYaxis()->SetTitleOffset(1.33);
            hv1Trkeven_and_absdiff_pt_40to80->GetYaxis()->SetLabelSize(0.06);
        }
        if (cbin >= 1 && cbin <= 3) {
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->CenterTitle();
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetTitleSize(0.07);
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetTitleOffset(1.00);
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetLabelSize(0.07);
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetLabelOffset(0.008);
        }
        if (cbin >= 4) {
            hv1Trkeven_and_absdiff_pt_40to80->SetYTitle("");
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetRangeUser(0, 12);
            hv1Trkeven_and_absdiff_pt_40to80->GetYaxis()->SetRangeUser(-0.15, 0.15);
        }
        if (cbin == 4) {
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->CenterTitle();
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetTitleSize(0.06);
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetTitleOffset(1.12);
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetLabelSize(0.06);
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetLabelOffset(0.018);
            hv1Trkeven_and_absdiff_pt_40to80->GetYaxis()->SetTitleSize(0.06);
            hv1Trkeven_and_absdiff_pt_40to80->GetYaxis()->SetTitleOffset(1.50);
            hv1Trkeven_and_absdiff_pt_40to80->GetYaxis()->SetLabelSize(0.05);
            hv1Trkeven_and_absdiff_pt_40to80->GetYaxis()->SetLabelOffset(0.010);
        }
        if (cbin >=5) {
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->CenterTitle();
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetTitleSize(0.07);
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetTitleOffset(1.00);
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetLabelSize(0.07);
            hv1Trkeven_and_absdiff_pt_40to80->GetXaxis()->SetLabelOffset(0.008);
        }
        hv1Trkeven_and_absdiff_pt_40to80->Draw();

        TPaveText * txv1Trkeven_and_absdiff_pt_40to80;
        if (cbin <= 3) {
            v1p_pt[anal][cbin+4]->Draw("same");
            v1m_pt[anal][cbin+4]->Draw("same");
            // v1odd_pt[anal][cbin+4]->Draw("same");
            if (cbin == 0) txv1Trkeven_and_absdiff_pt_40to80 = new TPaveText(0.75, 0.80, 0.97, 0.87,"NDC");
            else if (cbin >= 1 && cbin <= 3) txv1Trkeven_and_absdiff_pt_40to80 = new TPaveText(0.72, 0.80, 0.92, 0.87,"NDC");
            SetTPaveTxt(txv1Trkeven_and_absdiff_pt_40to80, 18);
            txv1Trkeven_and_absdiff_pt_40to80->AddText(Form("%d-%d%%",centBins[cbin+4],centBins[cbin+5]));
            txv1Trkeven_and_absdiff_pt_40to80->Draw();
        } else {
            absDiffv1_pt[anal][cbin]->Draw("same");

            TF1 * fit1 = new TF1("fit1", "pol0", 0, 12);
            fit1->SetLineColor(kBlue);
            absDiffv1_pt[anal][cbin]->Fit(fit1,"QR");
            double par0 = fit1->GetParameter(0);
            double par0E = fit1->GetParError(0);
            double par0Chi2 = fit1->GetChisquare();

            TPaveText * tx1Trkeven_and_absdiff_pt_fit;
            if (cbin == 4) tx1Trkeven_and_absdiff_pt_fit = new TPaveText(0.22, 0.21, 0.72, 0.35,"NDC");
            else tx1Trkeven_and_absdiff_pt_fit = new TPaveText(0.06, 0.21, 0.54, 0.35,"NDC");
            SetTPaveTxt(tx1Trkeven_and_absdiff_pt_fit, 16);
            tx1Trkeven_and_absdiff_pt_fit->AddText(Form("mean: %0.4f #pm %0.4f",par0,par0E));
            tx1Trkeven_and_absdiff_pt_fit->AddText(Form("#chi^{2}: %0.4f",par0Chi2));
            tx1Trkeven_and_absdiff_pt_fit->Draw();
        }

    }
    cv1Trkeven_and_absdiff_pt_40to80->cd(1);
    // TPaveText * txv1Trkeven_and_absdiff_pt_40to80_1 = new TPaveText(0.22, 0.68, 0.81, 0.87,"NDC");
    // SetTPaveTxt(txv1Trkeven_and_absdiff_pt_40to80_1, 18);
    // txv1Trkeven_and_absdiff_pt_40to80_1->AddText("PbPb #sqrt{s_{NN}} = 5.02 TeV");
    // txv1Trkeven_and_absdiff_pt_40to80_1->AddText("|#eta| < 2.4");
    // txv1Trkeven_and_absdiff_pt_40to80_1->Draw();

    TLegend * legv1Trkeven_and_absdiff_pt_40to80 = new TLegend(0.22, 0.71, 0.62, 0.87);
    SetLegend(legv1Trkeven_and_absdiff_pt_40to80, 18);
    legv1Trkeven_and_absdiff_pt_40to80->AddEntry(v1p_pt[anal][0]," 0 < #eta < 2.4","p");
    legv1Trkeven_and_absdiff_pt_40to80->AddEntry(v1m_pt[anal][0],"-2.4 < #eta < 0","p");
    legv1Trkeven_and_absdiff_pt_40to80->Draw();

    TPaveText * txv1TrkevenAbsDiff_pt_40to80_2 = new TPaveText(0.18, 0.93, 0.58, 1.0,"NDC");
    SetTPaveTxt(txv1TrkevenAbsDiff_pt_40to80_2, 18);
    txv1TrkevenAbsDiff_pt_40to80_2->AddText("#bf{CMS} #it{Preliminary}");
    txv1TrkevenAbsDiff_pt_40to80_2->Draw();

    cv1Trkeven_and_absdiff_pt_40to80->cd(5);
    TLegend * legv1TrkevenAbsDiff_pt_40to80 = new TLegend(0.24, 0.67, 0.51, 0.84);
    SetLegend(legv1TrkevenAbsDiff_pt_40to80, 18);
    legv1TrkevenAbsDiff_pt_40to80->SetHeader("|v_{1}(+#eta)| - |v_{1}(-#eta)|");
    legv1TrkevenAbsDiff_pt_40to80->AddEntry(fitleg," pol0 fit","l");
    legv1TrkevenAbsDiff_pt_40to80->Draw();

    cv1Trkeven_and_absdiff_pt_40to80->Print(Form("plots/intv1/intv1_pt/int%s/v1odd_with_absdiff_pt_40to80cent_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1Trkeven_and_absdiff_pt_40to80->Close();


}
