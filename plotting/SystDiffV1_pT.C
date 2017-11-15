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

TH1D * v1p_pt[nanals][ncentbins][netabins];
TH1D * v1m_pt[nanals][ncentbins][netabins];
TH1D * v1odd_pt[nanals][ncentbins][netabins];
TH1D * v1even_pt[nanals][ncentbins][netabins];

TH1D * v112p_pt[nanals][ncentbins][netabins];
TH1D * v112m_pt[nanals][ncentbins][netabins];
TH1D * v112odd_pt[nanals][ncentbins][netabins];
TH1D * v112even_pt[nanals][ncentbins][netabins];

TH1D * v123p_pt[nanals][ncentbins][netabins];
TH1D * v123m_pt[nanals][ncentbins][netabins];
TH1D * v123odd_pt[nanals][ncentbins][netabins];
TH1D * v123even_pt[nanals][ncentbins][netabins];

// ratios
TH1D * ratiov1pm_pT[nanals][ncentbins][netabins];
TH1D * ratiov1odd_pT[nanals][ncentbins][netabins/2];
TH1D * ratiov1even_pT[nanals][ncentbins][netabins/2];

TH1D * runParms[nanals];

TFile * tfin;

void SystDiffV1_pT()
{

    TH1::SetDefaultSumw2();

    tfin = new TFile("../outputs/final_outputs/v1Diff.root","read");

    //-- retrieve histograms from final output file
    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            for (int ebin = 0; ebin<netabins; ebin++) {
                string tag0 = Form("cent_%d-%d/eta_%0.1f-%0.1f",centBins[cbin],centBins[cbin+1],etabins[ebin],etabins[ebin+1]);
                string tag1 = Form("%s_c%d_e%d",AnalNames[i].data(),cbin,ebin);

                v1p_pt[i][cbin][ebin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%s/v1p_pt_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v1m_pt[i][cbin][ebin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%s/v1m_pt_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v1odd_pt[i][cbin][ebin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%s/v1odd_pt_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v1even_pt[i][cbin][ebin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%s/v1even_pt_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
            }

            for (int ebin = 0; ebin<netabins; ebin++) {
                ratiov1pm_pT[i][cbin][ebin] = (TH1D *) v1p_pt[i][cbin][ebin]->Clone(Form("ratiov1pm_pT_%s_%d_%d",AnalNames[i].data(),cbin,ebin));
                ratiov1pm_pT[i][cbin][ebin]->Divide(v1m_pt[i][cbin][netabins-ebin-1]);
            }
            for (int ebin = 0; ebin<netabins/2; ebin++) {
                ratiov1odd_pT[i][cbin][ebin] = (TH1D *) v1odd_pt[i][cbin][netabins/2+ebin]->Clone(Form("ratiov1odd_pT_%s_%d_%d",AnalNames[i].data(),cbin,ebin));
                ratiov1odd_pT[i][cbin][ebin]->Divide(v1odd_pt[i][cbin][netabins/2-ebin-1]);

                ratiov1even_pT[i][cbin][ebin] = (TH1D *) v1even_pt[i][cbin][netabins/2+ebin]->Clone(Form("ratiov1even_pT_%s_%d_%d",AnalNames[i].data(),cbin,ebin));
                ratiov1even_pT[i][cbin][ebin]->Divide(v1even_pt[i][cbin][netabins/2-ebin-1]);
            }
        }
    }

    //-- plotting options

    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            for (int ebin = 0; ebin<netabins; ebin++) {
                ratiov1pm_pT[i][cbin][ebin]->SetMarkerColor(kBlack);
                ratiov1pm_pT[i][cbin][ebin]->SetLineColor(kBlack);
                ratiov1pm_pT[i][cbin][ebin]->SetMarkerStyle(21);
                ratiov1pm_pT[i][cbin][ebin]->SetMarkerSize(1.1);
            }
            for (int ebin = 0; ebin<netabins/2; ebin++) {
                ratiov1odd_pT[i][cbin][ebin]->SetMarkerColor(kBlack);
                ratiov1odd_pT[i][cbin][ebin]->SetLineColor(kBlack);
                ratiov1odd_pT[i][cbin][ebin]->SetMarkerStyle(21);
                ratiov1odd_pT[i][cbin][ebin]->SetMarkerSize(1.1);

                ratiov1even_pT[i][cbin][ebin]->SetMarkerColor(kBlack);
                ratiov1even_pT[i][cbin][ebin]->SetLineColor(kBlack);
                ratiov1even_pT[i][cbin][ebin]->SetMarkerStyle(21);
                ratiov1even_pT[i][cbin][ebin]->SetMarkerSize(1.1);
            }
        }
    }


    //-- make plots
    if (!fopen("plots","r")) system("mkdir plots");
    if (!fopen("plots/diffv1","r")) system("mkdir plots/diffv1");
    if (!fopen("plots/diffv1/diffv1_pT","r")) system("mkdir plots/diffv1/diffv1_pT");

    int centMarkerStyle[] = {21, 24, 20, 25, 33, 27, 34, 28};
    float centMarkerSize[] = {1.1, 1.2, 1.2, 1.1, 1.7, 1.7, 1.5};

    int anal; // choice of analysis

    // differential v1(pT) using HF+/- for each centrality bin
    anal = 7;
    if (!fopen(Form("plots/diffv1/diffv1_pT/diff%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/diffv1/diffv1_pT/diff%s",AnalNames[anal].data()));

    TCanvas * cratiov1HFpm_pT[ncentbins];
    TH2D * hratiov1HFpm_pT_tmp = new TH2D("hratiov1HFpm_pT_tmp", "", 40, 0, 12, 40, -10, 10);
    hratiov1HFpm_pT_tmp->SetTitle("");
    hratiov1HFpm_pT_tmp->SetStats(0);
    hratiov1HFpm_pT_tmp->SetXTitle("p_{T} (GeV/c)");
    hratiov1HFpm_pT_tmp->SetYTitle("v_{1}^{HF+/-} ratio");
    hratiov1HFpm_pT_tmp->GetYaxis()->SetRangeUser(0.5,1.5);
    hratiov1HFpm_pT_tmp->SetNdivisions(509);

    TLine * lnpT_long = new TLine(0.0, 1.0, 12.0, 1.0);
    lnpT_long->SetLineColor(kBlack);
    lnpT_long->SetLineWidth(1);
    for (int cbin = 0; cbin<ncentbins; cbin++) {

        cratiov1HFpm_pT[cbin] = new TCanvas(Form("cratiov1HFpm_pT_cent%d-%d",centBins[cbin],centBins[cbin+1]),"cratiov1HFpm_pT",1100,850);
        cratiov1HFpm_pT[cbin]->Divide(4,3,0,0);

        for (int ebin = 0; ebin<netabins; ebin++) {
            TPad * padratiov1HFpm_pT = (TPad *) cratiov1HFpm_pT[cbin]->cd(ebin+1);
            if (gridlines) padratiov1HFpm_pT->SetGrid();
            if (ebin == 3 || ebin == 7 || ebin == 11) padratiov1HFpm_pT->SetRightMargin(0.02);
            TH2D * hratiov1HFpm_pT = (TH2D *) hratiov1HFpm_pT_tmp->Clone(Form("hratiov1HFpm_pT_%c_%d",cbin,ebin));
            if (ebin == 0 || ebin == 4) {
                hratiov1HFpm_pT->GetYaxis()->CenterTitle();
                hratiov1HFpm_pT->GetYaxis()->SetTitleSize(0.07);
                hratiov1HFpm_pT->GetYaxis()->SetTitleOffset(1.34);
                hratiov1HFpm_pT->GetYaxis()->SetLabelSize(0.06);
            } else if (ebin == 8) {
                hratiov1HFpm_pT->GetXaxis()->SetTitleSize(0.06);
                hratiov1HFpm_pT->GetXaxis()->SetTitleOffset(1.14);
                hratiov1HFpm_pT->GetYaxis()->CenterTitle();
                hratiov1HFpm_pT->GetYaxis()->SetTitleSize(0.06);
                hratiov1HFpm_pT->GetYaxis()->SetTitleOffset(1.48);
                hratiov1HFpm_pT->GetYaxis()->SetLabelSize(0.05);
            } else if (ebin>=9) {
                hratiov1HFpm_pT->GetXaxis()->SetTitleSize(0.07);
                hratiov1HFpm_pT->GetXaxis()->SetTitleOffset(1.00);
                hratiov1HFpm_pT->GetXaxis()->SetLabelSize(0.06);
                hratiov1HFpm_pT->GetXaxis()->SetLabelOffset(0.005);
            }
            hratiov1HFpm_pT->Draw();
            ratiov1pm_pT[anal][cbin][ebin]->Scale(-1);
            ratiov1pm_pT[anal][cbin][ebin]->Draw("same");
            lnpT_long->Draw();

            TPaveText * txratiov1HFpm_pT;
            if (ebin == 0 || ebin == 4) txratiov1HFpm_pT = new TPaveText(0.23, 0.05, 0.62, 0.15,"NDC");
            else if (ebin == 8) txratiov1HFpm_pT = new TPaveText(0.23, 0.19, 0.62, 0.29,"NDC");
            else if (ebin>=9) txratiov1HFpm_pT = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
            else txratiov1HFpm_pT = new TPaveText(0.05, 0.05, 0.50, 0.15,"NDC");
            SetTPaveTxt(txratiov1HFpm_pT, 18);
            txratiov1HFpm_pT->AddText(Form("%0.1f < #eta < %0.1f",etabins[ebin],etabins[ebin+1]));
            txratiov1HFpm_pT->Draw();
        }
        cratiov1HFpm_pT[cbin]->cd(1);
        TPaveText * txratiov1HFpm_pT_cent = new TPaveText(0.24, 0.84, 0.38, 0.95,"NDC");
        SetTPaveTxt(txratiov1HFpm_pT_cent, 18);
        txratiov1HFpm_pT_cent->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txratiov1HFpm_pT_cent->Draw();

        cratiov1HFpm_pT[cbin]->cd(5);
        TPaveText * txratiov1HFpm_pT_1 = new TPaveText(0.24, 0.72, 0.44, 0.92,"NDC");
        SetTPaveTxt(txratiov1HFpm_pT_1, 18);
        txratiov1HFpm_pT_1->AddText("#frac{|v_{1}^{HF+}(#eta)|}{|v_{1}^{HF-}(-#eta)|}");
        txratiov1HFpm_pT_1->Draw();

        cratiov1HFpm_pT[cbin]->Print(Form("plots/diffv1/diffv1_pT/diff%s/ratio_v1_pm_pT_%s_cent%d-%d.png",AnalNames[anal].data(),AnalNames[anal].data(),centBins[cbin],centBins[cbin+1]),"png");
        if (close_plots) cratiov1HFpm_pT[cbin]->Close();
    }



    // differential v1^odd(pT) for each centrality bin

    TCanvas * cratiov1HFodd_pT[ncentbins];
    TH2D * hratiov1HFodd_pT_tmp = new TH2D("hratiov1HFodd_pT_tmp", "", 40, 0, 12, 40, -10, 10);
    hratiov1HFodd_pT_tmp->SetTitle("");
    hratiov1HFodd_pT_tmp->SetStats(0);
    hratiov1HFodd_pT_tmp->SetXTitle("p_{T} (GeV/c)");
    hratiov1HFodd_pT_tmp->SetYTitle("v_{1}^{odd} ratio");
    hratiov1HFodd_pT_tmp->GetYaxis()->SetRangeUser(-1,3);
    hratiov1HFodd_pT_tmp->SetNdivisions(509);
    for (int cbin = 0; cbin<ncentbins; cbin++) {

        cratiov1HFodd_pT[cbin] = new TCanvas(Form("cratiov1HFodd_pT_cent%d-%d",centBins[cbin],centBins[cbin+1]),"cratiov1HFodd_pT",950,650);
        cratiov1HFodd_pT[cbin]->Divide(3,2,0,0);

        for (int ebin = 0; ebin<netabins/2; ebin++) {
            TPad * padratiov1HFodd_pT = (TPad *) cratiov1HFodd_pT[cbin]->cd(ebin+1);
            if (gridlines) padratiov1HFodd_pT->SetGrid();
            if (ebin == 2 || ebin == 5) padratiov1HFodd_pT->SetRightMargin(0.02);
            TH2D * hratiov1HFodd_pT = (TH2D *) hratiov1HFodd_pT_tmp->Clone(Form("hratiov1HFodd_pT_%c_%d",cbin,ebin));
            if (ebin == 0 || ebin == 4) {
                hratiov1HFodd_pT->GetYaxis()->CenterTitle();
                hratiov1HFodd_pT->GetYaxis()->SetTitleSize(0.07);
                hratiov1HFodd_pT->GetYaxis()->SetTitleOffset(1.34);
                hratiov1HFodd_pT->GetYaxis()->SetLabelSize(0.06);
            } else if (ebin == 8) {
                hratiov1HFodd_pT->GetXaxis()->SetTitleSize(0.06);
                hratiov1HFodd_pT->GetXaxis()->SetTitleOffset(1.14);
                hratiov1HFodd_pT->GetYaxis()->CenterTitle();
                hratiov1HFodd_pT->GetYaxis()->SetTitleSize(0.06);
                hratiov1HFodd_pT->GetYaxis()->SetTitleOffset(1.48);
                hratiov1HFodd_pT->GetYaxis()->SetLabelSize(0.05);
            } else if (ebin>=9) {
                hratiov1HFodd_pT->GetXaxis()->SetTitleSize(0.07);
                hratiov1HFodd_pT->GetXaxis()->SetTitleOffset(1.00);
                hratiov1HFodd_pT->GetXaxis()->SetLabelSize(0.06);
                hratiov1HFodd_pT->GetXaxis()->SetLabelOffset(0.005);
            }
            hratiov1HFodd_pT->Draw();
            ratiov1odd_pT[anal][cbin][ebin]->Scale(-1);
            ratiov1odd_pT[anal][cbin][ebin]->Draw("same");
            lnpT_long->Draw();

            TPaveText * txratiov1HFodd_pT;
            if (ebin == 0) txratiov1HFodd_pT = new TPaveText(0.25, 0.06, 0.59, 0.16,"NDC");
            else if (ebin == 1 || ebin == 2) txratiov1HFodd_pT = new TPaveText(0.05, 0.06, 0.49, 0.16,"NDC");
            else if (ebin == 3) txratiov1HFodd_pT = new TPaveText(0.25, 0.19, 0.59, 0.29,"NDC");
            else if (ebin>=4) txratiov1HFodd_pT = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
            else txratiov1HFodd_pT = new TPaveText(0.05, 0.19, 0.49, 0.29,"NDC");
            SetTPaveTxt(txratiov1HFodd_pT, 18);
            txratiov1HFodd_pT->AddText(Form("%0.1f < |#eta| < %0.1f",etabins[netabins/2+ebin],etabins[netabins/2-ebin-1]));
            txratiov1HFodd_pT->Draw();
        }
        cratiov1HFodd_pT[cbin]->cd(1);
        TPaveText * txratiov1HFodd_pT_cent = new TPaveText(0.24, 0.84, 0.38, 0.95,"NDC");
        SetTPaveTxt(txratiov1HFodd_pT_cent, 18);
        txratiov1HFodd_pT_cent->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txratiov1HFodd_pT_cent->Draw();

        cratiov1HFodd_pT[cbin]->cd(4);
        TPaveText * txratiov1HFodd_pT_1 = new TPaveText(0.24, 0.78, 0.44, 0.98,"NDC");
        SetTPaveTxt(txratiov1HFodd_pT_1, 18);
        txratiov1HFodd_pT_1->AddText("#frac{|v_{1}^{odd}(+#eta)|}{|v_{1}^{odd}(-#eta)|}");
        txratiov1HFodd_pT_1->Draw();

        cratiov1HFodd_pT[cbin]->Print(Form("plots/diffv1/diffv1_pT/diff%s/ratio_v1_odd_pT_%s_cent%d-%d.png",AnalNames[anal].data(),AnalNames[anal].data(),centBins[cbin],centBins[cbin+1]),"png");
        if (close_plots) cratiov1HFodd_pT[cbin]->Close();
    }



    // differential v1^even(pT) for each centrality bin

    TCanvas * cratiov1HFeven_pT[ncentbins];
    TH2D * hratiov1HFeven_pT_tmp = new TH2D("hratiov1HFeven_pT_tmp", "", 40, 0, 12, 40, -10, 10);
    hratiov1HFeven_pT_tmp->SetTitle("");
    hratiov1HFeven_pT_tmp->SetStats(0);
    hratiov1HFeven_pT_tmp->SetXTitle("p_{T} (GeV/c)");
    hratiov1HFeven_pT_tmp->SetYTitle("v_{1}^{even} ratio");
    hratiov1HFeven_pT_tmp->GetYaxis()->SetRangeUser(-1,3);
    hratiov1HFeven_pT_tmp->SetNdivisions(509);
    for (int cbin = 0; cbin<ncentbins; cbin++) {

        cratiov1HFeven_pT[cbin] = new TCanvas(Form("cratiov1HFeven_pT_cent%d-%d",centBins[cbin],centBins[cbin+1]),"cratiov1HFeven_pT",950,650);
        cratiov1HFeven_pT[cbin]->Divide(3,2,0,0);

        for (int ebin = 0; ebin<netabins/2; ebin++) {
            TPad * padratiov1HFeven_pT = (TPad *) cratiov1HFeven_pT[cbin]->cd(ebin+1);
            if (gridlines) padratiov1HFeven_pT->SetGrid();
            if (ebin == 2 || ebin == 5) padratiov1HFeven_pT->SetRightMargin(0.02);
            TH2D * hratiov1HFeven_pT = (TH2D *) hratiov1HFeven_pT_tmp->Clone(Form("hratiov1HFeven_pT_%c_%d",cbin,ebin));
            if (ebin == 0 || ebin == 4) {
                hratiov1HFeven_pT->GetYaxis()->CenterTitle();
                hratiov1HFeven_pT->GetYaxis()->SetTitleSize(0.07);
                hratiov1HFeven_pT->GetYaxis()->SetTitleOffset(1.34);
                hratiov1HFeven_pT->GetYaxis()->SetLabelSize(0.06);
            } else if (ebin == 8) {
                hratiov1HFeven_pT->GetXaxis()->SetTitleSize(0.06);
                hratiov1HFeven_pT->GetXaxis()->SetTitleOffset(1.14);
                hratiov1HFeven_pT->GetYaxis()->CenterTitle();
                hratiov1HFeven_pT->GetYaxis()->SetTitleSize(0.06);
                hratiov1HFeven_pT->GetYaxis()->SetTitleOffset(1.48);
                hratiov1HFeven_pT->GetYaxis()->SetLabelSize(0.05);
            } else if (ebin>=9) {
                hratiov1HFeven_pT->GetXaxis()->SetTitleSize(0.07);
                hratiov1HFeven_pT->GetXaxis()->SetTitleOffset(1.00);
                hratiov1HFeven_pT->GetXaxis()->SetLabelSize(0.06);
                hratiov1HFeven_pT->GetXaxis()->SetLabelOffset(0.005);
            }
            hratiov1HFeven_pT->Draw();
            ratiov1even_pT[anal][cbin][ebin]->Draw("same");
            lnpT_long->Draw();

            TPaveText * txratiov1HFeven_pT;
            if (ebin == 0) txratiov1HFeven_pT = new TPaveText(0.25, 0.06, 0.59, 0.16,"NDC");
            else if (ebin == 1 || ebin == 2) txratiov1HFeven_pT = new TPaveText(0.05, 0.06, 0.49, 0.16,"NDC");
            else if (ebin == 3) txratiov1HFeven_pT = new TPaveText(0.25, 0.19, 0.59, 0.29,"NDC");
            else if (ebin>=4) txratiov1HFeven_pT = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
            else txratiov1HFeven_pT = new TPaveText(0.05, 0.19, 0.49, 0.29,"NDC");
            SetTPaveTxt(txratiov1HFeven_pT, 18);
            txratiov1HFeven_pT->AddText(Form("%0.1f < |#eta| < %0.1f",etabins[netabins/2+ebin],etabins[netabins/2-ebin-1]));
            txratiov1HFeven_pT->Draw();
        }
        cratiov1HFeven_pT[cbin]->cd(1);
        TPaveText * txratiov1HFeven_pT_cent = new TPaveText(0.24, 0.84, 0.38, 0.95,"NDC");
        SetTPaveTxt(txratiov1HFeven_pT_cent, 18);
        txratiov1HFeven_pT_cent->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txratiov1HFeven_pT_cent->Draw();

        cratiov1HFeven_pT[cbin]->cd(4);
        TPaveText * txratiov1HFeven_pT_1 = new TPaveText(0.24, 0.78, 0.44, 0.98,"NDC");
        SetTPaveTxt(txratiov1HFeven_pT_1, 18);
        txratiov1HFeven_pT_1->AddText("#frac{|v_{1}^{even}(+#eta)|}{|v_{1}^{even}(-#eta)|}");
        txratiov1HFeven_pT_1->Draw();

        cratiov1HFeven_pT[cbin]->Print(Form("plots/diffv1/diffv1_pT/diff%s/ratio_v1_even_pT_%s_cent%d-%d.png",AnalNames[anal].data(),AnalNames[anal].data(),centBins[cbin],centBins[cbin+1]),"png");
        if (close_plots) cratiov1HFeven_pT[cbin]->Close();
    }



    // differential v1^even(pT) for each centrality bin
    anal = 15;
    if (!fopen(Form("plots/diffv1/diffv1_pT/diff%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/diffv1/diffv1_pT/diff%s",AnalNames[anal].data()));

    TCanvas * cratiov1Trackeven_pT[ncentbins];
    TH2D * hratiov1Trackeven_pT_tmp = new TH2D("hratiov1Trackeven_pT_tmp", "", 40, 0, 12, 40, -10, 10);
    hratiov1Trackeven_pT_tmp->SetTitle("");
    hratiov1Trackeven_pT_tmp->SetStats(0);
    hratiov1Trackeven_pT_tmp->SetXTitle("p_{T} (GeV/c)");
    hratiov1Trackeven_pT_tmp->SetYTitle("v_{1}^{even} ratio");
    hratiov1Trackeven_pT_tmp->GetYaxis()->SetRangeUser(-1,3);
    hratiov1Trackeven_pT_tmp->SetNdivisions(509);
    for (int cbin = 0; cbin<ncentbins; cbin++) {

        cratiov1Trackeven_pT[cbin] = new TCanvas(Form("cratiov1Trackeven_pT_cent%d-%d",centBins[cbin],centBins[cbin+1]),"cratiov1Trackeven_pT",950,650);
        cratiov1Trackeven_pT[cbin]->Divide(3,2,0,0);

        for (int ebin = 0; ebin<netabins/2; ebin++) {
            TPad * padratiov1Trackeven_pT = (TPad *) cratiov1Trackeven_pT[cbin]->cd(ebin+1);
            if (gridlines) padratiov1Trackeven_pT->SetGrid();
            if (ebin == 2 || ebin == 5) padratiov1Trackeven_pT->SetRightMargin(0.02);
            TH2D * hratiov1Trackeven_pT = (TH2D *) hratiov1Trackeven_pT_tmp->Clone(Form("hratiov1Trackeven_pT_%c_%d",cbin,ebin));
            if (ebin == 0 || ebin == 4) {
                hratiov1Trackeven_pT->GetYaxis()->CenterTitle();
                hratiov1Trackeven_pT->GetYaxis()->SetTitleSize(0.07);
                hratiov1Trackeven_pT->GetYaxis()->SetTitleOffset(1.34);
                hratiov1Trackeven_pT->GetYaxis()->SetLabelSize(0.06);
            } else if (ebin == 8) {
                hratiov1Trackeven_pT->GetXaxis()->SetTitleSize(0.06);
                hratiov1Trackeven_pT->GetXaxis()->SetTitleOffset(1.14);
                hratiov1Trackeven_pT->GetYaxis()->CenterTitle();
                hratiov1Trackeven_pT->GetYaxis()->SetTitleSize(0.06);
                hratiov1Trackeven_pT->GetYaxis()->SetTitleOffset(1.48);
                hratiov1Trackeven_pT->GetYaxis()->SetLabelSize(0.05);
            } else if (ebin>=9) {
                hratiov1Trackeven_pT->GetXaxis()->SetTitleSize(0.07);
                hratiov1Trackeven_pT->GetXaxis()->SetTitleOffset(1.00);
                hratiov1Trackeven_pT->GetXaxis()->SetLabelSize(0.06);
                hratiov1Trackeven_pT->GetXaxis()->SetLabelOffset(0.005);
            }
            hratiov1Trackeven_pT->Draw();
            ratiov1pm_pT[anal][cbin][ebin]->Draw("same");
            lnpT_long->Draw();

            TPaveText * txratiov1Trackeven_pT;
            if (ebin == 0) txratiov1Trackeven_pT = new TPaveText(0.25, 0.06, 0.59, 0.16,"NDC");
            else if (ebin == 1 || ebin == 2) txratiov1Trackeven_pT = new TPaveText(0.05, 0.06, 0.49, 0.16,"NDC");
            else if (ebin == 3) txratiov1Trackeven_pT = new TPaveText(0.25, 0.19, 0.59, 0.29,"NDC");
            else if (ebin>=4) txratiov1Trackeven_pT = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
            else txratiov1Trackeven_pT = new TPaveText(0.05, 0.19, 0.49, 0.29,"NDC");
            SetTPaveTxt(txratiov1Trackeven_pT, 18);
            txratiov1Trackeven_pT->AddText(Form("%0.1f < |#eta| < %0.1f",etabins[netabins/2+ebin],etabins[netabins/2-ebin-1]));
            txratiov1Trackeven_pT->Draw();
        }
        cratiov1Trackeven_pT[cbin]->cd(1);
        TPaveText * txratiov1Trackeven_pT_cent = new TPaveText(0.24, 0.84, 0.38, 0.95,"NDC");
        SetTPaveTxt(txratiov1Trackeven_pT_cent, 18);
        txratiov1Trackeven_pT_cent->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txratiov1Trackeven_pT_cent->Draw();

        cratiov1Trackeven_pT[cbin]->cd(4);
        TPaveText * txratiov1Trackeven_pT_1 = new TPaveText(0.24, 0.78, 0.44, 0.98,"NDC");
        SetTPaveTxt(txratiov1Trackeven_pT_1, 18);
        txratiov1Trackeven_pT_1->AddText("#frac{v_{1}^{even}(+#eta)}{v_{1}^{even}(-#eta)}");
        txratiov1Trackeven_pT_1->Draw();

        cratiov1Trackeven_pT[cbin]->Print(Form("plots/diffv1/diffv1_pT/diff%s/ratio_v1_even_pT_%s_cent%d-%d.png",AnalNames[anal].data(),AnalNames[anal].data(),centBins[cbin],centBins[cbin+1]),"png");
        if (close_plots) cratiov1Trackeven_pT[cbin]->Close();
    }

}
