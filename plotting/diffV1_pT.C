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

TH1D * runParms[nanals];

TFile * tfin;

void diffV1_pT()
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
                /*
                v112p_pt[i][cbin][ebin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%s/v112p_pt_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v112m_pt[i][cbin][ebin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%s/v112m_pt_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v112odd_pt[i][cbin][ebin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%s/v112odd_pt_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v112even_pt[i][cbin][ebin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%s/v112even_pt_%s",AnalNames[i].data(),tag0.data(),tag1.data()));

                v123p_pt[i][cbin][ebin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%s/v123p_pt_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v123m_pt[i][cbin][ebin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%s/v123m_pt_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v123odd_pt[i][cbin][ebin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%s/v123odd_pt_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                v123even_pt[i][cbin][ebin] = (TH1D *) tfin->Get(Form("%s/v1_pt/%s/v123even_pt_%s",AnalNames[i].data(),tag0.data(),tag1.data()));
                */
            }
        }
    }
    //tfin->Close();

    //-- plotting options

    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            for (int ebin = 0; ebin<netabins; ebin++) {
                v1p_pt[i][cbin][ebin]->SetMarkerColor(kRed);
                v1p_pt[i][cbin][ebin]->SetLineColor(kRed);
                v1p_pt[i][cbin][ebin]->SetMarkerStyle(21);
                v1p_pt[i][cbin][ebin]->SetMarkerSize(1.1);

                v1m_pt[i][cbin][ebin]->SetMarkerColor(kBlue);
                v1m_pt[i][cbin][ebin]->SetLineColor(kBlue);
                v1m_pt[i][cbin][ebin]->SetMarkerStyle(21);
                v1m_pt[i][cbin][ebin]->SetMarkerSize(1.1);

                v1odd_pt[i][cbin][ebin]->SetMarkerColor(kRed);
                v1odd_pt[i][cbin][ebin]->SetLineColor(kRed);
                v1odd_pt[i][cbin][ebin]->SetMarkerStyle(21);
                v1odd_pt[i][cbin][ebin]->SetMarkerSize(1.1);

                v1even_pt[i][cbin][ebin]->SetMarkerColor(kBlue);
                v1even_pt[i][cbin][ebin]->SetLineColor(kBlue);
                v1even_pt[i][cbin][ebin]->SetMarkerStyle(21);
                v1even_pt[i][cbin][ebin]->SetMarkerSize(1.1);
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

    TCanvas * cv1HFpm_pT[ncentbins];
    TH1D * hv1HFpm_pT_tmp = new TH1D("hv1HFpm_pT_tmp", "", 40, 0, 12);
    hv1HFpm_pT_tmp->SetTitle("");
    hv1HFpm_pT_tmp->SetStats(0);
    hv1HFpm_pT_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1HFpm_pT_tmp->SetYTitle("v_{1}");
    hv1HFpm_pT_tmp->GetYaxis()->SetRangeUser(-0.3, 0.3);
    hv1HFpm_pT_tmp->SetNdivisions(509);
    for (int cbin = 0; cbin<ncentbins; cbin++) {

        cv1HFpm_pT[cbin] = new TCanvas(Form("cv1HFpm_pT_cent%d-%d",centBins[cbin],centBins[cbin+1]),"cv1HFpm_pT",1100,850);
        cv1HFpm_pT[cbin]->Divide(4,3,0,0);

        for (int ebin = 0; ebin<netabins; ebin++) {
            TPad * padv1HFpm_pT = (TPad *) cv1HFpm_pT[cbin]->cd(ebin+1);
            if (gridlines) padv1HFpm_pT->SetGrid();
            if (ebin == 3 || ebin == 7 || ebin == 11) padv1HFpm_pT->SetRightMargin(0.02);
            TH1D * hv1HFpm_pT = (TH1D *) hv1HFpm_pT_tmp->Clone(Form("hv1HFpm_pT_%c_%d",cbin,ebin));
            // hv1HFpm_pT->GetYaxis()->SetRangeUser(-0.14-0.1*cbin,0.14+0.1*cbin);
            if (ebin == 0 || ebin == 4) {
                hv1HFpm_pT->GetYaxis()->CenterTitle();
                hv1HFpm_pT->GetYaxis()->SetTitleSize(0.07);
                hv1HFpm_pT->GetYaxis()->SetTitleOffset(1.34);
                hv1HFpm_pT->GetYaxis()->SetLabelSize(0.06);
            } else if (ebin == 8) {
                hv1HFpm_pT->GetXaxis()->SetTitleSize(0.06);
                hv1HFpm_pT->GetXaxis()->SetTitleOffset(1.14);
                hv1HFpm_pT->GetYaxis()->CenterTitle();
                hv1HFpm_pT->GetYaxis()->SetTitleSize(0.06);
                hv1HFpm_pT->GetYaxis()->SetTitleOffset(1.48);
                hv1HFpm_pT->GetYaxis()->SetLabelSize(0.05);
            } else if (ebin>=9) {
                hv1HFpm_pT->GetXaxis()->SetTitleSize(0.07);
                hv1HFpm_pT->GetXaxis()->SetTitleOffset(1.00);
                hv1HFpm_pT->GetXaxis()->SetLabelSize(0.06);
                hv1HFpm_pT->GetXaxis()->SetLabelOffset(0.005);
            }
            hv1HFpm_pT->Draw();
            v1p_pt[anal][cbin][ebin]->Draw("same");
            v1m_pt[anal][cbin][ebin]->Draw("same");

            TPaveText * txv1HFpm_pT;
            if (ebin == 0 || ebin == 4) txv1HFpm_pT = new TPaveText(0.23, 0.05, 0.62, 0.15,"NDC");
            else if (ebin == 8) txv1HFpm_pT = new TPaveText(0.23, 0.19, 0.62, 0.29,"NDC");
            else if (ebin>=9) txv1HFpm_pT = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
            else txv1HFpm_pT = new TPaveText(0.05, 0.05, 0.50, 0.15,"NDC");
            SetTPaveTxt(txv1HFpm_pT, 18);
            txv1HFpm_pT->AddText(Form("%0.1f < #eta < %0.1f",etabins[ebin],etabins[ebin+1]));
            txv1HFpm_pT->Draw();
        }
        cv1HFpm_pT[cbin]->cd(1);
        TPaveText * txv1HFpm_pT_cent = new TPaveText(0.24, 0.84, 0.38, 0.95,"NDC");
        SetTPaveTxt(txv1HFpm_pT_cent, 18);
        txv1HFpm_pT_cent->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1HFpm_pT_cent->Draw();

        cv1HFpm_pT[cbin]->cd(2);
        TLegend * legv1HFpm_pT = new TLegend(0.07, 0.74, 0.28, 0.94);
        SetLegend(legv1HFpm_pT, 18);
        legv1HFpm_pT->AddEntry(v1p_pt[anal][cbin][0],"HF+","p");
        legv1HFpm_pT->AddEntry(v1m_pt[anal][cbin][0],"HF-","p");
        legv1HFpm_pT->Draw();

        cv1HFpm_pT[cbin]->Print(Form("plots/diffv1/diffv1_pT/diff%s/v1_pm_pT_%s_cent%d-%d.png",AnalNames[anal].data(),AnalNames[anal].data(),centBins[cbin],centBins[cbin+1]),"png");
        if (close_plots) cv1HFpm_pT[cbin]->Close();
    }



    // differential v1(pT) using HF+/- with all centralities on the same plot

    TCanvas * cv1HFpm_pTCent;
    TH1D * hv1HFpm_pTCent_tmp = new TH1D("hv1HFpm_pTCent_tmp", "", 40, 0, 12);
    hv1HFpm_pTCent_tmp->SetTitle("");
    hv1HFpm_pTCent_tmp->SetStats(0);
    hv1HFpm_pTCent_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1HFpm_pTCent_tmp->SetYTitle("v_{1}");
    hv1HFpm_pTCent_tmp->GetYaxis()->SetRangeUser(-0.26,0.26);
    hv1HFpm_pTCent_tmp->SetNdivisions(509);

    cv1HFpm_pTCent = new TCanvas("cv1HFpm_pTCent","cv1HFpm_pTCent",1100,850);
    cv1HFpm_pTCent->Divide(4,3,0,0);

    TLegend * legv1HFpm_pTCent_0 = new TLegend(0.23, 0.65, 0.45, 0.94);
    SetLegend(legv1HFpm_pTCent_0, 16);

    for (int ebin = 0; ebin<netabins; ebin++) {
        TPad * padv1HFpm_pTCent = (TPad *) cv1HFpm_pTCent->cd(ebin+1);
        if (gridlines) padv1HFpm_pTCent->SetGrid();
        if (ebin == 3 || ebin == 7 || ebin == 11) padv1HFpm_pTCent->SetRightMargin(0.02);
        TH1D * hv1HFpm_pTCent = (TH1D *) hv1HFpm_pTCent_tmp->Clone(Form("hv1HFpm_pTCent_%d",ebin));
        if (ebin == 0 || ebin == 4) {
            hv1HFpm_pTCent->GetYaxis()->CenterTitle();
            hv1HFpm_pTCent->GetYaxis()->SetTitleSize(0.07);
            hv1HFpm_pTCent->GetYaxis()->SetTitleOffset(1.34);
            hv1HFpm_pTCent->GetYaxis()->SetLabelSize(0.06);
        } else if (ebin == 8) {
            hv1HFpm_pTCent->GetXaxis()->SetTitleSize(0.06);
            hv1HFpm_pTCent->GetXaxis()->SetTitleOffset(1.14);
            hv1HFpm_pTCent->GetYaxis()->CenterTitle();
            hv1HFpm_pTCent->GetYaxis()->SetTitleSize(0.06);
            hv1HFpm_pTCent->GetYaxis()->SetTitleOffset(1.48);
            hv1HFpm_pTCent->GetYaxis()->SetLabelSize(0.05);
        } else if (ebin>=9) {
            hv1HFpm_pTCent->GetXaxis()->SetTitleSize(0.07);
            hv1HFpm_pTCent->GetXaxis()->SetTitleOffset(1.00);
            hv1HFpm_pTCent->GetXaxis()->SetLabelSize(0.06);
            hv1HFpm_pTCent->GetXaxis()->SetLabelOffset(0.005);
        }
        hv1HFpm_pTCent->Draw();
        TH1D * v1p_pt_tmp[ncentbins];
        TH1D * v1m_pt_tmp[ncentbins];
        for (int cbin = 0; cbin<5; cbin++) {
            v1p_pt_tmp[cbin] = (TH1D *) v1p_pt[anal][cbin][ebin]->Clone(Form("v1p_pt_tmp%d",cbin));
            v1m_pt_tmp[cbin] = (TH1D *) v1m_pt[anal][cbin][ebin]->Clone(Form("v1m_pt_tmp%d",cbin));
            v1p_pt_tmp[cbin]->SetMarkerStyle(centMarkerStyle[cbin]);
            v1p_pt_tmp[cbin]->SetMarkerSize(centMarkerSize[cbin]);
            v1m_pt_tmp[cbin]->SetMarkerStyle(centMarkerStyle[cbin]);
            v1m_pt_tmp[cbin]->SetMarkerSize(centMarkerSize[cbin]);
            v1p_pt_tmp[cbin]->Draw("same");
            v1m_pt_tmp[cbin]->Draw("same");
            if (ebin == 0) {
                legv1HFpm_pTCent_0->AddEntry(v1m_pt_tmp[cbin],Form("%d-%d%%",centBins[cbin],centBins[cbin+1]),"p");
            }
        }

        TPaveText * txv1HFpm_pTCent;
        if (ebin == 0 || ebin == 4) txv1HFpm_pTCent = new TPaveText(0.23, 0.05, 0.62, 0.15,"NDC");
        else if (ebin == 8) txv1HFpm_pTCent = new TPaveText(0.23, 0.19, 0.62, 0.29,"NDC");
        else if (ebin>=9) txv1HFpm_pTCent = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
        else txv1HFpm_pTCent = new TPaveText(0.05, 0.05, 0.50, 0.15,"NDC");
        SetTPaveTxt(txv1HFpm_pTCent, 18);
        txv1HFpm_pTCent->AddText(Form("%0.1f < #eta < %0.1f",etabins[ebin],etabins[ebin+1]));
        txv1HFpm_pTCent->Draw();
    }

    cv1HFpm_pTCent->cd(1);
    legv1HFpm_pTCent_0->Draw();

    cv1HFpm_pTCent->cd(2);
    TLegend * legv1HFpm_pTCent_1 = new TLegend(0.07, 0.74, 0.28, 0.94);
    SetLegend(legv1HFpm_pTCent_1, 18);
    legv1HFpm_pTCent_1->AddEntry(v1p_pt[anal][0][0],"HF+","p");
    legv1HFpm_pTCent_1->AddEntry(v1m_pt[anal][0][0],"HF-","p");
    legv1HFpm_pTCent_1->Draw();

    cv1HFpm_pTCent->Print(Form("plots/diffv1/diffv1_pT/diff%s/v1_pm_pT_CentScan_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFpm_pTCent->Close();



    // differential v1(pT) using HF odd/even for each centrality bin
    anal = 7;

    TCanvas * cv1HFoddeven_pT[ncentbins];
    TH1D * hv1HFoddeven_pT_tmp = new TH1D("hv1HFoddeven_pT_tmp", "", 40, 0, 12);
    hv1HFoddeven_pT_tmp->SetTitle("");
    hv1HFoddeven_pT_tmp->SetStats(0);
    hv1HFoddeven_pT_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1HFoddeven_pT_tmp->SetYTitle("v_{1}");
    hv1HFoddeven_pT_tmp->GetYaxis()->SetRangeUser(-0.3, 0.3);
    hv1HFoddeven_pT_tmp->SetNdivisions(509);
    for (int cbin = 0; cbin<ncentbins; cbin++) {

        cv1HFoddeven_pT[cbin] = new TCanvas(Form("cv1HFoddeven_pT_cent%d-%d",centBins[cbin],centBins[cbin+1]),"cv1HFoddeven_pT",1100,850);
        cv1HFoddeven_pT[cbin]->Divide(4,3,0,0);

        for (int ebin = 0; ebin<netabins; ebin++) {
            TPad * padv1HFoddeven_pT = (TPad *) cv1HFoddeven_pT[cbin]->cd(ebin+1);
            if (gridlines) padv1HFoddeven_pT->SetGrid();
            if (ebin == 3 || ebin == 7 || ebin == 11) padv1HFoddeven_pT->SetRightMargin(0.02);
            TH1D * hv1HFoddeven_pT = (TH1D *) hv1HFoddeven_pT_tmp->Clone(Form("hv1HFoddeven_pT_%c_%d",cbin,ebin));
            // hv1HFoddeven_pT->GetYaxis()->SetRangeUser(-0.14-0.06*cbin,0.14+0.06*cbin);
            if (ebin == 0 || ebin == 4) {
                hv1HFoddeven_pT->GetYaxis()->CenterTitle();
                hv1HFoddeven_pT->GetYaxis()->SetTitleSize(0.07);
                hv1HFoddeven_pT->GetYaxis()->SetTitleOffset(1.34);
                hv1HFoddeven_pT->GetYaxis()->SetLabelSize(0.06);
            } else if (ebin == 8) {
                hv1HFoddeven_pT->GetXaxis()->SetTitleSize(0.06);
                hv1HFoddeven_pT->GetXaxis()->SetTitleOffset(1.14);
                hv1HFoddeven_pT->GetYaxis()->CenterTitle();
                hv1HFoddeven_pT->GetYaxis()->SetTitleSize(0.06);
                hv1HFoddeven_pT->GetYaxis()->SetTitleOffset(1.48);
                hv1HFoddeven_pT->GetYaxis()->SetLabelSize(0.05);
            } else if (ebin>=9) {
                hv1HFoddeven_pT->GetXaxis()->SetTitleSize(0.07);
                hv1HFoddeven_pT->GetXaxis()->SetTitleOffset(1.00);
                hv1HFoddeven_pT->GetXaxis()->SetLabelSize(0.06);
                hv1HFoddeven_pT->GetXaxis()->SetLabelOffset(0.005);
            }
            hv1HFoddeven_pT->Draw();
            v1odd_pt[anal][cbin][ebin]->Draw("same");
            v1even_pt[anal][cbin][ebin]->Draw("same");

            TPaveText * txv1HFoddeven_pT;
            if (ebin == 0 || ebin == 4) txv1HFoddeven_pT = new TPaveText(0.23, 0.05, 0.62, 0.15,"NDC");
            else if (ebin == 8) txv1HFoddeven_pT = new TPaveText(0.23, 0.19, 0.62, 0.29,"NDC");
            else if (ebin>=9) txv1HFoddeven_pT = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
            else txv1HFoddeven_pT = new TPaveText(0.05, 0.05, 0.50, 0.15,"NDC");
            SetTPaveTxt(txv1HFoddeven_pT, 18);
            txv1HFoddeven_pT->AddText(Form("%0.1f < #eta < %0.1f",etabins[ebin],etabins[ebin+1]));
            txv1HFoddeven_pT->Draw();
        }
        cv1HFoddeven_pT[cbin]->cd(1);
        TPaveText * txv1HFoddeven_pT_cent = new TPaveText(0.24, 0.84, 0.38, 0.95,"NDC");
        SetTPaveTxt(txv1HFoddeven_pT_cent, 18);
        txv1HFoddeven_pT_cent->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1HFoddeven_pT_cent->Draw();

        cv1HFoddeven_pT[cbin]->cd(2);
        TLegend * legv1HFoddeven_pT = new TLegend(0.07, 0.74, 0.28, 0.94);
        SetLegend(legv1HFoddeven_pT, 18);
        legv1HFoddeven_pT->AddEntry(v1odd_pt[anal][cbin][0],"v_{1}^{odd}","p");
        legv1HFoddeven_pT->AddEntry(v1even_pt[anal][cbin][0],"v_{1}^{even}","p");
        legv1HFoddeven_pT->Draw();

        cv1HFoddeven_pT[cbin]->Print(Form("plots/diffv1/diffv1_pT/diff%s/v1_oddeven_pT_%s_cent%d-%d.png",AnalNames[anal].data(),AnalNames[anal].data(),centBins[cbin],centBins[cbin+1]),"png");
        if (close_plots) cv1HFoddeven_pT[cbin]->Close();
    }



    // differential v1(pT) using HF odd/even with all centralities on the same plot

    TCanvas * cv1HFoddeven_pTCent;
    TH1D * hv1HFoddeven_pTCent_tmp = new TH1D("hv1HFoddeven_pTCent_tmp", "", 40, 0, 12);
    hv1HFoddeven_pTCent_tmp->SetTitle("");
    hv1HFoddeven_pTCent_tmp->SetStats(0);
    hv1HFoddeven_pTCent_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1HFoddeven_pTCent_tmp->SetYTitle("v_{1}");
    hv1HFoddeven_pTCent_tmp->GetYaxis()->SetRangeUser(-0.2,0.2);
    hv1HFoddeven_pTCent_tmp->SetNdivisions(509);

    cv1HFoddeven_pTCent = new TCanvas("cv1HFoddeven_pTCent","cv1HFoddeven_pTCent",1100,850);
    cv1HFoddeven_pTCent->Divide(4,3,0,0);

    TLegend * legv1HFoddeven_pTCent_0 = new TLegend(0.23, 0.65, 0.45, 0.94);
    SetLegend(legv1HFoddeven_pTCent_0, 16);

    for (int ebin = 0; ebin<netabins; ebin++) {
        TPad * padv1HFoddeven_pTCent = (TPad *) cv1HFoddeven_pTCent->cd(ebin+1);
        if (gridlines) padv1HFoddeven_pTCent->SetGrid();
        if (ebin == 3 || ebin == 7 || ebin == 11) padv1HFoddeven_pTCent->SetRightMargin(0.02);
        TH1D * hv1HFoddeven_pTCent = (TH1D *) hv1HFoddeven_pTCent_tmp->Clone(Form("hv1HFoddeven_pTCent_%d",ebin));
        if (ebin == 0 || ebin == 4) {
            hv1HFoddeven_pTCent->GetYaxis()->CenterTitle();
            hv1HFoddeven_pTCent->GetYaxis()->SetTitleSize(0.07);
            hv1HFoddeven_pTCent->GetYaxis()->SetTitleOffset(1.34);
            hv1HFoddeven_pTCent->GetYaxis()->SetLabelSize(0.06);
        } else if (ebin == 8) {
            hv1HFoddeven_pTCent->GetXaxis()->SetTitleSize(0.06);
            hv1HFoddeven_pTCent->GetXaxis()->SetTitleOffset(1.14);
            hv1HFoddeven_pTCent->GetYaxis()->CenterTitle();
            hv1HFoddeven_pTCent->GetYaxis()->SetTitleSize(0.06);
            hv1HFoddeven_pTCent->GetYaxis()->SetTitleOffset(1.48);
            hv1HFoddeven_pTCent->GetYaxis()->SetLabelSize(0.05);
        } else if (ebin>=9) {
            hv1HFoddeven_pTCent->GetXaxis()->SetTitleSize(0.07);
            hv1HFoddeven_pTCent->GetXaxis()->SetTitleOffset(1.00);
            hv1HFoddeven_pTCent->GetXaxis()->SetLabelSize(0.06);
            hv1HFoddeven_pTCent->GetXaxis()->SetLabelOffset(0.005);
        }
        hv1HFoddeven_pTCent->Draw();
        TH1D * v1odd_pt_tmp[ncentbins];
        TH1D * v1even_pt_tmp[ncentbins];
        for (int cbin = 0; cbin<5; cbin++) {
            v1odd_pt_tmp[cbin] = (TH1D *) v1odd_pt[anal][cbin][ebin]->Clone(Form("v1odd_pt_tmp%d",cbin));
            v1even_pt_tmp[cbin] = (TH1D *) v1even_pt[anal][cbin][ebin]->Clone(Form("v1even_pt_tmp%d",cbin));
            v1odd_pt_tmp[cbin]->SetMarkerStyle(centMarkerStyle[cbin]);
            v1odd_pt_tmp[cbin]->SetMarkerSize(centMarkerSize[cbin]);
            v1even_pt_tmp[cbin]->SetMarkerStyle(centMarkerStyle[cbin]);
            v1even_pt_tmp[cbin]->SetMarkerSize(centMarkerSize[cbin]);
            v1odd_pt_tmp[cbin]->Draw("same");
            v1even_pt_tmp[cbin]->Draw("same");
            if (ebin == 0) {
                legv1HFoddeven_pTCent_0->AddEntry(v1even_pt_tmp[cbin],Form("%d-%d%%",centBins[cbin],centBins[cbin+1]),"p");
            }
        }

        TPaveText * txv1HFoddeven_pTCent;
        if (ebin == 0 || ebin == 4) txv1HFoddeven_pTCent = new TPaveText(0.23, 0.05, 0.62, 0.15,"NDC");
        else if (ebin == 8) txv1HFoddeven_pTCent = new TPaveText(0.23, 0.19, 0.62, 0.29,"NDC");
        else if (ebin>=9) txv1HFoddeven_pTCent = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
        else txv1HFoddeven_pTCent = new TPaveText(0.05, 0.05, 0.50, 0.15,"NDC");
        SetTPaveTxt(txv1HFoddeven_pTCent, 18);
        txv1HFoddeven_pTCent->AddText(Form("%0.1f < #eta < %0.1f",etabins[ebin],etabins[ebin+1]));
        txv1HFoddeven_pTCent->Draw();
    }

    cv1HFoddeven_pTCent->cd(1);
    legv1HFoddeven_pTCent_0->Draw();

    cv1HFoddeven_pTCent->cd(2);
    TLegend * legv1HFoddeven_pTCent_1 = new TLegend(0.07, 0.74, 0.28, 0.94);
    SetLegend(legv1HFoddeven_pTCent_1, 18);
    legv1HFoddeven_pTCent_1->AddEntry(v1odd_pt[anal][0][0],"v_{1}^{odd}","p");
    legv1HFoddeven_pTCent_1->AddEntry(v1even_pt[anal][0][0],"v_{1}^{even}","p");
    legv1HFoddeven_pTCent_1->Draw();

    cv1HFoddeven_pTCent->Print(Form("plots/diffv1/diffv1_pT/diff%s/v1_oddeven_pT_CentScan_%s.png",AnalNames[anal].data(),AnalNames[anal].data()),"png");
    if (close_plots) cv1HFoddeven_pTCent->Close();



    // differential v1(pT) using tracker+/- for each centrality bin
    anal = 15;
    if (!fopen(Form("plots/diffv1/diffv1_pT/diff%s",AnalNames[anal].data()),"r")) system(Form("mkdir plots/diffv1/diffv1_pT/diff%s",AnalNames[anal].data()));

    TCanvas * cv1Trackpm_pT[ncentbins];
    TH1D * hv1Trackpm_pT_tmp = new TH1D("hv1Trackpm_pT_tmp", "", 40, 0, 12);
    hv1Trackpm_pT_tmp->SetTitle("");
    hv1Trackpm_pT_tmp->SetStats(0);
    hv1Trackpm_pT_tmp->SetXTitle("p_{T} (GeV/c)");
    hv1Trackpm_pT_tmp->SetYTitle("v_{1}");
    hv1Trackpm_pT_tmp->GetYaxis()->SetRangeUser(-0.075, 0.3);
    hv1Trackpm_pT_tmp->SetNdivisions(509);
    for (int cbin = 0; cbin<ncentbins; cbin++) {

        cv1Trackpm_pT[cbin] = new TCanvas(Form("cv1Trackpm_pT_cent%d-%d",centBins[cbin],centBins[cbin+1]),"cv1Trackpm_pT",1000,800);
        cv1Trackpm_pT[cbin]->Divide(4,3,0,0);

        for (int ebin = 0; ebin<netabins; ebin++) {
            TPad * padv1Trackpm_pT = (TPad *) cv1Trackpm_pT[cbin]->cd(ebin+1);
            if (gridlines) padv1Trackpm_pT->SetGrid();
            if (ebin == 3 || ebin == 7 || ebin == 11) padv1Trackpm_pT->SetRightMargin(0.02);
            TH1D * hv1Trackpm_pT = (TH1D *) hv1Trackpm_pT_tmp->Clone(Form("hv1Trackpm_pT_%c_%d",cbin,ebin));
            // hv1Trackpm_pT->GetYaxis()->SetRangeUser(-0.05-0.03, 0.15+0.03*cbin);
            if (ebin == 0 || ebin == 4) {
                hv1Trackpm_pT->GetYaxis()->CenterTitle();
                hv1Trackpm_pT->GetYaxis()->SetTitleSize(0.07);
                hv1Trackpm_pT->GetYaxis()->SetTitleOffset(1.34);
                hv1Trackpm_pT->GetYaxis()->SetLabelSize(0.06);
            } else if (ebin == 8) {
                hv1Trackpm_pT->GetXaxis()->SetTitleSize(0.06);
                hv1Trackpm_pT->GetXaxis()->SetTitleOffset(1.14);
                hv1Trackpm_pT->GetYaxis()->CenterTitle();
                hv1Trackpm_pT->GetYaxis()->SetTitleSize(0.06);
                hv1Trackpm_pT->GetYaxis()->SetTitleOffset(1.48);
                hv1Trackpm_pT->GetYaxis()->SetLabelSize(0.05);
            } else if (ebin>=9) {
                hv1Trackpm_pT->GetXaxis()->SetTitleSize(0.07);
                hv1Trackpm_pT->GetXaxis()->SetTitleOffset(1.00);
                hv1Trackpm_pT->GetXaxis()->SetLabelSize(0.06);
                hv1Trackpm_pT->GetXaxis()->SetLabelOffset(0.005);
            }
            hv1Trackpm_pT->Draw();
            v1p_pt[anal][cbin][ebin]->Draw("same");
            v1m_pt[anal][cbin][ebin]->Draw("same");

            TPaveText * txv1Trackpm_pT;
            if (ebin == 0 || ebin == 4) txv1Trackpm_pT = new TPaveText(0.23, 0.05, 0.62, 0.15,"NDC");
            else if (ebin == 8) txv1Trackpm_pT = new TPaveText(0.23, 0.19, 0.62, 0.29,"NDC");
            else if (ebin>=9) txv1Trackpm_pT = new TPaveText(0.05, 0.19, 0.50, 0.29,"NDC");
            else txv1Trackpm_pT = new TPaveText(0.05, 0.05, 0.50, 0.15,"NDC");
            SetTPaveTxt(txv1Trackpm_pT, 18);
            txv1Trackpm_pT->AddText(Form("%0.1f < #eta < %0.1f",etabins[ebin],etabins[ebin+1]));
            txv1Trackpm_pT->Draw();
        }
        cv1Trackpm_pT[cbin]->cd(1);
        TPaveText * txv1Trackpm_pT_cent = new TPaveText(0.24, 0.84, 0.38, 0.95,"NDC");
        SetTPaveTxt(txv1Trackpm_pT_cent, 18);
        txv1Trackpm_pT_cent->AddText(Form("%d-%d%%",centBins[cbin],centBins[cbin+1]));
        txv1Trackpm_pT_cent->Draw();

        cv1Trackpm_pT[cbin]->cd(2);
        TLegend * legv1Trackpm_pT = new TLegend(0.07, 0.74, 0.28, 0.94);
        SetLegend(legv1Trackpm_pT, 18);
        legv1Trackpm_pT->AddEntry(v1p_pt[anal][cbin][0],"Track+","p");
        legv1Trackpm_pT->AddEntry(v1m_pt[anal][cbin][0],"Track-","p");
        legv1Trackpm_pT->Draw();

        cv1Trackpm_pT[cbin]->Print(Form("plots/diffv1/diffv1_pT/diff%s/v1_pm_pT_%s_cent%d-%d.png",AnalNames[anal].data(),AnalNames[anal].data(),centBins[cbin],centBins[cbin+1]),"png");
        if (close_plots) cv1Trackpm_pT[cbin]->Close();
    }


}
