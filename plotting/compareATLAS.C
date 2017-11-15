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

Bool_t close_plots = kFALSE;
Bool_t gridlines = kFALSE;

# include "style.h"

# include "../../published_results/PhysRevC86_014907.h"

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

TH1D * runParms[nanals];

void compareATLAS()
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
        }
    }

    //-- retrieve ATLAS results
    TH1D * hATLASv1even[ncentbins];
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        hATLASv1even[cbin] = new TH1D(Form("ATLASv1even_c%d",cbin),Form("ATLASv1even_c%d",cbin), 85, 0.45, 9.05);
        hATLASv1even[cbin]->SetStats(kFALSE);
        hATLASv1even[cbin]->SetMarkerColor(kBlack);
        hATLASv1even[cbin]->SetMarkerStyle(27);
        hATLASv1even[cbin]->SetMarkerSize(0.2);
        hATLASv1even[cbin]->SetLineColor(kGray+1);
        hATLASv1even[cbin]->SetLineWidth(2);
        hATLASv1even[cbin]->SetFillColor(kGray+1);
        hATLASv1even[cbin]->SetFillStyle(1001);
    }
    for (int pbin = 1; pbin<=85; pbin++) {
        hATLASv1even[0]->SetBinContent(pbin,v1ATLAS_c00to05[pbin-1]);
        hATLASv1even[0]->SetBinError(pbin,v1ATLAS_c00to05_err[pbin-1]);
        hATLASv1even[1]->SetBinContent(pbin,v1ATLAS_c05to10[pbin-1]);
        hATLASv1even[1]->SetBinError(pbin,v1ATLAS_c05to10_err[pbin-1]);
        hATLASv1even[2]->SetBinContent(pbin,v1ATLAS_c10to20[pbin-1]);
        hATLASv1even[2]->SetBinError(pbin,v1ATLAS_c10to20_err[pbin-1]);
        hATLASv1even[3]->SetBinContent(pbin,v1ATLAS_c20to30[pbin-1]);
        hATLASv1even[3]->SetBinError(pbin,v1ATLAS_c20to30_err[pbin-1]);
        hATLASv1even[4]->SetBinContent(pbin,v1ATLAS_c30to40[pbin-1]);
        hATLASv1even[4]->SetBinError(pbin,v1ATLAS_c30to40_err[pbin-1]);
        hATLASv1even[5]->SetBinContent(pbin,v1ATLAS_c40to50[pbin-1]);
        hATLASv1even[5]->SetBinError(pbin,v1ATLAS_c40to50_err[pbin-1]);
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
        }
    }


    //-- make plots
    if (!fopen("plots","r")) system("mkdir plots");
    if (!fopen("plots/Comparisons","r")) system("mkdir plots/Comparisons");
    int anal = 15; // choice of analysis

    TCanvas * cATLAScompare = new TCanvas("cATLAScompare","cATLAScompare",1250,430);
    cATLAScompare->Divide(4,1,0,0);
    TH1D * hATLAScompare_tmp = new TH1D("hATLAScompare_tmp", "", 40, 0, 12);
    hATLAScompare_tmp->SetTitle("");
    hATLAScompare_tmp->SetStats(kFALSE);
    hATLAScompare_tmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hATLAScompare_tmp->GetYaxis()->SetTitle("v_{1}^{even}");
    hATLAScompare_tmp->GetXaxis()->CenterTitle(kTRUE);
    hATLAScompare_tmp->GetXaxis()->SetRangeUser(0, 9.2);
    hATLAScompare_tmp->GetYaxis()->SetRangeUser(-0.025, 0.2);
    hATLAScompare_tmp->GetYaxis()->SetDecimals(2);
    for (int cbin = 0; cbin<4; cbin++) {
        TPad * padATLAScompare = (TPad *) cATLAScompare->cd(cbin+1);
        if (gridlines) padATLAScompare->SetGrid();
        if (cbin == 3) padATLAScompare->SetRightMargin(0.02);
        padATLAScompare->SetTopMargin(0.07);
        TH1D * hATLAScompare = (TH1D *) hATLAScompare_tmp->Clone(Form("hATLAScompare_tmp_%d",cbin));
        if (cbin == 0) {
            hATLAScompare->GetXaxis()->SetTitleSize(0.05);
            hATLAScompare->GetXaxis()->SetTitleOffset(1.13);
            hATLAScompare->GetXaxis()->SetLabelSize(0.05);
            hATLAScompare->GetXaxis()->SetLabelOffset(0.016);

            hATLAScompare->GetYaxis()->SetTitleSize(0.07);
            hATLAScompare->GetYaxis()->SetTitleOffset(1.25);
            hATLAScompare->GetYaxis()->SetLabelSize(0.05);
            hATLAScompare->GetYaxis()->SetLabelOffset(0.010);
        } else {
            hATLAScompare->GetXaxis()->SetTitleSize(0.07);
            hATLAScompare->GetXaxis()->SetTitleOffset(0.96);
            hATLAScompare->GetXaxis()->SetLabelSize(0.06);
            hATLAScompare->GetXaxis()->SetLabelOffset(0.007);
        }
        hATLAScompare->Draw();
        hATLASv1even[cbin+2]->Draw("same E3");
        hATLASv1even[cbin+2]->Draw("same");
        v1odd_pt[anal][cbin]->Draw("same");

        TPaveText * txATLAScompare;
        if (cbin == 0) txATLAScompare = new TPaveText(0.23, 0.66, 0.41, 0.72,"NDC");
        else txATLAScompare = new TPaveText(0.06, 0.82, 0.26, 0.88,"NDC");
        SetTPaveTxt(txATLAScompare, 18);
        txATLAScompare->AddText(Form("%d-%d%%",centBins[cbin+1],centBins[cbin+2]));
        txATLAScompare->Draw();
    }
    cATLAScompare->cd(1);
    TLegend * legATLAScompare = new TLegend(0.23, 0.73, 0.58, 0.88);
    SetLegend(legATLAScompare, 18);
    legATLAScompare->AddEntry(v1odd_pt[anal][0],"CMS #sqrt{s_{NN}} = 5.02 TeV","p");
    legATLAScompare->AddEntry(hATLASv1even[0],"ATLAS #sqrt{s_{NN}} = 2.76 TeV","lp");
    legATLAScompare->Draw();

    TPaveText * txATLAScompare_1 = new TPaveText(0.18, 0.94, 0.58, 1.0,"NDC");
    SetTPaveTxt(txATLAScompare_1, 18);
    txATLAScompare_1->AddText("#bf{CMS} #it{Preliminary}");
    txATLAScompare_1->Draw();

    cATLAScompare->Print("plots/Comparisons/ATLAScompare.png","png");
    if (close_plots) cATLAScompare->Close();

}
