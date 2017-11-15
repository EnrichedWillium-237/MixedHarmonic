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

static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static const int netabins = 12;
static const double etabins[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,
                     0.4,  0.8,  1.2,  1.6, 2.0,  2.4};
static const int ncentbins = 8;
static const int centBins[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
static const int nanals = 16;
string AnalNames[] = {
    "v1SP",   "v1SP_mid",   "v1SP_102",   "v1SP_106",   "v1SP_110",   "v1SP_114",   "v1SP_118",   "v1SP_122",
    "v1SPmc", "v1SPmc_mid", "v1SPmc_102", "v1SPmc_106", "v1SPmc_110", "v1SPmc_114", "v1SPmc_118", "v1SPmc_122"
};

using namespace std;

static const bool subEvt2 = kTRUE;

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

TH1D * runParms[nanals];

# include "ErrCalcIntV1.C"

TGraphErrors * GetV1IntPt( string anal, string tag, int cbin, double etamin, double etamax, TGraphErrors * &g_p, TGraphErrors * &g_m );
TGraphErrors * GetV1IntEta( string anal, string tag, int cbin, double etamin, double etamax, TGraphErrors * &g_p, TGraphErrors * &g_m );

void GraphToHist( TGraphErrors * gin, TH1D * hout ) {
    int num = gin->GetN();
    Double_t x[100], y[100], yerr[100];
    for (int i = 0; i<num; i++) {
        gin->GetPoint(i, x[i], y[i]);
        yerr[i] = gin->GetErrorY(i);
        hout->SetBinContent(i+1, y[i]);
        hout->SetBinError(i+1, yerr[i]);
    }
}

void AnalyzeIntV1()
{

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    double etamin = -2.4;
    double etamax = 2.4;

    double ptmin = 0.3;
    double ptmax = 3.0;

    int minanal = 0;
    int maxanal = nanals;

    cout << "\nBeginning integral v1 analyzer...\n" << endl;

    for (int i = minanal; i<maxanal; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            v1p_pt[i][cbin] = new TH1D(Form("v1p_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);
            v1m_pt[i][cbin] = new TH1D(Form("v1m_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);
            v1odd_pt[i][cbin] = new TH1D(Form("v1odd_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);

            v112p_pt[i][cbin] = new TH1D(Form("v112p_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);
            v112m_pt[i][cbin] = new TH1D(Form("v112m_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);
            v112odd_pt[i][cbin] = new TH1D(Form("v112odd_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);

            v123p_pt[i][cbin] = new TH1D(Form("v123p_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);
            v123m_pt[i][cbin] = new TH1D(Form("v123m_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);
            v123odd_pt[i][cbin] = new TH1D(Form("v123odd_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);

            v1p_eta[i][cbin] = new TH1D(Form("v1p_eta_%s_%d",AnalNames[i].data(),cbin), "", netabins, etabins);
            v1m_eta[i][cbin] = new TH1D(Form("v1m_eta_%s_%d",AnalNames[i].data(),cbin), "", netabins, etabins);
            v1odd_eta[i][cbin] = new TH1D(Form("v1odd_eta_%s_%d",AnalNames[i].data(),cbin), "", netabins, etabins);

            v112p_eta[i][cbin] = new TH1D(Form("v112p_eta_%s_%d",AnalNames[i].data(),cbin), "", netabins, etabins);
            v112m_eta[i][cbin] = new TH1D(Form("v112m_eta_%s_%d",AnalNames[i].data(),cbin), "", netabins, etabins);
            v112odd_eta[i][cbin] = new TH1D(Form("v112odd_eta_%s_%d",AnalNames[i].data(),cbin), "", netabins, etabins);

            v123p_eta[i][cbin] = new TH1D(Form("v123p_eta_%s_%d",AnalNames[i].data(),cbin), "", netabins, etabins);
            v123m_eta[i][cbin] = new TH1D(Form("v123m_eta_%s_%d",AnalNames[i].data(),cbin), "", netabins, etabins);
            v123odd_eta[i][cbin] = new TH1D(Form("v123odd_eta_%s_%d",AnalNames[i].data(),cbin), "", netabins, etabins);
        }
        TFile * tfParms = new TFile(Form("../outputs/raw_outputs/results/%s.root",AnalNames[i].data()),"read");
        runParms[i] = (TH1D *) tfParms->Get(Form("%d-%d/runParms",centBins[0],centBins[1]));
    }


    //-- make histograms
    string tag = "3sub";
    if (subEvt2) tag = "2sub";
    for (int i = minanal; i<maxanal; i++) {
        cout << "  processing file: " << AnalNames[i].data() << endl;
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TGraphErrors * gA1;
            TGraphErrors * gB1;
            TGraphErrors * gAB1 = GetV1IntPt(Form("%s",AnalNames[i].data()), tag.data(), cbin, etamin, etamax, gA1, gB1);
            GraphToHist(gA1, v1p_pt[i][cbin]);
            GraphToHist(gB1, v1m_pt[i][cbin]);
            GraphToHist(gAB1, v1odd_pt[i][cbin]);
            v1even_pt[i][cbin] = (TH1D *) v1p_pt[i][cbin]->Clone(Form("v1even_pt_%s_%d",AnalNames[i].data(),cbin));
            v1even_pt[i][cbin]->Add(v1m_pt[i][cbin],-1);
            v1even_pt[i][cbin]->Scale(0.5);

            TGraphErrors * gA112;
            TGraphErrors * gB112;
            TGraphErrors * gAB112 = GetV1IntPt(Form("%s",AnalNames[i].data()), "112", cbin, etamin, etamax, gA112, gB112);
            GraphToHist(gA112, v112p_pt[i][cbin]);
            GraphToHist(gB112, v112m_pt[i][cbin]);
            GraphToHist(gAB112, v112odd_pt[i][cbin]);
            v112even_pt[i][cbin] = (TH1D *) v112p_pt[i][cbin]->Clone(Form("v112even_pt_%s_%d",AnalNames[i].data(),cbin));
            v112even_pt[i][cbin]->Add(v112m_pt[i][cbin],-1);
            v112even_pt[i][cbin]->Scale(0.5);

            TGraphErrors * gA123;
            TGraphErrors * gB123;
            TGraphErrors * gAB123 = GetV1IntPt(Form("%s",AnalNames[i].data()), "123", cbin, etamin, etamax, gA123, gB123);
            GraphToHist(gA123, v123p_pt[i][cbin]);
            GraphToHist(gB123, v123m_pt[i][cbin]);
            GraphToHist(gAB123, v123odd_pt[i][cbin]);
            v123even_pt[i][cbin] = (TH1D *) v123p_pt[i][cbin]->Clone(Form("v123even_pt_%s_%d",AnalNames[i].data(),cbin));
            v123even_pt[i][cbin]->Add(v123m_pt[i][cbin],-1);
            v123even_pt[i][cbin]->Scale(0.5);
        }
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TGraphErrors * gA1;
            TGraphErrors * gB1;
            TGraphErrors * gAB1 = GetV1IntEta(Form("%s",AnalNames[i].data()), tag.data(), cbin, ptmin, ptmax, gA1, gB1);
            GraphToHist(gA1, v1p_eta[i][cbin]);
            GraphToHist(gB1, v1m_eta[i][cbin]);
            GraphToHist(gAB1, v1odd_eta[i][cbin]);
            v1even_eta[i][cbin] = (TH1D *) v1p_eta[i][cbin]->Clone(Form("v1even_eta_%s_%d",AnalNames[i].data(),cbin));
            v1even_eta[i][cbin]->Add(v1m_eta[i][cbin],-1);
            v1even_eta[i][cbin]->Scale(0.5);

            TGraphErrors * gA112;
            TGraphErrors * gB112;
            TGraphErrors * gAB112 = GetV1IntEta(Form("%s",AnalNames[i].data()), "112", cbin, ptmin, ptmax, gA112, gB112);
            GraphToHist(gA112, v112p_eta[i][cbin]);
            GraphToHist(gB112, v112m_eta[i][cbin]);
            GraphToHist(gAB112, v112odd_eta[i][cbin]);
            v112even_eta[i][cbin] = (TH1D *) v112p_eta[i][cbin]->Clone(Form("v112even_eta_%s_%d",AnalNames[i].data(),cbin));
            v112even_eta[i][cbin]->Add(v112m_eta[i][cbin],-1);
            v112even_eta[i][cbin]->Scale(0.5);

            TGraphErrors * gA123;
            TGraphErrors * gB123;
            TGraphErrors * gAB123 = GetV1IntEta(Form("%s",AnalNames[i].data()), "123", cbin, ptmin, ptmax, gA123, gB123);
            GraphToHist(gA123, v123p_eta[i][cbin]);
            GraphToHist(gB123, v123m_eta[i][cbin]);
            GraphToHist(gAB123, v123odd_eta[i][cbin]);
            v123even_eta[i][cbin] = (TH1D *) v123p_eta[i][cbin]->Clone(Form("v123even_eta_%s_%d",AnalNames[i].data(),cbin));
            v123even_eta[i][cbin]->Add(v123m_eta[i][cbin],-1);
            v123even_eta[i][cbin]->Scale(0.5);
        }
    }
    if (!fopen("../outputs","r")) system("mkdir ../outputs");
    if (!fopen("../outputs/final_outputs","r")) system("mkdir ../outputs/final_outputs");
    TFile * tfout;
    if (subEvt2) tfout = new TFile("../outputs/final_outputs/v1Int_2sub.root","recreate");
    else tfout = new TFile("../outputs/final_outputs/v1Int.root","recreate");
    for (int i = minanal; i<maxanal; i++) {
        TDirectory * tdir = (TDirectory *) tfout->mkdir(Form("%s",AnalNames[i].data()));
        TDirectory * tdpt = (TDirectory *) tdir->mkdir("v1_pt");
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TDirectory * tdirCent = (TDirectory *) tdpt->mkdir(Form("%d-%d",centBins[cbin],centBins[cbin+1]));
            tdirCent->cd();
            v1p_pt[i][cbin]->Write();
            v1m_pt[i][cbin]->Write();
            v1odd_pt[i][cbin]->Write();
            v1even_pt[i][cbin]->Write();
            v112p_pt[i][cbin]->Write();
            v112m_pt[i][cbin]->Write();
            v112odd_pt[i][cbin]->Write();
            v112even_pt[i][cbin]->Write();
            v123p_pt[i][cbin]->Write();
            v123m_pt[i][cbin]->Write();
            v123odd_pt[i][cbin]->Write();
            v123even_pt[i][cbin]->Write();
        }
        TDirectory * tdeta = (TDirectory *) tdir->mkdir("v1_eta");
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TDirectory * tdirCent = (TDirectory *) tdeta->mkdir(Form("%d-%d",centBins[cbin],centBins[cbin+1]));
            tdirCent->cd();
            v1p_eta[i][cbin]->Write();
            v1m_eta[i][cbin]->Write();
            v1odd_eta[i][cbin]->Write();
            v1even_eta[i][cbin]->Write();
            v112p_eta[i][cbin]->Write();
            v112m_eta[i][cbin]->Write();
            v112odd_eta[i][cbin]->Write();
            v112even_eta[i][cbin]->Write();
            v123p_eta[i][cbin]->Write();
            v123m_eta[i][cbin]->Write();
            v123odd_eta[i][cbin]->Write();
            v123even_eta[i][cbin]->Write();
        }
        tdir->cd();
        runParms[i]->SetBinContent(3, ptmin);
        runParms[i]->SetBinContent(4, ptmax);
        runParms[i]->SetBinContent(5, etamin);
        runParms[i]->SetBinContent(6, etamax);
        runParms[i]->Write();
    }

    cout << "\n...Integral v1 results written out to ../outputs/final_outputs/v1Int.root \n" << endl;
    tfout->Close();

    ErrCalcIntV1();

}


TGraphErrors * GetV1IntPt( string anal, string tag, int cbin, double etamin, double etamax, TGraphErrors * &g_p, TGraphErrors * &g_m ) {

    bool sub2 = kFALSE;
    bool mix112 = kFALSE;
    bool mix123 = kFALSE;
    if (tag == "2sub") sub2 = kTRUE;
    else if (tag == "112") mix112 = kTRUE;
    else if (tag == "123") mix123 = kTRUE;
    TFile * tfin = new TFile(Form("../outputs/raw_outputs/results/%s.root",anal.data()),"read");

    TH1D * centbins = (TH1D *) tfin->Get(Form("%d-%d/centbins",centBins[cbin],centBins[cbin+1]));
    TH2D * ptav = (TH2D *) tfin->Get(Form("%d-%d/ptav_%d",centBins[cbin],centBins[cbin+1],cbin));
    TH2D * ptcnt = (TH2D *) tfin->Get(Form("%d-%d/ptcnt_%d",centBins[cbin],centBins[cbin+1],cbin));
    TH2D * badcnt = (TH2D *) tfin->Get(Form("%d-%d/badcnt_%d",centBins[cbin],centBins[cbin+1],cbin));
    TH1D * multTot = (TH1D *) tfin->Get(Form("%d-%d/multTot_%d",centBins[cbin],centBins[cbin+1],cbin));
    TH2D * q1_p;
    TH2D * q1_m;
    TH2D * w1_p;
    TH2D * w1_m;
    if (mix112) {
        q1_p = (TH2D *) tfin->Get(Form("%d-%d/q112_p_%d",centBins[cbin],centBins[cbin+1],cbin));
        q1_m = (TH2D *) tfin->Get(Form("%d-%d/q112_m_%d",centBins[cbin],centBins[cbin+1],cbin));
        w1_p = (TH2D *) tfin->Get(Form("%d-%d/w112_p_%d",centBins[cbin],centBins[cbin+1],cbin));
        w1_m = (TH2D *) tfin->Get(Form("%d-%d/w112_m_%d",centBins[cbin],centBins[cbin+1],cbin));
    } else if (mix123) {
        q1_p = (TH2D *) tfin->Get(Form("%d-%d/q123_p_%d",centBins[cbin],centBins[cbin+1],cbin));
        q1_m = (TH2D *) tfin->Get(Form("%d-%d/q123_m_%d",centBins[cbin],centBins[cbin+1],cbin));
        w1_p = (TH2D *) tfin->Get(Form("%d-%d/w123_p_%d",centBins[cbin],centBins[cbin+1],cbin));
        w1_m = (TH2D *) tfin->Get(Form("%d-%d/w123_m_%d",centBins[cbin],centBins[cbin+1],cbin));
    } else {
        q1_p = (TH2D *) tfin->Get(Form("%d-%d/q1_p_%d",centBins[cbin],centBins[cbin+1],cbin));
        q1_m = (TH2D *) tfin->Get(Form("%d-%d/q1_m_%d",centBins[cbin],centBins[cbin+1],cbin));
        w1_p = (TH2D *) tfin->Get(Form("%d-%d/w1_p_%d",centBins[cbin],centBins[cbin+1],cbin));
        w1_m = (TH2D *) tfin->Get(Form("%d-%d/w1_m_%d",centBins[cbin],centBins[cbin+1],cbin));
    }
    TH2D * qxav1 = (TH2D *) tfin->Get(Form("%d-%d/qxav1_%d",centBins[cbin],centBins[cbin+1],cbin));
    TH2D * qyav1 = (TH2D *) tfin->Get(Form("%d-%d/qyav1_%d",centBins[cbin],centBins[cbin+1],cbin));
    TH2D * qxycnt = (TH2D *) tfin->Get(Form("%d-%d/qxycnt_%d",centBins[cbin],centBins[cbin+1],cbin));

    double q1AB_p = ((TH2D *) tfin->Get(Form("%d-%d/q1AB_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1AC_p = ((TH2D *) tfin->Get(Form("%d-%d/q1AC_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1BC_p = ((TH2D *) tfin->Get(Form("%d-%d/q1BC_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1AB_m = ((TH2D *) tfin->Get(Form("%d-%d/q1AB_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1AC_m = ((TH2D *) tfin->Get(Form("%d-%d/q1AC_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1BC_m = ((TH2D *) tfin->Get(Form("%d-%d/q1BC_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1ABcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q1ABcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1ACcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q1ACcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1BCcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q1BCcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1ABcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q1ABcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1ACcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q1ACcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1BCcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q1BCcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);

    double q2AB_p = ((TH2D *) tfin->Get(Form("%d-%d/q2AB_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2AC_p = ((TH2D *) tfin->Get(Form("%d-%d/q2AC_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2BC_p = ((TH2D *) tfin->Get(Form("%d-%d/q2BC_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2AB_m = ((TH2D *) tfin->Get(Form("%d-%d/q2AB_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2AC_m = ((TH2D *) tfin->Get(Form("%d-%d/q2AC_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2BC_m = ((TH2D *) tfin->Get(Form("%d-%d/q2BC_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2ABcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q2ABcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2ACcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q2ACcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2BCcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q2BCcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2ABcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q2ABcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2ACcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q2ACcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2BCcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q2BCcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);

    double q3AB_p = ((TH2D *) tfin->Get(Form("%d-%d/q3AB_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3AC_p = ((TH2D *) tfin->Get(Form("%d-%d/q3AC_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3BC_p = ((TH2D *) tfin->Get(Form("%d-%d/q3BC_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3AB_m = ((TH2D *) tfin->Get(Form("%d-%d/q3AB_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3AC_m = ((TH2D *) tfin->Get(Form("%d-%d/q3AC_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3BC_m = ((TH2D *) tfin->Get(Form("%d-%d/q3BC_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3ABcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q3ABcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3ACcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q3ACcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3BCcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q3BCcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3ABcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q3ABcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3ACcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q3ACcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3BCcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q3BCcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);

    ptav->Divide(ptcnt);
    q1_p->Divide(w1_p);
    q1_m->Divide(w1_m);
    qxav1->Divide(qxycnt);
    qyav1->Divide(qxycnt);

    q1AB_p/=q1ABcnt_p;
    q1AC_p/=q1ACcnt_p;
    q1BC_p/=q1BCcnt_p;
    q1AB_m/=q1ABcnt_m;
    q1AC_m/=q1ACcnt_m;
    q1BC_m/=q1BCcnt_m;

    q2AB_p/=q2ABcnt_p;
    q2AC_p/=q2ACcnt_p;
    q2BC_p/=q2BCcnt_p;
    q2AB_m/=q2ABcnt_m;
    q2AC_m/=q2ACcnt_m;
    q2BC_m/=q2BCcnt_m;

    q3AB_p/=q3ABcnt_p;
    q3AC_p/=q3ACcnt_p;
    q3BC_p/=q3BCcnt_p;
    q3AB_m/=q3ABcnt_m;
    q3AC_m/=q3ACcnt_m;
    q3BC_m/=q3BCcnt_m;

    double res1_p, res1_m, res2_p, res2_m, res3_p, res3_m;

    if (sub2) {
        res1_p = sqrt( fabs(q1AB_p) );
        res1_m = sqrt( fabs(q1AB_m) );
    } else {
        res1_p = sqrt( fabs(q1AB_p*q1AC_p/q1BC_p) );
        res1_m = sqrt( fabs(q1AB_m*q1AC_m/q1BC_m) );
    }
    res2_p = sqrt( fabs(q2BC_p*q2AC_p/q2AB_p) );
    res2_m = sqrt( fabs(q2BC_m*q2AC_m/q2AB_m) );
    res3_p = sqrt( fabs(q3BC_p*q3AC_p/q3AB_p) );
    res3_m = sqrt( fabs(q3BC_m*q3AC_m/q3AB_m) );

    if (mix112) {
        q1_p->Scale(1./(res1_p*res2_p));
        q1_m->Scale(1./(res1_m*res2_m));
    } else if (mix123) {
        q1_p->Scale(1./(res2_p*res3_p));
        q1_m->Scale(1./(res2_m*res3_m));
    } else {
        q1_p->Scale(1./res1_p);
        q1_m->Scale(1./res1_m);
    }

    TH2D * q1_p_err[10];
    TH2D * q1_m_err[10];
    TH2D * w1_p_err[10];
    TH2D * w1_m_err[10];
    for (int k = 0; k<10; k++) {
        if (mix112) {
            q1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q112_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            q1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q112_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            w1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w112_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            w1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w112_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
        } else if (mix123) {
            q1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q123_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            q1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q123_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            w1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w123_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            w1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w123_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
        } else {
            q1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q1_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            q1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q1_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            w1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w1_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            w1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w1_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
        }
        q1_p_err[k]->Divide(w1_p_err[k]);
        q1_m_err[k]->Divide(w1_m_err[k]);

        double q1AB_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1AB_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1AC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1AC_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1BC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1BC_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1AB_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1AB_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1AC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1AC_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1BC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1BC_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1ABcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1ABcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1ACcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1ACcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1BCcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1BCcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1ABcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1ABcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1ACcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1ACcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1BCcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1BCcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);

        double q2AB_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2AB_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2AC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2AC_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2BC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2BC_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2AB_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2AB_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2AC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2AC_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2BC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2BC_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2ABcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2ABcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2ACcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2ACcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2BCcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2BCcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2ABcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2ABcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2ACcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2ACcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2BCcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2BCcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);

        double q3AB_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3AB_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3AC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3AC_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3BC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3BC_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3AB_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3AB_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3AC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3AC_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3BC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3BC_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3ABcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3ABcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3ACcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3ACcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3BCcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3BCcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3ABcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3ABcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3ACcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3ACcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3BCcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3BCcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);

        q1AB_p_err/=q1ABcnt_p_err;
        q1AC_p_err/=q1ACcnt_p_err;
        q1BC_p_err/=q1BCcnt_p_err;
        q1AB_m_err/=q1ABcnt_m_err;
        q1AC_m_err/=q1ACcnt_m_err;
        q1BC_m_err/=q1BCcnt_m_err;

        q2AB_p_err/=q2ABcnt_p_err;
        q2AC_p_err/=q2ACcnt_p_err;
        q2BC_p_err/=q2BCcnt_p_err;
        q2AB_m_err/=q2ABcnt_m_err;
        q2AC_m_err/=q2ACcnt_m_err;
        q2BC_m_err/=q2BCcnt_m_err;

        q3AB_p_err/=q3ABcnt_p_err;
        q3AC_p_err/=q3ACcnt_p_err;
        q3BC_p_err/=q3BCcnt_p_err;
        q3AB_m_err/=q3ABcnt_m_err;
        q3AC_m_err/=q3ACcnt_m_err;
        q3BC_m_err/=q3BCcnt_m_err;

        double res1_p_err, res1_m_err, res2_p_err, res2_m_err, res3_p_err, res3_m_err;

        if (sub2) {
            res1_p_err = sqrt( fabs(q1AB_p_err) );
            res1_m_err = sqrt( fabs(q1AB_m_err) );
        } else {
            res1_p_err = sqrt( fabs(q1AB_p_err*q1AC_p_err/q1BC_p_err) );
            res1_m_err = sqrt( fabs(q1AB_m_err*q1AC_m_err/q1BC_m_err) );
        }
        res2_p_err = sqrt( fabs(q2BC_p_err*q2AC_p_err/q2AB_p_err) );
        res2_m_err = sqrt( fabs(q2BC_m_err*q2AC_m_err/q2AB_m_err) );
        res3_p_err = sqrt( fabs(q3BC_p_err*q3AC_p_err/q3AB_p_err) );
        res3_m_err = sqrt( fabs(q3BC_m_err*q3AC_m_err/q3AB_m_err) );

        if (mix112) {
            q1_p_err[k]->Scale(1./(res1_p_err*res2_p_err));
            q1_m_err[k]->Scale(1./(res1_m_err*res2_m_err));
        } else if (mix123) {
            q1_p_err[k]->Scale(1./(res2_p_err*res3_p_err));
            q1_m_err[k]->Scale(1./(res2_m_err*res3_m_err));
        } else {
            q1_p_err[k]->Scale(1./res1_p_err);
            q1_m_err[k]->Scale(1./res1_m_err);
        }
    }

    int etamin1 = 0;
    int etamax1 = 0;
    int etamin2 = 0;
    int etamax2 = 0;

    TH1D * xpt = 0;
    TH1D * v1 = 0;
    TH1D * v1_p = 0;
    TH1D * v1_m = 0;
    TH1D * v12 = 0;
    TH1D * v1_p2 = 0;
    TH1D * v1_m2 = 0;
    if (etamin*etamax<0) {
        etamin1 = q1_p->GetYaxis()->FindBin(etamin);
        etamax1 = q1_p->GetYaxis()->FindBin(-0.001);
        etamin2 = q1_p->GetYaxis()->FindBin(0.001);
        etamax2 = q1_p->GetYaxis()->FindBin(etamax);

        xpt = (TH1D *) ptav->ProjectionX("xpt",etamin1,etamax2);
        double ebins_p = etamax1-etamin1+1;
        double ebins_m = etamax2-etamin2+1;
        xpt->Scale(1./(ebins_p+ebins_m));
        v1_p = (TH1D *) q1_p->ProjectionX("v1_p",etamin1,etamax1);
        v1_m = (TH1D *) q1_m->ProjectionX("v1_m",etamin2,etamax2);
        v1 = (TH1D *) v1_p->Clone("v1");
        v1->Add(v1_m);
        v1_p->Scale(1./ebins_p);
        v1_m->Scale(1./ebins_m);
        v1->Scale(1./(ebins_p+ebins_m));

        // calculate errors from subevents
        double v1tmp[20] = {0};
        double v1tmp_p[20] = {0};
        double v1tmp_m[20] = {0};
        double v1tmp2[20] = {0};
        double v1tmp_p2[20] = {0};
        double v1tmp_m2[20] = {0};

        TH1D * v1e_p;
        TH1D * v1e_m;
        TH1D * v1e;
        for (int i = 0; i<10; i++) {
            v1e_p = (TH1D *) q1_p_err[i]->ProjectionX(Form("v1e_p%d",i),etamin1,etamax1);
            v1e_m = (TH1D *) q1_m_err[i]->ProjectionX(Form("v1e_m%d",i),etamin2,etamax2);
            v1e = (TH1D *) v1e_p->Clone(Form("v1e%d",i));
            v1e->Add(v1e_m);
            v1e_p->Scale(1./ebins_p);
            v1e_m->Scale(1./ebins_m);
            v1e->Scale(1./(ebins_p+ebins_m));
            for (int j = 0; j<v1e->GetNbinsX(); j++) {
                v1tmp[j]+=v1e->GetBinContent(j+1);
                v1tmp_p[j]+=v1e_p->GetBinContent(j+1);
                v1tmp_m[j]+=v1e_m->GetBinContent(j+1);
                v1tmp2[j]+=pow(v1e->GetBinContent(j+1),2);
                v1tmp_p2[j]+=pow(v1e_p->GetBinContent(j+1),2);
                v1tmp_m2[j]+=pow(v1e_m->GetBinContent(j+1),2);
            }
        }
        for (int i = 0; i<v1->GetNbinsX(); i++) {
            v1tmp[i]/=10.;
            v1tmp_p[i]/=10.;
            v1tmp_m[i]/=10.;
            v1tmp2[i]/=10.;
            v1tmp_p2[i]/=10.;
            v1tmp_m2[i]/=10.;
            v1->SetBinError(  i+1, sqrt( (1./9.)*(v1tmp2[i]   - pow(v1tmp[i],2))   ));
            v1_p->SetBinError(i+1, sqrt( (1./9.)*(v1tmp_p2[i] - pow(v1tmp_p[i],2)) ));
            v1_m->SetBinError(i+1, sqrt( (1./9.)*(v1tmp_m2[i] - pow(v1tmp_m[i],2)) ));
        }
    } else {
        etamin1 = q1_p->GetYaxis()->FindBin(etamin);
        etamax1 = q1_p->GetYaxis()->FindBin(-0.001);
        xpt = (TH1D *) ptav->ProjectionX("xpt",etamin1,etamax1);
        double ebins = etamax1-etamin1+1;
        xpt->Scale(1./ebins);
        if (etamin<0) {
            v1 = (TH1D *) q1_p->ProjectionX("v1_p",etamin1,etamax1);
        } else {
            v1 = (TH1D *) q1_m->ProjectionX("v1_p",etamin1,etamax1);
        }
        v1->Scale(1./ebins);
    }

    TH1D * yld = ptcnt->ProjectionX("yld",etamin1,etamax2);
    double x[50];
    double y[50];
    double y_p[50];
    double y_m[50];
    double x_err[50];
    double y_err[50];
    double y_p_err[50];
    double y_m_err[50];
    int npt = 0;

    for (int pbin = 1; pbin<=xpt->GetNbinsX(); pbin++) {
        double pt = xpt->GetBinContent(pbin);
        if (pt>0.3 && pt<12.) {
            x[npt] = pt;
            y[npt] = v1->GetBinContent(pbin);
            y_p[npt] = v1_p->GetBinContent(pbin);
            y_m[npt] = v1_m->GetBinContent(pbin);
            x_err[npt] = 0;
            y_err[npt] = v1->GetBinError(pbin);
            y_p_err[npt] = v1_p->GetBinError(pbin);
            y_m_err[npt] = v1_m->GetBinError(pbin);
            ++npt;
        }
    }

    TGraphErrors * g = new TGraphErrors(npt, x, y, x_err, y_err);
    g_p = new TGraphErrors(npt, x, y_p, x_err, y_p_err);
    g_m = new TGraphErrors(npt, x, y_m, x_err, y_m_err);

    tfin->Close();

    return g;

}


TGraphErrors * GetV1IntEta( string anal, string tag, int cbin, double ptmin, double ptmax, TGraphErrors * &g_p, TGraphErrors * &g_m ) {

    bool sub2 = kFALSE;
    bool mix112 = kFALSE;
    bool mix123 = kFALSE;
    if (tag == "2sub") sub2 = kTRUE;
    else if (tag == "112") mix112 = kTRUE;
    else if (tag == "123") mix123 = kTRUE;
    TFile * tfin = new TFile(Form("../outputs/raw_outputs/results/%s.root",anal.data()),"read");

    TH1D * centbins = (TH1D *) tfin->Get(Form("%d-%d/centbins",centBins[cbin],centBins[cbin+1]));
    TH2D * ptav = (TH2D *) tfin->Get(Form("%d-%d/ptav_%d",centBins[cbin],centBins[cbin+1],cbin));
    TH2D * ptcnt = (TH2D *) tfin->Get(Form("%d-%d/ptcnt_%d",centBins[cbin],centBins[cbin+1],cbin));
    TH2D * badcnt = (TH2D *) tfin->Get(Form("%d-%d/badcnt_%d",centBins[cbin],centBins[cbin+1],cbin));
    TH1D * multTot = (TH1D *) tfin->Get(Form("%d-%d/multTot_%d",centBins[cbin],centBins[cbin+1],cbin));
    TH2D * q1_p;
    TH2D * q1_m;
    TH2D * w1_p;
    TH2D * w1_m;
    if (mix112) {
        q1_p = (TH2D *) tfin->Get(Form("%d-%d/q112_p_%d",centBins[cbin],centBins[cbin+1],cbin));
        q1_m = (TH2D *) tfin->Get(Form("%d-%d/q112_m_%d",centBins[cbin],centBins[cbin+1],cbin));
        w1_p = (TH2D *) tfin->Get(Form("%d-%d/w112_p_%d",centBins[cbin],centBins[cbin+1],cbin));
        w1_m = (TH2D *) tfin->Get(Form("%d-%d/w112_m_%d",centBins[cbin],centBins[cbin+1],cbin));
    } else if (mix123) {
        q1_p = (TH2D *) tfin->Get(Form("%d-%d/q123_p_%d",centBins[cbin],centBins[cbin+1],cbin));
        q1_m = (TH2D *) tfin->Get(Form("%d-%d/q123_m_%d",centBins[cbin],centBins[cbin+1],cbin));
        w1_p = (TH2D *) tfin->Get(Form("%d-%d/w123_p_%d",centBins[cbin],centBins[cbin+1],cbin));
        w1_m = (TH2D *) tfin->Get(Form("%d-%d/w123_m_%d",centBins[cbin],centBins[cbin+1],cbin));
    } else {
        q1_p = (TH2D *) tfin->Get(Form("%d-%d/q1_p_%d",centBins[cbin],centBins[cbin+1],cbin));
        q1_m = (TH2D *) tfin->Get(Form("%d-%d/q1_m_%d",centBins[cbin],centBins[cbin+1],cbin));
        w1_p = (TH2D *) tfin->Get(Form("%d-%d/w1_p_%d",centBins[cbin],centBins[cbin+1],cbin));
        w1_m = (TH2D *) tfin->Get(Form("%d-%d/w1_m_%d",centBins[cbin],centBins[cbin+1],cbin));
    }
    TH2D * qxav1 = (TH2D *) tfin->Get(Form("%d-%d/qxav1_%d",centBins[cbin],centBins[cbin+1],cbin));
    TH2D * qyav1 = (TH2D *) tfin->Get(Form("%d-%d/qyav1_%d",centBins[cbin],centBins[cbin+1],cbin));
    TH2D * qxycnt = (TH2D *) tfin->Get(Form("%d-%d/qxycnt_%d",centBins[cbin],centBins[cbin+1],cbin));

    double q1AB_p = ((TH2D *) tfin->Get(Form("%d-%d/q1AB_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1AC_p = ((TH2D *) tfin->Get(Form("%d-%d/q1AC_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1BC_p = ((TH2D *) tfin->Get(Form("%d-%d/q1BC_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1AB_m = ((TH2D *) tfin->Get(Form("%d-%d/q1AB_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1AC_m = ((TH2D *) tfin->Get(Form("%d-%d/q1AC_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1BC_m = ((TH2D *) tfin->Get(Form("%d-%d/q1BC_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1ABcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q1ABcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1ACcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q1ACcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1BCcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q1BCcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1ABcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q1ABcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1ACcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q1ACcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q1BCcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q1BCcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);

    double q2AB_p = ((TH2D *) tfin->Get(Form("%d-%d/q2AB_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2AC_p = ((TH2D *) tfin->Get(Form("%d-%d/q2AC_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2BC_p = ((TH2D *) tfin->Get(Form("%d-%d/q2BC_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2AB_m = ((TH2D *) tfin->Get(Form("%d-%d/q2AB_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2AC_m = ((TH2D *) tfin->Get(Form("%d-%d/q2AC_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2BC_m = ((TH2D *) tfin->Get(Form("%d-%d/q2BC_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2ABcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q2ABcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2ACcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q2ACcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2BCcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q2BCcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2ABcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q2ABcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2ACcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q2ACcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q2BCcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q2BCcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);

    double q3AB_p = ((TH2D *) tfin->Get(Form("%d-%d/q3AB_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3AC_p = ((TH2D *) tfin->Get(Form("%d-%d/q3AC_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3BC_p = ((TH2D *) tfin->Get(Form("%d-%d/q3BC_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3AB_m = ((TH2D *) tfin->Get(Form("%d-%d/q3AB_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3AC_m = ((TH2D *) tfin->Get(Form("%d-%d/q3AC_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3BC_m = ((TH2D *) tfin->Get(Form("%d-%d/q3BC_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3ABcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q3ABcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3ACcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q3ACcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3BCcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q3BCcnt_p_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3ABcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q3ABcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3ACcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q3ACcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);
    double q3BCcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q3BCcnt_m_%d",centBins[cbin],centBins[cbin+1],cbin)))->GetBinContent(1);

    ptav->Divide(ptcnt);
    q1_p->Divide(w1_p);
    q1_m->Divide(w1_m);
    qxav1->Divide(qxycnt);
    qyav1->Divide(qxycnt);

    q1AB_p/=q1ABcnt_p;
    q1AC_p/=q1ACcnt_p;
    q1BC_p/=q1BCcnt_p;
    q1AB_m/=q1ABcnt_m;
    q1AC_m/=q1ACcnt_m;
    q1BC_m/=q1BCcnt_m;

    q2AB_p/=q2ABcnt_p;
    q2AC_p/=q2ACcnt_p;
    q2BC_p/=q2BCcnt_p;
    q2AB_m/=q2ABcnt_m;
    q2AC_m/=q2ACcnt_m;
    q2BC_m/=q2BCcnt_m;

    q3AB_p/=q3ABcnt_p;
    q3AC_p/=q3ACcnt_p;
    q3BC_p/=q3BCcnt_p;
    q3AB_m/=q3ABcnt_m;
    q3AC_m/=q3ACcnt_m;
    q3BC_m/=q3BCcnt_m;

    double res1_p, res1_m, res2_p, res2_m, res3_p, res3_m;

    if (sub2) {
        res1_p = sqrt( fabs(q1AB_p) );
        res1_m = sqrt( fabs(q1AB_m) );
    } else {
        res1_p = sqrt( fabs(q1AB_p*q1AC_p/q1BC_p) );
        res1_m = sqrt( fabs(q1AB_m*q1AC_m/q1BC_m) );
    }
    res2_p = sqrt( fabs(q2BC_p*q2AC_p/q2AB_p) );
    res2_m = sqrt( fabs(q2BC_m*q2AC_m/q2AB_m) );
    res3_p = sqrt( fabs(q3BC_p*q3AC_p/q3AB_p) );
    res3_m = sqrt( fabs(q3BC_m*q3AC_m/q3AB_m) );

    if (mix112) {
        q1_p->Scale(1./(res1_p*res2_p));
        q1_m->Scale(1./(res1_m*res2_m));
    } else if (mix123) {
        q1_p->Scale(1./(res2_p*res3_p));
        q1_m->Scale(1./(res2_m*res3_m));
    } else {
        q1_p->Scale(1./res1_p);
        q1_m->Scale(1./res1_m);
    }

    TH2D * q1_p_err[10];
    TH2D * q1_m_err[10];
    TH2D * w1_p_err[10];
    TH2D * w1_m_err[10];
    for (int k = 0; k<10; k++) {
        if (mix112) {
            q1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q112_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            q1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q112_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            w1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w112_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            w1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w112_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
        } else if (mix123) {
            q1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q123_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            q1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q123_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            w1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w123_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            w1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w123_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
        } else {
            q1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q1_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            q1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q1_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            w1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w1_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
            w1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w1_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1));
        }
        q1_p_err[k]->Divide(w1_p_err[k]);
        q1_m_err[k]->Divide(w1_m_err[k]);

        double q1AB_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1AB_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1AC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1AC_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1BC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1BC_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1AB_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1AB_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1AC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1AC_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1BC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1BC_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1ABcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1ABcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1ACcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1ACcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1BCcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1BCcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1ABcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1ABcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1ACcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1ACcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q1BCcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1BCcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);

        double q2AB_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2AB_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2AC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2AC_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2BC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2BC_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2AB_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2AB_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2AC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2AC_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2BC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2BC_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2ABcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2ABcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2ACcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2ACcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2BCcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2BCcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2ABcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2ABcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2ACcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2ACcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q2BCcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2BCcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);

        double q3AB_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3AB_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3AC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3AC_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3BC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3BC_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3AB_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3AB_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3AC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3AC_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3BC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3BC_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3ABcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3ABcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3ACcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3ACcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3BCcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3BCcnt_p_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3ABcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3ABcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3ACcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3ACcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
        double q3BCcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3BCcnt_m_%d_%d",centBins[cbin],centBins[cbin+1],cbin,k+1)))->GetBinContent(1);

        q1AB_p_err/=q1ABcnt_p_err;
        q1AC_p_err/=q1ACcnt_p_err;
        q1BC_p_err/=q1BCcnt_p_err;
        q1AB_m_err/=q1ABcnt_m_err;
        q1AC_m_err/=q1ACcnt_m_err;
        q1BC_m_err/=q1BCcnt_m_err;

        q2AB_p_err/=q2ABcnt_p_err;
        q2AC_p_err/=q2ACcnt_p_err;
        q2BC_p_err/=q2BCcnt_p_err;
        q2AB_m_err/=q2ABcnt_m_err;
        q2AC_m_err/=q2ACcnt_m_err;
        q2BC_m_err/=q2BCcnt_m_err;

        q3AB_p_err/=q3ABcnt_p_err;
        q3AC_p_err/=q3ACcnt_p_err;
        q3BC_p_err/=q3BCcnt_p_err;
        q3AB_m_err/=q3ABcnt_m_err;
        q3AC_m_err/=q3ACcnt_m_err;
        q3BC_m_err/=q3BCcnt_m_err;

        double res1_p_err, res1_m_err, res2_p_err, res2_m_err, res3_p_err, res3_m_err;

        if (sub2) {
            res1_p_err = sqrt( fabs(q1AB_p_err) );
            res1_m_err = sqrt( fabs(q1AB_m_err) );
        } else {
            res1_p_err = sqrt( fabs(q1AB_p_err*q1AC_p_err/q1BC_p_err) );
            res1_m_err = sqrt( fabs(q1AB_m_err*q1AC_m_err/q1BC_m_err) );
        }
        res2_p_err = sqrt( fabs(q2BC_p_err*q2AC_p_err/q2AB_p_err) );
        res2_m_err = sqrt( fabs(q2BC_m_err*q2AC_m_err/q2AB_m_err) );
        res3_p_err = sqrt( fabs(q3BC_p_err*q3AC_p_err/q3AB_p_err) );
        res3_m_err = sqrt( fabs(q3BC_m_err*q3AC_m_err/q3AB_m_err) );

        if (mix112) {
            q1_p_err[k]->Scale(1./(res1_p_err*res2_p_err));
            q1_m_err[k]->Scale(1./(res1_m_err*res2_m_err));
        } else if (mix123) {
            q1_p_err[k]->Scale(1./(res2_p_err*res3_p_err));
            q1_m_err[k]->Scale(1./(res2_m_err*res3_m_err));
        } else {
            q1_p_err[k]->Scale(1./res1_p_err);
            q1_m_err[k]->Scale(1./res1_m_err);
        }
    }

    int ptmin1 = 0;
    int ptmax1 = 0;

    TH1D * xeta = 0;
    TH1D * v1 = 0;
    TH1D * v1_p = 0;
    TH1D * v1_m = 0;
    TH1D * v12 = 0;
    TH1D * v1_p2 = 0;
    TH1D * v1_m2 = 0;

    ptmin1 = q1_p->GetXaxis()->FindBin(ptmin);
    ptmax1 = q1_p->GetXaxis()->FindBin(ptmax);

    xeta = (TH1D *) ptav->ProjectionY("xeta",ptmin1,ptmax1);
    double pbins = ptmax1-ptmin1+1;
    xeta->Scale(1./pbins);
    v1_p = (TH1D *) q1_p->ProjectionY("v1_p",ptmin1,ptmax1);
    v1_m = (TH1D *) q1_m->ProjectionY("v1_m",ptmin1,ptmax1);
    v1 = (TH1D *) v1_p->Clone("v1");
    v1->Add(v1_m);
    v1_p->Scale(1./pbins);
    v1_m->Scale(1./pbins);
    v1->Scale(1./(2.*pbins));

    // calculate errors fromsubevents
    double v1tmp[50] = {0};
    double v1tmp_p[50] = {0};
    double v1tmp_m[50] = {0};
    double v1tmp2[50] = {0};
    double v1tmp_p2[50] = {0};
    double v1tmp_m2[50] = {0};

    TH1D * v1e_p;
    TH1D * v1e_m;
    TH1D * v1e;
    for (int i = 0; i<10; i++) {
        v1e_p = (TH1D *) q1_p_err[i]->ProjectionY(Form("v1e_p%d",i),ptmin1,ptmax1);
        v1e_m = (TH1D *) q1_m_err[i]->ProjectionY(Form("v1e_m%d",i),ptmin1,ptmax1);
        v1e = (TH1D *) v1e_p->Clone(Form("v1e%d",i));
        v1e->Add(v1e_m);
        v1e_p->Scale(1./pbins);
        v1e_m->Scale(1./pbins);
        v1e->Scale(1./(2*pbins));
        for (int j = 0; j<v1e->GetNbinsX(); j++) {
            v1tmp[j]+=v1e->GetBinContent(j+1);
            v1tmp_p[j]+=v1e_p->GetBinContent(j+1);
            v1tmp_m[j]+=v1e_m->GetBinContent(j+1);
            v1tmp2[j]+=pow(v1e->GetBinContent(j+1),2);
            v1tmp_p2[j]+=pow(v1e_p->GetBinContent(j+1),2);
            v1tmp_m2[j]+=pow(v1e_m->GetBinContent(j+1),2);
        }
    }
    for (int i = 0; i<v1->GetNbinsX(); i++) {
        v1tmp[i]/=10;
        v1tmp_p[i]/=10;
        v1tmp_m[i]/=10;
        v1tmp2[i]/=10;
        v1tmp_p2[i]/=10;
        v1tmp_m2[i]/=10;
        v1->SetBinError(  i+1, sqrt( (1./9.)*(v1tmp2[i]   - pow(v1tmp[i],2))   ));
        v1_p->SetBinError(i+1, sqrt( (1./9.)*(v1tmp_p2[i] - pow(v1tmp_p[i],2)) ));
        v1_m->SetBinError(i+1, sqrt( (1./9.)*(v1tmp_m2[i] - pow(v1tmp_m[i],2)) ));
    }

    TH1D * yld = ptcnt->ProjectionY("yld",ptmin1,ptmax1);
    double x[50];
    double y[50];
    double y_p[50];
    double y_m[50];
    double x_err[50];
    double y_err[50];
    double y_p_err[50];
    double y_m_err[50];
    int neta = 0;

    for (int ebin = 1; ebin<=xeta->GetNbinsX(); ebin++) {
        double eta = xeta->GetBinContent(ebin);
        if (eta>-2.4 && eta<2.4) {
            x[neta] = eta;
            y[neta] = v1->GetBinContent(ebin);
            y_p[neta] = v1_p->GetBinContent(ebin);
            y_m[neta] = v1_m->GetBinContent(ebin);
            x_err[neta] = 0;
            y_err[neta] = v1->GetBinError(ebin);
            y_p_err[neta] = v1_p->GetBinError(ebin);
            y_m_err[neta] = v1_m->GetBinError(ebin);
            ++neta;
        }
    }

    TGraphErrors * g = new TGraphErrors(neta, x, y, x_err, y_err);
    g_p = new TGraphErrors(neta, x, y_p, x_err, y_p_err);
    g_m = new TGraphErrors(neta, x, y_m, x_err, y_m_err);

    tfin->Close();

    return g;

}
