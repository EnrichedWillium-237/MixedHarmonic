# include "TAxis.h"
# include "TFile.h"
# include "TGraphErrors.h"
# include "TCanvas.h"
# include "TH1D.h"
# include "TH2D.h"
# include "TLegend.h"
# include "TMath.h"
# include "TPaveText.h"
# include "TROOT.h"
# include "TStopwatch.h"
# include "TStyle.h"
# include <fstream>
# include <iostream>

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

void AnalyzeDiffV1()
{
    int minanal = 0;
    int maxanal = nanals;

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    cout << "\nBeginning differential v1 analyzer...\n" << endl;
    TStopwatch * sw = new TStopwatch();
    sw->Start();

    for (int i = minanal; i<maxanal; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            for (int ebin = 0; ebin<netabins; ebin++) {
                v1p_pt[i][cbin][ebin] = new TH1D(Form("v1p_pt_%s_c%d_e%d",AnalNames[i].data(),cbin,ebin), "", nptbins, ptbins);
                v1m_pt[i][cbin][ebin] = new TH1D(Form("v1m_pt_%s_c%d_e%d",AnalNames[i].data(),cbin,ebin), "", nptbins, ptbins);
                v1odd_pt[i][cbin][ebin] = new TH1D(Form("v1odd_pt_%s_c%d_e%d",AnalNames[i].data(),cbin,ebin), "", nptbins, ptbins);
                v1even_pt[i][cbin][ebin] = new TH1D(Form("v1even_pt_%s_c%d_e%d",AnalNames[i].data(),cbin,ebin), "", nptbins, ptbins);

                v112p_pt[i][cbin][ebin] = new TH1D(Form("v112p_pt_%s_c%d_e%d",AnalNames[i].data(),cbin,ebin), "", nptbins, ptbins);
                v112m_pt[i][cbin][ebin] = new TH1D(Form("v112m_pt_%s_c%d_e%d",AnalNames[i].data(),cbin,ebin), "", nptbins, ptbins);
                v112odd_pt[i][cbin][ebin] = new TH1D(Form("v112odd_pt_%s_c%d_e%d",AnalNames[i].data(),cbin,ebin), "", nptbins, ptbins);
                v112even_pt[i][cbin][ebin] = new TH1D(Form("v112even_pt_%s_c%d_e%d",AnalNames[i].data(),cbin,ebin), "", nptbins, ptbins);

                v123p_pt[i][cbin][ebin] = new TH1D(Form("v123p_pt_%s_c%d_e%d",AnalNames[i].data(),cbin,ebin), "", nptbins, ptbins);
                v123m_pt[i][cbin][ebin] = new TH1D(Form("v123m_pt_%s_c%d_e%d",AnalNames[i].data(),cbin,ebin), "", nptbins, ptbins);
                v123odd_pt[i][cbin][ebin] = new TH1D(Form("v123odd_pt_%s_c%d_e%d",AnalNames[i].data(),cbin,ebin), "", nptbins, ptbins);
                v123even_pt[i][cbin][ebin] = new TH1D(Form("v123even_pt_%s_c%d_e%d",AnalNames[i].data(),cbin,ebin), "", nptbins, ptbins);

                v1p_pt[i][cbin][ebin]->SetXTitle(Form("%0.1f < #eta < %0.1f, %d to %d%%",etabins[ebin],etabins[ebin+1],centBins[cbin],centBins[cbin+1]));
                v1m_pt[i][cbin][ebin]->SetXTitle(Form("%0.1f < #eta < %0.1f, %d to %d%%",etabins[ebin],etabins[ebin+1],centBins[cbin],centBins[cbin+1]));
                v1odd_pt[i][cbin][ebin]->SetXTitle(Form("%0.1f < #eta < %0.1f, %d to %d%%",etabins[ebin],etabins[ebin+1],centBins[cbin],centBins[cbin+1]));
                v1even_pt[i][cbin][ebin]->SetXTitle(Form("%0.1f < #eta < %0.1f, %d to %d%%",etabins[ebin],etabins[ebin+1],centBins[cbin],centBins[cbin+1]));

                v112p_pt[i][cbin][ebin]->SetXTitle(Form("%0.1f < #eta < %0.1f, %d to %d%%",etabins[ebin],etabins[ebin+1],centBins[cbin],centBins[cbin+1]));
                v112m_pt[i][cbin][ebin]->SetXTitle(Form("%0.1f < #eta < %0.1f, %d to %d%%",etabins[ebin],etabins[ebin+1],centBins[cbin],centBins[cbin+1]));
                v112odd_pt[i][cbin][ebin]->SetXTitle(Form("%0.1f < #eta < %0.1f, %d to %d%%",etabins[ebin],etabins[ebin+1],centBins[cbin],centBins[cbin+1]));
                v112even_pt[i][cbin][ebin]->SetXTitle(Form("%0.1f < #eta < %0.1f, %d to %d%%",etabins[ebin],etabins[ebin+1],centBins[cbin],centBins[cbin+1]));

                v123p_pt[i][cbin][ebin]->SetXTitle(Form("%0.1f < #eta < %0.1f, %d to %d%%",etabins[ebin],etabins[ebin+1],centBins[cbin],centBins[cbin+1]));
                v123m_pt[i][cbin][ebin]->SetXTitle(Form("%0.1f < #eta < %0.1f, %d to %d%%",etabins[ebin],etabins[ebin+1],centBins[cbin],centBins[cbin+1]));
                v123odd_pt[i][cbin][ebin]->SetXTitle(Form("%0.1f < #eta < %0.1f, %d to %d%%",etabins[ebin],etabins[ebin+1],centBins[cbin],centBins[cbin+1]));
                v123even_pt[i][cbin][ebin]->SetXTitle(Form("%0.1f < #eta < %0.1f, %d to %d%%",etabins[ebin],etabins[ebin+1],centBins[cbin],centBins[cbin+1]));
            }
            for (int pbin = 0; pbin<nptbins; pbin++) {
                v1p_eta[i][cbin][pbin] = new TH1D(Form("v1p_eta_%s_c%d_p%d",AnalNames[i].data(),cbin,pbin), "", netabins, etabins);
                v1m_eta[i][cbin][pbin] = new TH1D(Form("v1m_eta_%s_c%d_p%d",AnalNames[i].data(),cbin,pbin), "", netabins, etabins);
                v1odd_eta[i][cbin][pbin] = new TH1D(Form("v1odd_eta_%s_c%d_p%d",AnalNames[i].data(),cbin,pbin), "", netabins, etabins);
                v1even_eta[i][cbin][pbin] = new TH1D(Form("v1even_eta_%s_c%d_p%d",AnalNames[i].data(),cbin,pbin), "", netabins, etabins);

                v112p_eta[i][cbin][pbin] = new TH1D(Form("v112p_eta_%s_c%d_p%d",AnalNames[i].data(),cbin,pbin), "", netabins, etabins);
                v112m_eta[i][cbin][pbin] = new TH1D(Form("v112m_eta_%s_c%d_p%d",AnalNames[i].data(),cbin,pbin), "", netabins, etabins);
                v112odd_eta[i][cbin][pbin] = new TH1D(Form("v112odd_eta_%s_c%d_p%d",AnalNames[i].data(),cbin,pbin), "", netabins, etabins);
                v112even_eta[i][cbin][pbin] = new TH1D(Form("v112even_eta_%s_c%d_p%d",AnalNames[i].data(),cbin,pbin), "", netabins, etabins);

                v123p_eta[i][cbin][pbin] = new TH1D(Form("v123p_eta_%s_c%d_p%d",AnalNames[i].data(),cbin,pbin), "", netabins, etabins);
                v123m_eta[i][cbin][pbin] = new TH1D(Form("v123m_eta_%s_c%d_p%d",AnalNames[i].data(),cbin,pbin), "", netabins, etabins);
                v123odd_eta[i][cbin][pbin] = new TH1D(Form("v123odd_eta_%s_c%d_p%d",AnalNames[i].data(),cbin,pbin), "", netabins, etabins);
                v123even_eta[i][cbin][pbin] = new TH1D(Form("v123even_eta_%s_c%d_p%d",AnalNames[i].data(),cbin,pbin), "", netabins, etabins);

                v1p_eta[i][cbin][pbin]->SetXTitle(Form("%0.2f < p_{T} < %0.2f (GeV/c), %d to %d%%",ptbins[pbin],ptbins[pbin+1],centBins[cbin],centBins[cbin+1]));
                v1m_eta[i][cbin][pbin]->SetXTitle(Form("%0.2f < p_{T} < %0.2f (GeV/c), %d to %d%%",ptbins[pbin],ptbins[pbin+1],centBins[cbin],centBins[cbin+1]));
                v1odd_eta[i][cbin][pbin]->SetXTitle(Form("%0.2f < p_{T} < %0.2f (GeV/c), %d to %d%%",ptbins[pbin],ptbins[pbin+1],centBins[cbin],centBins[cbin+1]));
                v1even_eta[i][cbin][pbin]->SetXTitle(Form("%0.2f < p_{T} < %0.2f (GeV/c), %d to %d%%",ptbins[pbin],ptbins[pbin+1],centBins[cbin],centBins[cbin+1]));

                v112p_eta[i][cbin][pbin]->SetXTitle(Form("%0.2f < p_{T} < %0.2f (GeV/c), %d to %d%%",ptbins[pbin],ptbins[pbin+1],centBins[cbin],centBins[cbin+1]));
                v112m_eta[i][cbin][pbin]->SetXTitle(Form("%0.2f < p_{T} < %0.2f (GeV/c), %d to %d%%",ptbins[pbin],ptbins[pbin+1],centBins[cbin],centBins[cbin+1]));
                v112odd_eta[i][cbin][pbin]->SetXTitle(Form("%0.2f < p_{T} < %0.2f (GeV/c), %d to %d%%",ptbins[pbin],ptbins[pbin+1],centBins[cbin],centBins[cbin+1]));
                v112even_eta[i][cbin][pbin]->SetXTitle(Form("%0.2f < p_{T} < %0.2f (GeV/c), %d to %d%%",ptbins[pbin],ptbins[pbin+1],centBins[cbin],centBins[cbin+1]));

                v123p_eta[i][cbin][pbin]->SetXTitle(Form("%0.2f < p_{T} < %0.2f (GeV/c), %d to %d%%",ptbins[pbin],ptbins[pbin+1],centBins[cbin],centBins[cbin+1]));
                v123m_eta[i][cbin][pbin]->SetXTitle(Form("%0.2f < p_{T} < %0.2f (GeV/c), %d to %d%%",ptbins[pbin],ptbins[pbin+1],centBins[cbin],centBins[cbin+1]));
                v123odd_eta[i][cbin][pbin]->SetXTitle(Form("%0.2f < p_{T} < %0.2f (GeV/c), %d to %d%%",ptbins[pbin],ptbins[pbin+1],centBins[cbin],centBins[cbin+1]));
                v123even_eta[i][cbin][pbin]->SetXTitle(Form("%0.2f < p_{T} < %0.2f (GeV/c), %d to %d%%",ptbins[pbin],ptbins[pbin+1],centBins[cbin],centBins[cbin+1]));
            }
        }
        TFile * tfParms = new TFile(Form("../outputs/raw_outputs/results/%s.root",AnalNames[i].data()),"read");
        runParms[i] = (TH1D *) tfParms->Get(Form("%d-%d/runParms",(int)centBins[0],(int)centBins[1]));
    }

    // retrieve raw histograms
    for (int anal = minanal; anal<maxanal; anal++) {

        bool sub2 = kFALSE;
        if (anal>7) sub2 = kTRUE;
        sw->Continue();
        double elapse = sw->RealTime();
        cout << "  processing file: " << AnalNames[anal].data() << "\ttime elapsed: " << elapse << " seconds" << endl;
        TFile * tfin = new TFile(Form("../outputs/raw_outputs/results/%s.root",AnalNames[anal].data()),"read");

        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TH1D * centbins = (TH1D *) tfin->Get(Form("%d-%d/centbins",(int)centBins[cbin],(int)centBins[cbin+1]));
            TH2D * ptav = (TH2D *) tfin->Get(Form("%d-%d/ptav_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * ptcnt = (TH2D *) tfin->Get(Form("%d-%d/ptcnt_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * badcnt = (TH2D *) tfin->Get(Form("%d-%d/badcnt_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH1D * multTot = (TH1D *) tfin->Get(Form("%d-%d/multTot_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * qxav1 = (TH2D *) tfin->Get(Form("%d-%d/qxav1_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * qyav1 = (TH2D *) tfin->Get(Form("%d-%d/qyav1_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * qxycnt = (TH2D *) tfin->Get(Form("%d-%d/qxycnt_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));

            TH2D * q1_p = (TH2D *) tfin->Get(Form("%d-%d/q1_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * q1_m = (TH2D *) tfin->Get(Form("%d-%d/q1_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * w1_p = (TH2D *) tfin->Get(Form("%d-%d/w1_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * w1_m = (TH2D *) tfin->Get(Form("%d-%d/w1_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));

            TH2D * q112_p = (TH2D *) tfin->Get(Form("%d-%d/q112_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * q112_m = (TH2D *) tfin->Get(Form("%d-%d/q112_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * w112_p = (TH2D *) tfin->Get(Form("%d-%d/w112_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * w112_m = (TH2D *) tfin->Get(Form("%d-%d/w112_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));

            TH2D * q123_p = (TH2D *) tfin->Get(Form("%d-%d/q123_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * q123_m = (TH2D *) tfin->Get(Form("%d-%d/q123_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * w123_p = (TH2D *) tfin->Get(Form("%d-%d/w123_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));
            TH2D * w123_m = (TH2D *) tfin->Get(Form("%d-%d/w123_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin));

            double q1AB_p = ((TH2D *) tfin->Get(Form("%d-%d/q1AB_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q1AC_p = ((TH2D *) tfin->Get(Form("%d-%d/q1AC_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q1BC_p = ((TH2D *) tfin->Get(Form("%d-%d/q1BC_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q1AB_m = ((TH2D *) tfin->Get(Form("%d-%d/q1AB_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q1AC_m = ((TH2D *) tfin->Get(Form("%d-%d/q1AC_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q1BC_m = ((TH2D *) tfin->Get(Form("%d-%d/q1BC_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q1ABcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q1ABcnt_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q1ACcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q1ACcnt_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q1BCcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q1BCcnt_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q1ABcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q1ABcnt_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q1ACcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q1ACcnt_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q1BCcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q1BCcnt_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);

            double q2AB_p = ((TH2D *) tfin->Get(Form("%d-%d/q2AB_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q2AC_p = ((TH2D *) tfin->Get(Form("%d-%d/q2AC_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q2BC_p = ((TH2D *) tfin->Get(Form("%d-%d/q2BC_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q2AB_m = ((TH2D *) tfin->Get(Form("%d-%d/q2AB_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q2AC_m = ((TH2D *) tfin->Get(Form("%d-%d/q2AC_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q2BC_m = ((TH2D *) tfin->Get(Form("%d-%d/q2BC_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q2ABcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q2ABcnt_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q2ACcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q2ACcnt_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q2BCcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q2BCcnt_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q2ABcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q2ABcnt_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q2ACcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q2ACcnt_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q2BCcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q2BCcnt_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);

            double q3AB_p = ((TH2D *) tfin->Get(Form("%d-%d/q3AB_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q3AC_p = ((TH2D *) tfin->Get(Form("%d-%d/q3AC_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q3BC_p = ((TH2D *) tfin->Get(Form("%d-%d/q3BC_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q3AB_m = ((TH2D *) tfin->Get(Form("%d-%d/q3AB_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q3AC_m = ((TH2D *) tfin->Get(Form("%d-%d/q3AC_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q3BC_m = ((TH2D *) tfin->Get(Form("%d-%d/q3BC_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q3ABcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q3ABcnt_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q3ACcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q3ACcnt_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q3BCcnt_p = ((TH2D *) tfin->Get(Form("%d-%d/q3BCcnt_p_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q3ABcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q3ABcnt_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q3ACcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q3ACcnt_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);
            double q3BCcnt_m = ((TH2D *) tfin->Get(Form("%d-%d/q3BCcnt_m_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin)))->GetBinContent(1);

            ptav->Divide(ptcnt);
            qxav1->Divide(qxycnt);
            qyav1->Divide(qxycnt);

            q1_p->Divide(w1_p);
            q1_m->Divide(w1_m);
            q112_p->Divide(w112_p);
            q112_m->Divide(w112_m);
            q123_p->Divide(w123_p);
            q123_m->Divide(w123_m);

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

            q1_p->Scale(1./res1_p);
            q1_m->Scale(1./res1_m);

            q112_p->Scale(1./(res1_p*res2_p));
            q112_m->Scale(1./(res1_m*res2_m));

            q123_p->Scale(1./(res2_p*res3_p));
            q123_m->Scale(1./(res2_m*res3_m));

            TH2D * q1_p_err[10];
            TH2D * q1_m_err[10];
            TH2D * w1_p_err[10];
            TH2D * w1_m_err[10];
            TH2D * q112_p_err[10];
            TH2D * q112_m_err[10];
            TH2D * w112_p_err[10];
            TH2D * w112_m_err[10];
            TH2D * q123_p_err[10];
            TH2D * q123_m_err[10];
            TH2D * w123_p_err[10];
            TH2D * w123_m_err[10];

            for (int k = 0; k<10; k++) {
                q1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q1_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1));
                q1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q1_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1));
                w1_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w1_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1));
                w1_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w1_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1));

                q112_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q112_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1));
                q112_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q112_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1));
                w112_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w112_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1));
                w112_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w112_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1));

                q123_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q123_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1));
                q123_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/q123_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1));
                w123_p_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w123_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1));
                w123_m_err[k] = (TH2D *) tfin->Get(Form("%d-%d/err_calc/w123_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1));

                q1_p_err[k]->Divide(w1_p_err[k]);
                q1_m_err[k]->Divide(w1_m_err[k]);

                q112_p_err[k]->Divide(w112_p_err[k]);
                q112_m_err[k]->Divide(w112_m_err[k]);

                q123_p_err[k]->Divide(w123_p_err[k]);
                q123_m_err[k]->Divide(w123_m_err[k]);

                double q1AB_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1AB_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q1AC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1AC_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q1BC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1BC_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q1AB_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1AB_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q1AC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1AC_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q1BC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1BC_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q1ABcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1ABcnt_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q1ACcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1ACcnt_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q1BCcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1BCcnt_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q1ABcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1ABcnt_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q1ACcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1ACcnt_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q1BCcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q1BCcnt_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);

                double q2AB_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2AB_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q2AC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2AC_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q2BC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2BC_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q2AB_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2AB_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q2AC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2AC_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q2BC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2BC_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q2ABcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2ABcnt_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q2ACcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2ACcnt_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q2BCcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2BCcnt_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q2ABcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2ABcnt_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q2ACcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2ACcnt_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q2BCcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q2BCcnt_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);

                double q3AB_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3AB_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q3AC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3AC_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q3BC_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3BC_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q3AB_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3AB_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q3AC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3AC_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q3BC_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3BC_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q3ABcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3ABcnt_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q3ACcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3ACcnt_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q3BCcnt_p_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3BCcnt_p_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q3ABcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3ABcnt_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q3ACcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3ACcnt_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);
                double q3BCcnt_m_err = ((TH2D *) tfin->Get(Form("%d-%d/err_calc/q3BCcnt_m_%d_%d",(int)centBins[cbin],(int)centBins[cbin+1],cbin,k+1)))->GetBinContent(1);

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

                q1_p_err[k]->Scale(1./res1_p_err);
                q1_m_err[k]->Scale(1./res1_m_err);

                q112_p_err[k]->Scale(1./(res1_p_err*res2_p_err));
                q112_m_err[k]->Scale(1./(res1_m_err*res2_m_err));

                q123_p_err[k]->Scale(1./(res2_p_err*res3_p_err));
                q123_m_err[k]->Scale(1./(res2_m_err*res3_m_err));
            }


            //-- get 1D histos for differential v1

            //-- v1(pt)
            for (int ebin = 0; ebin<netabins; ebin++) {
                TH1D * xpt = 0;
                TH1D * v1 = 0;
                TH1D * v1_p = 0;
                TH1D * v1_m = 0;
                TH1D * v112 = 0;
                TH1D * v112_p = 0;
                TH1D * v112_m = 0;
                TH1D * v123 = 0;
                TH1D * v123_p = 0;
                TH1D * v123_m = 0;

                xpt = (TH1D *) ptav->ProjectionX("xpt",ebin,(int)(ebin+0.001));
                v1_p = new TH1D("v1_p", "", nptbins, ptbins);
                v1_m = new TH1D("v1_m", "", nptbins, ptbins);
                for (int pbin = 0; pbin<nptbins; pbin++) {
                    v1_p->SetBinContent(pbin+1, q1_p->GetBinContent(pbin+1, ebin+1));
                    v1_m->SetBinContent(pbin+1, q1_m->GetBinContent(pbin+1, ebin+1));
                }
                v1 = (TH1D *) v1_p->Clone("v1");
                v1->Add(v1_m);
                v1->Scale(0.5);

                v112_p = new TH1D("v112_p", "", nptbins, ptbins);
                v112_m = new TH1D("v112_m", "", nptbins, ptbins);
                for (int pbin = 0; pbin<nptbins; pbin++) {
                    v112_p->SetBinContent(pbin+1, q112_p->GetBinContent(pbin+1, ebin+1));
                    v112_m->SetBinContent(pbin+1, q112_m->GetBinContent(pbin+1, ebin+1));
                }
                v112 = (TH1D *) v112_p->Clone("v112");
                v112->Add(v112_m);
                v112->Scale(0.5);

                v123_p = new TH1D("v123_p", "", nptbins, ptbins);
                v123_m = new TH1D("v123_m", "", nptbins, ptbins);
                for (int pbin = 0; pbin<nptbins; pbin++) {
                    v123_p->SetBinContent(pbin+1, q123_p->GetBinContent(pbin+1, ebin+1));
                    v123_m->SetBinContent(pbin+1, q123_m->GetBinContent(pbin+1, ebin+1));
                }
                v123 = (TH1D *) v123_p->Clone("v1");
                v123->Add(v123_m);
                v123->Scale(0.5);

                // calculate errors from subevents
                double v1tmp[20] = {0};
                double v1tmp_p[20] = {0};
                double v1tmp_m[20] = {0};
                double v1tmp2[20] = {0};
                double v1tmp_p2[20] = {0};
                double v1tmp_m2[20] = {0};

                double v112tmp[20] = {0};
                double v112tmp_p[20] = {0};
                double v112tmp_m[20] = {0};
                double v112tmp2[20] = {0};
                double v112tmp_p2[20] = {0};
                double v112tmp_m2[20] = {0};

                double v123tmp[20] = {0};
                double v123tmp_p[20] = {0};
                double v123tmp_m[20] = {0};
                double v123tmp2[20] = {0};
                double v123tmp_p2[20] = {0};
                double v123tmp_m2[20] = {0};

                TH1D * v1e_p;
                TH1D * v1e_m;
                TH1D * v1e;
                TH1D * v112e_p;
                TH1D * v112e_m;
                TH1D * v112e;
                TH1D * v123e_p;
                TH1D * v123e_m;
                TH1D * v123e;
                for (int i = 0; i<10; i++) {
                    v1e_p = new TH1D(Form("v1e_p%d_%d_%d_%d",anal,cbin,ebin,i), "", nptbins, ptbins);
                    v1e_m = new TH1D(Form("v1e_m%d_%d_%d_%d",anal,cbin,ebin,i), "", nptbins, ptbins);
                    for (int pbin = 0; pbin<nptbins; pbin++) {
                        v1e_p->SetBinContent(pbin+1, q1_p_err[i]->GetBinContent(pbin+1, ebin+1));
                        v1e_m->SetBinContent(pbin+1, q1_m_err[i]->GetBinContent(pbin+1, ebin+1));
                    }
                    v1e = (TH1D *) v1e_p->Clone(Form("v1e%d",i));
                    v1e->Scale(0.5);

                    v112e_p = new TH1D(Form("v112e_p%d_%d_%d_%d",anal,cbin,ebin,i), "", nptbins, ptbins);
                    v112e_m = new TH1D(Form("v112e_m%d_%d_%d_%d",anal,cbin,ebin,i), "", nptbins, ptbins);
                    for (int pbin = 0; pbin<nptbins; pbin++) {
                        v112e_p->SetBinContent(pbin+1, q112_p_err[i]->GetBinContent(pbin+1, ebin+1));
                        v112e_m->SetBinContent(pbin+1, q112_m_err[i]->GetBinContent(pbin+1, ebin+1));
                    }
                    v112e = (TH1D *) v112e_p->Clone(Form("v1e%d",i));
                    v112e->Scale(0.5);

                    v123e_p = new TH1D(Form("v123e_p%d_%d_%d_%d",anal,cbin,ebin,i), "", nptbins, ptbins);
                    v123e_m = new TH1D(Form("v123e_m%d_%d_%d_%d",anal,cbin,ebin,i), "", nptbins, ptbins);
                    for (int pbin = 0; pbin<nptbins; pbin++) {
                        v123e_p->SetBinContent(pbin+1, q123_p_err[i]->GetBinContent(pbin+1, ebin+1));
                        v123e_m->SetBinContent(pbin+1, q123_m_err[i]->GetBinContent(pbin+1, ebin+1));
                    }
                    v123e = (TH1D *) v123e_p->Clone(Form("v123e%d",i));
                    v123e->Scale(0.5);
                    for (int j = 0; j<v1e->GetNbinsX(); j++) {
                        v1tmp[j]+=v1e->GetBinContent(j+1);
                        v1tmp_p[j]+=v1e_p->GetBinContent(j+1);
                        v1tmp_m[j]+=v1e_m->GetBinContent(j+1);
                        v1tmp2[j]+=pow(v1e->GetBinContent(j+1),2);
                        v1tmp_p2[j]+=pow(v1e_p->GetBinContent(j+1),2);
                        v1tmp_m2[j]+=pow(v1e_m->GetBinContent(j+1),2);
                    }
                    for (int j = 0; j<v112e->GetNbinsX(); j++) {
                        v112tmp[j]+=v112e->GetBinContent(j+1);
                        v112tmp_p[j]+=v112e_p->GetBinContent(j+1);
                        v112tmp_m[j]+=v112e_m->GetBinContent(j+1);
                        v112tmp2[j]+=pow(v112e->GetBinContent(j+1),2);
                        v112tmp_p2[j]+=pow(v112e_p->GetBinContent(j+1),2);
                        v112tmp_m2[j]+=pow(v112e_m->GetBinContent(j+1),2);
                    }
                    for (int j = 0; j<v123e->GetNbinsX(); j++) {
                        v123tmp[j]+=v123e->GetBinContent(j+1);
                        v123tmp_p[j]+=v123e_p->GetBinContent(j+1);
                        v123tmp_m[j]+=v123e_m->GetBinContent(j+1);
                        v123tmp2[j]+=pow(v123e->GetBinContent(j+1),2);
                        v123tmp_p2[j]+=pow(v123e_p->GetBinContent(j+1),2);
                        v123tmp_m2[j]+=pow(v123e_m->GetBinContent(j+1),2);
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

                    v1p_pt[anal][cbin][ebin]->SetBinContent(i+1, v1_p->GetBinContent(i+1));
                    v1p_pt[anal][cbin][ebin]->SetBinError(i+1, v1_p->GetBinError(i+1));
                    v1m_pt[anal][cbin][ebin]->SetBinContent(i+1, v1_m->GetBinContent(i+1));
                    v1m_pt[anal][cbin][ebin]->SetBinError(i+1, v1_m->GetBinError(i+1));
                    v1odd_pt[anal][cbin][ebin]->SetBinContent(i+1, v1->GetBinContent(i+1));
                    v1odd_pt[anal][cbin][ebin]->SetBinError(i+1, v1->GetBinError(i+1));
                }
                v1even_pt[anal][cbin][ebin] = (TH1D *) v1_p->Clone(Form("v1even_pt_%s_c%d_e%d",AnalNames[anal].data(),cbin,ebin));
                v1even_pt[anal][cbin][ebin]->Add(v1_m,-1);
                v1even_pt[anal][cbin][ebin]->Scale(0.5);
                for (int i = 0; i<v112->GetNbinsX(); i++) {
                    v112tmp[i]/=10.;
                    v112tmp_p[i]/=10.;
                    v112tmp_m[i]/=10.;
                    v112tmp2[i]/=10.;
                    v112tmp_p2[i]/=10.;
                    v112tmp_m2[i]/=10.;
                    v112->SetBinError(  i+1, sqrt( (1./9.)*(v112tmp2[i]   - pow(v112tmp[i],2))   ));
                    v112_p->SetBinError(i+1, sqrt( (1./9.)*(v112tmp_p2[i] - pow(v112tmp_p[i],2)) ));
                    v112_m->SetBinError(i+1, sqrt( (1./9.)*(v112tmp_m2[i] - pow(v112tmp_m[i],2)) ));

                    v112p_pt[anal][cbin][ebin]->SetBinContent(i+1, v112_p->GetBinContent(i+1));
                    v112p_pt[anal][cbin][ebin]->SetBinError(i+1, v112_p->GetBinError(i+1));
                    v112m_pt[anal][cbin][ebin]->SetBinContent(i+1, v112_m->GetBinContent(i+1));
                    v112m_pt[anal][cbin][ebin]->SetBinError(i+1, v112_m->GetBinError(i+1));
                    v112odd_pt[anal][cbin][ebin]->SetBinContent(i+1, v112->GetBinContent(i+1));
                    v112odd_pt[anal][cbin][ebin]->SetBinError(i+1, v112->GetBinError(i+1));
                }
                v112even_pt[anal][cbin][ebin] = (TH1D *) v112_p->Clone(Form("v112odd_pt_%s_c%d_e%d",AnalNames[anal].data(),cbin,ebin));
                v112even_pt[anal][cbin][ebin]->Add(v112_m,-1);
                v112even_pt[anal][cbin][ebin]->Scale(0.5);
                for (int i = 0; i<v123->GetNbinsX(); i++) {
                    v123tmp[i]/=10.;
                    v123tmp_p[i]/=10.;
                    v123tmp_m[i]/=10.;
                    v123tmp2[i]/=10.;
                    v123tmp_p2[i]/=10.;
                    v123tmp_m2[i]/=10.;
                    v123->SetBinError(  i+1, sqrt( (1./9.)*(v123tmp2[i]   - pow(v123tmp[i],2))   ));
                    v123_p->SetBinError(i+1, sqrt( (1./9.)*(v123tmp_p2[i] - pow(v123tmp_p[i],2)) ));
                    v123_m->SetBinError(i+1, sqrt( (1./9.)*(v123tmp_m2[i] - pow(v123tmp_m[i],2)) ));

                    v123p_pt[anal][cbin][ebin]->SetBinContent(i+1, v123_p->GetBinContent(i+1));
                    v123p_pt[anal][cbin][ebin]->SetBinError(i+1, v123_p->GetBinError(i+1));
                    v123m_pt[anal][cbin][ebin]->SetBinContent(i+1, v123_m->GetBinContent(i+1));
                    v123m_pt[anal][cbin][ebin]->SetBinError(i+1, v123_m->GetBinError(i+1));
                    v123odd_pt[anal][cbin][ebin]->SetBinContent(i+1, v123->GetBinContent(i+1));
                    v123odd_pt[anal][cbin][ebin]->SetBinError(i+1, v123->GetBinError(i+1));
                }
                v123even_pt[anal][cbin][ebin] = (TH1D *) v123_p->Clone(Form("v123odd_pt_%s_c%d_e%d",AnalNames[anal].data(),cbin,ebin));
                v123even_pt[anal][cbin][ebin]->Add(v123_m,-1);
                v123even_pt[anal][cbin][ebin]->Scale(0.5);

                xpt->Delete();
                v1->Delete();
                v1_p->Delete();
                v1_m->Delete();
                v1e->Delete();
                v1e_p->Delete();
                v1e_m->Delete();

                v112->Delete();
                v112_p->Delete();
                v112_m->Delete();
                v112e->Delete();
                v112e_p->Delete();
                v112e_m->Delete();

                v123->Delete();
                v123_p->Delete();
                v123_m->Delete();
                v123e->Delete();
                v123e_p->Delete();
                v123e_m->Delete();
            }  // end v1(pT) loop


            //-- v1(eta)
            for (int pbin = 0; pbin<nptbins; pbin++) {
                TH1D * xeta = 0;
                TH1D * v1 = 0;
                TH1D * v1_p = 0;
                TH1D * v1_m = 0;
                TH1D * v112 = 0;
                TH1D * v112_p = 0;
                TH1D * v112_m = 0;
                TH1D * v123 = 0;
                TH1D * v123_p = 0;
                TH1D * v123_m = 0;

                xeta = (TH1D *) ptav->ProjectionY("xeta",pbin,(int)(pbin+0.001));
                v1_p = new TH1D("v1_p", "", netabins, etabins);
                v1_m = new TH1D("v1_m", "", netabins, etabins);
                for (int ebin = 0; ebin<netabins; ebin++) {
                    v1_p->SetBinContent(ebin+1, q1_p->GetBinContent(pbin+1, ebin+1));
                    v1_m->SetBinContent(ebin+1, q1_m->GetBinContent(pbin+1, ebin+1));
                }
                v1 = (TH1D *) v1_p->Clone("v1");
                v1->Add(v1_m);
                v1->Scale(0.5);

                v112_p = new TH1D("v112_p", "", netabins, etabins);
                v112_m = new TH1D("v112_m", "", netabins, etabins);
                for (int ebin = 0; ebin<netabins; ebin++) {
                    v112_p->SetBinContent(ebin+1, q112_p->GetBinContent(pbin+1, ebin+1));
                    v112_m->SetBinContent(ebin+1, q112_m->GetBinContent(pbin+1, ebin+1));
                }
                v112 = (TH1D *) v112_p->Clone("v112");
                v112->Add(v112_m);
                v112->Scale(0.5);

                v123_p = new TH1D("v123_p", "", netabins, etabins);
                v123_m = new TH1D("v123_m", "", netabins, etabins);
                for (int ebin = 0; ebin<netabins; ebin++) {
                    v123_p->SetBinContent(ebin+1, q123_p->GetBinContent(pbin+1, ebin+1));
                    v123_m->SetBinContent(ebin+1, q123_m->GetBinContent(pbin+1, ebin+1));
                }
                v123 = (TH1D *) v123_p->Clone("v123");
                v123->Add(v123_m);
                v123->Scale(0.5);

                // calculate errors from subevents
                double v1tmp[20] = {0};
                double v1tmp_p[20] = {0};
                double v1tmp_m[20] = {0};
                double v1tmp2[20] = {0};
                double v1tmp_p2[20] = {0};
                double v1tmp_m2[20] = {0};

                double v112tmp[20] = {0};
                double v112tmp_p[20] = {0};
                double v112tmp_m[20] = {0};
                double v112tmp2[20] = {0};
                double v112tmp_p2[20] = {0};
                double v112tmp_m2[20] = {0};

                double v123tmp[20] = {0};
                double v123tmp_p[20] = {0};
                double v123tmp_m[20] = {0};
                double v123tmp2[20] = {0};
                double v123tmp_p2[20] = {0};
                double v123tmp_m2[20] = {0};

                TH1D * v1e_p;
                TH1D * v1e_m;
                TH1D * v1e;
                TH1D * v112e_p;
                TH1D * v112e_m;
                TH1D * v112e;
                TH1D * v123e_p;
                TH1D * v123e_m;
                TH1D * v123e;
                for (int i = 0; i<10; i++) {
                    v1e_p = new TH1D(Form("v1e_p_eta%d_%d_%d_%d",anal,cbin,pbin,i), "", netabins, etabins);
                    v1e_m = new TH1D(Form("v1e_m_eta%d_%d_%d_%d",anal,cbin,pbin,i), "", netabins, etabins);
                    for (int ebin = 0; ebin<netabins; ebin++) {
                        v1e_p->SetBinContent(ebin+1, q1_p_err[i]->GetBinContent(pbin+1, ebin+1));
                        v1e_m->SetBinContent(ebin+1, q1_m_err[i]->GetBinContent(pbin+1, ebin+1));
                    }
                    v1e = (TH1D *) v1e_p->Clone(Form("v1e%d",i));
                    v1e->Scale(0.5);

                    v112e_p = new TH1D(Form("v112e_p_eta%d_%d_%d_%d",anal,cbin,pbin,i), "", netabins, etabins);
                    v112e_m = new TH1D(Form("v112e_m_eta%d_%d_%d_%d",anal,cbin,pbin,i), "", netabins, etabins);
                    for (int ebin = 0; ebin<netabins; ebin++) {
                        v112e_p->SetBinContent(ebin+1, q112_p_err[i]->GetBinContent(pbin+1, ebin+1));
                        v112e_m->SetBinContent(ebin+1, q112_m_err[i]->GetBinContent(pbin+1, ebin+1));
                    }
                    v112e = (TH1D *) v112e_p->Clone(Form("v112e%d",i));
                    v112e->Scale(0.5);

                    v123e_p = new TH1D(Form("v123e_p_eta%d_%d_%d_%d",anal,cbin,pbin,i), "", netabins, etabins);
                    v123e_m = new TH1D(Form("v123e_m_eta%d_%d_%d_%d",anal,cbin,pbin,i), "", netabins, etabins);
                    for (int ebin = 0; ebin<netabins; ebin++) {
                        v123e_p->SetBinContent(ebin+1, q123_p_err[i]->GetBinContent(pbin+1, ebin+1));
                        v123e_m->SetBinContent(ebin+1, q123_m_err[i]->GetBinContent(pbin+1, ebin+1));
                    }
                    v123e = (TH1D *) v123e_p->Clone(Form("v123e%d",i));
                    v123e->Scale(0.5);
                    for (int j = 0; j<v1e->GetNbinsX(); j++) {
                        v1tmp[j]+=v1e->GetBinContent(j+1);
                        v1tmp_p[j]+=v1e_p->GetBinContent(j+1);
                        v1tmp_m[j]+=v1e_m->GetBinContent(j+1);
                        v1tmp2[j]+=pow(v1e->GetBinContent(j+1),2);
                        v1tmp_p2[j]+=pow(v1e_p->GetBinContent(j+1),2);
                        v1tmp_m2[j]+=pow(v1e_m->GetBinContent(j+1),2);
                    }
                    for (int j = 0; j<v112e->GetNbinsX(); j++) {
                        v112tmp[j]+=v112e->GetBinContent(j+1);
                        v112tmp_p[j]+=v112e_p->GetBinContent(j+1);
                        v112tmp_m[j]+=v112e_m->GetBinContent(j+1);
                        v112tmp2[j]+=pow(v112e->GetBinContent(j+1),2);
                        v112tmp_p2[j]+=pow(v112e_p->GetBinContent(j+1),2);
                        v112tmp_m2[j]+=pow(v112e_m->GetBinContent(j+1),2);
                    }
                    for (int j = 0; j<v123e->GetNbinsX(); j++) {
                        v123tmp[j]+=v123e->GetBinContent(j+1);
                        v123tmp_p[j]+=v123e_p->GetBinContent(j+1);
                        v123tmp_m[j]+=v123e_m->GetBinContent(j+1);
                        v123tmp2[j]+=pow(v123e->GetBinContent(j+1),2);
                        v123tmp_p2[j]+=pow(v123e_p->GetBinContent(j+1),2);
                        v123tmp_m2[j]+=pow(v123e_m->GetBinContent(j+1),2);
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

                    v1p_eta[anal][cbin][pbin]->SetBinContent(i+1, v1_p->GetBinContent(i+1));
                    v1p_eta[anal][cbin][pbin]->SetBinError(i+1, v1_p->GetBinError(i+1));
                    v1m_eta[anal][cbin][pbin]->SetBinContent(i+1, v1_m->GetBinContent(i+1));
                    v1m_eta[anal][cbin][pbin]->SetBinError(i+1, v1_m->GetBinError(i+1));
                    v1odd_eta[anal][cbin][pbin]->SetBinContent(i+1, v1->GetBinContent(i+1));
                    v1odd_eta[anal][cbin][pbin]->SetBinError(i+1, v1->GetBinError(i+1));
                }
                v1even_eta[anal][cbin][pbin] = (TH1D *) v1_p->Clone(Form("v1even_eta_%s_c%d_p%d",AnalNames[anal].data(),cbin,pbin));
                v1even_eta[anal][cbin][pbin]->Add(v1_m,-1);
                v1even_eta[anal][cbin][pbin]->Scale(0.5);
                for (int i = 0; i<v112->GetNbinsX(); i++) {
                    v112tmp[i]/=10.;
                    v112tmp_p[i]/=10.;
                    v112tmp_m[i]/=10.;
                    v112tmp2[i]/=10.;
                    v112tmp_p2[i]/=10.;
                    v112tmp_m2[i]/=10.;
                    v112->SetBinError(  i+1, sqrt( (1./9.)*(v112tmp2[i]   - pow(v112tmp[i],2))   ));
                    v112_p->SetBinError(i+1, sqrt( (1./9.)*(v112tmp_p2[i] - pow(v112tmp_p[i],2)) ));
                    v112_m->SetBinError(i+1, sqrt( (1./9.)*(v112tmp_m2[i] - pow(v112tmp_m[i],2)) ));

                    v112p_eta[anal][cbin][pbin]->SetBinContent(i+1, v112_p->GetBinContent(i+1));
                    v112p_eta[anal][cbin][pbin]->SetBinError(i+1, v112_p->GetBinError(i+1));
                    v112m_eta[anal][cbin][pbin]->SetBinContent(i+1, v112_m->GetBinContent(i+1));
                    v112m_eta[anal][cbin][pbin]->SetBinError(i+1, v112_m->GetBinError(i+1));
                    v112odd_eta[anal][cbin][pbin]->SetBinContent(i+1, v112->GetBinContent(i+1));
                    v112odd_eta[anal][cbin][pbin]->SetBinError(i+1, v112->GetBinError(i+1));
                }
                v112even_eta[anal][cbin][pbin] = (TH1D *) v112_p->Clone(Form("v112even_eta_%s_c%d_p%d",AnalNames[anal].data(),cbin,pbin));
                v112even_eta[anal][cbin][pbin]->Add(v112_m,-1);
                v112even_eta[anal][cbin][pbin]->Scale(0.5);
                for (int i = 0; i<v123->GetNbinsX(); i++) {
                    v123tmp[i]/=10.;
                    v123tmp_p[i]/=10.;
                    v123tmp_m[i]/=10.;
                    v123tmp2[i]/=10.;
                    v123tmp_p2[i]/=10.;
                    v123tmp_m2[i]/=10.;
                    v123->SetBinError(  i+1, sqrt( (1./9.)*(v123tmp2[i]   - pow(v123tmp[i],2))   ));
                    v123_p->SetBinError(i+1, sqrt( (1./9.)*(v123tmp_p2[i] - pow(v123tmp_p[i],2)) ));
                    v123_m->SetBinError(i+1, sqrt( (1./9.)*(v123tmp_m2[i] - pow(v123tmp_m[i],2)) ));

                    v123p_eta[anal][cbin][pbin]->SetBinContent(i+1, v123_p->GetBinContent(i+1));
                    v123p_eta[anal][cbin][pbin]->SetBinError(i+1, v123_p->GetBinError(i+1));
                    v123m_eta[anal][cbin][pbin]->SetBinContent(i+1, v123_m->GetBinContent(i+1));
                    v123m_eta[anal][cbin][pbin]->SetBinError(i+1, v123_m->GetBinError(i+1));
                    v123odd_eta[anal][cbin][pbin]->SetBinContent(i+1, v123->GetBinContent(i+1));
                    v123odd_eta[anal][cbin][pbin]->SetBinError(i+1, v123->GetBinError(i+1));
                }
                v123even_eta[anal][cbin][pbin] = (TH1D *) v123_p->Clone(Form("v123even_eta_%s_c%d_p%d",AnalNames[anal].data(),cbin,pbin));
                v123even_eta[anal][cbin][pbin]->Add(v123_m,-1);
                v123even_eta[anal][cbin][pbin]->Scale(0.5);

                v1->Delete();
                v1_p->Delete();
                v1_m->Delete();
                v1e->Delete();
                v1e_p->Delete();
                v1e_m->Delete();
                v112->Delete();
                v112_p->Delete();
                v112_m->Delete();
                v112e->Delete();
                v112e_p->Delete();
                v112e_m->Delete();
                v123->Delete();
                v123_p->Delete();
                v123_m->Delete();
                v123e->Delete();
                v123e_p->Delete();
                v123e_m->Delete();
            }  // end v1(eta) loop
        }  // end cent loop
        gROOT->GetListOfFiles()->Remove(tfin);
    }  // end analysis loop


    // write histograms to output file
    if (!fopen("../outputs","r")) system("mkdir ../outputs");
    if (!fopen("../outputs/final_outputs","r")) system("mkdir ../outputs/final_outputs");
    TFile * tfout = new TFile("../outputs/final_outputs/v1Diff.root","recreate");

    for (int i = minanal; i<maxanal; i++) {
        TDirectory * tdAnal = (TDirectory *) tfout->mkdir(Form("%s",AnalNames[i].data()));
        TDirectory * tdPt = (TDirectory *) tdAnal->mkdir("v1_pt");
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TDirectory * tdPtCent = (TDirectory *) tdPt->mkdir(Form("cent_%d-%d",centBins[cbin],centBins[cbin+1]));
            for (int ebin = 0; ebin<netabins; ebin++) {
                TDirectory * tdPtCentEtaRange = (TDirectory *) tdPtCent->mkdir(Form("eta_%0.1f-%0.1f",etabins[ebin],etabins[ebin+1]));
                tdPtCentEtaRange->cd();
                v1p_pt[i][cbin][ebin]->Write();
                v1m_pt[i][cbin][ebin]->Write();
                v1odd_pt[i][cbin][ebin]->Write();
                v1even_pt[i][cbin][ebin]->Write();
                v112p_pt[i][cbin][ebin]->Write();
                v112m_pt[i][cbin][ebin]->Write();
                v112odd_pt[i][cbin][ebin]->Write();
                v112even_pt[i][cbin][ebin]->Write();
                v123p_pt[i][cbin][ebin]->Write();
                v123m_pt[i][cbin][ebin]->Write();
                v123odd_pt[i][cbin][ebin]->Write();
                v123even_pt[i][cbin][ebin]->Write();
            }
        }
        TDirectory * tdEta = (TDirectory *) tdAnal->mkdir("v1_eta");
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TDirectory * tdEtaCent = (TDirectory *) tdEta->mkdir(Form("cent_%d-%d",centBins[cbin],centBins[cbin+1]));
            for (int pbin = 0; pbin<nptbins; pbin++) {
                TDirectory * tdEtaCentPtRange = (TDirectory *) tdEtaCent->mkdir(Form("pt_%0.1f-%0.1f",ptbins[pbin],ptbins[pbin+1]));
                tdEtaCentPtRange->cd();
                v1p_eta[i][cbin][pbin]->Write();
                v1m_eta[i][cbin][pbin]->Write();
                v1odd_eta[i][cbin][pbin]->Write();
                v1even_eta[i][cbin][pbin]->Write();
                v112p_eta[i][cbin][pbin]->Write();
                v112m_eta[i][cbin][pbin]->Write();
                v112odd_eta[i][cbin][pbin]->Write();
                v112even_eta[i][cbin][pbin]->Write();
                v123p_eta[i][cbin][pbin]->Write();
                v123m_eta[i][cbin][pbin]->Write();
                v123odd_eta[i][cbin][pbin]->Write();
                v123even_eta[i][cbin][pbin]->Write();
            }
        }
        tdAnal->cd();
        runParms[i]->Write();
    }

    cout << "\n ...Differential v1 results written out to ../outputs/final_outputs/v1Diff.root \n" << endl;

}
