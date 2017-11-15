# include "TArrayD.h"
# include "TDirectory.h"
# include "TFile.h"
# include "TH1.h"
# include "TH2.h"
# include "TList.h"
# include "TMath.h"
# include "TRandom3.h"
# include "TStopwatch.h"
# include "TString.h"
# include "TTree.h"
# include <cmath>
# include <iostream>
# include <string>
# include <vector>

# include "src/HiEvtPlaneList.h"

using namespace hi;

TRandom3 * ran;

static const double MaxCent = 80;
static const int MaxEvents = -1;
int MaxFiles = 3000;
static const double VtxCut = 15.;
static const int ncentbins = 8;
static const double centBins[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
static const int nanals = 16;
string AnalNames[] = {
    "v1SP",   "v1SP_mid",   "v1SP_102",   "v1SP_106",   "v1SP_110",   "v1SP_114",   "v1SP_118",   "v1SP_122",
    "v1SPmc", "v1SPmc_mid", "v1SPmc_102", "v1SPmc_106", "v1SPmc_110", "v1SPmc_114", "v1SPmc_118", "v1SPmc_122"
};

string rpnames[hi::NumEPNames];
int epA1p, epB1p, epC1p, epA1m, epB1m, epC1m;
int epA2p, epB2p, epC2p, epA2m, epB2m, epC2m;
int epA3p, epB3p, epC3p, epA3m, epB3m, epC3m;

//----------------------------------
// Tree variables
//
double centval;
int ntrks;
double vtx;
Double_t epang[hi::NumEPNames];
Double_t qx[hi::NumEPNames];
Double_t qy[hi::NumEPNames];
Double_t epmult[hi::NumEPNames];
Double_t sumw[hi::NumEPNames];
TH2D * qxtrk1_;
TH2D * qytrk1_;
TH2D * qcnt_;
TH2D * qw_;
TH2D * avpt_;

Int_t NumEvents;
Int_t TotNumEvents;
//----------------------------------

TH1D * centbins;
TH2D * ptav[ncentbins];
TH2D * ptcnt[ncentbins];
TH2D * badcnt[ncentbins];
TH1D * multTot[ncentbins][11];
TH2D * q1_p[ncentbins][11];
TH2D * q1_m[ncentbins][11];
TH2D * q1_pm[ncentbins][11];
TH2D * q112_p[ncentbins][11];
TH2D * q112_m[ncentbins][11];
TH2D * q112_pm[ncentbins][11];
TH2D * q123_p[ncentbins][11];
TH2D * q123_m[ncentbins][11];
TH2D * q123_pm[ncentbins][11];
TH2D * w1_p[ncentbins][11];
TH2D * w1_m[ncentbins][11];
TH2D * w1_pm[ncentbins][11];
TH2D * w112_p[ncentbins][11];
TH2D * w112_m[ncentbins][11];
TH2D * w112_pm[ncentbins][11];
TH2D * w123_p[ncentbins][11];
TH2D * w123_m[ncentbins][11];
TH2D * w123_pm[ncentbins][11];
TH2D * qxav1[ncentbins];
TH2D * qyav1[ncentbins];
TH2D * qxycnt[ncentbins];

TH1D * q1AB_p[ncentbins][11];
TH1D * q1AC_p[ncentbins][11];
TH1D * q1BC_p[ncentbins][11];
TH1D * q1AB_m[ncentbins][11];
TH1D * q1AC_m[ncentbins][11];
TH1D * q1BC_m[ncentbins][11];
TH1D * q1ABcnt_p[ncentbins][11];
TH1D * q1ACcnt_p[ncentbins][11];
TH1D * q1BCcnt_p[ncentbins][11];
TH1D * q1ABcnt_m[ncentbins][11];
TH1D * q1ACcnt_m[ncentbins][11];
TH1D * q1BCcnt_m[ncentbins][11];

TH1D * q2AB_p[ncentbins][11];
TH1D * q2AC_p[ncentbins][11];
TH1D * q2BC_p[ncentbins][11];
TH1D * q2AB_m[ncentbins][11];
TH1D * q2AC_m[ncentbins][11];
TH1D * q2BC_m[ncentbins][11];
TH1D * q2ABcnt_p[ncentbins][11];
TH1D * q2ACcnt_p[ncentbins][11];
TH1D * q2BCcnt_p[ncentbins][11];
TH1D * q2ABcnt_m[ncentbins][11];
TH1D * q2ACcnt_m[ncentbins][11];
TH1D * q2BCcnt_m[ncentbins][11];

TH1D * q3AB_p[ncentbins][11];
TH1D * q3AC_p[ncentbins][11];
TH1D * q3BC_p[ncentbins][11];
TH1D * q3AB_m[ncentbins][11];
TH1D * q3AC_m[ncentbins][11];
TH1D * q3BC_m[ncentbins][11];
TH1D * q3ABcnt_p[ncentbins][11];
TH1D * q3ACcnt_p[ncentbins][11];
TH1D * q3BCcnt_p[ncentbins][11];
TH1D * q3ABcnt_m[ncentbins][11];
TH1D * q3ACcnt_m[ncentbins][11];
TH1D * q3BCcnt_m[ncentbins][11];

TH1D * runParms;

# include "src/GetEventInfo.h"
# include "src/Qcalc.h"
# include "src/EPSetup.h"

void ReadTree( GetEventInfo * info, string anal, string inlist );
void GetNumEvents( string inlist );
void GenerateV1( string anal = "", string inlist = "", int FileLimit = 100 ) {

    MaxFiles = FileLimit;
    cout << "\nStarting analysis for " << anal << "\n" << endl;
    FILE * flist;
    flist = fopen(inlist.data(),"r");
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    ran = new TRandom3();

    //-- locate information about data structure
    GetEventInfo * info = 0;
    char buf[120];
    while (fgets(buf,120,flist)!=NULL) {
        buf[strlen(buf) - 1] = 0;
        TString inFile = buf;
        FILE * ftest = fopen(inFile.Data(),"r");
        if (ftest == NULL) continue;
        fclose(ftest);
        cout << "Found: " << inFile.Data() << endl;
        info = new GetEventInfo(inFile);
        if (info->status == 0) {cout << inFile.Data() << " not found or has an error" << endl; return;}

        EPSetup( anal );

        centbins = new TH1D("centbins", "centbins", ncentbins, centBins);
        centbins->SetDirectory(0);

        for (int cbin = 0; cbin<ncentbins; cbin++) {
            ptav[cbin] = (TH2D *) info->getTemplate()->Clone(Form("ptav_%d",cbin));
            ptav[cbin]->SetDirectory(0);
            ptcnt[cbin] = (TH2D *) info->getTemplate()->Clone(Form("ptcnt_%d",cbin));
            ptcnt[cbin]->SetDirectory(0);
            badcnt[cbin] = (TH2D *) info->getTemplate()->Clone(Form("badcnt_%d",cbin));
            multTot[cbin][0] = new TH1D(Form("multTot_%d",cbin), "", 200, 0, 4000);
            multTot[cbin][0]->SetDirectory(0);
            q1_p[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("q1_p_%d",cbin));
            q1_p[cbin][0]->SetDirectory(0);
            q1_m[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("q1_m_%d",cbin));
            q1_m[cbin][0]->SetDirectory(0);
            q1_pm[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("q1_pm_%d",cbin));
            q1_pm[cbin][0]->SetDirectory(0);
            q112_p[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("q112_p_%d",cbin));
            q112_p[cbin][0]->SetDirectory(0);
            q112_m[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("q112_m_%d",cbin));
            q112_m[cbin][0]->SetDirectory(0);
            q112_pm[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("q112_pm_%d",cbin));
            q112_pm[cbin][0]->SetDirectory(0);
            q123_p[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("q123_p_%d",cbin));
            q123_p[cbin][0]->SetDirectory(0);
            q123_m[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("q123_m_%d",cbin));
            q123_m[cbin][0]->SetDirectory(0);
            q123_pm[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("q123_pm_%d",cbin));
            q123_pm[cbin][0]->SetDirectory(0);
            w1_p[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("w1_p_%d",cbin));
            w1_p[cbin][0]->SetDirectory(0);
            w1_m[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("w1_m_%d",cbin));
            w1_m[cbin][0]->SetDirectory(0);
            w1_pm[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("w1_pm_%d",cbin));
            w1_pm[cbin][0]->SetDirectory(0);
            w112_p[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("w112_p_%d",cbin));
            w112_p[cbin][0]->SetDirectory(0);
            w112_m[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("w112_m_%d",cbin));
            w112_m[cbin][0]->SetDirectory(0);
            w112_pm[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("w112_pm_%d",cbin));
            w112_pm[cbin][0]->SetDirectory(0);
            w123_p[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("w123_p_%d",cbin));
            w123_p[cbin][0]->SetDirectory(0);
            w123_m[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("w123_m_%d",cbin));
            w123_m[cbin][0]->SetDirectory(0);
            w123_pm[cbin][0] = (TH2D *) info->getTemplate()->Clone(Form("w123_pm_%d",cbin));
            w123_pm[cbin][0]->SetDirectory(0);
            qxav1[cbin] = (TH2D *) info->getTemplate()->Clone(Form("qxav1_%d",cbin));
            qxav1[cbin]->SetDirectory(0);
            qyav1[cbin] = (TH2D *) info->getTemplate()->Clone(Form("qyav1_%d",cbin));
            qyav1[cbin]->SetDirectory(0);
            qxycnt[cbin] = (TH2D *) info->getTemplate()->Clone(Form("qxycnt_%d",cbin));
            qxycnt[cbin]->SetDirectory(0);

            q1AB_p[cbin][0] = new TH1D(Form("q1AB_p_%d",cbin), "", 1, 0, 1);
            q1AB_p[cbin][0]->SetDirectory(0);
            q1AC_p[cbin][0] = new TH1D(Form("q1AC_p_%d",cbin), "", 1, 0, 1);
            q1AC_p[cbin][0]->SetDirectory(0);
            q1BC_p[cbin][0] = new TH1D(Form("q1BC_p_%d",cbin), "", 1, 0, 1);
            q1BC_p[cbin][0]->SetDirectory(0);
            q1AB_m[cbin][0] = new TH1D(Form("q1AB_m_%d",cbin), "", 1, 0, 1);
            q1AB_m[cbin][0]->SetDirectory(0);
            q1AC_m[cbin][0] = new TH1D(Form("q1AC_m_%d",cbin), "", 1, 0, 1);
            q1AC_m[cbin][0]->SetDirectory(0);
            q1BC_m[cbin][0] = new TH1D(Form("q1BC_m_%d",cbin), "", 1, 0, 1);
            q1BC_m[cbin][0]->SetDirectory(0);
            q1ABcnt_p[cbin][0] = new TH1D(Form("q1ABcnt_p_%d",cbin), "", 1, 0, 1);
            q1ABcnt_p[cbin][0]->SetDirectory(0);
            q1ACcnt_p[cbin][0] = new TH1D(Form("q1ACcnt_p_%d",cbin), "", 1, 0, 1);
            q1ACcnt_p[cbin][0]->SetDirectory(0);
            q1BCcnt_p[cbin][0] = new TH1D(Form("q1BCcnt_p_%d",cbin), "", 1, 0, 1);
            q1BCcnt_p[cbin][0]->SetDirectory(0);
            q1ABcnt_m[cbin][0] = new TH1D(Form("q1ABcnt_m_%d",cbin), "", 1, 0, 1);
            q1ABcnt_m[cbin][0]->SetDirectory(0);
            q1ACcnt_m[cbin][0] = new TH1D(Form("q1ACcnt_m_%d",cbin), "", 1, 0, 1);
            q1ACcnt_m[cbin][0]->SetDirectory(0);
            q1BCcnt_m[cbin][0] = new TH1D(Form("q1BCcnt_m_%d",cbin), "", 1, 0, 1);
            q1BCcnt_m[cbin][0]->SetDirectory(0);

            q2AB_p[cbin][0] = new TH1D(Form("q2AB_p_%d",cbin), "", 1, 0, 1);
            q2AB_p[cbin][0]->SetDirectory(0);
            q2AC_p[cbin][0] = new TH1D(Form("q2AC_p_%d",cbin), "", 1, 0, 1);
            q2AC_p[cbin][0]->SetDirectory(0);
            q2BC_p[cbin][0] = new TH1D(Form("q2BC_p_%d",cbin), "", 1, 0, 1);
            q2BC_p[cbin][0]->SetDirectory(0);
            q2AB_m[cbin][0] = new TH1D(Form("q2AB_m_%d",cbin), "", 1, 0, 1);
            q2AB_m[cbin][0]->SetDirectory(0);
            q2AC_m[cbin][0] = new TH1D(Form("q2AC_m_%d",cbin), "", 1, 0, 1);
            q2AC_m[cbin][0]->SetDirectory(0);
            q2BC_m[cbin][0] = new TH1D(Form("q2BC_m_%d",cbin), "", 1, 0, 1);
            q2BC_m[cbin][0]->SetDirectory(0);
            q2ABcnt_p[cbin][0] = new TH1D(Form("q2ABcnt_p_%d",cbin), "", 1, 0, 1);
            q2ABcnt_p[cbin][0]->SetDirectory(0);
            q2ACcnt_p[cbin][0] = new TH1D(Form("q2ACcnt_p_%d",cbin), "", 1, 0, 1);
            q2ACcnt_p[cbin][0]->SetDirectory(0);
            q2BCcnt_p[cbin][0] = new TH1D(Form("q2BCcnt_p_%d",cbin), "", 1, 0, 1);
            q2BCcnt_p[cbin][0]->SetDirectory(0);
            q2ABcnt_m[cbin][0] = new TH1D(Form("q2ABcnt_m_%d",cbin), "", 1, 0, 1);
            q2ABcnt_m[cbin][0]->SetDirectory(0);
            q2ACcnt_m[cbin][0] = new TH1D(Form("q2ACcnt_m_%d",cbin), "", 1, 0, 1);
            q2ACcnt_m[cbin][0]->SetDirectory(0);
            q2BCcnt_m[cbin][0] = new TH1D(Form("q2BCcnt_m_%d",cbin), "", 1, 0, 1);
            q2BCcnt_m[cbin][0]->SetDirectory(0);

            q3AB_p[cbin][0] = new TH1D(Form("q3AB_p_%d",cbin), "", 1, 0, 1);
            q3AB_p[cbin][0]->SetDirectory(0);
            q3AC_p[cbin][0] = new TH1D(Form("q3AC_p_%d",cbin), "", 1, 0, 1);
            q3AC_p[cbin][0]->SetDirectory(0);
            q3BC_p[cbin][0] = new TH1D(Form("q3BC_p_%d",cbin), "", 1, 0, 1);
            q3BC_p[cbin][0]->SetDirectory(0);
            q3AB_m[cbin][0] = new TH1D(Form("q3AB_m_%d",cbin), "", 1, 0, 1);
            q3AB_m[cbin][0]->SetDirectory(0);
            q3AC_m[cbin][0] = new TH1D(Form("q3AC_m_%d",cbin), "", 1, 0, 1);
            q3AC_m[cbin][0]->SetDirectory(0);
            q3BC_m[cbin][0] = new TH1D(Form("q3BC_m_%d",cbin), "", 1, 0, 1);
            q3BC_m[cbin][0]->SetDirectory(0);
            q3ABcnt_p[cbin][0] = new TH1D(Form("q3ABcnt_p_%d",cbin), "", 1, 0, 1);
            q3ABcnt_p[cbin][0]->SetDirectory(0);
            q3ACcnt_p[cbin][0] = new TH1D(Form("q3ACcnt_p_%d",cbin), "", 1, 0, 1);
            q3ACcnt_p[cbin][0]->SetDirectory(0);
            q3BCcnt_p[cbin][0] = new TH1D(Form("q3BCcnt_p_%d",cbin), "", 1, 0, 1);
            q3BCcnt_p[cbin][0]->SetDirectory(0);
            q3ABcnt_m[cbin][0] = new TH1D(Form("q3ABcnt_m_%d",cbin), "", 1, 0, 1);
            q3ABcnt_m[cbin][0]->SetDirectory(0);
            q3ACcnt_m[cbin][0] = new TH1D(Form("q3ACcnt_m_%d",cbin), "", 1, 0, 1);
            q3ACcnt_m[cbin][0]->SetDirectory(0);
            q3BCcnt_m[cbin][0] = new TH1D(Form("q3BCcnt_m_%d",cbin), "", 1, 0, 1);
            q3BCcnt_m[cbin][0]->SetDirectory(0);

            for (int k = 1; k<=10; k++) {
                multTot[cbin][k] = new TH1D(Form("multTot_%d_%d",cbin,k), "", 200, 0, 4000);
                multTot[cbin][k]->SetDirectory(0);
                q1_p[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("q1_p_%d_%d",cbin,k));
                q1_p[cbin][k]->SetDirectory(0);
                q1_m[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("q1_m_%d_%d",cbin,k));
                q1_m[cbin][k]->SetDirectory(0);
                q1_pm[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("q1_pm_%d_%d",cbin,k));
                q1_pm[cbin][k]->SetDirectory(0);
                q112_p[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("q112_p_%d_%d",cbin,k));
                q112_p[cbin][k]->SetDirectory(0);
                q112_m[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("q112_m_%d_%d",cbin,k));
                q112_m[cbin][k]->SetDirectory(0);
                q112_pm[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("q112_pm_%d_%d",cbin,k));
                q112_pm[cbin][k]->SetDirectory(0);
                q123_p[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("q123_p_%d_%d",cbin,k));
                q123_p[cbin][k]->SetDirectory(0);
                q123_m[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("q123_m_%d_%d",cbin,k));
                q123_m[cbin][k]->SetDirectory(0);
                q123_pm[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("q123_pm_%d_%d",cbin,k));
                q123_pm[cbin][k]->SetDirectory(0);
                w1_p[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("w1_p_%d_%d",cbin,k));
                w1_p[cbin][k]->SetDirectory(0);
                w1_m[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("w1_m_%d_%d",cbin,k));
                w1_m[cbin][k]->SetDirectory(0);
                w1_pm[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("w1_pm_%d_%d",cbin,k));
                w1_pm[cbin][k]->SetDirectory(0);
                w112_p[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("w112_p_%d_%d",cbin,k));
                w112_p[cbin][k]->SetDirectory(0);
                w112_m[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("w112_m_%d_%d",cbin,k));
                w112_m[cbin][k]->SetDirectory(0);
                w112_pm[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("w112_pm_%d_%d",cbin,k));
                w112_pm[cbin][k]->SetDirectory(0);
                w123_p[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("w123_p_%d_%d",cbin,k));
                w123_p[cbin][k]->SetDirectory(0);
                w123_m[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("w123_m_%d_%d",cbin,k));
                w123_m[cbin][k]->SetDirectory(0);
                w123_pm[cbin][k] = (TH2D *) info->getTemplate()->Clone(Form("w123_pm_%d_%d",cbin,k));
                w123_pm[cbin][k]->SetDirectory(0);

                q1AB_p[cbin][k] = new TH1D(Form("q1AB_p_%d_%d",cbin,k), "", 1, 0, 1);
                q1AB_p[cbin][k]->SetDirectory(0);
                q1AC_p[cbin][k] = new TH1D(Form("q1AC_p_%d_%d",cbin,k), "", 1, 0, 1);
                q1AC_p[cbin][k]->SetDirectory(0);
                q1BC_p[cbin][k] = new TH1D(Form("q1BC_p_%d_%d",cbin,k), "", 1, 0, 1);
                q1BC_p[cbin][k]->SetDirectory(0);
                q1AB_m[cbin][k] = new TH1D(Form("q1AB_m_%d_%d",cbin,k), "", 1, 0, 1);
                q1AB_m[cbin][k]->SetDirectory(0);
                q1AC_m[cbin][k] = new TH1D(Form("q1AC_m_%d_%d",cbin,k), "", 1, 0, 1);
                q1AC_m[cbin][k]->SetDirectory(0);
                q1BC_m[cbin][k] = new TH1D(Form("q1BC_m_%d_%d",cbin,k), "", 1, 0, 1);
                q1BC_m[cbin][k]->SetDirectory(0);
                q1ABcnt_p[cbin][k] = new TH1D(Form("q1ABcnt_p_%d_%d",cbin,k), "", 1, 0, 1);
                q1ABcnt_p[cbin][k]->SetDirectory(0);
                q1ACcnt_p[cbin][k] = new TH1D(Form("q1ACcnt_p_%d_%d",cbin,k), "", 1, 0, 1);
                q1ACcnt_p[cbin][k]->SetDirectory(0);
                q1BCcnt_p[cbin][k] = new TH1D(Form("q1BCcnt_p_%d_%d",cbin,k), "", 1, 0, 1);
                q1BCcnt_p[cbin][k]->SetDirectory(0);
                q1ABcnt_m[cbin][k] = new TH1D(Form("q1ABcnt_m_%d_%d",cbin,k), "", 1, 0, 1);
                q1ABcnt_m[cbin][k]->SetDirectory(0);
                q1ACcnt_m[cbin][k] = new TH1D(Form("q1ACcnt_m_%d_%d",cbin,k), "", 1, 0, 1);
                q1ACcnt_m[cbin][k]->SetDirectory(0);
                q1BCcnt_m[cbin][k] = new TH1D(Form("q1BCcnt_m_%d_%d",cbin,k), "", 1, 0, 1);
                q1BCcnt_m[cbin][k]->SetDirectory(0);

                q2AB_p[cbin][k] = new TH1D(Form("q2AB_p_%d_%d",cbin,k), "", 1, 0, 1);
                q2AB_p[cbin][k]->SetDirectory(0);
                q2AC_p[cbin][k] = new TH1D(Form("q2AC_p_%d_%d",cbin,k), "", 1, 0, 1);
                q2AC_p[cbin][k]->SetDirectory(0);
                q2BC_p[cbin][k] = new TH1D(Form("q2BC_p_%d_%d",cbin,k), "", 1, 0, 1);
                q2BC_p[cbin][k]->SetDirectory(0);
                q2AB_m[cbin][k] = new TH1D(Form("q2AB_m_%d_%d",cbin,k), "", 1, 0, 1);
                q2AB_m[cbin][k]->SetDirectory(0);
                q2AC_m[cbin][k] = new TH1D(Form("q2AC_m_%d_%d",cbin,k), "", 1, 0, 1);
                q2AC_m[cbin][k]->SetDirectory(0);
                q2BC_m[cbin][k] = new TH1D(Form("q2BC_m_%d_%d",cbin,k), "", 1, 0, 1);
                q2BC_m[cbin][k]->SetDirectory(0);
                q2ABcnt_p[cbin][k] = new TH1D(Form("q2ABcnt_p_%d_%d",cbin,k), "", 1, 0, 1);
                q2ABcnt_p[cbin][k]->SetDirectory(0);
                q2ACcnt_p[cbin][k] = new TH1D(Form("q2ACcnt_p_%d_%d",cbin,k), "", 1, 0, 1);
                q2ACcnt_p[cbin][k]->SetDirectory(0);
                q2BCcnt_p[cbin][k] = new TH1D(Form("q2BCcnt_p_%d_%d",cbin,k), "", 1, 0, 1);
                q2BCcnt_p[cbin][k]->SetDirectory(0);
                q2ABcnt_m[cbin][k] = new TH1D(Form("q2ABcnt_m_%d_%d",cbin,k), "", 1, 0, 1);
                q2ABcnt_m[cbin][k]->SetDirectory(0);
                q2ACcnt_m[cbin][k] = new TH1D(Form("q2ACcnt_m_%d_%d",cbin,k), "", 1, 0, 1);
                q2ACcnt_m[cbin][k]->SetDirectory(0);
                q2BCcnt_m[cbin][k] = new TH1D(Form("q2BCcnt_m_%d_%d",cbin,k), "", 1, 0, 1);
                q2BCcnt_m[cbin][k]->SetDirectory(0);

                q3AB_p[cbin][k] = new TH1D(Form("q3AB_p_%d_%d",cbin,k), "", 1, 0, 1);
                q3AB_p[cbin][k]->SetDirectory(0);
                q3AC_p[cbin][k] = new TH1D(Form("q3AC_p_%d_%d",cbin,k), "", 1, 0, 1);
                q3AC_p[cbin][k]->SetDirectory(0);
                q3BC_p[cbin][k] = new TH1D(Form("q3BC_p_%d_%d",cbin,k), "", 1, 0, 1);
                q3BC_p[cbin][k]->SetDirectory(0);
                q3AB_m[cbin][k] = new TH1D(Form("q3AB_m_%d_%d",cbin,k), "", 1, 0, 1);
                q3AB_m[cbin][k]->SetDirectory(0);
                q3AC_m[cbin][k] = new TH1D(Form("q3AC_m_%d_%d",cbin,k), "", 1, 0, 1);
                q3AC_m[cbin][k]->SetDirectory(0);
                q3BC_m[cbin][k] = new TH1D(Form("q3BC_m_%d_%d",cbin,k), "", 1, 0, 1);
                q3BC_m[cbin][k]->SetDirectory(0);
                q3ABcnt_p[cbin][k] = new TH1D(Form("q3ABcnt_p_%d_%d",cbin,k), "", 1, 0, 1);
                q3ABcnt_p[cbin][k]->SetDirectory(0);
                q3ACcnt_p[cbin][k] = new TH1D(Form("q3ACcnt_p_%d_%d",cbin,k), "", 1, 0, 1);
                q3ACcnt_p[cbin][k]->SetDirectory(0);
                q3BCcnt_p[cbin][k] = new TH1D(Form("q3BCcnt_p_%d_%d",cbin,k), "", 1, 0, 1);
                q3BCcnt_p[cbin][k]->SetDirectory(0);
                q3ABcnt_m[cbin][k] = new TH1D(Form("q3ABcnt_m_%d_%d",cbin,k), "", 1, 0, 1);
                q3ABcnt_m[cbin][k]->SetDirectory(0);
                q3ACcnt_m[cbin][k] = new TH1D(Form("q3ACcnt_m_%d_%d",cbin,k), "", 1, 0, 1);
                q3ACcnt_m[cbin][k]->SetDirectory(0);
                q3BCcnt_m[cbin][k] = new TH1D(Form("q3BCcnt_m_%d_%d",cbin,k), "", 1, 0, 1);
                q3BCcnt_m[cbin][k]->SetDirectory(0);
            }
        }
        fclose(flist);
        break;
    }

    GetNumEvents( inlist );

    runParms = new TH1D("runParms", "runParms", 6, 0, 6);
    const char * ParmsLabel[] = {"NumEvnts", "VtxCut","p_{T}^{min}","p_{T}^{max}","#eta^{min}","#eta^{max}"};
    // add more labels as I think of them
    for (int i = 1; i<=6; i++) runParms->GetXaxis()->SetBinLabel(i, ParmsLabel[i-1]);

    cout<<"\n\nAnalysis: "<<anal<<endl;
    cout<<"\n1st-order subevents... "<<endl;
    cout<<"Pos side epA1: "<<hi::EPNames[epA1p]<<"\tepB1: "<<hi::EPNames[epB1p]<<"\tepC1: "<<hi::EPNames[epC1p]<<endl;
    cout<<"Neg side epA1: "<<hi::EPNames[epA1m]<<"\tepB1: "<<hi::EPNames[epB1m]<<"\tepC1: "<<hi::EPNames[epC1m]<<endl;
    cout<<"\n2nd-order subevents... "<<endl;
    cout<<"Pos side epA2: "<<hi::EPNames[epA2p]<<"\tepB2: "<<hi::EPNames[epB2p]<<"\tepC2: "<<hi::EPNames[epC2p]<<endl;
    cout<<"Neg side epA2: "<<hi::EPNames[epA2m]<<"\tepB2: "<<hi::EPNames[epB2m]<<"\tepC2: "<<hi::EPNames[epC2m]<<endl;
    cout<<"\n3rd-order subevents... "<<endl;
    cout<<"Pos side epA3: "<<hi::EPNames[epA3p]<<"\tepB3: "<<hi::EPNames[epB3p]<<"\tepC3: "<<hi::EPNames[epC3p]<<endl;
    cout<<"Neg side epA3: "<<hi::EPNames[epA3m]<<"\tepB3: "<<hi::EPNames[epB3m]<<"\tepC3: "<<hi::EPNames[epC3m]<<endl;
    cout<<"\n"<<endl;

    ReadTree( info, anal, inlist );

}

void ReadTree( GetEventInfo * info, string anal, string inlist ) {

    TStopwatch * sw = new TStopwatch();
    centbins->Reset();
    TFile * tfin;
    TFile * tfout;
    TDirectory * tdcent[ncentbins];
    int filecnt = 0;
    int NumEvnts = 0;
    int NEvt = TotNumEvents*(MaxCent/100.);
    int nbins = ncentbins;
    FILE * flist = fopen(inlist.data(),"r");
    char buf[120];
    string outFile = "results/"+anal;
    string outFext = outFile+".root";
    tfout = new TFile(outFext.data(),"RECREATE");
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        ptav[cbin]->Reset();
        ptcnt[cbin]->Reset();
        badcnt[cbin]->Reset();
    }
    sw->Start();

    while (fgets(buf,120,flist)!=NULL) {
        if (filecnt>=MaxFiles) continue;
        buf[strlen(buf) - 1] = 0;
        TString inFile = buf;
        FILE * ftest = fopen(inFile.Data(),"r");
        if (ftest == NULL) continue;
        fclose(ftest);
        string inf = inFile.Data();
        tfin = new TFile(inFile.Data(),"read");
        if (tfin->IsZombie()) continue;
        if (tfin->TestBit(TFile::kRecovered)) continue;
        tfin->ResetErrno();
        ++filecnt;

        TTree * tree = (TTree *) tfin->Get("vnanalyzer/tree");
        tree->SetBranchAddress("Cent",      &centval);
        tree->SetBranchAddress("ntrkflat",  &ntrks);
        tree->SetBranchAddress("Vtx",       &vtx);
        tree->SetBranchAddress("sumw",      &sumw);
        tree->SetBranchAddress("qx",        &qx);
        tree->SetBranchAddress("qy",        &qy);
        tree->SetBranchAddress("mult",      &epmult);
        tree->SetBranchAddress("qxtrk1",    &qxtrk1_);
        tree->SetBranchAddress("qytrk1",    &qytrk1_);
        tree->SetBranchAddress("qcnt",      &qcnt_);
        tree->SetBranchAddress("avpt",      &avpt_);
        qxtrk1_->SetOption("colz");
        qytrk1_->SetOption("colz");

        for (int ievent = 0; ievent<tree->GetEntries(); ievent++) {
            if (MaxEvents>0 && NumEvnts>=MaxEvents) break;

            tree->GetEntry(ievent);
            if (fabs(vtx)>VtxCut) continue;
            if (centval>MaxCent) continue;
            int cbin = centbins->FindBin(centval)-1;
            if (cbin>=ncentbins) continue;
            if (cbin<0) continue;

            if (MaxEvents<0) {
                if ((int)fmod(NumEvnts, TotNumEvents/20) == 0) {
                    sw->Continue();
                    double elapse = sw->RealTime();
                    cout<<(int)(100*(NumEvnts/(double)NEvt)+0.5)<<" Elapsed: "<<elapse<<"\t Time per event: "<<elapse/(double)NumEvnts<<endl;
                }
            } else {
                if ((int)fmod(NumEvnts, MaxEvents/20) == 0) {
                    sw->Continue();
                    double elapse = sw->RealTime();
                    cout<<(int)(100*(NumEvnts/(double)MaxEvents)+0.5)<<" Elapsed: "<<elapse<<"\t Time per event: "<<elapse/(double)NumEvnts<<endl;
                }
            }

            centbins->Fill(centval);
            ptav[cbin]->Add(avpt_);
            ptcnt[cbin]->Add(qcnt_);
            qxav1[cbin]->Add(qxtrk1_);
            qyav1[cbin]->Add(qytrk1_);
            qxycnt[cbin]->Add(qcnt_);

            v1SP( cbin, qxtrk1_, qytrk1_, qcnt_, qx, qy, sumw, ntrks );

            ++NumEvnts;
        }

        tfin->Close();
    }

    runParms->SetBinContent(1, NumEvnts);
    runParms->SetBinContent(2, VtxCut);
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        ptav[cbin]->GetXaxis()->SetRangeUser(0.301, 11.999);
        ptcnt[cbin]->GetXaxis()->SetRangeUser(0.301, 11.999);
        qxav1[cbin]->GetXaxis()->SetRangeUser(0.301, 11.999);
        qyav1[cbin]->GetXaxis()->SetRangeUser(0.301, 11.999);
        qxycnt[cbin]->GetXaxis()->SetRangeUser(0.301, 11.999);
        for (int k = 0; k<=10; k++) {
            q1_p[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            q1_m[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            q1_pm[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            q112_p[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            q112_m[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            q112_pm[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            q123_p[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            q123_m[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            q123_pm[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            w1_p[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            w1_m[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            w1_pm[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            w112_p[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            w112_m[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            w112_pm[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            w123_p[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            w123_m[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
            w123_pm[cbin][k]->GetXaxis()->SetRangeUser(0.301, 11.999);
        }

        tdcent[cbin] = tfout->mkdir(Form("%d-%d",(int)centBins[cbin],(int)centBins[cbin+1]));
        tdcent[cbin]->cd();

        q1_p[cbin][0]->Write();
        q1_m[cbin][0]->Write();
        q1_pm[cbin][0]->Write();
        q112_p[cbin][0]->Write();
        q112_m[cbin][0]->Write();
        q112_pm[cbin][0]->Write();
        q123_p[cbin][0]->Write();
        q123_m[cbin][0]->Write();
        q123_pm[cbin][0]->Write();
        w1_p[cbin][0]->Write();
        w1_m[cbin][0]->Write();
        w1_pm[cbin][0]->Write();
        w112_p[cbin][0]->Write();
        w112_m[cbin][0]->Write();
        w112_pm[cbin][0]->Write();
        w123_p[cbin][0]->Write();
        w123_m[cbin][0]->Write();
        w123_pm[cbin][0]->Write();

        q1AB_p[cbin][0]->Write();
        q1AC_p[cbin][0]->Write();
        q1BC_p[cbin][0]->Write();
        q1AB_m[cbin][0]->Write();
        q1AC_m[cbin][0]->Write();
        q1BC_m[cbin][0]->Write();
        q1ABcnt_p[cbin][0]->Write();
        q1ACcnt_p[cbin][0]->Write();
        q1BCcnt_p[cbin][0]->Write();
        q1ABcnt_m[cbin][0]->Write();
        q1ACcnt_m[cbin][0]->Write();
        q1BCcnt_m[cbin][0]->Write();

        q2AB_p[cbin][0]->Write();
        q2AC_p[cbin][0]->Write();
        q2BC_p[cbin][0]->Write();
        q2AB_m[cbin][0]->Write();
        q2AC_m[cbin][0]->Write();
        q2BC_m[cbin][0]->Write();
        q2ABcnt_p[cbin][0]->Write();
        q2ACcnt_p[cbin][0]->Write();
        q2BCcnt_p[cbin][0]->Write();
        q2ABcnt_m[cbin][0]->Write();
        q2ACcnt_m[cbin][0]->Write();
        q2BCcnt_m[cbin][0]->Write();

        q3AB_p[cbin][0]->Write();
        q3AC_p[cbin][0]->Write();
        q3BC_p[cbin][0]->Write();
        q3AB_m[cbin][0]->Write();
        q3AC_m[cbin][0]->Write();
        q3BC_m[cbin][0]->Write();
        q3ABcnt_p[cbin][0]->Write();
        q3ACcnt_p[cbin][0]->Write();
        q3BCcnt_p[cbin][0]->Write();
        q3ABcnt_m[cbin][0]->Write();
        q3ACcnt_m[cbin][0]->Write();
        q3BCcnt_m[cbin][0]->Write();

        runParms->Write();
        centbins->Write();
        ptav[cbin]->Write();
        ptcnt[cbin]->Write();
        badcnt[cbin]->Write();
        multTot[cbin][0]->Write();
        qxav1[cbin]->Write();
        qyav1[cbin]->Write();
        qxycnt[cbin]->Write();

        TDirectory * tderr = (TDirectory *) tdcent[cbin]->mkdir("err_calc");
        tderr->cd();
        for (int k = 1; k<=10; k++) {
            q1_p[cbin][k]->Write();
            q1_m[cbin][k]->Write();
            q1_pm[cbin][k]->Write();
            q112_p[cbin][k]->Write();
            q112_m[cbin][k]->Write();
            q112_pm[cbin][k]->Write();
            q123_p[cbin][k]->Write();
            q123_m[cbin][k]->Write();
            q123_pm[cbin][k]->Write();
            w1_p[cbin][k]->Write();
            w1_m[cbin][k]->Write();
            w1_pm[cbin][k]->Write();
            w112_p[cbin][k]->Write();
            w112_m[cbin][k]->Write();
            w112_pm[cbin][k]->Write();
            w123_p[cbin][k]->Write();
            w123_m[cbin][k]->Write();
            w123_pm[cbin][k]->Write();

            q1AB_p[cbin][k]->Write();
            q1AC_p[cbin][k]->Write();
            q1BC_p[cbin][k]->Write();
            q1AB_m[cbin][k]->Write();
            q1AC_m[cbin][k]->Write();
            q1BC_m[cbin][k]->Write();
            q1ABcnt_p[cbin][k]->Write();
            q1ACcnt_p[cbin][k]->Write();
            q1BCcnt_p[cbin][k]->Write();
            q1ABcnt_m[cbin][k]->Write();
            q1ACcnt_m[cbin][k]->Write();
            q1BCcnt_m[cbin][k]->Write();

            q2AB_p[cbin][k]->Write();
            q2AC_p[cbin][k]->Write();
            q2BC_p[cbin][k]->Write();
            q2AB_m[cbin][k]->Write();
            q2AC_m[cbin][k]->Write();
            q2BC_m[cbin][k]->Write();
            q2ABcnt_p[cbin][k]->Write();
            q2ACcnt_p[cbin][k]->Write();
            q2BCcnt_p[cbin][k]->Write();
            q2ABcnt_m[cbin][k]->Write();
            q2ACcnt_m[cbin][k]->Write();
            q2BCcnt_m[cbin][k]->Write();

            q3AB_p[cbin][k]->Write();
            q3AC_p[cbin][k]->Write();
            q3BC_p[cbin][k]->Write();
            q3AB_m[cbin][k]->Write();
            q3AC_m[cbin][k]->Write();
            q3BC_m[cbin][k]->Write();
            q3ABcnt_p[cbin][k]->Write();
            q3ACcnt_p[cbin][k]->Write();
            q3BCcnt_p[cbin][k]->Write();
            q3ABcnt_m[cbin][k]->Write();
            q3ACcnt_m[cbin][k]->Write();
            q3BCcnt_m[cbin][k]->Write();

            multTot[cbin][k]->Write();
        }
    }
    tfout->Close();

    cout << "\n...leaving main event loop" << endl;
    cout << " Total number of events processed: " << NumEvents << endl;
    cout << " Centrality range: " << (int)centBins[0] << "-" << (int)centBins[nbins] << "%" << endl;
    cout << " Vertex cut: " << VtxCut << " cm" << endl;
    cout << " MaxEvents: " << MaxEvents << endl;
    cout << " Number of events accepted: " << NumEvnts << "\n" << endl;

}

void GetNumEvents( string inlist ) {
    cout << "Determining number of events... " << endl;
    TString tr = "";
    TotNumEvents = 0;
    NumEvents = 0;
    TFile * tfin ;
    int filecnt = 0;
    char buf[120];
    FILE * flist = fopen(inlist.data(),"r");

    while (fgets(buf,120,flist)!=NULL) {
        if (filecnt>=MaxFiles) continue;
        buf[strlen(buf) - 1] = 0;
        TString inFile = buf;
        FILE * ftest = fopen(inFile.Data(),"r");
        if (ftest == NULL) continue;
        fclose(ftest);
        tfin = new TFile(inFile.Data(),"read");
        if (tfin->IsZombie()) continue;
        if (tfin->TestBit(TFile::kRecovered)) continue;
        tfin->ResetErrno();
        ++filecnt;

        TTree * tree = (TTree *) tfin->Get("vnanalyzer/tree");
        NumEvents += tree->GetEntries();
        TotNumEvents += tree->GetEntries();
        cout << " Reading file " << inFile.Data() << endl;
        cout << " Filecnt: " << filecnt << "\tNumber of events in file: " << tree->GetEntries()
        << "\t\tRunning total: " << NumEvents << endl;
        tfin->Close();
    }
    fclose(flist);
    cout << " Total number of events found: " << TotNumEvents << endl;
    if (MaxEvents>0) cout << " Number of events to be analyzed: " << MaxEvents << endl;
    else cout << "Analyzing all events" << endl;
    return;
}
