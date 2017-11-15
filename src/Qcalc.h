void v1SP( int cbin, TH2D * qxtrk1_, TH2D * qytrk1_, TH2D * qcnt_, double * qx, double * qy, double * sumw, int ntrks ) {

    int Ntrk = ntrks;
    multTot[cbin][0]->Fill(Ntrk);

    double A1x_p = qx[epA1p];
    double A1y_p = qy[epA1p];
    double wA1_p = sumw[epA1p];
    double A1x_m = qx[epA1m];
    double A1y_m = qy[epA1m];
    double wA1_m = sumw[epA1m];

    double A2x_p = qx[epA2p];
    double A2y_p = qy[epA2p];
    double wA2_p = sumw[epA2p];
    double A2x_m = qx[epA2m];
    double A2y_m = qy[epA2m];
    double wA2_m = sumw[epA2m];

    double A3x_p = qx[epA3p];
    double A3y_p = qy[epA3p];
    double wA3_p = sumw[epA3p];
    double A3x_m = qx[epA3m];
    double A3y_m = qy[epA3m];
    double wA3_m = sumw[epA3m];

    double B1x_p = qx[epB1p];
    double B1y_p = qy[epB1p];
    double wB1_p = sumw[epB1p];
    double B1x_m = qx[epB1m];
    double B1y_m = qy[epB1m];
    double wB1_m = sumw[epB1m];

    double B2x_p = qx[epB2p];
    double B2y_p = qy[epB2p];
    double wB2_p = sumw[epB2p];
    double B2x_m = qx[epB2m];
    double B2y_m = qy[epB2m];
    double wB2_m = sumw[epB2m];

    double B3x_p = qx[epB3p];
    double B3y_p = qy[epB3p];
    double wB3_p = sumw[epB3p];
    double B3x_m = qx[epB3m];
    double B3y_m = qy[epB3m];
    double wB3_m = sumw[epB3m];

    double C1x_p = qx[epC1p];
    double C1y_p = qy[epC1p];
    double wC1_p = sumw[epC1p];
    double C1x_m = qx[epC1m];
    double C1y_m = qy[epC1m];
    double wC1_m = sumw[epC1m];

    double C2x_p = qx[epC2p];
    double C2y_p = qy[epC2p];
    double wC2_p = sumw[epC2p];
    double C2x_m = qx[epC2m];
    double C2y_m = qy[epC2m];
    double wC2_m = sumw[epC2m];

    double C3x_p = qx[epC3p];
    double C3y_p = qy[epC3p];
    double wC3_p = sumw[epC3p];
    double C3x_m = qx[epC3m];
    double C3y_m = qy[epC3m];
    double wC3_m = sumw[epC3m];

    q1_p[cbin][0]->Add(qxtrk1_, A1x_p);
    q1_p[cbin][0]->Add(qytrk1_, A1y_p);
    q1_m[cbin][0]->Add(qxtrk1_, A1x_m);
    q1_m[cbin][0]->Add(qytrk1_, A1y_m);
    q1_pm[cbin][0]->Add(q1_p[cbin][0]);
    q1_pm[cbin][0]->Add(q1_m[cbin][0]);
    q112_p[cbin][0]->Add(qxtrk1_, A1x_p*C2x_p + A1y_p*C2y_p);
    q112_p[cbin][0]->Add(qytrk1_, A1x_p*C2y_p - A1y_p*C2x_p);
    q112_m[cbin][0]->Add(qxtrk1_, A1x_m*C2x_m + A1y_m*C2y_m);
    q112_m[cbin][0]->Add(qytrk1_, A1x_m*C2y_m - A1y_m*C2x_m);
    q112_pm[cbin][0]->Add(q112_p[cbin][0]);
    q112_pm[cbin][0]->Add(q112_m[cbin][0]);
    q123_p[cbin][0]->Add(qxtrk1_, A2x_p*C3x_p + A2y_p*C3y_p);
    q123_p[cbin][0]->Add(qytrk1_, A2x_p*C3y_p - A2y_p*C3x_p);
    q123_m[cbin][0]->Add(qxtrk1_, A2x_m*C3x_m + A2y_m*C3y_m);
    q123_m[cbin][0]->Add(qytrk1_, A2x_m*C3y_m - A2y_m*C3x_m);
    q123_pm[cbin][0]->Add(q123_p[cbin][0]);
    q123_pm[cbin][0]->Add(q123_m[cbin][0]);
    w1_p[cbin][0]->Add(qcnt_, wA1_p);
    w1_m[cbin][0]->Add(qcnt_, wA1_m);
    w1_pm[cbin][0]->Add(w1_p[cbin][0]);
    w1_pm[cbin][0]->Add(w1_m[cbin][0]);
    w112_p[cbin][0]->Add(qcnt_, wA1_p*wC2_p);
    w112_m[cbin][0]->Add(qcnt_, wA1_m*wC2_m);
    w112_pm[cbin][0]->Add(w112_p[cbin][0]);
    w112_pm[cbin][0]->Add(w112_m[cbin][0]);
    w123_p[cbin][0]->Add(qcnt_, wA2_p*wC3_p);
    w123_m[cbin][0]->Add(qcnt_, wA2_m*wC3_m);
    w123_pm[cbin][0]->Add(w123_p[cbin][0]);
    w123_pm[cbin][0]->Add(w123_m[cbin][0]);

    q1AB_p[cbin][0]->Fill(0., A1x_p*B1x_p + A1y_p*B1y_p);
    q1AC_p[cbin][0]->Fill(0., A1x_p*C1x_p + A1y_p*C1y_p);
    q1BC_p[cbin][0]->Fill(0., B1x_p*C1x_p + B1y_p*C1y_p);
    q1AB_m[cbin][0]->Fill(0., A1x_m*B1x_m + A1y_m*B1y_m);
    q1AC_m[cbin][0]->Fill(0., A1x_m*C1x_m + A1y_m*C1y_m);
    q1BC_m[cbin][0]->Fill(0., B1x_m*C1x_m + B1y_m*C1y_m);
    q1ABcnt_p[cbin][0]->Fill(0., wA1_p*wB1_p);
    q1ACcnt_p[cbin][0]->Fill(0., wA1_p*wC1_p);
    q1BCcnt_p[cbin][0]->Fill(0., wB1_p*wC1_p);
    q1ABcnt_m[cbin][0]->Fill(0., wA1_m*wB1_m);
    q1ACcnt_m[cbin][0]->Fill(0., wA1_m*wC1_m);
    q1BCcnt_m[cbin][0]->Fill(0., wB1_m*wC1_m);

    q2AB_p[cbin][0]->Fill(0., A2x_p*B2x_p + A2y_p*B2y_p);
    q2AC_p[cbin][0]->Fill(0., A2x_p*C2x_p + A2y_p*C2y_p);
    q2BC_p[cbin][0]->Fill(0., B2x_p*C2x_p + B2y_p*C2y_p);
    q2AB_m[cbin][0]->Fill(0., A2x_m*B2x_m + A2y_m*B2y_m);
    q2AC_m[cbin][0]->Fill(0., A2x_m*C2x_m + A2y_m*C2y_m);
    q2BC_m[cbin][0]->Fill(0., B2x_m*C2x_m + B2y_m*C2y_m);
    q2ABcnt_p[cbin][0]->Fill(0., wA2_p*wB2_p);
    q2ACcnt_p[cbin][0]->Fill(0., wA2_p*wC2_p);
    q2BCcnt_p[cbin][0]->Fill(0., wB2_p*wC2_p);
    q2ABcnt_m[cbin][0]->Fill(0., wA2_m*wB2_m);
    q2ACcnt_m[cbin][0]->Fill(0., wA2_m*wC2_m);
    q2BCcnt_m[cbin][0]->Fill(0., wB2_m*wC2_m);

    q3AB_p[cbin][0]->Fill(0., A3x_p*B3x_p + A3y_p*B3y_p);
    q3AC_p[cbin][0]->Fill(0., A3x_p*C3x_p + A3y_p*C3y_p);
    q3BC_p[cbin][0]->Fill(0., B3x_p*C3x_p + B3y_p*C3y_p);
    q3AB_m[cbin][0]->Fill(0., A3x_m*B3x_m + A3y_m*B3y_m);
    q3AC_m[cbin][0]->Fill(0., A3x_m*C3x_m + A3y_m*C3y_m);
    q3BC_m[cbin][0]->Fill(0., B3x_m*C3x_m + B3y_m*C3y_m);
    q3ABcnt_p[cbin][0]->Fill(0., wA3_p*wB3_p);
    q3ACcnt_p[cbin][0]->Fill(0., wA3_p*wC3_p);
    q3BCcnt_p[cbin][0]->Fill(0., wB3_p*wC3_p);
    q3ABcnt_m[cbin][0]->Fill(0., wA3_m*wB3_m);
    q3ACcnt_m[cbin][0]->Fill(0., wA3_m*wC3_m);
    q3BCcnt_m[cbin][0]->Fill(0., wB3_m*wC3_m);

    int k = (int)(ran->Uniform(0,9.999)) + 1;

    multTot[cbin][k]->Fill(Ntrk);

    q1_p[cbin][k]->Add(qxtrk1_, A1x_p);
    q1_p[cbin][k]->Add(qytrk1_, A1y_p);
    q1_m[cbin][k]->Add(qxtrk1_, A1x_m);
    q1_m[cbin][k]->Add(qytrk1_, A1y_m);
    q1_pm[cbin][k]->Add(q1_p[cbin][k]);
    q1_pm[cbin][k]->Add(q1_m[cbin][k]);
    q112_p[cbin][k]->Add(qxtrk1_, A1x_p*C2x_p + A1y_p*C2y_p);
    q112_p[cbin][k]->Add(qytrk1_, A1x_p*C2y_p - A1y_p*C2x_p);
    q112_m[cbin][k]->Add(qxtrk1_, A1x_m*C2x_m + A1y_m*C2y_m);
    q112_m[cbin][k]->Add(qytrk1_, A1x_m*C2y_m - A1y_m*C2x_m);
    q112_pm[cbin][k]->Add(q112_p[cbin][k]);
    q112_pm[cbin][k]->Add(q112_m[cbin][k]);
    q123_p[cbin][k]->Add(qxtrk1_, A2x_p*C3x_p + A2y_p*C3y_p);
    q123_p[cbin][k]->Add(qytrk1_, A2x_p*C3y_p - A2y_p*C3x_p);
    q123_m[cbin][k]->Add(qxtrk1_, A2x_m*C3x_m + A2y_m*C3y_m);
    q123_m[cbin][k]->Add(qytrk1_, A2x_m*C3y_m - A2y_m*C3x_m);
    q123_pm[cbin][k]->Add(q123_p[cbin][k]);
    q123_pm[cbin][k]->Add(q123_m[cbin][k]);
    w1_p[cbin][k]->Add(qcnt_,wA1_p);
    w1_m[cbin][k]->Add(qcnt_,wA1_m);
    w1_pm[cbin][k]->Add(w1_p[cbin][k]);
    w1_pm[cbin][k]->Add(w1_m[cbin][k]);
    w112_p[cbin][k]->Add(qcnt_, wA1_p*wC2_p);
    w112_m[cbin][k]->Add(qcnt_, wA1_m*wC2_m);
    w112_pm[cbin][k]->Add(w112_p[cbin][k]);
    w112_pm[cbin][k]->Add(w112_m[cbin][k]);
    w123_p[cbin][k]->Add(qcnt_, wA2_p*wC3_p);
    w123_m[cbin][k]->Add(qcnt_, wA2_m*wC3_m);
    w123_pm[cbin][k]->Add(w123_p[cbin][k]);
    w123_pm[cbin][k]->Add(w123_m[cbin][k]);

    q1AB_p[cbin][k]->Fill(0., A1x_p*B1x_p + A1y_p*B1y_p);
    q1AC_p[cbin][k]->Fill(0., A1x_p*C1x_p + A1y_p*C1y_p);
    q1BC_p[cbin][k]->Fill(0., B1x_p*C1x_p + B1y_p*C1y_p);
    q1AB_m[cbin][k]->Fill(0., A1x_m*B1x_m + A1y_m*B1y_m);
    q1AC_m[cbin][k]->Fill(0., A1x_m*C1x_m + A1y_m*C1y_m);
    q1BC_m[cbin][k]->Fill(0., B1x_m*C1x_m + B1y_m*C1y_m);
    q1ABcnt_p[cbin][k]->Fill(0., wA1_p*wB1_p);
    q1ACcnt_p[cbin][k]->Fill(0., wA1_p*wC1_p);
    q1BCcnt_p[cbin][k]->Fill(0., wB1_p*wC1_p);
    q1ABcnt_m[cbin][k]->Fill(0., wA1_m*wB1_m);
    q1ACcnt_m[cbin][k]->Fill(0., wA1_m*wC1_m);
    q1BCcnt_m[cbin][k]->Fill(0., wB1_m*wC1_m);

    q2AB_p[cbin][k]->Fill(0., A2x_p*B2x_p + A2y_p*B2y_p);
    q2AC_p[cbin][k]->Fill(0., A2x_p*C2x_p + A2y_p*C2y_p);
    q2BC_p[cbin][k]->Fill(0., B2x_p*C2x_p + B2y_p*C2y_p);
    q2AB_m[cbin][k]->Fill(0., A2x_m*B2x_m + A2y_m*B2y_m);
    q2AC_m[cbin][k]->Fill(0., A2x_m*C2x_m + A2y_m*C2y_m);
    q2BC_m[cbin][k]->Fill(0., B2x_m*C2x_m + B2y_m*C2y_m);
    q2ABcnt_p[cbin][k]->Fill(0., wA2_p*wB2_p);
    q2ACcnt_p[cbin][k]->Fill(0., wA2_p*wC2_p);
    q2BCcnt_p[cbin][k]->Fill(0., wB2_p*wC2_p);
    q2ABcnt_m[cbin][k]->Fill(0., wA2_m*wB2_m);
    q2ACcnt_m[cbin][k]->Fill(0., wA2_m*wC2_m);
    q2BCcnt_m[cbin][k]->Fill(0., wB2_m*wC2_m);

    q3AB_p[cbin][k]->Fill(0., A3x_p*B3x_p + A3y_p*B3y_p);
    q3AC_p[cbin][k]->Fill(0., A3x_p*C3x_p + A3y_p*C3y_p);
    q3BC_p[cbin][k]->Fill(0., B3x_p*C3x_p + B3y_p*C3y_p);
    q3AB_m[cbin][k]->Fill(0., A3x_m*B3x_m + A3y_m*B3y_m);
    q3AC_m[cbin][k]->Fill(0., A3x_m*C3x_m + A3y_m*C3y_m);
    q3BC_m[cbin][k]->Fill(0., B3x_m*C3x_m + B3y_m*C3y_m);
    q3ABcnt_p[cbin][k]->Fill(0., wA3_p*wB3_p);
    q3ACcnt_p[cbin][k]->Fill(0., wA3_p*wC3_p);
    q3BCcnt_p[cbin][k]->Fill(0., wB3_p*wC3_p);
    q3ABcnt_m[cbin][k]->Fill(0., wA3_m*wB3_m);
    q3ACcnt_m[cbin][k]->Fill(0., wA3_m*wC3_m);
    q3BCcnt_m[cbin][k]->Fill(0., wB3_m*wC3_m);

}
