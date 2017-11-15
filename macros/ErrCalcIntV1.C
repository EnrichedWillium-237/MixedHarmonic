
TH1D * v1Stat_pt[nanals][ncentbins];
TH1D * v1Syst_pt[nanals][ncentbins];
TH1D * RelErrv1_pt[nanals][ncentbins];

TH1D * v1Stat_eta[nanals][ncentbins];
TH1D * v1Syst_eta[nanals][ncentbins];
TH1D * RelErrv1_eta[nanals][ncentbins];

TF1 * fitRel_pt[nanals][ncentbins];
TF1 * fitRel_eta[nanals][ncentbins];

void ErrCalcIntV1()
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

            v1Stat_pt[i][cbin] = new TH1D(Form("v1Stat_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);
            v1Syst_pt[i][cbin] = new TH1D(Form("v1Syst_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);
            RelErrv1_pt[i][cbin] = new TH1D(Form("RelErrv1_pt_%s_%d",AnalNames[i].data(),cbin), "", nptbins, ptbins);

            v1p_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v1p_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v1m_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v1m_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v1odd_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v1odd_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v1even_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v1even_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));

            v112p_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v112p_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v112m_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v112m_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v112odd_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v112odd_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v112even_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v112even_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));

            v123p_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v123p_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v123m_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v123m_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v123odd_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v123odd_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));
            v123even_eta[i][cbin] = (TH1D *) tfin->Get(Form("%s/v1_eta/%d-%d/v123even_eta_%s_%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1],AnalNames[i].data(),cbin));

            v1Stat_eta[i][cbin] = new TH1D(Form("v1Stat_eta_%s_%d",AnalNames[i].data(),cbin), "", netabins, etabins);
            v1Syst_eta[i][cbin] = new TH1D(Form("v1Syst_eta_%s_%d",AnalNames[i].data(),cbin), "", netabins, etabins);
            RelErrv1_eta[i][cbin] = new TH1D(Form("RelErrv1_eta_%s_%d",AnalNames[i].data(),cbin), "", netabins, etabins);
        }
    }

    //if (!fopen("data","r")) system("mkdir data");

    //-- symmetrize v1 (v1even for tracker incorrectly called "odd" in output file)
    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            for (int pbin = 0; pbin<nptbins; pbin++) {
                double v1_pos = v1p_pt[i][cbin]->GetBinContent(pbin+1);
                double v1_neg = v1m_pt[i][cbin]->GetBinContent(pbin+1);
                double v1_pos_err = v1p_pt[i][cbin]->GetBinError(pbin+1);
                double v1_neg_err = v1m_pt[i][cbin]->GetBinError(pbin+1);

                double v1stat;
                if (i<=7) v1stat = 0.5*(v1_pos - v1_neg);
                else v1stat = 0.5*(v1_pos + v1_neg);
                double v1stat_err = 0.5*sqrt( v1_pos_err*v1_pos_err + v1_neg_err*v1_neg_err );
                v1Stat_pt[i][cbin]->SetBinContent(pbin+1, v1stat);
                v1Stat_pt[i][cbin]->SetBinError(pbin+1, v1stat_err);

                double v1rel = (fabs(v1_pos) - fabs(v1_neg)) / (fabs(v1_pos) + fabs(v1_neg));
                double v1rel_err = v1rel * sqrt( pow(v1stat_err/(fabs(v1_pos) - fabs(v1_neg)),2) + pow(v1stat_err/(fabs(v1_pos) + fabs(v1_neg)),2) );
                RelErrv1_pt[i][cbin]->SetBinContent(pbin+1, v1rel);
                RelErrv1_pt[i][cbin]->SetBinError(pbin+1, v1rel_err);
            }
            for (int ebin = 0; ebin<netabins; ebin++) {
                double v1_pos = v1odd_eta[i][cbin]->GetBinContent(netabins-ebin);
                double v1_neg = v1odd_eta[i][cbin]->GetBinContent(ebin+1);
                double v1_pos_err = v1odd_eta[i][cbin]->GetBinError(netabins-ebin);
                double v1_neg_err = v1odd_eta[i][cbin]->GetBinError(ebin+1);

                double v1stat;
                if (cbin<=7) v1stat = 0.5*(v1_pos - v1_neg);
                else v1stat = 0.5*(v1_pos + v1_neg);
                double v1stat_err = 0.5*sqrt( v1_pos_err*v1_pos_err + v1_neg_err*v1_neg_err );
                v1Stat_eta[i][cbin]->SetBinContent(netabins-ebin, v1stat);
                v1Stat_eta[i][cbin]->SetBinError(netabins-ebin, v1stat_err);

                double v1rel = (fabs(v1_pos) - fabs(v1_neg)) / (fabs(v1_pos) + fabs(v1_neg));
                double v1rel_err = v1rel * sqrt( pow(v1stat_err/(fabs(v1_pos) - fabs(v1_neg)),2) + pow(v1stat_err/(fabs(v1_pos) + fabs(v1_neg)),2) );
                RelErrv1_eta[i][cbin]->SetBinContent(netabins-ebin, v1rel);
                RelErrv1_eta[i][cbin]->SetBinError(netabins-ebin, v1rel_err);
            }
        }
    }

    //-- calculate percent error
    for (int i = 0; i<nanals; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            fitRel_pt[i][cbin] = new TF1(Form("fitRel_pt_%s_%d-%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1]), "pol0", 0, 12);
            RelErrv1_pt[i][cbin]->Fit(fitRel_pt[i][cbin],"QR");
            double relerr_pt = fitRel_pt[i][cbin]->GetParameter(0);
            for (int pbin = 0; pbin<nptbins; pbin++) {
                double v1val = v1Stat_pt[i][cbin]->GetBinContent(pbin+1);
                v1Syst_pt[i][cbin]->SetBinContent(pbin+1, v1val);
                v1Syst_pt[i][cbin]->SetBinError(pbin+1, relerr_pt*v1val);
            }
            fitRel_eta[i][cbin] = new TF1(Form("fitRel_eta_%s_%d-%d",AnalNames[i].data(),centBins[cbin],centBins[cbin+1]), "pol0", 0, 2.4);
            RelErrv1_eta[i][cbin]->Fit(fitRel_eta[i][cbin],"QR");
            double relerr_eta = fitRel_eta[i][cbin]->GetParameter(0);
            for (int ebin = 0; ebin<netabins; ebin++) {
                double v1val = v1Stat_eta[i][cbin]->GetBinContent(ebin+1);
                v1Syst_eta[i][cbin]->SetBinContent(ebin+1, v1val);
                v1Syst_eta[i][cbin]->SetBinError(ebin+1, relerr_eta*v1val);
            }
        }
    }


    TFile * tfout;
    if (subEvt2) tfout = new TFile("../outputs/final_outputs/v1Int_2sub_ErrCalc.root","recreate");
    else tfout = new TFile("../outputs/final_outputs/v1Int_ErrCalc.root","recreate");
    for (int i = 0; i<nanals; i++) {
        TDirectory * tdAnal = (TDirectory *) tfout->mkdir(Form("%s",AnalNames[i].data()));
        TDirectory * tdPt = (TDirectory *) tdAnal->mkdir("v1_pt");
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TDirectory * tdPtCent = (TDirectory *) tdPt->mkdir(Form("%d-%d",centBins[cbin],centBins[cbin+1]));
            tdPtCent->cd();
            v1Stat_pt[i][cbin]->Write();
            v1Syst_pt[i][cbin]->Write();
            RelErrv1_pt[i][cbin]->Write();
            fitRel_pt[i][cbin]->Write();
        }
        TDirectory * tdEta = (TDirectory *) tdAnal->mkdir("v1_eta");
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TDirectory * tdPtEta = (TDirectory *) tdEta->mkdir(Form("%d-%d",centBins[cbin],centBins[cbin+1]));
            tdPtEta->cd();
            v1Stat_eta[i][cbin]->Write();
            v1Syst_eta[i][cbin]->Write();
            RelErrv1_eta[i][cbin]->Write();
            fitRel_eta[i][cbin]->Write();
        }
    }

    tfout->Close();

}
