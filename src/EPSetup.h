// Positive and negative side event plane selection
//
// Note: When using the HF event plane, suffix refers to subevent C
// When using the tracker event plane, suffix refers to subevent A

void EPSetup( string anal ) {

    TString AnalysisType = anal;

    epA1p = hi::HFp1;
    epB1p = hi::HFm1;
    epC1p = hi::trackp1;
    epA1m = hi::HFm1;
    epB1m = hi::HFp1;
    epC1m = hi::trackm1;

    epA2p = hi::HFp2;
    epB2p = hi::HFm2;
    epC2p = hi::trackp222;
    epA2m = hi::HFm2;
    epB2m = hi::HFp2;
    epC2m = hi::trackm222;

    epA3p = hi::HFp3;
    epB3p = hi::HFm3;
    epC3p = hi::trackp322;
    epA3m = hi::HFm3;
    epB3m = hi::HFp3;
    epC3m = hi::trackm322;

    if (AnalysisType.Contains("mid")) {
        epC1p = hi::trackmid1;
        epC1m = hi::trackmid1;

        epC2p = hi::trackmid2;
        epC2m = hi::trackmid2;

        epC3p = hi::trackmid3;
        epC3m = hi::trackmid3;
    }

    if (AnalysisType.Contains("_102")) {
        epC1p = hi::trackp102;
        epC1m = hi::trackm102;

        epC2p = hi::trackp202;
        epC2m = hi::trackm202;

        epC3p = hi::trackp302;
        epC3m = hi::trackm302;
    }

    if (AnalysisType.Contains("_106")) {
        epC1p = hi::trackp106;
        epC1m = hi::trackm106;

        epC2p = hi::trackp206;
        epC2m = hi::trackm206;

        epC3p = hi::trackp306;
        epC3m = hi::trackm306;
    }

    if (AnalysisType.Contains("_110")) {
        epC1p = hi::trackp110;
        epC1m = hi::trackm110;

        epC2p = hi::trackp210;
        epC2m = hi::trackm210;

        epC3p = hi::trackp310;
        epC3m = hi::trackm310;
    }

    if (AnalysisType.Contains("_114")) {
        epC1p = hi::trackp114;
        epC1m = hi::trackm114;

        epC2p = hi::trackp214;
        epC2m = hi::trackm214;

        epC3p = hi::trackp314;
        epC3m = hi::trackm314;
    }

    if (AnalysisType.Contains("_118")) {
        epC1p = hi::trackp118;
        epC1m = hi::trackm118;

        epC2p = hi::trackp218;
        epC2m = hi::trackm218;

        epC3p = hi::trackp318;
        epC3m = hi::trackm318;
    }

    if (AnalysisType.Contains("_122")) {
        epC1p = hi::trackp122;
        epC1m = hi::trackm122;

        epC2p = hi::trackp222;
        epC2m = hi::trackm222;

        epC3p = hi::trackp322;
        epC3m = hi::trackm322;
    }

    if (AnalysisType.Contains("mc")) {
        epA1p = hi::trackp1mc;
        epB1p = hi::trackm122mc;
        epC1p = hi::trackp122mc;

        epA1m = hi::trackm1mc;
        epB1m = hi::trackp122mc;
        epC1m = hi::trackm122mc;
    }

    if (AnalysisType.Contains("mc") && AnalysisType.Contains("mid")) {
        epA1p = hi::trackmid1mc;
        epB1p = hi::trackp122mc;
        epC1p = hi::trackm122mc;

        epA1m = hi::trackmid1mc;
        epB1m = hi::trackm122mc;
        epC1m = hi::trackp122mc;
    }

    if (AnalysisType.Contains("mc") && AnalysisType.Contains("_102")) {
        epA1p = hi::trackp102mc;
        epB1p = hi::trackm110mc;
        epC1p = hi::trackm122mc;

        epA1m = hi::trackm102mc;
        epB1m = hi::trackp110mc;
        epC1m = hi::trackp122mc;
    }

    if (AnalysisType.Contains("mc") && AnalysisType.Contains("_106")) {
        epA1p = hi::trackp106mc;
        epB1p = hi::trackm110mc;
        epC1p = hi::trackm122mc;

        epA1m = hi::trackm106mc;
        epB1m = hi::trackp110mc;
        epC1m = hi::trackm122mc;
    }

    if (AnalysisType.Contains("mc") && AnalysisType.Contains("_110")) {
        epA1p = hi::trackp110mc;
        epB1p = hi::trackm1mc;
        epC1p = hi::trackm122mc;

        epA1m = hi::trackm110mc;
        epB1m = hi::trackp1mc;
        epC1m = hi::trackp122mc;
    }

    if (AnalysisType.Contains("mc") && AnalysisType.Contains("_114")) {
        epA1p = hi::trackp114mc;
        epB1p = hi::trackm1mc;
        epC1p = hi::trackm122mc;

        epA1m = hi::trackm114mc;
        epB1m = hi::trackp1mc;
        epC1m = hi::trackp122mc;
    }

    if (AnalysisType.Contains("mc") && AnalysisType.Contains("_118")) {
        epA1p = hi::trackp118mc;
        epB1p = hi::trackm1mc;
        epC1p = hi::trackm122mc;

        epA1m = hi::trackm118mc;
        epB1m = hi::trackp1mc;
        epC1m = hi::trackp122mc;
    }

    if (AnalysisType.Contains("mc") && AnalysisType.Contains("_122")) {
        epA1p = hi::trackp122mc;
        epB1p = hi::trackm1mc;
        epC1p = hi::trackm122mc;

        epA1m = hi::trackm122mc;
        epB1m = hi::trackp1mc;
        epC1m = hi::trackp122mc;
    }

}
