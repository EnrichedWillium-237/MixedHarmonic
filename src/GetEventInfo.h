#ifndef GETEVENTINFO
#define GETEVENTINFO

class GetEventInfo {
public:
    GetEventInfo(TString inFile);
    ~GetEventInfo() {tfin->Delete();}
    TTree * getTree() {return tree;}
    int status;
    int getNumEtaBins() {return netabins;}
    int getNumPtBins() {return nptbins;}
    TH1D * getEtaHist() {return heta;}
    TH1D * getPtHist() {return hpt;}
    TH2D * getTemplate() {return hTemplate;}
    TString getFileName() {return filename;}
private:
    TTree * tree;
    int netabins;
    int nptbins;
    TH1D * heta;
    TH1D * hpt;
    TH2D * hTemplate;
    TFile * tfin;
    TString filename;
};

GetEventInfo::GetEventInfo( TString inFile ) {
    status = 0;
    filename = inFile;
    hTemplate = 0;
    TFile * tfin = new TFile(inFile.Data(),"read");
    tfin    = new TFile(inFile.Data(),"read");
    if(tfin->IsZombie()) return;
    if(tfin->TestBit(TFile::kRecovered)) return;
    tfin->ResetErrno();

    const TArrayD * etabins;
    const TArrayD * ptbins;
    tree = (TTree *) tfin->Get("vnanalyzer/tree");
    TH2D * h1 = (TH2D *) tfin->Get("vnanalyzer/qxtrk1");
    if (h1 == 0) return;
    netabins = h1->GetNbinsY();
    nptbins = h1->GetNbinsX();
    etabins = h1->GetYaxis()->GetXbins();
    ptbins = h1->GetXaxis()->GetXbins();
    if (!hTemplate) hTemplate = (TH2D *) h1->Clone("hTemplate");
    hTemplate->SetXTitle("p_{T} (GeV/c)");
    hTemplate->SetYTitle("#eta");
    hTemplate->SetOption("colz");
    if (hTemplate) hTemplate->Reset();
    double pttmp[100];
    double etatmp[100];

    pttmp[0] = 0.3; // hard-coded to min pT (otherwise defaults to 0 GeV/c)
    for (int i = 1; i<=nptbins; i++) {
        pttmp[i] = ptbins->At(i);
    }
    for (int i = 0; i<=netabins; i++) {
        etatmp[i] = etabins->At(i);
    }
    hpt = new TH1D(Form("hpt_%s",inFile.Data()), "hpt", nptbins, pttmp);
    heta = new TH1D(Form("heta_%s",inFile.Data()), "heta", netabins, etatmp);

    // cout << "Event info retrieved... " << endl;
    // cout << "  number of pt bins: " << nptbins << endl;
    // for (int i = 0; i<nptbins; i++) cout << "  " << pttmp[i];
    // cout << "  " << pttmp[nptbins] << endl;
    // cout << "  number of eta bins: " << netabins << endl;
    // for (int i = 0; i<netabins; i++) cout << "  " << etatmp[i];
    // cout << "  " << etatmp[netabins] << endl;

    status = 1;
    return;
}

#endif
