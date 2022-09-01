double BarlowLLH(double mc, double w2, double data)
{
        // Solving for the quadratic equation,
        // beta^2 + (mu * sigma^2 - 1)beta - data * sigma^2) = 0
        // where sigma^2 is the relative variance.
        double rel_var = w2 / (mc * mc);
        double b       = (mc * rel_var) - 1;
        double c       = 4 * data * rel_var;

        double beta   = (-b + std::sqrt(b * b + c)) / 2.0;
        double mc_hat = mc * beta;

        // Calculate the following LLH:
        //-2lnL = 2 * beta*mc - data + data * ln(data / (beta*mc)) + (beta-1)^2 / sigma^2
        // where sigma^2 is the same as above.
        double chi2 = 0.0;
        //if(data <= 0.0)
        //{
        //    chi2 = 2 * mc_hat;
        //    chi2 += (beta - 1) * (beta - 1) / rel_var;
        //}
        if(mc_hat > 0.0)
        {
            chi2 = 2 * (mc_hat - data);
            if(data > 0.0)
                chi2 += 2 * data * std::log(data / mc_hat);

            chi2 += (beta - 1) * (beta - 1) / rel_var;
        }

        return (chi2 >= 0.0) ? chi2 : 0.0;
}

//#############################################################################
void SetStyleVariables(TStyle *t2kStyle){

  // use plain black on white colors
  t2kStyle->SetFrameBorderMode(0);
  t2kStyle->SetCanvasBorderMode(0);
  t2kStyle->SetPadBorderMode(0);
  t2kStyle->SetPadColor(0);
  t2kStyle->SetCanvasColor(0);
  t2kStyle->SetStatColor(0);
  //t2kStyle->SetFillColor(0);
  t2kStyle->SetLegendBorderSize(1);

  // set the paper & margin sizes
  t2kStyle->SetPaperSize(20,26);
  t2kStyle->SetPadTopMargin(0.1);
  t2kStyle->SetPadRightMargin(0.05); //0.15 
  t2kStyle->SetPadBottomMargin(0.15);
  t2kStyle->SetPadLeftMargin(0.15);

  t2kStyle->SetCanvasDefH(550);

  // use large Times-Roman fonts
  t2kStyle->SetTextFont(132);
  t2kStyle->SetTextSize(0.08);
  t2kStyle->SetLabelFont(132,"x");
  t2kStyle->SetLabelFont(132,"y");
  t2kStyle->SetLabelFont(132,"z");
  t2kStyle->SetLabelSize(0.06,"x");
  t2kStyle->SetTitleSize(0.072,"x");
  t2kStyle->SetTitleOffset(0.88,"x");
  t2kStyle->SetLabelSize(0.05,"y");
  t2kStyle->SetTitleSize(0.08,"y");
  t2kStyle->SetTitleOffset(0.7,"y");
  t2kStyle->SetLabelSize(0.06,"z");
  t2kStyle->SetTitleSize(0.06,"z");
  t2kStyle->SetLabelFont(132,"t");
  t2kStyle->SetTitleFont(132,"x");
  t2kStyle->SetTitleFont(132,"y");
  t2kStyle->SetTitleFont(132,"z");
  t2kStyle->SetTitleFont(132,"t");
  t2kStyle->SetTitleFillColor(0);
  t2kStyle->SetTitleX(0.25);
  t2kStyle->SetTitleFontSize(0.08);
  t2kStyle->SetTitleFont(132,"pad");

  //t2kStyle->SetPadGridX(true);
  //t2kStyle->SetPadGridY(true);

  // use bold lines and markers
  //  t2kStyle->SetMarkerStyle(20);
  t2kStyle->SetHistLineWidth(1.85);
  t2kStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  //  t2kStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  t2kStyle->SetOptTitle(0);
  t2kStyle->SetOptStat(0);
  t2kStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  //t2kStyle->SetPadTickX(1);
  t2kStyle->SetPadTickY(1);

  t2kStyle->SetPalette(1,0);  // use the nice red->blue palette
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  //Double_t green[NRGBs] = { 0.00, 0.00, 0.00, 0.00, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                   NCont);
  t2kStyle->SetNumberContours(NCont);

  // End of definition of t2kStyle
}

void makePlots_kinematics(int TKI=0){
    TStyle *t2kstyle = new TStyle("T2K","T2K approved plots style");
    SetStyleVariables(t2kstyle);
    gROOT->SetStyle("T2K");

    const int nDpttbins = 5;
    const double Dpttbins[nDpttbins+1] = {-700,-300,-100,100,300,700};
    const int nPnbins = 4;
    const double Pnbins[nPnbins+1] = {0,120,240,600,1500};
    const int nDatbins = 3;
    const double Datbins[nDatbins+1] = {0,60,120,180};
    double Bins[100];
    int nBins = -1;
    std::string mcfile;
    std::string datafile;
    double tkimax = -1;
    std::string ylabel;
    std::string xlabel;
    std::string varName;
    switch (TKI) {
        case 0:
            mcfile = "../inputs/neut_prod6T_dptt_geomfixed.root";
            datafile = "../inputs/data_dptt_prod6T.root";
            nBins = nDpttbins;
            for (int i=0;i<nDpttbins+1;i++) Bins[i]=Dpttbins[i];
            tkimax = 700;
            xlabel = "Reconstructed #it{#deltap_{TT}} (MeV/c)";
            varName = "dptt";
            break;
        default:
            printf("***Warning: Not a valid TKI variable***\n");
            break;
    }

    TFile* fMC = new TFile(mcfile.c_str());
    TTree* selectedEvents = (TTree*)fMC->Get("selectedEvents");
    TFile* fData = new TFile(datafile.c_str());
    TTree* selectedEventsData = (TTree*)fData->Get("selectedEvents");

    const int nPmombins = 6;
    const double Pmombins[nPmombins+1] = {405,575,700,825,950,1075,1320};
    const int nTopo = 12;
    const int nSamples = 5;
    TH1D* hMC_Prefit[nSamples][nTopo];
    const int nKinematics = 7;
    TH1D* hMC_Prefit_kinematics[nKinematics][nSamples][nTopo]; // muon, pion, proton kinematic plots + proton multiplicity
    TH2D* hMC_2D[nSamples][3];
    TH1D* hMC_Prefit_pmom[nTopo];
    for (int i=0;i<nTopo;i++) 
        for (int j=0;j<nSamples;j++) {
            hMC_Prefit[j][i] = new TH1D("","",nBins,Bins);
            hMC_Prefit_kinematics[0][j][i] = new TH1D("","",10,225,7700);
            hMC_Prefit_kinematics[1][j][i] = new TH1D("","",10,0.342,1);
            hMC_Prefit_kinematics[2][j][i] = new TH1D("","",10,405,1320);
            hMC_Prefit_kinematics[3][j][i] = new TH1D("","",10,0.342,1);
            hMC_Prefit_kinematics[4][j][i] = new TH1D("","",10,135,1320);
            hMC_Prefit_kinematics[5][j][i] = new TH1D("","",10,0.342,1);
            hMC_Prefit_kinematics[6][j][i] = new TH1D("","",3,1,4);
    }
    TH1D* hData[nSamples];
    TH1D* hData_kinematics[nKinematics][nSamples];
    for (int j=0;j<nSamples;j++){
        hData_kinematics[0][j] = new TH1D("","",10,225,7700);
        hData_kinematics[1][j] = new TH1D("","",10,0.342,1);
        hData_kinematics[2][j] = new TH1D("","",10,405,1320);
        hData_kinematics[3][j] = new TH1D("","",10,0.342,1);
        hData_kinematics[4][j] = new TH1D("","",10,135,1320);
        hData_kinematics[5][j] = new TH1D("","",10,0.342,1);
        hData_kinematics[6][j] = new TH1D("","",3,1,4);
    }
    TH1D* hData_pmom = new TH1D("","",nPmombins,Pmombins);

    Int_t sample;
    Float_t D1true;
    Float_t D2true;
    Int_t nutype;
    Int_t topology;
    Int_t reaction;
    Int_t target;
    Float_t D1Reco;
    Float_t D2Reco;
    Float_t weight;
    Float_t weightNom;
    Float_t weightMC;
    Float_t muMomRec, muCosThetaRec, pMomRec, pCosThetaRec, piMomRec, piCosThetaRec;
    Int_t NTPCproton;

    selectedEvents->SetBranchAddress("reaction", &reaction);
    selectedEvents->SetBranchAddress("target", &target);
    selectedEvents->SetBranchAddress("cutBranch", &sample);
    selectedEvents->SetBranchAddress("topology", &topology);
    selectedEvents->SetBranchAddress("D1True", &D1true);
    selectedEvents->SetBranchAddress("D1Rec", &D1Reco);
    selectedEvents->SetBranchAddress("D2True", &D2true);
    selectedEvents->SetBranchAddress("D2Rec", &D2Reco);
    selectedEvents->SetBranchAddress("weight", &weight);
    selectedEvents->SetBranchAddress("muMomRec", &muMomRec);
    selectedEvents->SetBranchAddress("muCosThetaRec", &muCosThetaRec);
    selectedEvents->SetBranchAddress("pMomRec", &pMomRec);
    selectedEvents->SetBranchAddress("pCosThetaRec", &pCosThetaRec);
    selectedEvents->SetBranchAddress("piMomRec", &piMomRec);
    selectedEvents->SetBranchAddress("piCosThetaRec", &piCosThetaRec);
    selectedEvents->SetBranchAddress("NTPCproton", &NTPCproton);

    Float_t nMCEvents[nSamples][nTopo]={0};
    Float_t nMCEventsSum[nSamples]={0};
    for (int i=0;i<selectedEvents->GetEntries();i++) { 
        selectedEvents->GetEntry(i);
        if(topology<0&&topology>11) continue;
        if (fabs(D1Reco)<tkimax) { // exclude the OOPS events
            hMC_Prefit_kinematics[0][sample][topology]->Fill(muMomRec,weight);
            hMC_Prefit_kinematics[1][sample][topology]->Fill(muCosThetaRec,weight);
            hMC_Prefit_kinematics[2][sample][topology]->Fill(pMomRec,weight);
            hMC_Prefit_kinematics[3][sample][topology]->Fill(pCosThetaRec,weight);
            hMC_Prefit_kinematics[4][sample][topology]->Fill(piMomRec,weight);
            hMC_Prefit_kinematics[5][sample][topology]->Fill(piCosThetaRec,weight);
            hMC_Prefit_kinematics[6][sample][topology]->Fill(NTPCproton,weight);
        }
        if (abs(D1Reco)<tkimax) {
            nMCEvents[sample][topology]+=weight;
            nMCEventsSum[sample]+=weight;
        }
    }

    THStack *hs_evhist_prefit[nSamples];
    THStack *hs_evhist_prefit_kinematics[nKinematics][nSamples];
    for (int i=0;i<nSamples;i++) {
        hs_evhist_prefit[i]= new THStack("","");
        for (int j=0;j<nKinematics;j++) {
            hs_evhist_prefit_kinematics[j][i]= new THStack("","");
        }
    }

    THStack *hs_evhist_prefit_pmom = new THStack("","");
    std::string topo_multipi_types[nTopo]={"CC0#pi","CC1#pi^{+}0p","CC1#pi^{+}1p","CC1#pi^{+}Np","","CC1#pi^{+}1#pi^{-}","CC1#pi^{+}X#pi^{0}","CC-other-X#pi^{0}","CC-other-0#pi^{0}","NC, #bar{#nu}_{#mu}, #nu_{e}, #bar{#nu}_{e}","OOFV","OOPS"};
    std::string samName[nSamples]={"Signal sample","CC1#pi^{+}1#pi^{-} enriched","CC1#pi^{+}X#pi^{0} enriched","CC-other-X#pi^{0} enriched","CC-other-0#pi^{0} enriched"};
    int topology_colors[nTopo]= {2        , 3              , 4              , 7              , 0, 31,  51                   , 5                 ,48           ,6 , 1,8};
    for (int i=0;i<nTopo;i++) {
        for (int j=0;j<nSamples;j++){
            for (int k=0;k<nKinematics;k++) {
                hMC_Prefit_kinematics[k][j][nTopo-1-i]->SetLineColor(topology_colors[nTopo-1-i]);
                hMC_Prefit_kinematics[k][j][nTopo-1-i]->SetFillColor(topology_colors[nTopo-1-i]);
                hMC_Prefit_kinematics[k][j][nTopo-1-i]->SetFillStyle(3154);
                hs_evhist_prefit_kinematics[k][j]->Add(hMC_Prefit_kinematics[k][j][nTopo-1-i]);
            }
        }
    }

    selectedEventsData->SetBranchAddress("cutBranch", &sample);
    selectedEventsData->SetBranchAddress("D1Rec", &D1Reco);
    selectedEventsData->SetBranchAddress("D2Rec", &D2Reco);
    selectedEventsData->SetBranchAddress("weight", &weight);
    selectedEventsData->SetBranchAddress("muMomRec", &muMomRec);
    selectedEventsData->SetBranchAddress("muCosThetaRec", &muCosThetaRec);
    selectedEventsData->SetBranchAddress("pMomRec", &pMomRec);
    selectedEventsData->SetBranchAddress("pCosThetaRec", &pCosThetaRec);
    selectedEventsData->SetBranchAddress("piMomRec", &piMomRec);
    selectedEventsData->SetBranchAddress("piCosThetaRec", &piCosThetaRec);
    selectedEventsData->SetBranchAddress("NTPCproton", &NTPCproton);
    for (int i=0;i<selectedEventsData->GetEntries();i++) { 
        selectedEventsData->GetEntry(i);
        weight = 1;
        if (abs(D1Reco)<tkimax) {
            hData_kinematics[0][sample]->Fill(muMomRec,weight);
            hData_kinematics[1][sample]->Fill(muCosThetaRec,weight);
            hData_kinematics[2][sample]->Fill(pMomRec,weight);
            hData_kinematics[3][sample]->Fill(pCosThetaRec,weight);
            hData_kinematics[4][sample]->Fill(piMomRec,weight);
            hData_kinematics[5][sample]->Fill(piCosThetaRec,weight);
            hData_kinematics[6][sample]->Fill(NTPCproton,weight);
        }
    }
    for (int i=0;i<nSamples;i++) for (int j=0;j<nKinematics;j++) {
        hData_kinematics[j][i]->SetLineColor(kBlack);
        hData_kinematics[j][i]->SetMarkerStyle(kFullCircle);
    }
    
    std::string xlabel_kinematics[]={
        "Reconstructed #it{p_{#mu}} (MeV/c)",
        "Reconstructed cos#it{#theta_{#mu}} (MeV/c)",
        "Reconstructed #it{p}_{p} (MeV/c)",
        "Reconstructed cos#it{#theta}_{p} (MeV/c)",
        "Reconstructed #it{p_{#pi}} (MeV/c)",
        "Reconstructed cos#it{#theta_{#pi}} (MeV/c)",
        "Reconstructed proton multiplicity"
    };
    std::string varName_kinematics[]={
        "pmu","cthmu","pp","cthp","ppi","cthpi","np"
    };
    for (int j=0;j<nSamples;j++) for (int k=0;k<nKinematics;k++) {
        if (j>0&&k<nKinematics-1) continue;
        //if (j>0) continue;
        //if (j==0) c1 = new TCanvas();
        //else c1 = new TCanvas("","",500,500);
        TCanvas* c1 = new TCanvas();
        std::cout<<"Sample "<<j<<", "<<hData_kinematics[k][j]->Integral()<<" events."<<std::endl;
        double ymax = 1.2*std::max(hs_evhist_prefit_kinematics[k][j]->GetMaximum(),hData_kinematics[k][j]->GetMaximum());
        hs_evhist_prefit_kinematics[k][j]->SetMaximum(ymax);
        hs_evhist_prefit_kinematics[k][j]->Draw("hist");
        hs_evhist_prefit_kinematics[k][j]->GetXaxis()->SetTitle(xlabel_kinematics[k].c_str());
        hs_evhist_prefit_kinematics[k][j]->GetYaxis()->SetTitle("Events per bin");
        //hs_evhist_prefit[j]->GetYaxis()->SetRangeUser(0,ymax);
        //c1->Update();
        hData_kinematics[k][j]->SetMarkerSize(1.3);
        hData_kinematics[k][j]->Draw("E1 X0 same");
        TPaveLabel *title;
        if (j>0) title = new TPaveLabel(.05,.81,.8,1.1,samName[j].c_str(),"NDC"); 
        else title = new TPaveLabel(.05,.9,.8,1.,samName[j].c_str(),"NDC"); 
        title->SetBorderSize(0);
        title->SetFillStyle(0);
        if (j>0||k==nKinematics-1) title->Draw();
        TLegend* legendTopo;
        //if (!(TKI==1&&j==0)&&!(TKI==2&&j>0)) continue;
        if (TKI!=2) legendTopo = new TLegend(0.5,0.34,0.9,0.89);
        //else if (TKI==2) legendTopo = new TLegend(0.18,0.39,0.58,0.89);
        else legendTopo = new TLegend(0.,0.,1,1);
        //if (j>0) legendTopo->AddEntry((TObject*)0, samName[j].c_str(), "");
        legendTopo->AddEntry(hData_kinematics[k][j],"Data","lep");
        for (int i=0;i<nTopo;i++){
            if (i==4) continue;
            double frac = nMCEvents[j][i]/nMCEventsSum[j]*100.;
            //std::cout<<"frac = "<<nMCEvents[j][i]<<"/"<<nMCEventsSum[j]<<"="<<frac<<std::endl;
            legendTopo->AddEntry(hMC_Prefit_kinematics[k][j][i],Form("%s (%4.1f%%)",topo_multipi_types[i].c_str(),frac),"f");
        }
        legendTopo->SetBorderSize(0);
        if (TKI==2) legendTopo->SetFillStyle(0);
        if (j==0&&k==0) legendTopo->Draw();
        if (k==nKinematics-1) legendTopo->Draw();
        if (j==0) c1->SaveAs(Form("%s_sig.pdf",varName_kinematics[k].c_str()));
        if (j>0) c1->SaveAs(Form("%s_CR%i.pdf",varName_kinematics[k].c_str(),j));
    }
    if (TKI!=0) return;
    
}


void makePlots_TKI(int TKI){
    TStyle *t2kstyle = new TStyle("T2K","T2K approved plots style");
    SetStyleVariables(t2kstyle);
    gROOT->SetStyle("T2K");

    const int nDpttbins = 5;
    const double Dpttbins[nDpttbins+1] = {-700,-300,-100,100,300,700};
    const int nPnbins = 4;
    const double Pnbins[nPnbins+1] = {0,120,240,600,1500};
    const int nDatbins = 3;
    const double Datbins[nDatbins+1] = {0,60,120,180};
    double Bins[100];
    int nBins = -1;
    std::string mcfile;
    std::string datafile;
    double tkimax = -1;
    std::string ylabel;
    std::string xlabel;
    std::string varName;
    switch (TKI) {
        case 0:
            mcfile = "test_dptt_prod6T_geomfixed.root";
            datafile = "data_dptt_prod6T.root";
            nBins = nDpttbins;
            for (int i=0;i<nDpttbins+1;i++) Bins[i]=Dpttbins[i];
            tkimax = 700;
            xlabel = "Reconstructed #it{#deltap_{TT}} (MeV/c)";
            varName = "dptt";
            break;
        case 1:
            mcfile = "test_pN_prod6T_geomfixed.root";
            datafile = "data_pN_prod6T.root";
            nBins = nPnbins;
            for (int i=0;i<nPnbins+1;i++) Bins[i]=Pnbins[i];
            tkimax = 1500;
            xlabel = "Reconstructed #it{p_{N}} (MeV/c)";
            varName = "pN";
            break;
        case 2:
            mcfile = "test_dat_prod6T_geomfixed.root";
            datafile = "data_dat_prod6T.root";
            nBins = nDatbins;
            for (int i=0;i<nDatbins+1;i++) Bins[i]=Datbins[i];
            tkimax = 180;
            xlabel = "Reconstructed #it{#delta#alpha_{T}} (deg)";
            varName = "daT";
            break;
        default:
            printf("***Warning: Not a valid TKI variable***\n");
            break;
    }

    TFile* fMC = new TFile(mcfile.c_str());
    TTree* selectedEvents = (TTree*)fMC->Get("selectedEvents");
    TFile* fData = new TFile(datafile.c_str());
    TTree* selectedEventsData = (TTree*)fData->Get("selectedEvents");

    const int nPmombins = 6;
    const double Pmombins[nPmombins+1] = {405,575,700,825,950,1075,1320};
    const int nTopo = 12;
    const int nSamples = 5;
    TH1D* hMC_Prefit[nSamples][nTopo];
    TH2D* hMC_2D[nSamples][3];
    TH1D* hMC_Prefit_pmom[nTopo];
    for (int i=0;i<nTopo;i++) {
        hMC_Prefit_pmom[i] = new TH1D("","",nPmombins,Pmombins);
        for (int j=0;j<nSamples;j++)
            hMC_Prefit[j][i] = new TH1D("","",nBins,Bins);
    }
    TH1D* hData[nSamples];
    for (int j=0;j<nSamples;j++){
        hData[j] = new TH1D("","",nBins,Bins);
        for (int i=0;i<3;i++) hMC_2D[j][i] = new TH2D("","",nBins,Bins,nPmombins,Pmombins);
    }
    TH1D* hData_pmom = new TH1D("","",nPmombins,Pmombins);

    Int_t sample;
    Float_t D1true;
    Float_t D2true;
    Int_t nutype;
    Int_t topology;
    Int_t reaction;
    Int_t target;
    Float_t D1Reco;
    Float_t D2Reco;
    Float_t weight;
    Float_t weightNom;
    Float_t weightMC;

    selectedEvents->SetBranchAddress("nutype", &nutype);
    selectedEvents->SetBranchAddress("reaction", &reaction);
    selectedEvents->SetBranchAddress("target", &target);
    selectedEvents->SetBranchAddress("sample", &sample);
    selectedEvents->SetBranchAddress("topology", &topology);
    selectedEvents->SetBranchAddress("D1True", &D1true);
    selectedEvents->SetBranchAddress("D1Reco", &D1Reco);
    selectedEvents->SetBranchAddress("D2True", &D2true);
    selectedEvents->SetBranchAddress("D2Reco", &D2Reco);
    selectedEvents->SetBranchAddress("weight", &weight);
    selectedEvents->SetBranchAddress("weightNom", &weightNom);
    selectedEvents->SetBranchAddress("weightMC", &weightMC);

    Float_t nMCEvents[nSamples][nTopo]={0};
    Float_t nMCEventsSum[nSamples]={0};
    for (int i=0;i<selectedEvents->GetEntries();i++) { 
        selectedEvents->GetEntry(i);
        if(topology<0&&topology>11) continue;
        hMC_Prefit[sample][topology]->Fill(D1Reco,weightMC);
        hMC_2D[sample][0]->Fill(D1Reco,D2Reco,weightMC);
        hMC_2D[sample][1]->Fill(D1Reco,D2Reco,weight);
        if (abs(D1Reco)<tkimax) {
            nMCEvents[sample][topology]+=weightMC;
            nMCEventsSum[sample]+=weightMC;
        }
        if (sample==0) {
            if (abs(D1Reco)<tkimax) hMC_Prefit_pmom[topology]->Fill(D2Reco,weightMC);
        }
    }
    THStack *hs_evhist_prefit[nSamples];
    for (int i=0;i<nSamples;i++) hs_evhist_prefit[i]= new THStack("","");
    THStack *hs_evhist_prefit_pmom = new THStack("","");
    std::string topo_multipi_types[nTopo]={"CC0#pi","CC1#pi^{+}0p","CC1#pi^{+}1p","CC1#pi^{+}Np","","CC1#pi^{+}1#pi^{-}","CC1#pi^{+}X#pi^{0}","CC-other-X#pi^{0}","CC-other-0#pi^{0}","NC, #bar{#nu}_{#mu}, #nu_{e}, #bar{#nu}_{e}","OOFV","OOPS"};
    std::string samName[nSamples]={"Signal","CC1#pi^{+}1#pi^{-} enriched","CC1#pi^{+}X#pi^{0} enriched","CC-other-X#pi^{0} enriched","CC-other-0#pi^{0} enriched"};
    int topology_colors[nTopo]= {2        , 3              , 4              , 7              , 0, 31,  51                   , 5                 ,48           ,6 , 1,8};
    for (int i=0;i<nTopo;i++) {
        for (int j=0;j<nSamples;j++){
            hMC_Prefit[j][nTopo-1-i]->SetLineColor(topology_colors[nTopo-1-i]);
            hMC_Prefit[j][nTopo-1-i]->SetFillColor(topology_colors[nTopo-1-i]);
            hMC_Prefit[j][nTopo-1-i]->SetFillStyle(3154);
            hs_evhist_prefit[j]->Add(hMC_Prefit[j][nTopo-1-i]);
        }
        hMC_Prefit_pmom[nTopo-1-i]->SetLineColor(topology_colors[nTopo-1-i]);
        hMC_Prefit_pmom[nTopo-1-i]->SetFillColor(topology_colors[nTopo-1-i]);
        hMC_Prefit_pmom[nTopo-1-i]->SetFillStyle(3154);
        hs_evhist_prefit_pmom->Add(hMC_Prefit_pmom[nTopo-1-i]);
    }

    selectedEventsData->SetBranchAddress("cutBranch", &sample);
    selectedEventsData->SetBranchAddress("D1Rec", &D1Reco);
    selectedEventsData->SetBranchAddress("D2Rec", &D2Reco);
    selectedEventsData->SetBranchAddress("weight", &weight);
    for (int i=0;i<selectedEventsData->GetEntries();i++) { 
        selectedEventsData->GetEntry(i);
        hData[sample]->Fill(D1Reco,weight);
        hMC_2D[sample][2]->Fill(D1Reco,D2Reco,weight);
        if (sample==0) {
            if (abs(D1Reco)<tkimax) hData_pmom->Fill(D2Reco,weight);
        }
    }
    for (int i=0;i<nSamples;i++) {
        hData[i]->SetLineColor(kBlack);
        hData[i]->SetMarkerStyle(kFullCircle);
    }
    hData_pmom->SetLineColor(kBlack);
    hData_pmom->SetMarkerStyle(kFullCircle);

    TH1D* hPrefit1D[nPmombins];
    TH1D* hPostfit1D[nPmombins];
    TH1D* hData1D[nPmombins];
    for (int i=0;i<nPmombins;i++) {
        hPrefit1D[i] = hMC_2D[0][0]->ProjectionX(Form("Prefit_Bin%i",i+1),i+1,i+1);
        hPostfit1D[i] = hMC_2D[0][1]->ProjectionX(Form("Postfit_Bin%i",i+1),i+1,i+1);
        hData1D[i] = hMC_2D[0][2]->ProjectionX(Form("Data_Bin%i",i+1),i+1,i+1);
        double* exp_w  = hPostfit1D[i]->GetArray();
        double* exp_w2 = hPostfit1D[i]->GetSumw2()->GetArray();
        double* exp_w_prefit  = hPrefit1D[i]->GetArray();
        double* exp_w2_prefit = hPrefit1D[i]->GetSumw2()->GetArray();
        double* data   = hData1D[i]->GetArray();
        double chi2 = 0.0;
        double chi2_prefit = 0.0;
        for (int j=1;j<=hData1D[i]->GetNbinsX();j++){
            chi2 += BarlowLLH(exp_w[j], exp_w2[j], data[j]);
            chi2_prefit += BarlowLLH(exp_w_prefit[j], exp_w2_prefit[j], data[j]);
        }
        //std::cout<<"chi2 = "<<chi2<<std::endl;
        TCanvas* c = new TCanvas("","",700,600);
        c->SetFrameBorderMode(0);
        hPrefit1D[i]->GetXaxis()->SetTitle(xlabel.c_str());
        hPrefit1D[i]->GetYaxis()->SetTitle("Events per bin");
        hPrefit1D[i]->SetLineColor(kBlack);
        hPrefit1D[i]->SetLineStyle(9);
        hPrefit1D[i]->SetMinimum(0);
        hPostfit1D[i]->SetLineColor(kBlue);
        hData1D[i]->SetLineColor(kBlack);
        hData1D[i]->SetMarkerStyle(kFullCircle);
        hData1D[i]->SetMarkerSize(1.3);
        double maxy = std::max(hPrefit1D[i]->GetMaximum(),hPostfit1D[i]->GetMaximum());
        maxy = std::max(hData1D[i]->GetMaximum(),maxy);
        maxy*=1.3;
        if (TKI==1) maxy*=1.1;
        hPrefit1D[i]->SetMaximum(maxy);
        hPrefit1D[i]->Draw("hist");
        hPostfit1D[i]->Draw("hist same");
        hData1D[i]->Draw("E01 X0 same");
        TPaveLabel *title = new TPaveLabel(.1,.9,1.,.99,Form("Signal sample, %iMeV/c<p_{p}<%iMeV/c",(int)Pmombins[i],(int)Pmombins[i+1]),"NDC"); 
        title->SetBorderSize(0);
        title->SetFillStyle(0);
        title->SetTextSize(0.7);
        title->Draw(); 
        TPaveText *chi2label;
        if (TKI==0||TKI==1) chi2label = new TPaveText(.65,.6,0.9,.85,"blNDC");
        else chi2label = new TPaveText(.41,.6,0.66,.85,"blNDC");
        chi2label->AddText(Form("Pre-fit #chi^{2}_{stat}=%4.1f",chi2_prefit)); 
        chi2label->AddText(Form("Post-fit #chi^{2}_{stat}=%4.1f",chi2)); 
        chi2label->SetBorderSize(1);
        chi2label->SetFillStyle(0);
        //chi2label->SetTextSize(0.7);
        chi2label->Draw(); 
        //if (i>0) continue;
        TLegend* legendFit;
        if (TKI==1) legendFit = new TLegend(0.4,0.61,0.7,0.91);
        else if (TKI==0) legendFit = new TLegend(0.20,0.53,0.5,0.83);
        else if (TKI==2) legendFit = new TLegend(0.17,0.48,0.47,0.78);
        legendFit->AddEntry(hData1D[i],"Data","lep");
        legendFit->AddEntry(hPrefit1D[i],"Pre-fit","l");
        legendFit->AddEntry(hPostfit1D[i],"Post-fit","l");
        legendFit->SetBorderSize(0);
        legendFit->SetFillStyle(0);
        if (i==0) legendFit->Draw();
        c->SaveAs(Form("%s_sig_postfit_p%i.pdf",varName.c_str(),i+1));
    }

    for (int i=1;i<nSamples;i++) {
        hPrefit1D[i] = hMC_2D[i][0]->ProjectionX(Form("Prefit_CR%i",i),1,nPmombins);
        hPostfit1D[i] = hMC_2D[i][1]->ProjectionX(Form("Postfit_CR%i",i),1,nPmombins);
        hData1D[i] = hMC_2D[i][2]->ProjectionX(Form("Data_CR%i",i),1,nPmombins);
        double* exp_w  = hPostfit1D[i]->GetArray();
        double* exp_w2 = hPostfit1D[i]->GetSumw2()->GetArray();
        double* data   = hData1D[i]->GetArray();
        double* exp_w_prefit  = hPrefit1D[i]->GetArray();
        double* exp_w2_prefit = hPrefit1D[i]->GetSumw2()->GetArray();
        double chi2 = 0.0;
        double chi2_prefit = 0.0;
        for (int j=1;j<=hData1D[i]->GetNbinsX();j++){
            chi2 += BarlowLLH(exp_w[j], exp_w2[j], data[j]);
            chi2_prefit += BarlowLLH(exp_w_prefit[j], exp_w2_prefit[j], data[j]);
        }
        TCanvas* c = new TCanvas("","",700,600);
        c->SetFrameBorderMode(0);
        hData1D[i]->GetXaxis()->SetTitle(xlabel.c_str());
        hData1D[i]->GetYaxis()->SetTitle("Events per bin");
        hPrefit1D[i]->SetLineColor(kBlack);
        hPrefit1D[i]->SetLineStyle(9);
        hData1D[i]->SetMinimum(0);
        hPostfit1D[i]->SetLineColor(kBlue);
        hData1D[i]->SetLineColor(kBlack);
        hData1D[i]->SetMarkerStyle(kFullCircle);
        hData1D[i]->SetMarkerSize(1.3);
        double maxy = std::max(hPrefit1D[i]->GetMaximum(),hPostfit1D[i]->GetMaximum());
        maxy = std::max(hData1D[i]->GetMaximum(),maxy);
        maxy*=1.3;
        hData1D[i]->SetMaximum(maxy);
        hData1D[i]->Draw("E01 X0 same");
        hPrefit1D[i]->Draw("hist same");
        hPostfit1D[i]->Draw("hist same");
        TPaveLabel *title = new TPaveLabel(.3,.9,.7,1.,Form("%s sample",samName[i].c_str()),"NDC"); 
        title->SetBorderSize(0);
        title->SetFillStyle(0);
        title->SetTextSize(0.7);
        title->Draw(); 
        TPaveText *chi2label;
        if (TKI==0) chi2label = new TPaveText(.7,.6,0.92,.85,"blNDC");
        else if (TKI==1) chi2label = new TPaveText(.65,.2,0.9,.45,"blNDC");
        else chi2label = new TPaveText(.2,.6,0.45,.85,"blNDC");
        chi2label->AddText(Form("Pre-fit #chi^{2}_{stat}=%4.1f",chi2_prefit)); 
        chi2label->AddText(Form("Post-fit #chi^{2}_{stat}=%4.1f",chi2)); 
        chi2label->SetBorderSize(1);
        chi2label->SetFillStyle(0);
        chi2label->Draw(); 
        c->SaveAs(Form("%s_cr%i_postfit.pdf",varName.c_str(),i));
    }

    for (int j=0;j<nSamples;j++) {
        TCanvas* c1;
        if (j==0) c1 = new TCanvas();
        else c1 = new TCanvas("","",500,500);
        std::cout<<"Sample "<<j<<", "<<hData[j]->Integral()<<" events."<<std::endl;
        double ymax = 1.2*std::max(hs_evhist_prefit[j]->GetMaximum(),hData[j]->GetMaximum());
        hs_evhist_prefit[j]->SetMaximum(ymax);
        hs_evhist_prefit[j]->Draw("hist");
        hs_evhist_prefit[j]->GetXaxis()->SetTitle(xlabel.c_str());
        hs_evhist_prefit[j]->GetYaxis()->SetTitle("Events per bin");
        //hs_evhist_prefit[j]->GetYaxis()->SetRangeUser(0,ymax);
        //c1->Update();
        hData[j]->SetMarkerSize(1.3);
        hData[j]->Draw("E1 X0 same");
        TPaveLabel *title = new TPaveLabel(.05,.81,.8,1.1,samName[j].c_str(),"NDC"); 
        title->SetBorderSize(0);
        title->SetFillStyle(0);
        if (j>0) title->Draw();
        TLegend* legendTopo;
        //if (!(TKI==1&&j==0)&&!(TKI==2&&j>0)) continue;
        if (TKI!=2) legendTopo = new TLegend(0.5,0.45,0.9,0.99);
        //else if (TKI==2) legendTopo = new TLegend(0.18,0.39,0.58,0.89);
        else legendTopo = new TLegend(0.,0.,1,1);
        //if (j>0) legendTopo->AddEntry((TObject*)0, samName[j].c_str(), "");
        legendTopo->AddEntry(hData[j],"Data","lep");
        for (int i=0;i<nTopo;i++){
            if (i==4) continue;
            double frac = nMCEvents[j][i]/nMCEventsSum[j]*100.;
            //std::cout<<"frac = "<<nMCEvents[j][i]<<"/"<<nMCEventsSum[j]<<"="<<frac<<std::endl;
            legendTopo->AddEntry(hMC_Prefit[j][i],Form("%s (%4.1f%%)",topo_multipi_types[i].c_str(),frac),"f");
        }
        legendTopo->SetBorderSize(0);
        if (TKI==2) legendTopo->SetFillStyle(0);
        if (j==0&&TKI==1) legendTopo->Draw();
        if (j==0) c1->SaveAs(Form("%s_sig.pdf",varName.c_str()));
        else c1->SaveAs(Form("%s_cr%i.pdf",varName.c_str(),j));
        if (j>0&&TKI==2){
            TCanvas* clegend = new TCanvas(Form("cl%i",j),Form("cl%i",j),235,320);
            legendTopo->Draw();
            clegend->SaveAs(Form("legend_cr%i.pdf",j));
        }
    }
    if (TKI!=0) return;
    TCanvas* c2 = new TCanvas();
    double ymax = 1.2*std::max(hs_evhist_prefit_pmom->GetMaximum(),hData_pmom->GetMaximum());
    hs_evhist_prefit_pmom->SetMaximum(ymax);
    hs_evhist_prefit_pmom->Draw("hist");
    hs_evhist_prefit_pmom->GetXaxis()->SetTitle("Reconstructed #it{p}_{p} (MeV/c)");
    hs_evhist_prefit_pmom->GetYaxis()->SetTitle("Events per bin");
    hData_pmom->SetMarkerSize(1.3);
    hData_pmom->Draw("E1 X0 same");
    c2->SaveAs("pmom_sig.pdf");
}

void makePlots_Samples(){
    TH1::SetDefaultSumw2(true);
    makePlots_kinematics();
    //makePlots_TKI(0);
    //makePlots_TKI(1);
    //makePlots_TKI(2);
}