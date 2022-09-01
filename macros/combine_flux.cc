void combine_flux()
{
    const int nFiles = 5;
    std::string filename[]={
        "/user/kmtsui/Downloads/tuned13av7p1/run2/nd5_tuned13av7p1_13anom_run2_numode.root",
        "/user/kmtsui/Downloads/tuned13av7p1/run3b/nd5_tuned13av7p1_13anom_run3b_numode.root",
        "/user/kmtsui/Downloads/tuned13av7p1/run3c/nd5_tuned13av7p1_13anom_run3c_numode.root",
        "/user/kmtsui/Downloads/tuned13av7p1/run4/nd5_tuned13av7p1_13anom_run4_numode.root",
        "/user/kmtsui/Downloads/tuned13av7p1/run8/nd5_tuned13av7p1_13anom_run8_numode.root",
    };
    double norm[]={
        0.79/11.53,
        0.2137667/11.53,
        1.371495/11.53,
        3.42/11.53,
        5.73/11.53
    };

    TH1D* fluxhist[nFiles];
    for (int i=0;i<nFiles;i++){
        TFile* f = new TFile(filename[i].c_str());
        fluxhist[i]=(TH1D*)f->Get("enu_nd5_tuned13a_numu");
    }

    TFile* out_file = TFile::Open("flux_combine.root", "RECREATE");
    TH1D* fluxhist_combine = (TH1D*)fluxhist[0]->Clone();
    fluxhist_combine->Reset();
    for (int i=0;i<nFiles;i++)
        fluxhist_combine->Add(fluxhist[i],norm[i]);
    out_file->cd();
    fluxhist_combine->Write("enu_nd5_tuned13a_numu");
    out_file->Close();
   
}
