#include "AnaSample.hh"
using xsllh::FitBin;

// ctor
AnaSample::AnaSample(int sample_id, const std::string& name, const std::string& detector,
                     const std::string& binning, TTree* t_data, int llh_func)
    : m_sample_id(sample_id)
    , m_name(name)
    , m_detector(detector)
    , m_binning(binning)
    , m_data_tree(t_data)
    , m_norm(1.0)
    , m_llh_func(llh_func)
{
    
    TH1::SetDefaultSumw2(true);
    SetBinning(m_binning);

    bin_manager.SetBinning("/hepstore/kmtsui/T2K/work/xsLLhFitter_super/inputs/p_pi_mu_mom.txt");
    m_nDetBins = bin_manager.GetNbins();

    std::cout << TAG << m_name << ", ID " << m_sample_id << std::endl
              << TAG << "Detector: " << m_detector << std::endl
              << TAG << "Bin edges: " << std::endl;

    for(const auto& bin : m_bin_edges)
    {
        std::cout << bin.D2low << " " << bin.D2high << " " << bin.D1low << " " << bin.D1high
                  << std::endl;
    }

    m_hpred    = nullptr;
    m_hmc      = nullptr;
    m_hmc_true = nullptr;
    m_hsig     = nullptr;
    m_hdata    = nullptr;

    m_hpred_pmom    = nullptr;
    m_hdata_pmom    = nullptr;
    m_hpred_ptheta    = nullptr;
    m_hdata_ptheta    = nullptr;
    m_hpred_pimom    = nullptr;
    m_hdata_pimom    = nullptr;
    m_hpred_pitheta    = nullptr;
    m_hdata_pitheta    = nullptr;

  for (int i=0;i<100;i++) {
    m_hpred_p_ptheta[i]    = nullptr;
    m_hdata_p_ptheta[i]    = nullptr;
  }

    m_hpred_pi_ptheta    = nullptr;
    m_hdata_pi_ptheta    = nullptr;

    m_hpred_2d = nullptr;
    m_hdata_2d = nullptr;

    MakeHistos(); // with default binning

    std::cout << TAG << "MakeHistos called." << std::endl;
}

AnaSample::~AnaSample()
{
    if(m_hpred != nullptr)
        delete m_hpred;
    if(m_hmc != nullptr)
        delete m_hmc;
    if(m_hmc_true != nullptr)
        delete m_hmc_true;
    if(m_hdata != nullptr)
        delete m_hdata;
    if(m_hsig != nullptr)
        delete m_hsig;

    if(m_hpred_pmom != nullptr)
        delete m_hpred_pmom;
    if(m_hdata_pmom != nullptr)
        delete m_hdata_pmom;
    if(m_hpred_ptheta != nullptr)
        delete m_hpred_ptheta;
    if(m_hdata_ptheta != nullptr)
        delete m_hdata_ptheta;
    if(m_hpred_pimom != nullptr)
        delete m_hpred_pimom;
    if(m_hdata_pimom != nullptr)
        delete m_hdata_pimom;
    if(m_hpred_pitheta != nullptr)
        delete m_hpred_pitheta;
    if(m_hdata_pitheta != nullptr)
        delete m_hdata_pitheta;

  for (int i=0;i<m_nbins;i++) {
    if(m_hpred_p_ptheta[i] != nullptr)
        delete m_hpred_p_ptheta[i];
    if(m_hdata_p_ptheta[i] != nullptr)
        delete m_hdata_p_ptheta[i];
  }
    if(m_hpred_pi_ptheta != nullptr)
        delete m_hpred_pi_ptheta;
    if(m_hdata_pi_ptheta != nullptr)
        delete m_hdata_pi_ptheta;

    if(m_hpred_2d != nullptr)
        delete m_hpred_2d;
    if(m_hdata_2d != nullptr)
        delete m_hdata_2d;

}

void AnaSample::SetBinning(const std::string& binning)
{
    m_binning = binning;
    m_nbins   = 0;

    std::ifstream fin(m_binning, std::ios::in);
    if(!fin.is_open())
    {
        std::cerr << ERR << "In AnaSample::SetBinning().\n"
                  << ERR << "Failed to open binning file: " << m_binning << std::endl;
    }
    else
    {
        std::string line;
        while(std::getline(fin, line))
        {
            std::stringstream ss(line);
            double D1_1, D1_2, D2_1, D2_2;
            if(!(ss >> D2_1 >> D2_2 >> D1_1 >> D1_2))
            {
                std::cerr << TAG << "Bad line format: " << line << std::endl;
                continue;
            }
            m_bin_edges.emplace_back(FitBin(D1_1, D1_2, D2_1, D2_2));
        }
        m_nbins = m_bin_edges.size();
    }
}

// ClearEvents -- clears all events from event vector
void AnaSample::ClearEvents() { m_events.clear(); }

// GetN -- get number of events stored
int AnaSample::GetN() const { return (int)m_events.size(); }

AnaEvent* AnaSample::GetEvent(int evnum)
{
    if(m_events.empty())
    {
        std::cerr << "[ERROR]: In AnaSample::GetEvent()" << std::endl;
        std::cerr << "[ERROR]: No events are found in " << m_name << " sample." << std::endl;
        return nullptr;
    }
    else if(evnum >= m_events.size())
    {
        std::cerr << "[ERROR]: In AnaSample::GetEvent()" << std::endl;
        std::cerr << "[ERROR]: Event number out of bounds in " << m_name << " sample." << std::endl;
        return nullptr;
    }

    return &m_events.at(evnum);
}

void AnaSample::AddEvent(const AnaEvent& event) { m_events.push_back(event); }

void AnaSample::ResetWeights()
{
    for(auto& event : m_events)
        event.SetEvWght(1.0);
}

void AnaSample::PrintStats() const
{
    double mem_kb = sizeof(m_events) * m_events.size() / 1000.0;
    std::cout << TAG << "Sample " << m_name << " ID = " << m_sample_id << std::endl;
    std::cout << TAG << "Num of events = " << m_events.size() << std::endl;
    std::cout << TAG << "Memory used   = " << mem_kb << " kB." << std::endl;
}

void AnaSample::MakeHistos()
{
    if(m_hpred != nullptr)
        delete m_hpred;
    m_hpred = new TH1D(Form("%s_pred_recD1D2", m_name.c_str()),
                       Form("%s_pred_recD1D2", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hpred->SetDirectory(0);

    if(m_hmc != nullptr)
        delete m_hmc;
    m_hmc = new TH1D(Form("%s_mc_recD1D2", m_name.c_str()), Form("%s_mc_recD1D2", m_name.c_str()),
                     m_nbins, 0, m_nbins);
    m_hmc->SetDirectory(0);

    if(m_hmc_true != nullptr)
        delete m_hmc_true;
    m_hmc_true = new TH1D(Form("%s_mc_trueD1D2", m_name.c_str()),
                          Form("%s_mc_trueD1D2", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hmc_true->SetDirectory(0);

    if(m_hsig != nullptr)
        delete m_hsig;
    m_hsig = new TH1D(Form("%s_mc_trueSignal", m_name.c_str()),
                      Form("%s_mc_trueSignal", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hsig->SetDirectory(0);


    const int nPMombins = 3;
    const double PMombins[nPMombins+1] = {450,700,1000,1200};
    const double PRecMombins[nPMombins+1] = {405,700,1000,1320};
    const int nPiMombins = 3;
    const double PiMombins[nPiMombins+1] = {135,400,600,1320};
    const int nPThetabins = 3;
    const double PThetabins[nPThetabins+1] = {0.342,0.6,0.8,1.};

    if(m_hpred_pmom != nullptr)
        delete m_hpred_pmom;
    m_hpred_pmom = new TH1D(Form("%s_pred_pmom", m_name.c_str()),
//                       Form("%s_pred_pmom", m_name.c_str()), nPMombins, PRecMombins);
                       Form("%s_pred_pmom", m_name.c_str()), 6, 405, 1320);
    m_hpred_pmom->SetDirectory(0);

    if(m_hdata_pmom != nullptr)
        delete m_hdata_pmom;
    m_hdata_pmom = new TH1D(Form("%s_data_pmom", m_name.c_str()),
//                       Form("%s_data_pmom", m_name.c_str()), nPMombins, PRecMombins);
                       Form("%s_data_pmom", m_name.c_str()), 6, 405, 1320);
    m_hdata_pmom->SetDirectory(0);

    if(m_hpred_ptheta != nullptr)
        delete m_hpred_ptheta;
    m_hpred_ptheta = new TH1D(Form("%s_pred_ptheta", m_name.c_str()),
//                       Form("%s_pred_ptheta", m_name.c_str()), nPThetabins, PThetabins);
                       Form("%s_pred_ptheta", m_name.c_str()), 6, 0.342, 1.);
    m_hpred_ptheta->SetDirectory(0);

    if(m_hdata_ptheta != nullptr)
        delete m_hdata_ptheta;
    m_hdata_ptheta = new TH1D(Form("%s_data_ptheta", m_name.c_str()),
//                       Form("%s_data_ptheta", m_name.c_str()), nPThetabins, PThetabins);
                       Form("%s_data_ptheta", m_name.c_str()), 6, 0.342, 1.);
    m_hdata_ptheta->SetDirectory(0);

    if(m_hpred_pimom != nullptr)
        delete m_hpred_pimom;
    m_hpred_pimom = new TH1D(Form("%s_pred_pimom", m_name.c_str()),
                       Form("%s_pred_pimom", m_name.c_str()), 6, 135, 1320);
    m_hpred_pimom->SetDirectory(0);

    if(m_hdata_pimom != nullptr)
        delete m_hdata_pimom;
    m_hdata_pimom = new TH1D(Form("%s_data_pimom", m_name.c_str()),
                       Form("%s_data_pimom", m_name.c_str()), 6, 135, 1320);
    m_hdata_pimom->SetDirectory(0);

    if(m_hpred_pitheta != nullptr)
        delete m_hpred_pitheta;
    m_hpred_pitheta = new TH1D(Form("%s_pred_pitheta", m_name.c_str()),
                       Form("%s_pred_pitheta", m_name.c_str()), 6, 0.342, 1.);
    m_hpred_pitheta->SetDirectory(0);

    if(m_hdata_pitheta != nullptr)
        delete m_hdata_pitheta;
    m_hdata_pitheta = new TH1D(Form("%s_data_pitheta", m_name.c_str()),
                       Form("%s_data_pitheta", m_name.c_str()), 6, 0.342, 1.);
    m_hdata_pitheta->SetDirectory(0);

  for (int i=0;i<m_nbins;i++) {
    if(m_hpred_p_ptheta[i] != nullptr)
        delete m_hpred_p_ptheta[i];
    m_hpred_p_ptheta[i] = new TH2D(Form("%s_pred_p_ptheta_bin%d", m_name.c_str(), i),
                       Form("%s_pred_p_ptheta_bin%d", m_name.c_str(), i), nPMombins, PRecMombins, nPThetabins, PThetabins);
    m_hpred_p_ptheta[i]->SetDirectory(0);

    if(m_hdata_p_ptheta[i] != nullptr)
        delete m_hdata_p_ptheta[i];
    m_hdata_p_ptheta[i] = new TH2D(Form("%s_data_p_ptheta_bin%d", m_name.c_str(), i),
                       Form("%s_data_p_ptheta_bin%d", m_name.c_str(), i), nPMombins, PRecMombins, nPThetabins, PThetabins);
    m_hdata_p_ptheta[i]->SetDirectory(0);
  }
    if(m_hpred_pi_ptheta != nullptr)
        delete m_hpred_pi_ptheta;
    m_hpred_pi_ptheta = new TH2D(Form("%s_pred_pi_ptheta", m_name.c_str()),
                       Form("%s_pred_pi_ptheta", m_name.c_str()), nPiMombins, PiMombins, nPThetabins, PThetabins);
    m_hpred_pi_ptheta->SetDirectory(0);

    if(m_hdata_pi_ptheta != nullptr)
        delete m_hdata_pi_ptheta;
    m_hdata_pi_ptheta = new TH2D(Form("%s_data_pi_ptheta", m_name.c_str()),
                       Form("%s_data_pi_ptheta", m_name.c_str()), nPiMombins, PiMombins, nPThetabins, PThetabins);
    m_hdata_pi_ptheta->SetDirectory(0);

    if(m_hpred_2d != nullptr)
        delete m_hpred_2d;
    m_hpred_2d = new TH2D(Form("%s_pred_2d", m_name.c_str()),
                       Form("%s_pred_2d", m_name.c_str()), m_nbins, 0, m_nbins, m_nDetBins, 0, m_nDetBins);
    m_hpred_2d->SetDirectory(0);

    if(m_hdata_2d != nullptr)
        delete m_hdata_2d;
    m_hdata_2d = new TH2D(Form("%s_data_2d", m_name.c_str()),
                       Form("%s_data_2d", m_name.c_str()), m_nbins, 0, m_nbins, m_nDetBins, 0, m_nDetBins);
    m_hdata_2d->SetDirectory(0);

    std::cout << TAG << m_nbins << " bins inside MakeHistos()." << std::endl;
}

void AnaSample::SetData(TObject* hdata)
{
    if(m_hdata != nullptr)
        delete m_hdata;
    m_hdata = (TH1D*)hdata->Clone(Form("%s_data", m_name.c_str()));
    m_hdata->SetDirectory(0);
}

int AnaSample::GetBinIndex(const double D1, const double D2) const
{
    for(int i = 0; i < m_bin_edges.size(); ++i)
    {
        if(D1 >= m_bin_edges[i].D1low && D1 < m_bin_edges[i].D1high && D2 >= m_bin_edges[i].D2low
           && D2 < m_bin_edges[i].D2high)
        {
            return i;
        }
    }
    return -1;
}

void AnaSample::FillEventHist(int datatype, bool stat_fluc)
{
    if(m_hpred != nullptr)
        m_hpred->Reset();
    if(m_hmc != nullptr)
        m_hmc->Reset();
    if(m_hmc_true != nullptr)
        m_hmc_true->Reset();
    if(m_hsig != nullptr)
        m_hsig->Reset();
    if(m_hpred_pmom != nullptr)
        m_hpred_pmom->Reset();
    if(m_hpred_ptheta != nullptr)
        m_hpred_ptheta->Reset();
    if(m_hpred_pimom != nullptr)
        m_hpred_pimom->Reset();
    if(m_hpred_pitheta != nullptr)
        m_hpred_pitheta->Reset();

    if(m_hpred_2d != nullptr)
        m_hpred_2d->Reset();

  for (int i=0;i<m_nbins;i++)
    if(m_hpred_p_ptheta[i] != nullptr)
        m_hpred_p_ptheta[i]->Reset();


    if(m_hpred_pi_ptheta != nullptr)
        m_hpred_pi_ptheta->Reset();

    for(std::size_t i = 0; i < m_events.size(); ++i)
    {
        double D1_rec  = m_events[i].GetRecoD1();
        double D2_rec  = m_events[i].GetRecoD2();
        double D1_true = m_events[i].GetTrueD1();
        double D2_true = m_events[i].GetTrueD2();
        double wght    = datatype >= 0 ? m_events[i].GetEvWght() : m_events[i].GetEvWghtMC();

        int anybin_index_rec  = GetBinIndex(D1_rec, D2_rec);
        int anybin_index_true = GetBinIndex(D1_true, D2_true);

        m_hpred->Fill(anybin_index_rec + 0.5, wght);
        m_hmc->Fill(anybin_index_rec + 0.5, wght);
        m_hmc_true->Fill(anybin_index_true + 0.5, wght);

        if(m_events[i].isSignalEvent())
            m_hsig->Fill(anybin_index_true + 0.5, wght);

        /*double pmom  = m_events[i].GetPrMom();
        double ptheta  = m_events[i].GetPrCostheta();
        m_hpred_pmom->Fill(pmom,wght);
        m_hpred_ptheta->Fill(ptheta,wght);
        double pimom  = m_events[i].GetPiMom();
        double pitheta  = m_events[i].GetPiCostheta();
        m_hpred_pimom->Fill(pimom,wght);
        m_hpred_pitheta->Fill(pitheta,wght);

        m_hpred_p_ptheta[anybin_index_rec]->Fill(pmom,ptheta,wght);
        m_hpred_pi_ptheta->Fill(pimom,pitheta,wght);

        int anybin_index_binmanager = bin_manager.GetBinIndex(std::vector<double>{m_events[i].GetMuMom(),m_events[i].GetMuCostheta(),pimom,pitheta,pmom,ptheta});
        m_hpred_2d->Fill(anybin_index_rec,anybin_index_binmanager,wght);*/
    }

    m_hpred->Scale(m_norm);
    m_hmc->Scale(m_norm);
    m_hsig->Scale(m_norm);

    m_hpred_pmom->Scale(m_norm);
    m_hpred_ptheta->Scale(m_norm);
    m_hpred_pimom->Scale(m_norm);
    m_hpred_pitheta->Scale(m_norm);

    m_hpred_2d->Scale(m_norm);

    for (int i=0;i<m_nbins;i++)	m_hpred_p_ptheta[i]->Scale(m_norm);

    m_hpred_pi_ptheta->Scale(m_norm);

    if(datatype == 0 || datatype == -1)
        return;

    else if(datatype == 1)
    {
        SetData(m_hpred);
        m_hdata->Reset();
        m_hdata_pmom->Reset();
        m_hdata_ptheta->Reset();
        m_hdata_pimom->Reset();
        m_hdata_pitheta->Reset();

      for (int i=0;i<m_nbins;i++)
        m_hdata_p_ptheta[i]->Reset();

        m_hdata_pi_ptheta->Reset();

        m_hdata_2d->Reset();

        if(stat_fluc)
            std::cout << TAG << "Applying statistical fluctuations..." << std::endl;

        for(int j = 1; j <= m_hpred->GetNbinsX(); ++j)
        {
            double val = m_hpred->GetBinContent(j);
            if(stat_fluc)
                val = gRandom->Poisson(val);
#ifndef NDEBUG
            if(val == 0.0)
            {
                std::cout << "[WARNING] In AnaSample::FillEventHist()\n"
                          << "[WARNING] " << m_name << " bin " << j
                          << " has 0 entries. This may cause a problem with chi2 computations."
                          << std::endl;
                continue;
            }
#endif
            m_hdata->SetBinContent(j, val);
        }

/*        for(int j = 1; j <= m_hpred_pmom->GetNbinsX(); ++j)
        {
            double val = m_hpred_pmom->GetBinContent(j);
            if(stat_fluc)
                val = gRandom->Poisson(val);
#ifndef NDEBUG
            if(val == 0.0)
            {
                std::cout << "[WARNING] In AnaSample::FillEventHist()\n"
                          << "[WARNING] " << "Proton_mom" << " bin " << j
                          << " has 0 entries. This may cause a problem with chi2 computations."
                          << std::endl;
                continue;
            }
#endif
            m_hdata_pmom->SetBinContent(j, val);
        }

        for(int j = 1; j <= m_hpred_ptheta->GetNbinsX(); ++j)
        {
            double val = m_hpred_ptheta->GetBinContent(j);
            if(stat_fluc)
                val = gRandom->Poisson(val);
#ifndef NDEBUG
            if(val == 0.0)
            {
                std::cout << "[WARNING] In AnaSample::FillEventHist()\n"
                          << "[WARNING] " << "Proton_theta" << " bin " << j
                          << " has 0 entries. This may cause a problem with chi2 computations."
                          << std::endl;
                continue;
            }
#endif
            m_hdata_ptheta->SetBinContent(j, val);
        }

        for(int j = 1; j <= m_hpred_pimom->GetNbinsX(); ++j)
        {
            double val = m_hpred_pimom->GetBinContent(j);
            if(stat_fluc)
                val = gRandom->Poisson(val);
#ifndef NDEBUG
            if(val == 0.0)
            {
                std::cout << "[WARNING] In AnaSample::FillEventHist()\n"
                          << "[WARNING] " << "Pi_mom" << " bin " << j
                          << " has 0 entries. This may cause a problem with chi2 computations."
                          << std::endl;
                continue;
            }
#endif
            m_hdata_pimom->SetBinContent(j, val);
        }

        for(int j = 1; j <= m_hpred_pitheta->GetNbinsX(); ++j)
        {
            double val = m_hpred_pitheta->GetBinContent(j);
            if(stat_fluc)
                val = gRandom->Poisson(val);
#ifndef NDEBUG
            if(val == 0.0)
            {
                std::cout << "[WARNING] In AnaSample::FillEventHist()\n"
                          << "[WARNING] " << "Pi_theta" << " bin " << j
                          << " has 0 entries. This may cause a problem with chi2 computations."
                          << std::endl;
                continue;
            }
#endif
            m_hdata_pitheta->SetBinContent(j, val);
        }

      for (int i=0;i<m_nbins;i++)
        for(int j = 1; j <= m_hpred_p_ptheta[i]->GetNbinsX(); ++j) for(int k = 1; k <= m_hpred_p_ptheta[i]->GetNbinsY(); ++k)
            m_hdata_p_ptheta[i]->SetBinContent(j, k, m_hpred_p_ptheta[i]->GetBinContent(j,k));

        for(int j = 1; j <= m_hpred_pi_ptheta->GetNbinsX(); ++j) for(int k = 1; k <= m_hpred_pi_ptheta->GetNbinsY(); ++k)
            m_hdata_pi_ptheta->SetBinContent(j, k, m_hpred_pi_ptheta->GetBinContent(j,k));

        for(int j = 1; j <= m_hpred_2d->GetNbinsX(); ++j) for(int k = 1; k <= m_hpred_2d->GetNbinsY(); ++k)
            m_hdata_2d->SetBinContent(j, k, m_hpred_2d->GetBinContent(j,k));*/

      for(std::size_t i = 0; i < m_events.size(); ++i)
      {
        double D1_rec  = m_events[i].GetRecoD1();
        double D2_rec  = m_events[i].GetRecoD2();
        double D1_true = m_events[i].GetTrueD1();
        double D2_true = m_events[i].GetTrueD2();
        double wght    = datatype >= 0 ? m_events[i].GetEvWght() : m_events[i].GetEvWghtMC();

        int anybin_index_rec  = GetBinIndex(D1_rec, D2_rec);
        int anybin_index_true = GetBinIndex(D1_true, D2_true);

        double pmom  = m_events[i].GetPrMom();
        double ptheta  = m_events[i].GetPrCostheta();
        double pimom  = m_events[i].GetPiMom();
        double pitheta  = m_events[i].GetPiCostheta();

        int anybin_index_binmanager = bin_manager.GetBinIndex(std::vector<double>{m_events[i].GetMuMom(),m_events[i].GetMuCostheta(),pimom,pitheta,pmom,ptheta});
        m_hdata_2d->Fill(anybin_index_rec,anybin_index_binmanager,wght);
      }

    }

    else if(datatype == 2 || datatype == 3)
    {
        SetData(m_hpred);
        m_hdata->Reset();
        m_hdata_pmom->Reset();
        m_hdata_ptheta->Reset();
        m_hdata_pimom->Reset();
        m_hdata_pitheta->Reset();
        m_hdata_2d->Reset();

      for (int i=0;i<m_nbins;i++)
        m_hdata_p_ptheta[i]->Reset();

        m_hdata_pi_ptheta->Reset();

        float D1_rec_tree, D2_rec_tree, wght;
        int cut_branch;

        m_data_tree->SetBranchAddress("cutBranch", &cut_branch);
        m_data_tree->SetBranchAddress("weight", &wght);
        m_data_tree->SetBranchAddress("D1Rec", &D1_rec_tree);
        m_data_tree->SetBranchAddress("D2Rec", &D2_rec_tree);

        float pMomRec, pCosThetaRec;
        m_data_tree->SetBranchAddress("pMomRec", &pMomRec);
        m_data_tree->SetBranchAddress("pCosThetaRec", &pCosThetaRec);

        float piMomRec, piCosThetaRec;
        m_data_tree->SetBranchAddress("piMomRec", &piMomRec);
        m_data_tree->SetBranchAddress("piCosThetaRec", &piCosThetaRec);

        float muMomRec, muCosThetaRec;
        m_data_tree->SetBranchAddress("muMomRec", &muMomRec);
        m_data_tree->SetBranchAddress("muCosThetaRec", &muCosThetaRec);

        long int n_entries = m_data_tree->GetEntries();
        for(std::size_t i = 0; i < n_entries; ++i)
        {
            m_data_tree->GetEntry(i);
            if(cut_branch != m_sample_id)
                continue;

            int anybin_index = GetBinIndex(D1_rec_tree, D2_rec_tree);
            if(anybin_index != -1)
            {
                m_hdata->Fill(anybin_index + 0.5, wght);
            }
#ifndef NDEBUG
            else
            {
                std::cout << "[WARNING] In AnaSample::FillEventHist()\n"
                          << "[WARNING] No bin for current data event.\n"
                          << "[WARNING] D1 Reco: " << D1_rec_tree << std::endl
                          << "[WARNING] D2 Reco: " << D2_rec_tree << std::endl;
            }
#endif
            m_hdata_pmom->Fill(pMomRec,wght);
            m_hdata_ptheta->Fill(pCosThetaRec,wght);
            m_hdata_pimom->Fill(piMomRec,wght);
            m_hdata_pitheta->Fill(piCosThetaRec,wght);
            if(anybin_index != -1) m_hdata_p_ptheta[anybin_index]->Fill(pMomRec,pCosThetaRec,wght);
            m_hdata_pi_ptheta->Fill(piMomRec,piCosThetaRec,wght);

            int anybin_index_binmanager = bin_manager.GetBinIndex(std::vector<double>{muMomRec,muCosThetaRec,piMomRec,piCosThetaRec,pMomRec,pCosThetaRec});
            m_hdata_2d->Fill(anybin_index,anybin_index_binmanager,wght);
        }

        if(stat_fluc && datatype == 2)
        {
            if(stat_fluc)
                std::cout << TAG << "Applying statistical fluctuations..." << std::endl;

            for(unsigned int i = 1; i <= m_hdata->GetNbinsX(); ++i)
            {
                double val = gRandom->Poisson(m_hdata->GetBinContent(i));
                m_hdata->SetBinContent(i, val);
            }
        }

#ifndef NDEBUG
        std::cout << TAG << "Data histogram filled: " << std::endl;
        m_hdata->Print();
#endif
    }

    else
    {
        std::cout << "[WARNING]: In AnaSample::FillEventHist()\n"
                  << "[WARNING]: Invalid data type to fill histograms!\n";
    }
}

double AnaSample::GetProtonMomChi2()
{
/*
    if(m_hdata_pmom == nullptr)
    {
        std::cerr << "[ERROR]: In AnaSample::GetProtonMomChi2()\n"
                  << "[ERROR]: Need to define pmom data histogram." << std::endl;
        return 0.0;
    }
    if(m_hdata_ptheta == nullptr)
    {
        std::cerr << "[ERROR]: In AnaSample::GetProtonMomChi2()\n"
                  << "[ERROR]: Need to define ptheta data histogram." << std::endl;
        return 0.0;
    }

    int nbins = m_hpred_pmom->GetNbinsX();
    if(nbins != m_hdata_pmom->GetNbinsX())
    {
        std::cerr << "[ERROR]: In AnaSample::GetProtonMomChi2()\n"
                  << "[ERROR]: pmom Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbins << ", Data bins: " << m_hdata_pmom->GetNbinsX()
                  << std::endl;
        return 0.0;
    }

    double chi2 = 0.0;
    for(int j = 1; j <= nbins; ++j)
    {
        double obs = m_hdata_pmom->GetBinContent(j);
        double exp = m_hpred_pmom->GetBinContent(j);
        if(exp > 0.0)
        {
            // added when external fake datasets (you cannot reweight when simply 0)
            // this didn't happen when all from same MC since if exp=0 then obs =0

            chi2 += 2 * (exp - obs);
            if(obs > 0.0)
                chi2 += 2 * obs * TMath::Log(obs / exp);

            if(chi2 < 0.0)
            {
#ifndef NDEBUG
                std::cerr << "[WARNING]: In AnaSample::GetProtonMomChi2()\n"
                          << "[WARNING]: Stat chi2 is less than 0: " << chi2 << ", setting to 0."
                          << std::endl;
                std::cerr << "[WARNING]: exp and obs is: " << exp << " and " << obs << "."
                          << std::endl;
#endif
                chi2 = 0.0;
            }
        }
    }

    nbins = m_hpred_ptheta->GetNbinsX();
    if(nbins != m_hdata_ptheta->GetNbinsX())
    {
        std::cerr << "[ERROR]: In AnaSample::GetProtonMomChi2()\n"
                  << "[ERROR]: ptheta Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbins << ", Data bins: " << m_hdata_ptheta->GetNbinsX()
                  << std::endl;
        return 0.0;
    }

    for(int j = 1; j <= nbins; ++j)
    {
        double obs = m_hdata_ptheta->GetBinContent(j);
        double exp = m_hpred_ptheta->GetBinContent(j);
        if(exp > 0.0)
        {
            // added when external fake datasets (you cannot reweight when simply 0)
            // this didn't happen when all from same MC since if exp=0 then obs =0

            chi2 += 2 * (exp - obs);
            if(obs > 0.0)
                chi2 += 2 * obs * TMath::Log(obs / exp);

            if(chi2 < 0.0)
            {
#ifndef NDEBUG
                std::cerr << "[WARNING]: In AnaSample::GetProtonMomChi2()\n"
                          << "[WARNING]: Stat chi2 is less than 0: " << chi2 << ", setting to 0."
                          << std::endl;
                std::cerr << "[WARNING]: exp and obs is: " << exp << " and " << obs << "."
                          << std::endl;
#endif
                chi2 = 0.0;
            }
        }
    }
*/
    double chi2 = 0.0;

  for (int i=0;i<m_nbins-1;i++) {
    if(m_hdata_p_ptheta[i] == nullptr)
    {
        std::cerr << "[ERROR]: In AnaSample::GetProtonMomChi2()\n"
                  << "[ERROR]: Need to define p_ptheta data histogram." << std::endl;
        return 0.0;
    }

    int nbinsX = m_hpred_p_ptheta[i]->GetNbinsX();
    if(nbinsX != m_hdata_p_ptheta[i]->GetNbinsX())
    {
        std::cerr << "[ERROR]: In AnaSample::GetProtonMomChi2()\n"
                  << "[ERROR]: p mom Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbinsX << ", Data bins: " << m_hdata_p_ptheta[i]->GetNbinsX()
                  << std::endl;
        return 0.0;
    }

    int nbinsY = m_hpred_p_ptheta[i]->GetNbinsY();
    if(nbinsY != m_hdata_p_ptheta[i]->GetNbinsY())
    {
        std::cerr << "[ERROR]: In AnaSample::GetProtonMomChi2()\n"
                  << "[ERROR]: p theta Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbinsY << ", Data bins: " << m_hdata_p_ptheta[i]->GetNbinsY()
                  << std::endl;
        return 0.0;
    }

    for(int j = 1; j <= nbinsX; ++j) for(int k = 1; k <= nbinsY; ++k)
    {
        double obs = m_hdata_p_ptheta[i]->GetBinContent(j,k);
        double exp = m_hpred_p_ptheta[i]->GetBinContent(j,k);
        if(exp > 0.0)
        {
            // added when external fake datasets (you cannot reweight when simply 0)
            // this didn't happen when all from same MC since if exp=0 then obs =0

            chi2 += 2 * (exp - obs);
            if(obs > 0.0)
                chi2 += 2 * obs * TMath::Log(obs / exp);

            if(chi2 < 0.0)
            {
#ifndef NDEBUG
                std::cerr << "[WARNING]: In AnaSample::GetProtonMomChi2()\n"
                          << "[WARNING]: Stat chi2 is less than 0: " << chi2 << ", setting to 0."
                          << std::endl;
                std::cerr << "[WARNING]: exp and obs is: " << exp << " and " << obs << "."
                          << std::endl;
#endif
                chi2 = 0.0;
            }
        }
    }
  }

    double OOPS_obs = m_hdata->GetBinContent(m_nbins);
    double OOPS_exp = m_hpred->GetBinContent(m_nbins);

        if(OOPS_exp > 0.0)
        {
            // added when external fake datasets (you cannot reweight when simply 0)
            // this didn't happen when all from same MC since if exp=0 then obs =0

            chi2 += 2 * (OOPS_exp - OOPS_obs);
            if(OOPS_obs > 0.0)
                chi2 += 2 * OOPS_obs * TMath::Log(OOPS_obs / OOPS_exp);

            if(chi2 < 0.0)
            {
#ifndef NDEBUG
                std::cerr << "[WARNING]: In AnaSample::GetProtonMomChi2()\n"
                          << "[WARNING]: Stat chi2 is less than 0: " << chi2 << ", setting to 0."
                          << std::endl;
                std::cerr << "[WARNING]: exp and obs is: " << OOPS_exp << " and " << OOPS_obs << "."
                          << std::endl;
#endif
                chi2 = 0.0;
            }
        }

    if(chi2 != chi2)
    {
        std::cerr << "[WARNING]: In AnaSample::GetProtonMomChi2()\n"
                  << "[WARNING]: Stat chi2 is nan, setting to 0." << std::endl;
        chi2 = 0.0;
    }

    return chi2;

}

double AnaSample::GetPionMomChi2()
{
/*
    if(m_hdata_pimom == nullptr)
    {
        std::cerr << "[ERROR]: In AnaSample::GetPionMomChi2()\n"
                  << "[ERROR]: Need to define pimom data histogram." << std::endl;
        return 0.0;
    }
    if(m_hdata_pitheta == nullptr)
    {
        std::cerr << "[ERROR]: In AnaSample::GetPionMomChi2()\n"
                  << "[ERROR]: Need to define pitheta data histogram." << std::endl;
        return 0.0;
    }

    int nbins = m_hpred_pimom->GetNbinsX();
    if(nbins != m_hdata_pimom->GetNbinsX())
    {
        std::cerr << "[ERROR]: In AnaSample::GetPionMomChi2()\n"
                  << "[ERROR]: pimom Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbins << ", Data bins: " << m_hdata_pimom->GetNbinsX()
                  << std::endl;
        return 0.0;
    }

    double chi2 = 0.0;
    for(int j = 1; j <= nbins; ++j)
    {
        double obs = m_hdata_pimom->GetBinContent(j);
        double exp = m_hpred_pimom->GetBinContent(j);
        if(exp > 0.0)
        {
            // added when external fake datasets (you cannot reweight when simply 0)
            // this didn't happen when all from same MC since if exp=0 then obs =0

            chi2 += 2 * (exp - obs);
            if(obs > 0.0)
                chi2 += 2 * obs * TMath::Log(obs / exp);

            if(chi2 < 0.0)
            {
#ifndef NDEBUG
                std::cerr << "[WARNING]: In AnaSample::GetPionMomChi2()\n"
                          << "[WARNING]: Stat chi2 is less than 0: " << chi2 << ", setting to 0."
                          << std::endl;
                std::cerr << "[WARNING]: exp and obs is: " << exp << " and " << obs << "."
                          << std::endl;
#endif
                chi2 = 0.0;
            }
        }
    }

    nbins = m_hpred_pitheta->GetNbinsX();
    if(nbins != m_hdata_pitheta->GetNbinsX())
    {
        std::cerr << "[ERROR]: In AnaSample::GetPionMomChi2()\n"
                  << "[ERROR]: pitheta Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbins << ", Data bins: " << m_hdata_pitheta->GetNbinsX()
                  << std::endl;
        return 0.0;
    }

    for(int j = 1; j <= nbins; ++j)
    {
        double obs = m_hdata_pitheta->GetBinContent(j);
        double exp = m_hpred_pitheta->GetBinContent(j);
        if(exp > 0.0)
        {
            // added when external fake datasets (you cannot reweight when simply 0)
            // this didn't happen when all from same MC since if exp=0 then obs =0

            chi2 += 2 * (exp - obs);
            if(obs > 0.0)
                chi2 += 2 * obs * TMath::Log(obs / exp);

            if(chi2 < 0.0)
            {
#ifndef NDEBUG
                std::cerr << "[WARNING]: In AnaSample::GetPionMomChi2()\n"
                          << "[WARNING]: Stat chi2 is less than 0: " << chi2 << ", setting to 0."
                          << std::endl;
                std::cerr << "[WARNING]: exp and obs is: " << exp << " and " << obs << "."
                          << std::endl;
#endif
                chi2 = 0.0;
            }
        }
    }
*/

    if(m_hdata_pi_ptheta == nullptr)
    {
        std::cerr << "[ERROR]: In AnaSample::GetPionMomChi2()\n"
                  << "[ERROR]: Need to define pi_ptheta data histogram." << std::endl;
        return 0.0;
    }

    int nbinsX = m_hpred_pi_ptheta->GetNbinsX();
    if(nbinsX != m_hdata_pi_ptheta->GetNbinsX())
    {
        std::cerr << "[ERROR]: In AnaSample::GetPionMomChi2()\n"
                  << "[ERROR]: pi mom Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbinsX << ", Data bins: " << m_hdata_pi_ptheta->GetNbinsX()
                  << std::endl;
        return 0.0;
    }

    int nbinsY = m_hpred_pi_ptheta->GetNbinsY();
    if(nbinsY != m_hdata_pi_ptheta->GetNbinsY())
    {
        std::cerr << "[ERROR]: In AnaSample::GetPionMomChi2()\n"
                  << "[ERROR]: pi theta Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbinsY << ", Data bins: " << m_hdata_pi_ptheta->GetNbinsY()
                  << std::endl;
        return 0.0;
    }

    double chi2 = 0.0;
    for(int j = 1; j <= nbinsX; ++j) for(int k = 1; k <= nbinsY; ++k)
    {
        double obs = m_hdata_pi_ptheta->GetBinContent(j,k);
        double exp = m_hpred_pi_ptheta->GetBinContent(j,k);
        if(exp > 0.0)
        {
            // added when external fake datasets (you cannot reweight when simply 0)
            // this didn't happen when all from same MC since if exp=0 then obs =0

            chi2 += 2 * (exp - obs);
            if(obs > 0.0)
                chi2 += 2 * obs * TMath::Log(obs / exp);

            if(chi2 < 0.0)
            {
#ifndef NDEBUG
                std::cerr << "[WARNING]: In AnaSample::GetProtonMomChi2()\n"
                          << "[WARNING]: Stat chi2 is less than 0: " << chi2 << ", setting to 0."
                          << std::endl;
                std::cerr << "[WARNING]: exp and obs is: " << exp << " and " << obs << "."
                          << std::endl;
#endif
                chi2 = 0.0;
            }
        }
    }


    if(chi2 != chi2)
    {
        std::cerr << "[WARNING]: In AnaSample::GetPionMomChi2()\n"
                  << "[WARNING]: Stat chi2 is nan, setting to 0." << std::endl;
        chi2 = 0.0;
    }

    return chi2;

}

double AnaSample::Calc2DChi2()
{
    if(m_hdata_2d == nullptr)
    {
        std::cerr << "[ERROR]: In AnaSample::Calc2DChi2()\n"
                  << "[ERROR]: Need to define data histogram." << std::endl;
        return 0.0;
    }

    if(m_hpred_2d != nullptr)
        m_hpred_2d->Reset();

      for(std::size_t i = 0; i < m_events.size(); ++i)
      {
        double D1_rec  = m_events[i].GetRecoD1();
        double D2_rec  = m_events[i].GetRecoD2();
        double D1_true = m_events[i].GetTrueD1();
        double D2_true = m_events[i].GetTrueD2();
        double wght    = m_events[i].GetEvWght();

        int anybin_index_rec  = GetBinIndex(D1_rec, D2_rec);
        int anybin_index_true = GetBinIndex(D1_true, D2_true);

        double pmom  = m_events[i].GetPrMom();
        double ptheta  = m_events[i].GetPrCostheta();
        double pimom  = m_events[i].GetPiMom();
        double pitheta  = m_events[i].GetPiCostheta();

        int anybin_index_binmanager = bin_manager.GetBinIndex(std::vector<double>{m_events[i].GetMuMom(),m_events[i].GetMuCostheta(),pimom,pitheta,pmom,ptheta});
        m_hpred_2d->Fill(anybin_index_rec,anybin_index_binmanager,wght);
      }

    m_hpred_2d->Scale(m_norm);

    int nbinsX = m_hpred_2d->GetNbinsX();
    if(nbinsX != m_hdata_2d->GetNbinsX())
    {
        std::cerr << "[ERROR]: In AnaSample::GetPionMomChi2()\n"
                  << "[ERROR]: pi mom Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbinsX << ", Data bins: " << m_hdata_2d->GetNbinsX()
                  << std::endl;
        return 0.0;
    }

    int nbinsY = m_hpred_2d->GetNbinsY();
    if(nbinsY != m_hdata_2d->GetNbinsY())
    {
        std::cerr << "[ERROR]: In AnaSample::GetPionMomChi2()\n"
                  << "[ERROR]: pi theta Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbinsY << ", Data bins: " << m_hdata_2d->GetNbinsY()
                  << std::endl;
        return 0.0;
    }

    double chi2 = 0.0;
    for(int j = 1; j <= nbinsX-1; ++j) for(int k = 1; k <= nbinsY; ++k)
    {
        double obs = m_hdata_2d->GetBinContent(j,k);
        double exp = m_hpred_2d->GetBinContent(j,k);
        if(exp > 0.0)
        {
            // added when external fake datasets (you cannot reweight when simply 0)
            // this didn't happen when all from same MC since if exp=0 then obs =0

            chi2 += 2 * (exp - obs);
            if(obs > 0.0)
                chi2 += 2 * obs * TMath::Log(obs / exp);

            if(chi2 < 0.0)
            {
#ifndef NDEBUG
                std::cerr << "[WARNING]: In AnaSample::CalcChi2()\n"
                          << "[WARNING]: Stat chi2 is less than 0: " << chi2 << ", setting to 0."
                          << std::endl;
                std::cerr << "[WARNING]: exp and obs is: " << exp << " and " << obs << "."
                          << std::endl;
#endif
                chi2 = 0.0;
            }
        }
    }

    double OOPS_obs = m_hdata->GetBinContent(m_nbins);
    double OOPS_exp = m_hpred->GetBinContent(m_nbins);

        if(OOPS_exp > 0.0)
        {
            // added when external fake datasets (you cannot reweight when simply 0)
            // this didn't happen when all from same MC since if exp=0 then obs =0

            chi2 += 2 * (OOPS_exp - OOPS_obs);
            if(OOPS_obs > 0.0)
                chi2 += 2 * OOPS_obs * TMath::Log(OOPS_obs / OOPS_exp);

            if(chi2 < 0.0)
            {
#ifndef NDEBUG
                std::cerr << "[WARNING]: In AnaSample::GetProtonMomChi2()\n"
                          << "[WARNING]: Stat chi2 is less than 0: " << chi2 << ", setting to 0."
                          << std::endl;
                std::cerr << "[WARNING]: exp and obs is: " << OOPS_exp << " and " << OOPS_obs << "."
                          << std::endl;
#endif
                chi2 = 0.0;
            }
        }

    if(chi2 != chi2)
    {
        std::cerr << "[WARNING]: In AnaSample::CalcChi2()\n"
                  << "[WARNING]: Stat chi2 is nan, setting to 0." << std::endl;
        chi2 = 0.0;
    }

    return chi2;
}


double AnaSample::CalcChi2() const
{
    if(m_hdata == nullptr)
    {
        std::cerr << "[ERROR]: In AnaSample::CalcChi2()\n"
                  << "[ERROR]: Need to define data histogram." << std::endl;
        return 0.0;
    }

    int nbins = m_hpred->GetNbinsX();
    if(nbins != m_hdata->GetNbinsX())
    {
        std::cerr << "[ERROR]: In AnaSample::CalcChi2()\n"
                  << "[ERROR]: Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbins << ", Data bins: " << m_hdata->GetNbinsX()
                  << std::endl;
        return 0.0;
    }

    double* exp_w  = m_hpred->GetArray();
    double* exp_w2 = m_hpred->GetSumw2()->GetArray();
    double* data   = m_hdata->GetArray();

    double chi2 = 0.0;
    for(int j = 1; j <= nbins; ++j)
    {
        if (m_llh_func==0) chi2 += PoissonLLH(exp_w[j], exp_w2[j], data[j]);
        else if (m_llh_func==1) chi2 += BarlowLLH(exp_w[j], exp_w2[j], data[j]);
        else if (m_llh_func==2) chi2 += EffLLH(exp_w[j], exp_w2[j], data[j]);
        else std::cerr << "[ERROR]: In AnaSample::CalcChi2()\n"
                       << "[ERROR]: Undefined likelihood function: m_llh_func = "<< m_llh_func << std::endl;

        /*double obs = m_hdata->GetBinContent(j);
        double exp = m_hpred->GetBinContent(j);
        if(exp > 0.0)
        {
            // added when external fake datasets (you cannot reweight when simply 0)
            // this didn't happen when all from same MC since if exp=0 then obs =0

            chi2 += 2 * (exp - obs);
            if(obs > 0.0)
                chi2 += 2 * obs * TMath::Log(obs / exp);

            if(chi2 < 0.0)
            {
#ifndef NDEBUG
                std::cerr << "[WARNING]: In AnaSample::CalcChi2()\n"
                          << "[WARNING]: Stat chi2 is less than 0: " << chi2 << ", setting to 0."
                          << std::endl;
                std::cerr << "[WARNING]: exp and obs is: " << exp << " and " << obs << "."
                          << std::endl;
#endif
                chi2 = 0.0;
            }
        }*/
    }

    if(chi2 != chi2)
    {
        std::cerr << "[WARNING]: In AnaSample::CalcChi2()\n"
                  << "[WARNING]: Stat chi2 is nan, setting to 0." << std::endl;
        chi2 = 0.0;
    }

    return chi2;
}

void AnaSample::Write(TDirectory* dirout, const std::string& bsname, int fititer)
{
    dirout->cd();
    m_hpred->Write(Form("%s_pred", bsname.c_str()));
    m_hmc_true->Write(Form("%s_true", bsname.c_str()));

    m_hpred_pmom->Write(Form("%s_pmom_pred", bsname.c_str()));
    m_hdata_pmom->Write(Form("%s_pmom_data", bsname.c_str()));
    m_hpred_ptheta->Write(Form("%s_ptheta_pred", bsname.c_str()));
    m_hdata_ptheta->Write(Form("%s_ptheta_data", bsname.c_str()));

  for (int i=0;i<m_nbins;i++) {
    m_hpred_p_ptheta[i]->Write(Form("%s_p_ptheta_pred_bin%d", bsname.c_str(), i));
    m_hdata_p_ptheta[i]->Write(Form("%s_p_ptheta_data_bin%d", bsname.c_str(), i));
  }

    m_hpred_pi_ptheta->Write(Form("%s_pi_ptheta_pred", bsname.c_str()));
    m_hdata_pi_ptheta->Write(Form("%s_pi_ptheta_data", bsname.c_str()));

    m_hpred_2d->Write(Form("%s_pred_2d", bsname.c_str()));
    m_hdata_2d->Write(Form("%s_data_2d", bsname.c_str()));

    if(fititer == 0)
    {
        m_hmc->Write(Form("%s_mc", bsname.c_str()));
        if(m_hdata != nullptr)
            m_hdata->Write(Form("%s_data", bsname.c_str()));
    }
}

void AnaSample::GetSampleBreakdown(TDirectory* dirout, const std::string& tag,
                                   const std::vector<std::string>& topology, bool save)
{
    const int ntopology = topology.size();
    int compos[ntopology];
    std::vector<TH1D> hAnybin_rec;
    std::vector<TH1D> hAnybin_true;

    for(int i = 0; i < ntopology; ++i)
    {
        compos[i] = 0;
        hAnybin_rec.emplace_back(
            TH1D(Form("%s_Anybins_rec_%s_%s", m_name.c_str(), topology[i].c_str(), tag.c_str()),
                 Form("%s_Anybins_rec_%s_%s", m_name.c_str(), topology[i].c_str(), tag.c_str()),
                 m_nbins, 0, m_nbins));
        hAnybin_rec[i].SetDirectory(0);
        hAnybin_rec[i].GetXaxis()->SetTitle("Bin Index");

        hAnybin_true.emplace_back(
            TH1D(Form("%s_Anybins_true_%s_%s", m_name.c_str(), topology[i].c_str(), tag.c_str()),
                 Form("%s_Anybins_true_%s_%s", m_name.c_str(), topology[i].c_str(), tag.c_str()),
                 m_nbins, 0, m_nbins));
        hAnybin_true[i].SetDirectory(0);
        hAnybin_true[i].GetXaxis()->SetTitle("Bin Index");
    }

    int Ntot = GetN();
    for(std::size_t i = 0; i < m_events.size(); ++i)
    {
        double D1_rec    = m_events[i].GetRecoD1();
        double D2_rec    = m_events[i].GetRecoD2();
        double D1_true   = m_events[i].GetTrueD1();
        double D2_true   = m_events[i].GetTrueD2();
        double wght      = m_events[i].GetEvWght();
        int evt_topology = m_events[i].GetTopology();

        compos[evt_topology]++;
        int anybin_index_rec  = GetBinIndex(D1_rec, D2_rec);
        int anybin_index_true = GetBinIndex(D1_true, D2_true);
        hAnybin_rec[evt_topology].Fill(anybin_index_rec + 0.5, wght);
        hAnybin_true[evt_topology].Fill(anybin_index_true + 0.5, wght);
    }

    dirout->cd();
    for(int i = 0; i < ntopology; ++i)
    {
        hAnybin_true[i].Scale(m_norm);
        hAnybin_rec[i].Scale(m_norm);

        if(save == true)
        {
            hAnybin_true[i].Write();
            hAnybin_rec[i].Write();
        }
    }

    std::cout << TAG << "GetSampleBreakdown()\n"
              << "============ Sample " << m_name << " ===========" << std::endl;

    for(int j = 0; j < ntopology; ++j)
    {
        std::cout << std::setw(10) << topology[j] << std::setw(5) << j << std::setw(5) << compos[j]
                  << std::setw(10) << ((1.0 * compos[j]) / Ntot) * 100.0 << "%" << std::endl;
    }

    std::cout << std::setw(10) << "Total" << std::setw(5) << " " << std::setw(5) << Ntot
              << std::setw(10) << "100.00%" << std::endl;
}

double AnaSample::PoissonLLH(double mc, double w2, double data) const
{
        // Standard Poisson LLH.
        double chi2 = 0.0;
        if(mc > 0.0)
        {
            chi2 = 2 * (mc - data);
            if(data > 0.0)
                chi2 += 2 * data * TMath::Log(data / mc);
        }

        return (chi2 >= 0.0) ? chi2 : 0.0;
}

double AnaSample::EffLLH(double mc, double w2, double data) const
{
        // Effective LLH based on Tianlu's paper.
        if(mc <= 0.0)
            return 0.0;

        const double b = mc / w2;
        const double a = (mc * b) + 1.0;
        const double k = data;

        return -2 * (a * std::log(b) + std::lgamma(k + a) - std::lgamma(k + 1)
               - ((k + a) * std::log1p(b)) - std::lgamma(a));
}

double AnaSample::BarlowLLH(double mc, double w2, double data) const
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
