#include "DetParameters.hh"
using xsllh::FitBin;

DetParameters::DetParameters(const std::string& name)
{
    m_name = name;
}

DetParameters::~DetParameters() { ; }

bool DetParameters::SetBinning(const std::string& file_name, std::vector<FitBin>& bins)
{
    std::ifstream fin(file_name, std::ios::in);
    if(!fin.is_open())
    {
        std::cerr << ERR << "In DetParameters::SetBinning()\n"
                  << ERR << "Failed to open binning file: " << file_name << std::endl;
        return false;
    }

    else
    {
        std::string line;
        while(getline(fin, line))
        {
            std::stringstream ss(line);
            double D1_1, D1_2, D2_1, D2_2;
            if(!(ss>>D2_1>>D2_2>>D1_1>>D1_2))
            {
                std::cout << WAR << "In DetParameters::SetBinning()\n"
                          << WAR << "Bad line format: " << line << std::endl;
                continue;
            }
            bins.emplace_back(FitBin(D1_1, D1_2, D2_1, D2_2));
        }
        fin.close();

        std::cout << TAG << "Fit binning: \n";
        for(std::size_t i = 0; i < bins.size(); ++i)
        {
            std::cout << std::setw(3) << i
                      << std::setw(5) << bins[i].D2low
                      << std::setw(5) << bins[i].D2high
                      << std::setw(5) << bins[i].D1low
                      << std::setw(5) << bins[i].D1high << std::endl;
        }

        return true;
    }
}

int DetParameters::GetBinIndex(const int sam, double D1, double D2) const
{
    int bin = BADBIN;
    const std::vector<FitBin> &temp_bins = m_sample_bins.at(sam);

    for(int i = 0; i < temp_bins.size(); ++i)
    {
        if(D1 >= temp_bins[i].D1low && D1 < temp_bins[i].D1high &&
           D2 >= temp_bins[i].D2low && D2 < temp_bins[i].D2high)
        {
            bin = i;
            break;
        }
    }
    return bin;
}

void DetParameters::InitEventMap(std::vector<AnaSample*>& sample, int mode)
{
    InitParameters();
    m_evmap.clear();

    if(mode == 2)
        std::cout << TAG << "Not using detector reweighting." << std::endl;

    for(std::size_t s = 0; s < sample.size(); ++s)
    {
        /*if (sample[s]->GetSampleID()==-1) {
          std::cout<< TAG << "Not loading detector reweighting for truth events." << std::endl;continue;
        }*/
        std::vector<int> sample_map;
        for(int i = 0; i < sample[s]->GetN(); ++i)
        {
            AnaEvent* ev = sample[s]->GetEvent(i);
            double D1 = ev->GetRecoD1();
            double D2 = ev->GetRecoD2();
            int bin   = GetBinIndex(sample[s]->GetSampleID(), D1, D2);
#ifndef NDEBUG
            if(bin == BADBIN && sample[s]->GetSampleID()!=-1)
            {
                std::cout << WAR << m_name << ", Event: " << i << std::endl
                          << WAR << "D1 = " << D1 << ", D2 = " << D2 << ", falls outside bin ranges." << std::endl
                          << WAR << "This event will be ignored in the analysis." << std::endl;
            }
#endif
            if (useDetbinning) {
                double muMomRec = ev->GetMuMom();
                double muCosThetaRec = ev->GetMuCostheta();
                double pMomRec = ev->GetPrMom();
                double pCosThetaRec = ev->GetPrCostheta();
                double piMomRec = ev->GetPiMom();
                double piCosThetaRec = ev->GetPiCostheta();
                std::vector<double> vars;
                vars.push_back(muMomRec);vars.push_back(muCosThetaRec);
                vars.push_back(piMomRec);vars.push_back(piCosThetaRec);
                vars.push_back(pMomRec);vars.push_back(pCosThetaRec);
                int idx = m_bin_manager[s].GetBinIndex(vars);
                if (s==0 && bin==m_sample_bins.at(s).size()-1) bin = m_bin_manager[s].GetNbins();
                else bin = idx;
            }
            // If event is signal let the c_i params handle the reweighting:
            if(mode == 1 && ev->isSignalEvent())
                bin = PASSEVENT;
            else if(mode == 2)
                bin = PASSEVENT;
            sample_map.push_back(bin);
        }
        m_evmap.push_back(sample_map);
    }
}

void DetParameters::ReWeight(AnaEvent* event, const std::string& det, int nsample, int nevent, std::vector<double>& params)
{
    if(m_evmap.empty()) // need to build an event map first
    {
        std::cerr << ERR << "In DetParameters::ReWeight()\n"
                  << ERR << "Need to build event map index for " << m_name << std::endl;
        return;
    }

    const int bin = m_evmap[nsample][nevent];

    if(bin == PASSEVENT)
        return;
    if(bin == BADBIN)
        event->AddEvWght(0.0);
    else
    {
        if(bin > params.size())
        {
            std::cout << WAR << "In DetParameters::ReWeight()\n"
                      << WAR << "Number of bins in " << m_name << " does not match num of parameters.\n"
                      << WAR << "Setting event weight to zero." << std::endl;
            event->AddEvWght(0.0);
        }

        event->AddEvWght(params[bin + m_sample_offset.at(event->GetSampleType())]);
    }
}

void DetParameters::InitParameters()
{
    unsigned int offset = 0;
    for(const auto& sam : v_samples)
    {
        m_sample_offset.emplace(std::make_pair(sam, offset));
        const int nbins = useDetbinning ?   (sam==0 ? m_bin_manager.at(sam).GetNbins()+1:m_bin_manager.at(sam).GetNbins())
                                          : m_sample_bins.at(sam).size();
        for(int i = 0; i < nbins; ++i)
        {
            pars_name.push_back(Form("%s_sam%d_%d", m_name.c_str(), sam, i));
            pars_prior.push_back(1.0);
            pars_step.push_back(0.05);
            pars_limlow.push_back(0.0);
            pars_limhigh.push_back(2.0);
            pars_fixed.push_back(false);
        }

        std::cout << TAG << "Total " << nbins << " parameters at "
                  << offset << " for sample ID " << sam << std::endl;
        offset += nbins;
    }

    Npar = pars_name.size();
    pars_original = pars_prior;

    if(m_decompose)
    {
        pars_prior = eigen_decomp -> GetDecompParameters(pars_prior);
        pars_limlow = std::vector<double>(Npar, -100);
        pars_limhigh = std::vector<double>(Npar, 100);

        const int idx = eigen_decomp -> GetInfoFraction(m_info_frac);
        for(int i = idx; i < Npar; ++i)
            pars_fixed[i] = true;

        std::cout << TAG << "Decomposed parameters.\n"
                  << TAG << "Keeping the " << idx << " largest eigen values.\n"
                  << TAG << "Corresponds to " << m_info_frac * 100.0
                  << "\% total variance.\n";
    }
}

void DetParameters::AddDetector(const std::string& det, std::vector<AnaSample*>& v_sample, bool match_bins, const std::string& detbinning)
{
    std::cout << TAG << "Adding detector " << det << " for " << m_name << std::endl;

    for(const auto& sample : v_sample)
    {
        if(sample->GetDetector() != det)
            continue;

        const int sample_id = sample->GetSampleID();
        v_samples.emplace_back(sample_id);

        std::cout << TAG << "Adding sample " << sample->GetName()
                  << " with ID " << sample_id << " to fit." << std::endl;

        if(match_bins)
            m_sample_bins.emplace(std::make_pair(sample_id, sample->GetBinEdges()));
        else
        {
            std::vector<FitBin> temp_vector;
            if(SetBinning(sample->GetDetBinning(), temp_vector))
            {
                m_sample_bins.emplace(std::make_pair(sample_id, temp_vector));
            }
            else
                std::cout << WAR << "Adding sample binning for DetParameters failed." << std::endl;

            if(detbinning!="") m_bin_manager.push_back(BinManager(detbinning));
        }
    }
}
