#include "XsecParameters.hh"

XsecParameters::XsecParameters(const std::string& name)
{
    m_name = name;
    Npar = 0;
}

XsecParameters::~XsecParameters() { ; }

void XsecParameters::InitEventMap(std::vector<AnaSample*>& sample, int mode)
{
    InitParameters();
    if(Npar == 0)
    {
        std::cerr << ERR << "In XsecParameters::InitEventMap\n"
                  << ERR << "No parameters delcared. Not building event map."
                  << std::endl;
    }
    m_dial_evtmap.clear();
    m_norm_evtmap.clear();

    for(std::size_t s = 0; s < sample.size(); ++s)
    {
        std::vector<std::vector<int>> sample_map;
        std::vector<std::vector<int>> sample_map_norm;
        for(int i = 0; i < sample[s] -> GetN(); i++)
        {
            AnaEvent* ev = sample[s] -> GetEvent(i);
            std::vector<int> dial_index_map;

            std::vector<XsecDial> &v_dials = m_dials.at(sample[s] -> GetDetector());
            int num_dials = v_dials.size();

            for(int d = 0; d < num_dials; ++d)
            {
                double q2 = ev -> GetQ2True() / 1.0E6; //MeV to GeV conversion.
                //int idx = v_dials.at(d).GetSplineIndex(ev -> GetTopology(), ev -> GetReaction(), q2);

                //int idx = v_dials.at(d).GetSplineIndex(std::vector<int>{ev -> GetTopology(), ev -> GetReaction()},
                //                                       std::vector<double>{q2});
                //int idx = v_dials.at(d).GetSplineIndex(std::vector<int>{ev->GetTopology(), ev->GetReaction()},
                //                                       std::vector<double>{ev->GetTrueD2(), ev->GetTrueD1()});
                int idx = v_dials.at(d).GetSplineIndex(std::vector<int>{ev->GetTopology(), ev->GetSampleType()},
                                                       std::vector<double>{ev->GetRecoD1(), ev->GetTrueD1()},
                                                       std::vector<double>{ev->GetPrMomTrue(), ev->GetPrMom()},ev->GetTarget());
                if (ev->GetSampleType()==-1) idx=-1;
                if (ev->GetSampleType()==-1) {//truth tree
                   idx = v_dials.at(d).GetSplineIndex(std::vector<int>{ev->GetTopology(), 5},
                                                      std::vector<double>{0.1, ev->GetTrueD1()},
                                                      std::vector<double>{ev->GetPrMomTrue(), 500.},ev->GetTarget());
                }
                if(idx == BADBIN)
                {
                    std::cout << WAR << "Event falls outside spline range.\n"
                              << WAR << "This event will be ignored in the analysis."
                              << std::endl;
                    ev -> AddEvWght(0.0);
                }

                if(mode == 1 && ev -> isSignalEvent())
                    idx = PASSEVENT;

                dial_index_map.push_back(idx);
            }

            sample_map.emplace_back(dial_index_map);

            std::vector<int> norm_index_map;
            std::vector<XsecNorm> &v_norms = m_norms.at(sample[s] -> GetDetector());
            int num_norms = v_norms.size();

            for(int d = 0; d < num_norms; ++d)
            {
                int nDef = v_norms.at(d).GetNDef();
                bool useNorm = true;
                for (int def = 0; def<nDef; def++) {
                  if (!v_norms.at(d).IsApplied(def,ev->GetTopology(),ev->GetEventVarValue(v_norms.at(d).GetDef(def)))) {useNorm=false;break;}
                }

                int idx = useNorm;

                if(mode == 1 && ev -> isSignalEvent())
                    idx = PASSEVENT;

                norm_index_map.push_back(idx);
            }

            sample_map_norm.emplace_back(norm_index_map);
        }

        m_dial_evtmap.emplace_back(sample_map);
        m_norm_evtmap.emplace_back(sample_map_norm);
    }
}

void XsecParameters::InitParameters()
{
    unsigned int offset = 0;
    for(const auto& det : v_detectors)
    {
        m_offset.insert(std::make_pair(det, offset));
        for(const auto& d : m_dials.at(det))
        {
            pars_name.push_back(Form("%s_%s", det.c_str(), d.GetName().c_str()));
            pars_prior.push_back(d.GetNominal());
            pars_step.push_back(d.GetStep());
            pars_limlow.push_back(d.GetLimitLow());
            pars_limhigh.push_back(d.GetLimitHigh());
            pars_fixed.push_back(false);

            std::cout << TAG << "Added " << det << "_" << d.GetName()
                      << std::endl;
        }

        for(const auto& d : m_norms.at(det))
        {
            pars_name.push_back(Form("%s_%s", det.c_str(), d.GetName().c_str()));
            pars_prior.push_back(d.GetNominal());
            pars_step.push_back(d.GetStep());
            pars_limlow.push_back(d.GetLimitLow());
            pars_limhigh.push_back(d.GetLimitHigh());
            pars_fixed.push_back(d.IsFixed());

            std::cout << TAG << "Added " << det << "_" << d.GetName()
                      << std::endl;
        }

        std::cout << TAG << "Total " << m_dials.at(det).size() + m_norms.at(det).size() << " parameters at "
                  << offset << " for " << det << std::endl;

        offset += m_dials.at(det).size()+m_norms.at(det).size();
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

void XsecParameters::ReWeight(AnaEvent* event, const std::string& det, int nsample, int nevent, std::vector<double>& params)
{
    if(m_dial_evtmap.empty()) // need to build an event map first
    {
        std::cerr << ERR << "In XsecParameters::ReWeight()\n"
                  << ERR << "Need to build event map index for " << m_name << std::endl;
        return;
    }

    std::vector<XsecDial> &v_dials = m_dials.at(det);
    int num_dials = v_dials.size();
    double weight = 1.0;

    for(int d = 0; d < num_dials; ++d)
    {
        int idx = m_dial_evtmap[nsample][nevent][d];
        double dial_weight = v_dials[d].GetSplineValue(idx, params[d + m_offset.at(det)]);
        if (dial_weight<0||dial_weight>10) dial_weight=1.;
        weight *= dial_weight;

        /*
        if(dial_weight != 1.0 && det == "INGRID")
        {
            std::cout << "--------------" << std::endl;
            std::cout << "Ev T: " << event -> GetTopology() << std::endl
                      << "Ev R: " << event -> GetReaction() << std::endl
                      << "Ev Q: " << event -> GetQ2() << std::endl;
            std::cout << "Ev I: " << idx << std::endl;
            std::cout << "Ev W: " << dial_weight << std::endl;
            std::cout << "Dl V: " << params[d + m_offset.at(det)] << std::endl;
            std::cout << "Dl N: " << det << "_" << v_dials[d].GetName() << std::endl;
            std::cout << "Sp N: " << v_dials[d].GetSplineName(idx) << std::endl;
        }
        */
    }

    std::vector<XsecNorm> &v_norms = m_norms.at(det);
    int num_norms = v_norms.size();
    for(int d = 0; d < num_norms; ++d)
    {
        int idx = m_norm_evtmap[nsample][nevent][d];
        if (idx==1)  weight *= params[d + m_offset.at(det) + num_dials];
    }

    event -> AddEvWght(weight);
}

void XsecParameters::AddDetector(const std::string& det, const std::string& config)
{
    std::cout << TAG << "Adding detector " << det << " for " << m_name << std::endl;
    std::fstream f;
    f.open(config, std::ios::in);

    json j;
    f >> j;

    std::string input_dir = std::string(std::getenv("XSLLHFITTER"))
                            + j["input_dir"].get<std::string>();
    std::cout << TAG << "Adding the following dials." << std::endl;

    std::vector<int> global_dimensions = j["dimensions"].get<std::vector<int>>();
    std::vector<XsecDial> v_dials;
    for(const auto& dial : j["dials"])
    {
        if(dial["use"] == true)
        {
            std::string fname_binning = input_dir + dial["binning"].get<std::string>();
            std::string fname_splines = input_dir + dial["splines"].get<std::string>();
            std::vector<int> dimensions = dial.value("dimensions", global_dimensions);

            XsecDial x(dial["name"], fname_binning, fname_splines);
            x.SetVars(dial["nominal"], dial["step"], dial["limit_lo"], dial["limit_hi"]);
            x.SetDimensions(dimensions);
            x.Print(false);
            v_dials.emplace_back(x);
        }
    }

    std::vector<XsecNorm> v_norms;
    for(const auto& norm : j["norms"])
    {
        if(norm["use"] == true)
        {
            std::vector<int> topologies = norm["topologies"].get<std::vector<int>>();
            std::vector<std::string> definition=norm["definition"].get<std::vector<std::string>>();
            std::vector<std::vector<double>> def_value =norm["def_value"].get<std::vector<std::vector<double>>>();

            XsecNorm x(norm["name"], topologies, definition, def_value, norm["fixed"]);
            x.SetVars(norm["nominal"], norm["step"], norm["limit_lo"], norm["limit_hi"]);
            x.Print(false);
            v_norms.emplace_back(x);
        }
    }

    v_detectors.emplace_back(det);
    m_dials.insert(std::make_pair(det, v_dials));
    m_norms.insert(std::make_pair(det, v_norms));
}
