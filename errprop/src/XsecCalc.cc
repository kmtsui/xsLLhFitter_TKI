#include "XsecCalc.hh"
using json = nlohmann::json;

XsecCalc::XsecCalc(const std::string& json_config)
    : num_toys(0)
    , rng_seed(0)
    , num_signals(0)
    , total_signal_bins(0)
    , postfit_cov(nullptr)
    , postfit_cor(nullptr)
    , toy_thrower(nullptr)
    , eff_variation(nullptr)
{
    std::cout << TAG << "Reading error propagation options." << std::endl;
    std::fstream f;
    f.open(json_config, std::ios::in);

    json j;
    f >> j;

    std::string input_dir
        = std::string(std::getenv("XSLLHFITTER")) + j["input_dir"].get<std::string>();

    input_file = input_dir + j["input_fit_file"].get<std::string>();
    output_file = j["output_file"].get<std::string>();

    extra_hists = j.value("extra_hists", "");
    if(!extra_hists.empty())
        extra_hists = input_dir + extra_hists;

    num_toys = j["num_toys"];
    rng_seed = j["rng_seed"];

    do_incompl_chol = j["decomposition"].value("incomplete_chol", false);
    dropout_tol = j["decomposition"].value("drop_tolerance", 1.0E-3);
    do_force_posdef = j["decomposition"].value("do_force_posdef", false);
    force_padd = j["decomposition"].value("force_posdef_val", 1.0E-9);

    bool use_det_binning = j["use_det_binning"];

    std::string sel_json_config = input_dir + j["sel_config"].get<std::string>();
    std::string tru_json_config = input_dir + j["tru_config"].get<std::string>();

    std::cout << TAG << "Input file from fit: " << input_file << std::endl
              << TAG << "Output xsec file: " << output_file << std::endl
              << TAG << "Num. toys: " << num_toys << std::endl
              << TAG << "RNG  seed: " << rng_seed << std::endl
              << TAG << "Selected events config: " << sel_json_config << std::endl
              << TAG << "True events config: " << tru_json_config << std::endl;

    std::cout << TAG << "Reading post-fit file..." << std::endl;
    TH1::AddDirectory(false);
    ReadFitFile(input_file);

    std::cout << TAG << "Initializing fit objects..." << std::endl;
    selected_events = new FitObj(sel_json_config, "selectedEvents", false, use_det_binning);
    true_events = new FitObj(tru_json_config, "trueEvents", true);
    total_signal_bins = selected_events->GetNumSignalBins();

    InitNormalization(j["sig_norm"], input_dir);
    std::cout << TAG << "Finished initialization." << std::endl;

    std::string eff_file_string = j["eff_file"].get<std::string>();
    std::string eff_hist = j["eff_hist"].get<std::string>();
    TFile* eff_file = TFile::Open(eff_file_string.c_str(), "READ");
    if (eff_file && eff_file->Get(eff_hist.c_str())) {
      std::cout << TAG << "Getting efficiency variation from : " << eff_file_string << std::endl
                << TAG << "Efficiency histogram: " << eff_hist << std::endl;
      eff_variation = (TH1D*)eff_file->Get("eff_best_fit");
      RNG = new TRandom3(rng_seed);
      eff_file->Close();
    }
}

XsecCalc::~XsecCalc()
{
    delete toy_thrower;
    delete selected_events;
    delete true_events;

    delete postfit_cov;
    delete postfit_cor;
}

void XsecCalc::ReadFitFile(const std::string& file)
{
    if(postfit_cov != nullptr)
        delete postfit_cov;
    if(postfit_cor != nullptr)
        delete postfit_cor;
    postfit_param.clear();
    prefit_param.clear();
    true_param.clear();

    std::cout << TAG << "Opening " << file << std::endl;
    input_file = file;

    TFile* postfit_file = TFile::Open(file.c_str(), "READ");
    postfit_cov = (TMatrixDSym*)postfit_file->Get("res_cov_matrix");
    postfit_cor = (TMatrixDSym*)postfit_file->Get("res_cor_matrix");

    TVectorD* postfit_param_root = (TVectorD*)postfit_file->Get("res_vector");
    for(int i = 0; i < postfit_param_root->GetNoElements(); ++i)
        postfit_param.emplace_back((*postfit_param_root)[i]);

    TVectorD* true_param_root;
    if (postfit_file->Get("vec_par_all_iter0")){
        true_param_root = (TVectorD*)postfit_file->Get("vec_par_all_iter0");
        for(int i = 0; i < true_param_root->GetNoElements(); ++i)
            true_param.emplace_back((*true_param_root)[i]);
    }

    postfit_file->Close();
    use_prefit_cov = false;

    std::cout << TAG << "Successfully read fit file." << std::endl;
    InitToyThrower();
}

void XsecCalc::UsePrefitCov()
{
    if(selected_events == nullptr)
        std::cout << TAG << "FitObj not initialized for prefit covariance." << std::endl;
    else
    {
        use_prefit_cov = true;
        InitToyThrower();
    }
}

void XsecCalc::InitToyThrower()
{
    std::cout << TAG << "Initializing toy-thrower..." << std::endl;
    if(toy_thrower != nullptr)
        delete toy_thrower;

    TMatrixDSym cov_mat;
    if(use_prefit_cov)
    {
        std::cout << TAG << "Using prefit covariance matrix." << std::endl;
        cov_mat.ResizeTo(selected_events->GetNpar(), selected_events->GetNpar());
        cov_mat = selected_events->GetPrefitCov();
    }
    else
    {
        std::cout << TAG << "Using postfit covariance matrix." << std::endl;
        cov_mat.ResizeTo(postfit_cov->GetNrows(), postfit_cov->GetNrows());
        cov_mat = *postfit_cov;
    }

    toy_thrower = new ToyThrower(cov_mat, rng_seed, false, 1E-48);
    if(do_force_posdef)
    {
        if(!toy_thrower->ForcePosDef(force_padd, 1E-48))
        {
            std::cout << ERR << "Covariance matrix could not be made positive definite.\n"
                << "Exiting." << std::endl;
            exit(1);
        }
    }

    if(do_incompl_chol)
    {
        std::cout << TAG << "Performing incomplete Cholesky decomposition." << std::endl;
        toy_thrower->IncompCholDecomp(dropout_tol, true);
    }
    else
    {
        std::cout << TAG << "Performing ROOT Cholesky decomposition." << std::endl;
        toy_thrower->SetupDecomp(1E-48);
    }
}

void XsecCalc::InitNormalization(const nlohmann::json& j, const std::string input_dir)
{

    for(const auto& sig_def : selected_events->GetSignalDef())
    {
        if(sig_def.use_signal == true)
        {
            std::cout << TAG << "Adding normalization parameters for " << sig_def.name << " signal."
                      << std::endl;

            json s;
            try
            {
                s = j.at(sig_def.name);
            }
            catch(json::exception& e)
            {
                std::cout << ERR << "Signal " << sig_def.name
                          << " not found in error propagation config file." << std::endl;
                exit(1);
            }

            SigNorm n;
            n.name = sig_def.name;
            n.detector = sig_def.detector;
            n.flux_file = s["flux_file"];
            n.flux_name = s["flux_hist"];
            n.flux_int = s["flux_int"];
            n.flux_err = s["flux_err"];
            n.use_flux_fit = s["use_flux_fit"];
            n.num_targets_val = s["num_targets_val"];
            n.num_targets_err = s["num_targets_err"];
            n.is_rel_err = s["relative_err"];
            n.potvalue=s["POT"];
            if(n.is_rel_err)
            {
                n.flux_err = n.flux_int * n.flux_err;
                n.num_targets_err = n.num_targets_val * n.num_targets_err;
            }

            BinManager bm(sig_def.binning);
            n.nbins = bm.GetNbins();

            std::cout << TAG << "Num. bins: " << n.nbins << std::endl
                      << TAG << "Flux file: " << n.flux_file << std::endl
                      << TAG << "Flux hist: " << n.flux_name << std::endl
                      << TAG << "Flux integral: " << n.flux_int << std::endl
                      << TAG << "Flux error: " << n.flux_err << std::endl
                      << TAG << "Use flux fit: " << std::boolalpha << n.use_flux_fit << std::endl
                      << TAG << "Num. targets: " << n.num_targets_val << std::endl
                      << TAG << "Num. targets err: " << n.num_targets_err << std::endl
                      << TAG << "Relative err: " << std::boolalpha << n.is_rel_err << std::endl;

            std::string temp_fname = input_dir + n.flux_file;
            TFile* temp_file = TFile::Open(temp_fname.c_str(), "READ");
            temp_file->cd();
            TH1D* temp_hist = (TH1D*)temp_file->Get(n.flux_name.c_str());
            n.flux_hist = *temp_hist;
            temp_file->Close();

            unsigned int nbins = 100;
            temp_hist
                = new TH1D("", "", nbins, n.flux_int - 5 * n.flux_err, n.flux_int + 5 * n.flux_err);
            n.flux_throws = *temp_hist;
            temp_hist = new TH1D("", "", nbins, n.num_targets_val - 5 * n.num_targets_err,
                                 n.num_targets_val + 5 * n.num_targets_err);
            n.target_throws = *temp_hist;
            v_normalization.push_back(n);
        }
    }
    num_signals = v_normalization.size();
}

void XsecCalc::ReweightNominal()
{
    selected_events->ReweightNominal();
    true_events->ReweightNominal();
}

void XsecCalc::CalcTrueXsec()
{
    if (true_param.size()!=postfit_param.size()){
        std::cout << TAG << "In XsecCalc::CalcTrueXsec(), true_param not initialzied" << std::endl;
        return;
    }

    ReweightParam(true_param);

    auto sel_hists = selected_events->GetSignalHist();
    auto tru_hists = true_events->GetSignalHist();

    ApplyEff(sel_hists, tru_hists, false);

    ApplyNorm(sel_hists, true_param, false);

    sel_true_fit = ConcatHist(sel_hists, "sel_best_fit");
    tru_true_fit = ConcatHist(tru_hists, "tru_best_fit");
}


void XsecCalc::ReweightBestFit()
{
    ReweightParam(postfit_param);

    auto sel_hists = selected_events->GetSignalHist();
    auto tru_hists = true_events->GetSignalHist();
    /*for (int i =0;i<sel_hists.size();i++) {
       std::cout<<"sel_hists["<<i<<"].Integral()="<<sel_hists[i].Integral()<<std::endl;
       for(int j = 1; j <= sel_hists[i].GetNbinsX(); ++j) std::cout<<"Bin["<<j<<"]="<<sel_hists[i].GetBinContent(j)<<std::endl;
    }*/

    sel_best_fit_before_apply=ConcatHist(sel_hists, "sel_best_fit_before_apply");

    /*for (int i =0;i<tru_hists.size();i++) {
      std::cout<<"tru_hists["<<i<<"].Integral()="<<tru_hists[i].Integral()<<std::endl;
       for(int j = 1; j <= tru_hists[i].GetNbinsX(); ++j) {
         std::cout<<"Bin["<<j<<"]="<<tru_hists[i].GetBinContent(j)<<std::endl;
         if (tru_hists[i].GetBinContent(j)==0) {
           std::cout<<"Maually set bin content to non-zero!Remember to change!"<<std::endl;
           tru_hists[i].SetBinContent(j,sel_hists[i].GetBinContent(j)*10);
         }
       }
    }*/

    ApplyEff(sel_hists, tru_hists, false);

    sel_best_fit_after_eff=ConcatHist(sel_hists, "sel_best_fit_after_eff");

    ApplyNorm(sel_hists, postfit_param, false);

    all_bin_width=GetAllBinWidth(sel_hists,perMeV);

    sel_best_fit = ConcatHist(sel_hists, "sel_best_fit");
    tru_best_fit = ConcatHist(tru_hists, "tru_best_fit");
    signal_best_fit = std::move(sel_hists);
    sel_norm_best_fit = GetTotalXsec(sel_best_fit);
}

double XsecCalc::GetTotalXsec(TH1D vec_hist)
{
    double totalXsec = 0;
    for (int i=0;i<vec_hist.GetNbinsX()-1;i++){
        totalXsec += vec_hist.GetBinContent(i+1)*all_bin_width[i];
    }
    return totalXsec;
}

void XsecCalc::ReweightParam(const std::vector<double>& param)
{
    selected_events->ReweightEvents(param);
    true_events->ReweightEvents(param);
}

void XsecCalc::GenerateToys(bool throw_eff, int fluxrw) { GenerateToys(num_toys,throw_eff, fluxrw); }

void XsecCalc::GenerateToys(const int ntoys,bool throw_eff, int fluxrw)
{
    num_toys = ntoys;
    std::cout << TAG << "Throwing " << num_toys << " toys..." << std::endl;

    int fluxparpos[2];
    selected_events->GetFluxParPos(fluxparpos);
    if (fluxrw==1) {
        std::cout << TAG << "Fixing at best-fit flux..." << std::endl
                  << TAG << fluxparpos[1]<<" Flux par starts at "<<fluxparpos[0]<<  std::endl;
    }

    ProgressBar pbar(60, "#");
    pbar.SetRainbow();
    pbar.SetPrefix(std::string(TAG + "Throwing "));

    toys_sel_norm.clear();


    for(int i = 0; i < ntoys; ++i)
    {
        if(i % 10 == 0 || i == (ntoys - 1))
            pbar.Print(i, ntoys - 1);

        const unsigned int npar = postfit_param.size();
        std::vector<double> toy(npar, 0.0);
        toy_thrower->Throw(toy);

        if (fluxrw==1) {
            for (int ip=fluxparpos[0]; ip<fluxparpos[0]+fluxparpos[1]; ++ip)
                toy[ip] = 0;
        }

        std::transform(toy.begin(), toy.end(), postfit_param.begin(), toy.begin(),
                       std::plus<double>());
        for(int ip = 0; ip < npar; ++ip)
        {
            if(toy[ip] < 0.0)
                toy[ip] = 0.01;
        }

        selected_events->ReweightEvents(toy);
        true_events->ReweightEvents(toy);

        auto sel_hists = selected_events->GetSignalHist();
        auto tru_hists = true_events->GetSignalHist();

    /*for (int i =0;i<tru_hists.size();i++) {
      std::cout<<"tru_hists["<<i<<"].Integral()="<<tru_hists[i].Integral()<<std::endl;
       for(int j = 1; j <= tru_hists[i].GetNbinsX(); ++j) {
         std::cout<<"Bin["<<j<<"]="<<tru_hists[i].GetBinContent(j)<<std::endl;
         if (tru_hists[i].GetBinContent(j)==0) {
           std::cout<<"Maually set bin content to non-zero!Remember to change!"<<std::endl;
           tru_hists[i].SetBinContent(j,sel_hists[i].GetBinContent(j)*10);
         }
       }
    }*/

        ApplyEff(sel_hists, tru_hists, true, throw_eff);
        ApplyNorm(sel_hists, toy, true);

        TH1D toy_hist = ConcatHist(sel_hists, ("sel_signal_toy" + std::to_string(i)));
        toys_sel_events.emplace_back(toy_hist);
        toys_tru_events.emplace_back(ConcatHist(tru_hists, ("tru_signal_toy" + std::to_string(i))));

        toys_sel_norm.emplace_back(GetTotalXsec(toy_hist));

        /*
        total_signal_bins = npar;
        std::string temp = "toy" + std::to_string(i);
        TH1D h_toy(temp.c_str(), temp.c_str(), total_signal_bins, 0, total_signal_bins);
        for(int i = 0; i < total_signal_bins; ++i)
            h_toy.SetBinContent(i+1, toy[i]);
        toys_sel_events.emplace_back(h_toy);
        */
    }
}

void XsecCalc::ApplyEff(std::vector<TH1D>& sel_hist, std::vector<TH1D>& tru_hist, bool is_toy, bool throw_eff)
{
    std::vector<TH1D> eff_hist;
    for(int i = 0; i < num_signals; ++i)
    {
        TH1D eff = sel_hist[i];
        eff.Divide(&sel_hist[i], &tru_hist[i]);

        if (throw_eff && eff_variation != nullptr) {
          //std::cout<<"Throwing efficiency"<<std::endl;
          for (int j = 1; j <= sel_hist[i].GetNbinsX(); ++j) {
            eff.SetBinContent(j,eff_variation->GetBinContent(j)+RNG->Gaus()*(eff_variation->GetBinError(j)+0.01));
          }
        }

        eff_hist.emplace_back(eff);

        for(int j = 1; j <= sel_hist[i].GetNbinsX(); ++j)
        {
            double bin_eff = eff.GetBinContent(j);
            double bin_val = sel_hist[i].GetBinContent(j);
            sel_hist[i].SetBinContent(j, bin_val / bin_eff);
        }
    }

    if(is_toy)
    {
        std::string eff_name = "eff_combined_toy" + std::to_string(toys_eff.size());
        toys_eff.emplace_back(ConcatHist(eff_hist, eff_name));
    }
    else
    {
        eff_best_fit = ConcatHist(eff_hist, "eff_best_fit");
    }
}

std::vector<double> XsecCalc::GetAllBinWidth(std::vector<TH1D>& sel_hist, const double unit_scale)
{
    std::vector<double> binWidth;
    for(unsigned int i = 0; i < num_signals; ++i)
    {
        BinManager bm = selected_events->GetBinManager(i);
        for(int j = 0; j < sel_hist[i].GetNbinsX(); ++j)
        {
            binWidth.push_back(bm.GetBinWidth(j) / unit_scale);
        }
    }
    return binWidth;
}

void XsecCalc::ApplyNorm(std::vector<TH1D>& vec_hist, const std::vector<double>& param, bool is_toy)
{
    for(unsigned int i = 0; i < num_signals; ++i)
    {
        ApplyTargets(i, vec_hist[i], is_toy);
        ApplyFlux(i, vec_hist[i], param, is_toy);
        ApplyBinWidth(i, vec_hist[i], perMeV);
    }
}

void XsecCalc::ApplyTargets(const unsigned int signal_id, TH1D& hist, bool is_toy)
{
    double num_targets = 1.0;
    if(is_toy)
    {
        num_targets = toy_thrower->ThrowSinglePar(v_normalization[signal_id].num_targets_val,
                                                  v_normalization[signal_id].num_targets_err);
        v_normalization[signal_id].target_throws.Fill(num_targets);
    }
    else
        num_targets = v_normalization[signal_id].num_targets_val;

    //std::cout<<"num_targets = "<<num_targets<<std::endl;

    hist.Scale(1.0 / num_targets);
}

void XsecCalc::ApplyFlux(const unsigned int signal_id, TH1D& hist, const std::vector<double>& param,
                         bool is_toy)
{
    if(v_normalization[signal_id].use_flux_fit)
    {
        double flux_int = selected_events->ReweightFluxHist(
            param, v_normalization[signal_id].flux_hist, v_normalization[signal_id].detector);
        flux_int *= v_normalization[signal_id].potvalue/1.E21/0.05; //flux unit in flux file
        hist.Scale(1.0 / flux_int);
        //std::cout<<"flux_int = "<<flux_int<<std::endl;
        if(is_toy)
            v_normalization[signal_id].flux_throws.Fill(flux_int);
    }
    else
    {
        double flux_int = 1.0;
        if(is_toy)
        {
            flux_int = toy_thrower->ThrowSinglePar(v_normalization[signal_id].flux_int,
                                                   v_normalization[signal_id].flux_err);
            v_normalization[signal_id].flux_throws.Fill(flux_int);
        }
        else
            flux_int = v_normalization[signal_id].flux_int;

        //std::cout<<"flux_int = "<<flux_int<<std::endl;

        hist.Scale(1.0 / flux_int);
    }
}

void XsecCalc::ApplyBinWidth(const unsigned int signal_id, TH1D& hist, const double unit_scale)
{
    BinManager bm = selected_events->GetBinManager(signal_id);
    for(int i = 0; i < hist.GetNbinsX(); ++i)
    {
        const double bin_width = bm.GetBinWidth(i) / unit_scale;
        const double bin_value = hist.GetBinContent(i + 1);

        //std::cout<<"bin_width = "<<bin_width<<std::endl;

        hist.SetBinContent(i + 1, bin_value / bin_width);
    }
}

TH1D XsecCalc::ConcatHist(const std::vector<TH1D>& vec_hist, const std::string& hist_name)
{
    TH1D hist_combined(hist_name.c_str(), hist_name.c_str(), total_signal_bins, 0,
                       total_signal_bins);

    unsigned int bin = 1;
    for(const auto& hist : vec_hist)
    {
        for(int i = 1; i <= hist.GetNbinsX(); ++i)
            hist_combined.SetBinContent(bin++, hist.GetBinContent(i));
    }

    return hist_combined;
}

void XsecCalc::CalcCovariance(bool use_best_fit, bool fit_hydrogen, int idxH, int nPbins, double binwidth)
{
    std::cout << TAG << "Calculating covariance matrix..." << std::endl;
    std::cout << TAG << "Using " << num_toys << " toys." << std::endl;

    TH1D h_cov;
    TH1D h_eff_cov;
    if(use_best_fit)
    {
        ReweightBestFit();
        h_cov = sel_best_fit;
        h_eff_cov = eff_best_fit;
        std::cout << TAG << "Using best fit for covariance." << std::endl;
    }
    else
    {
        TH1D h_mean("", "", total_signal_bins, 0, total_signal_bins);
        TH1D h_eff_mean("", "", total_signal_bins, 0, total_signal_bins);
        for(const auto& hist : toys_sel_events)
        {
            for(int i = 0; i < total_signal_bins; ++i)
                h_mean.Fill(i + 0.5, hist.GetBinContent(i + 1));
        }
        h_mean.Scale(1.0 / (1.0 * num_toys));
        h_cov = h_mean;

        for(const auto& hist : toys_eff)
        {
            for(int i = 0; i < total_signal_bins; ++i)
                h_eff_mean.Fill(i + 0.5, hist.GetBinContent(i + 1));
        }
        h_eff_mean.Scale(1.0 / (1.0 * num_toys));
        h_eff_cov = h_mean;

        for(int i = 1; i <= total_signal_bins; ++i)
            std::cout << "Bin " << i << ": " << h_cov.GetBinContent(i) << std::endl;

        std::cout << TAG << "Using mean of toys for covariance." << std::endl;
    }

    xsec_cov.ResizeTo(total_signal_bins, total_signal_bins);
    xsec_cov.Zero();

    xsec_cor.ResizeTo(total_signal_bins, total_signal_bins);
    xsec_cor.Zero();

    xsec_shape_cov.ResizeTo(total_signal_bins, total_signal_bins);
    xsec_shape_cov.Zero();

    eff_cov.ResizeTo(total_signal_bins, total_signal_bins);
    eff_cov.Zero();

    eff_cor.ResizeTo(total_signal_bins, total_signal_bins);
    eff_cor.Zero();

    for(const auto& hist : toys_sel_events)
    {
      if (nPbins<0)
        for(int i = 0; i < total_signal_bins; ++i)
        {
            for(int j = 0; j < total_signal_bins; ++j)
            {
                double x = hist.GetBinContent(i + 1) - h_cov.GetBinContent(i + 1);
                double y = hist.GetBinContent(j + 1) - h_cov.GetBinContent(j + 1);
                if (fit_hydrogen) {
                    if (idxH<0 || idxH==i) x += hist.GetBinContent(total_signal_bins-1) - h_cov.GetBinContent(total_signal_bins-1);
                    if (idxH<0 || idxH==j) y += hist.GetBinContent(total_signal_bins-1) - h_cov.GetBinContent(total_signal_bins-1);
                }
                xsec_cov(i, j) += x * y / (1.0 * num_toys);

            }
        }
      else if (nPbins>0)
        for(int i = 0; i < total_signal_bins/nPbins; ++i)
        {
            for(int j = 0; j < total_signal_bins/nPbins; ++j)
            {
                //double x = hist.GetBinContent(i + 1) - h_cov.GetBinContent(i + 1);
                //double y = hist.GetBinContent(j + 1) - h_cov.GetBinContent(j + 1);
                int ndptt = total_signal_bins/nPbins;
                double x = 0;
                double y = 0;
                for (int k=0;k<nPbins;k++) x += hist.GetBinContent(i + k*ndptt+ 1) - h_cov.GetBinContent(i + k*ndptt+ + 1);
                for (int k=0;k<nPbins;k++) y += hist.GetBinContent(j + k*ndptt+ 1) - h_cov.GetBinContent(j + k*ndptt+ + 1);
                if (binwidth>0) {x*=binwidth;y*=binwidth;}
                xsec_cov(i, j) += x * y / (1.0 * num_toys);
            }
        }
    }

    for(const auto& hist : toys_eff)
    {
        for(int i = 0; i < total_signal_bins; ++i)
        {
            for(int j = 0; j < total_signal_bins; ++j)
            {
                double x = hist.GetBinContent(i + 1) - h_eff_cov.GetBinContent(i + 1);
                double y = hist.GetBinContent(j + 1) - h_eff_cov.GetBinContent(j + 1);
                eff_cov(i, j) += x * y / (1.0 * num_toys);
            }
        }
    }

    for(int itoy=0;itoy<toys_sel_events.size();itoy++){
        for(int i = 0; i < total_signal_bins; ++i){
            for(int j = 0; j < total_signal_bins; ++j){
                double x = toys_sel_events[itoy].GetBinContent(i + 1)/toys_sel_norm[itoy] - h_cov.GetBinContent(i + 1)/sel_norm_best_fit;
                double y = toys_sel_events[itoy].GetBinContent(j + 1)/toys_sel_norm[itoy] - h_cov.GetBinContent(j + 1)/sel_norm_best_fit;
                if (fit_hydrogen) {
                    if (idxH<0 || idxH==i) x += toys_sel_events[itoy].GetBinContent(total_signal_bins-1)/toys_sel_norm[itoy] - h_cov.GetBinContent(total_signal_bins-1)/sel_norm_best_fit;
                    if (idxH<0 || idxH==j) y += toys_sel_events[itoy].GetBinContent(total_signal_bins-1)/toys_sel_norm[itoy] - h_cov.GetBinContent(total_signal_bins-1)/sel_norm_best_fit;
                }
                xsec_shape_cov(i, j) += x * y / (1.0 * num_toys) * all_bin_width[i] * all_bin_width[j];
            }
        }
    }

    for(int i = 0; i < total_signal_bins; ++i)
    {
        for(int j = 0; j < total_signal_bins; ++j)
        {
            const double x = xsec_cov(i, i);
            const double y = xsec_cov(j, j);
            const double z = xsec_cov(i, j);
            xsec_cor(i, j) = z / (sqrt(x * y));

            if(std::isnan(xsec_cor(i, j)))
                xsec_cor(i, j) = 0.0;
        }
    }

    for(int i = 0; i < total_signal_bins; ++i)
    {
        for(int j = 0; j < total_signal_bins; ++j)
        {
            const double x = eff_cov(i, i);
            const double y = eff_cov(j, j);
            const double z = eff_cov(i, j);
            eff_cor(i, j) = z / (sqrt(x * y));

            if(std::isnan(eff_cor(i, j)))
                eff_cor(i, j) = 0.0;
        }
    }

    if (nPbins>0) {
        for(int i = 0; i < total_signal_bins/nPbins; ++i) {
           int ndptt = total_signal_bins/nPbins;
           double x = 0;
           for (int k=0;k<nPbins;k++) x += sel_best_fit.GetBinContent(i + k*ndptt+ 1);
           if (binwidth>0) {x*=binwidth;}
           sel_best_fit.SetBinContent(i + 1, x);
        }
    }

    if (fit_hydrogen) 
        for(int i = 0; i < total_signal_bins-2; ++i) 
            if (idxH<0 || idxH==i) {
                sel_best_fit.SetBinContent(i + 1, sel_best_fit.GetBinContent(i+1)+sel_best_fit.GetBinContent(total_signal_bins-1));
                if (sel_true_fit.GetNbinsX()==sel_best_fit.GetNbinsX())
                    sel_true_fit.SetBinContent(i + 1, sel_true_fit.GetBinContent(i+1)+sel_true_fit.GetBinContent(total_signal_bins-1));
            }

    for(int i = 0; i < total_signal_bins; ++i) 
        sel_best_fit.SetBinError(i + 1, sqrt(xsec_cov(i, i)));
    
    unsigned int idx = 0;
    for(int n = 0; n < signal_best_fit.size(); ++n)
    {
        unsigned int nbins = v_normalization.at(n).nbins;
        for(int i = 0; i < nbins; ++i)
            signal_best_fit.at(n).SetBinError(i + 1, sqrt(xsec_cov(i+idx,i+idx)));

        idx += nbins;
    }

    std::cout << TAG << "Covariance and correlation matrix calculated." << std::endl;
    std::cout << TAG << "Errors applied to histograms." << std::endl;
}

void XsecCalc::SaveOutput(bool save_toys)
{
    TFile* file = TFile::Open(output_file.c_str(), "RECREATE");
    std::cout << TAG << "Saving output to " << output_file << std::endl;

    file->cd();
    if(save_toys)
    {
        std::cout << TAG << "Saving toys to file." << std::endl;
        for(int i = 0; i < num_toys; ++i)
        {
            toys_sel_events.at(i).Write();
            toys_tru_events.at(i).Write();
            toys_eff.at(i).Write();
        }
    }

    TH1D* eff_hists[100];
    for (int i=0;i<100;i++) eff_hists[i]=new TH1D("","",1000,0,1);
    for (int i=0;i<num_toys;i++) 
    {
      for (int j = 0; j < toys_eff[i].GetNbinsX(); ++j) eff_hists[j]->Fill(toys_eff[i].GetBinContent(j+1));
    }
    for (int j=0;j<eff_best_fit.GetNbinsX();j++) eff_best_fit.SetBinError(j+1,eff_hists[j]->GetRMS());
    

    sel_best_fit_before_apply.Write("sel_best_fit_before_apply");
    sel_best_fit_after_eff.Write("sel_best_fit_after_eff");
    sel_best_fit.Write("sel_best_fit");
    tru_best_fit.Write("tru_best_fit");
    eff_best_fit.Write("eff_best_fit");

    sel_true_fit.Write("sel_true_fit");
    tru_true_fit.Write("tru_true_fit");

    xsec_cov.Write("xsec_cov");
    xsec_cor.Write("xsec_cor");
    xsec_shape_cov.Write("xsec_shape_cov");

    eff_cov.Write("eff_cov");
    eff_cor.Write("eff_cor");

    //if(use_prefit_cov)
        selected_events->GetPrefitCov().Write("prefit_cov");

    postfit_cov->Write("postfit_cov");
    postfit_cor->Write("postfit_cor");

    TVectorD postfit_param_root(postfit_param.size(), postfit_param.data());
    postfit_param_root.Write("postfit_param");

    SaveSignalHist(file);

    for(const auto& n : v_normalization)
    {
        std::string name = n.name + "_flux_nominal";
        n.flux_hist.Write(name.c_str());

        name = n.name + "_flux_throws";
        n.flux_throws.Write(name.c_str());

        name = n.name + "_target_throws";
        n.target_throws.Write(name.c_str());
    }

    SaveExtra(file);
    file->Close();
}

void XsecCalc::SaveSignalHist(TFile* file)
{
    file->cd();
    for(int id = 0; id < num_signals; ++id)
    {
        signal_best_fit.at(id).Write();

        BinManager bm = selected_events->GetBinManager(id);
        auto cos_edges = bm.GetEdgeVector(0);
        auto pmu_edges = bm.GetEdgeVector(1);

        std::vector<std::vector<double>> bin_edges;
        bin_edges.emplace_back(std::vector<double>());
        bin_edges.at(bin_edges.size()-1).emplace_back(pmu_edges.at(0).first);
        bin_edges.at(bin_edges.size()-1).emplace_back(pmu_edges.at(0).second);

        for(int m = 1; m < cos_edges.size(); ++m)
        {
            if(cos_edges[m] != cos_edges[m-1])
                bin_edges.emplace_back(std::vector<double>());

            bin_edges.at(bin_edges.size()-1).emplace_back(pmu_edges.at(m).first);
            bin_edges.at(bin_edges.size()-1).emplace_back(pmu_edges.at(m).second);
        }

        for(int n = 0; n < bin_edges.size(); ++n)
        {
            std::sort(bin_edges.at(n).begin(), bin_edges.at(n).end());
            auto iter = std::unique(bin_edges.at(n).begin(), bin_edges.at(n).end());
            bin_edges.at(n).erase(iter, bin_edges.at(n).end());
        }

        unsigned int offset = 0;
        for(int k = 0; k < bin_edges.size(); ++k)
        {
            std::string name = v_normalization.at(id).name + "_cos_bin" + std::to_string(k);
            TH1D temp(name.c_str(), name.c_str(), bin_edges.at(k).size()-1, bin_edges.at(k).data());

            for(int l = 1; l <= temp.GetNbinsX(); ++l)
            {
                temp.SetBinContent(l, signal_best_fit.at(id).GetBinContent(l+offset));
                temp.SetBinError(l, signal_best_fit.at(id).GetBinError(l+offset));
            }
            offset += temp.GetNbinsX();
            temp.GetXaxis()->SetRange(1,temp.GetNbinsX()-1);
            temp.Write();
        }
    }
}

void XsecCalc::SaveExtra(TFile* output)
{
    if(extra_hists.empty())
        return;

    std::ifstream fin(extra_hists, std::ios::in);
    if(!fin.is_open())
    {
        std::cout << ERR << "Failed to open file: " << extra_hists << std::endl;
        return;
    }
    else
    {
        std::cout << TAG << "Saving extra histograms from fit output." << std::endl;
        TFile* file = TFile::Open(input_file.c_str(), "READ");
        output->cd();

        std::string line;
        while(std::getline(fin, line))
        {
            if(line.front() != COMMENT_CHAR)
            {
                if (!file->Get(line.c_str())) continue;
                TH1D* temp = (TH1D*)file->Get(line.c_str());
                temp->Write(line.c_str());
            }
        }
        file->Close();
    }
}
