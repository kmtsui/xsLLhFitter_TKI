#ifndef XSECCALC_HH
#define XSECCALC_HH

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TMatrixT.h>
#include <TVectorT.h>
#include "TRandom3.h"

#include "BinManager.hh"
#include "ColorOutput.hh"
#include "FitObj.hh"
#include "OptParser.hh"
#include "ProgressBar.hh"
#include "ToyThrower.hh"

using TMatrixD    = TMatrixT<double>;
using TMatrixDSym = TMatrixTSym<double>;
using TVectorD    = TVectorT<double>;

struct SigNorm
{
    TH1D flux_hist;
    TH1D flux_throws;
    TH1D target_throws;
    std::string name;
    std::string detector;
    unsigned int nbins;

    std::string flux_file;
    std::string flux_name;
    double flux_int;
    double flux_err;
    bool use_flux_fit;

    double num_targets_val;
    double num_targets_err;
    bool is_rel_err;

    double potvalue;
};

class XsecCalc
{
public:
    XsecCalc(const std::string& json_config);
    ~XsecCalc();

    void ReadFitFile(const std::string& file);
    void UsePrefitCov();
    void CalcTrueXsec();

    void ReweightParam(const std::vector<double>& param);
    void ReweightBestFit();
    void ReweightNominal();
    void GenerateToys(bool throw_eff=false, int fluxrw = -1);
    void GenerateToys(const int ntoys, bool throw_eff=false, int fluxrw = -1);

    void CalcCovariance(bool use_best_fit, bool fit_hydrogen=false, int idxH = -1, int nPbins = -1, double binwidth = -1);

    void ApplyEff(std::vector<TH1D>& sel_hist, std::vector<TH1D>& tru_hist, bool is_toy, bool throw_eff=false);
    void ApplyNorm(std::vector<TH1D>& vec_hist, const std::vector<double>& param, bool is_toy);
    void ApplyTargets(const unsigned int signal_id, TH1D& hist, bool is_toy);
    void ApplyFlux(const unsigned int signal_id, TH1D& hist, const std::vector<double>& param,
                   bool is_toy);
    void ApplyBinWidth(const unsigned int signal_id, TH1D& hist, const double unit_scale);

    std::vector<double> GetAllBinWidth(std::vector<TH1D>& sel_hist, const double unit_scale);
    double GetTotalXsec(TH1D vec_hist);

    TH1D ConcatHist(const std::vector<TH1D>& vec_hists, const std::string& hist_name = "");

    std::vector<TH1D> GetSelSignal() { return selected_events->GetSignalHist(); };
    TH1D GetSelSignal(const int signal_id) { return selected_events->GetSignalHist(signal_id); };

    std::vector<TH1D> GetTruSignal() { return true_events->GetSignalHist(); };
    TH1D GetTruSignal(const int signal_id) { return true_events->GetSignalHist(signal_id); };

    void SaveOutput(bool save_toys = false);
    void SaveSignalHist(TFile* file);
    void SaveExtra(TFile* file);
    void SetOutputFile(const std::string& override_file) { output_file = override_file; };

    std::string GetOutputFileName() const { return output_file; };
    std::string GetInputFileName() const { return input_file; };
    unsigned int GetNumToys() const { return num_toys; };
    void SetNumToys(const int n) { num_toys = n; };

private:
    void InitToyThrower();
    void InitFluxHist();
    void InitNormalization(const nlohmann::json& j, const std::string input_dir);

    FitObj* selected_events;
    FitObj* true_events;

    ToyThrower* toy_thrower;

    TMatrixDSym* postfit_cov;
    TMatrixDSym* postfit_cor;
    std::vector<double> postfit_param;
    std::vector<double> true_param;
    std::vector<double> prefit_param;

    TMatrixDSym xsec_cov;
    TMatrixDSym xsec_cor;
    TMatrixDSym xsec_shape_cov;

    TMatrixDSym eff_cov;
    TMatrixDSym eff_cor;

    TH1D sel_best_fit_before_apply;
    TH1D sel_best_fit_after_eff;
    TH1D sel_shape_best_fit;
    std::vector<double> all_bin_width;

    TH1D sel_best_fit;
    TH1D tru_best_fit;
    TH1D eff_best_fit;
    std::vector<TH1D> signal_best_fit;
    std::vector<TH1D> toys_sel_events;
    std::vector<TH1D> toys_tru_events;
    std::vector<TH1D> toys_eff;
    std::vector<SigNorm> v_normalization;
    double sel_norm_best_fit;
    std::vector<double> toys_sel_norm;

    TH1D sel_true_fit;
    TH1D tru_true_fit;

    TH1D* eff_variation;
    TRandom3* RNG;

    std::string input_file;
    std::string output_file;
    std::string extra_hists;

    unsigned int num_toys;
    unsigned int rng_seed;
    unsigned int num_signals;
    unsigned int total_signal_bins;

    bool do_incompl_chol;
    bool do_force_posdef;
    double dropout_tol;
    double force_padd;

    bool use_prefit_cov;
    const double perMeV = 1.0;
    const double perGeV = 1000.0;

    const std::string TAG = color::YELLOW_STR + "[XsecExtract]: " + color::RESET_STR;
    const std::string ERR = color::RED_STR + "[ERROR]: " + color::RESET_STR;
    const std::string COMMENT_CHAR = "#";
};

#endif
