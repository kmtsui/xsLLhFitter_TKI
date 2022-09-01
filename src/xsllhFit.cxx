#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

#include "AnaSample.hh"
#include "AnaTreeMC.hh"
#include "ColorOutput.hh"
#include "DetParameters.hh"
#include "FitParameters.hh"
#include "FluxParameters.hh"
#include "OptParser.hh"
#include "XsecFitter.hh"
#include "XsecParameters.hh"

int main(int argc, char** argv)
{
    const std::string TAG = color::CYAN_STR + "[xsFit]: " + color::RESET_STR;
    const std::string ERR = color::RED_STR + color::BOLD_STR
                            + "[ERROR]: " + color::RESET_STR;

    std::cout << "------------------------------------------------\n"
              << TAG << color::RainbowText("Welcome to the Super-xsLLhFitter.\n")
              << TAG << color::RainbowText("Initializing the fit machinery...") << std::endl;

    const std::string xslf_env = std::getenv("XSLLHFITTER");
    if(xslf_env.empty())
    {
        std::cerr << ERR << "Environment variable \"XSLLHFITTER\" not set." << std::endl
                  << ERR << "Cannot determine source tree location." << std::endl;
        return 1;
    }

    std::string json_file;
    std::string outname;
    int randomSeed=-1;
    double regStrength=-1;
    bool dry_run = false;
    char option;
    bool useDetBinning = false;
    bool minosErr = false;
    bool tkiFit = false;
    bool tkiConstraint = false;
    bool throwPrefit = false;
    bool fixFitParameters = false;
    int llh_func = 0;
    while((option = getopt(argc, argv, "j:o:s:r:l:nhdmtcpf")) != -1)
    {
        switch(option)
        {
            case 'j':
                json_file = optarg;
                break;
            case 'l':
                llh_func = atoi(optarg);
                break;
            case 'n':
                dry_run = true;
                break;
            case 'm':
                minosErr = true;
                break;
            case 't':
                tkiFit = true;
                break;
            case 'c':
                tkiConstraint = true;
                break;
            case 'f':
                fixFitParameters = true;
                break;
            case 'p':
                throwPrefit = true;
                break;
            case 'd':
                useDetBinning = true;
                break;
            case 'o':
                outname = optarg;
                break;
            case 's':
                randomSeed = atoi(optarg);
                break;
            case 'r':
                regStrength = atof(optarg);
                break;
            case 'h':
                std::cout << "USAGE: "
                          << argv[0] << "\nOPTIONS:\n"
                          << "-j : JSON input\n"
                          << "-n : Dry run - Set up but do not run fit.\n";
            default:
                return 0;
        }
    }

    OptParser parser;
    if(!parser.ParseJSON(json_file))
    {
        std::cerr << ERR << "JSON parsing failed. Exiting.\n";
        return 1;
    }
    parser.PrintOptions();

    std::string input_dir = parser.input_dir;
    std::string fname_data = parser.fname_data;
    std::string fname_mc   = parser.fname_mc;
    std::string fname_output = parser.fname_output;
    if (outname.size()>0) fname_output=outname;
    std::cout<<"fname_output = "<<fname_output<<std::endl;;
    std::cout<<"regStrength = "<<regStrength<<std::endl;;
    std::vector<std::string> topology = parser.sample_topology;
    std::string det_binning = parser.det_binning;

    const double potD  = parser.data_POT;
    const double potMC = parser.mc_POT;
    int seed = parser.rng_seed;
    if (randomSeed>-1) seed=randomSeed;
    std::cout<<"seed = "<<seed<<std::endl;;
    int threads = parser.num_threads;

    //Setup data trees
    TFile* fdata = TFile::Open(fname_data.c_str(), "READ");
    TTree* tdata = (TTree*)(fdata->Get("selectedEvents"));

    std::cout << TAG << "Configure file parsing finished." << std::endl;
    std::cout << TAG << "Opening " << fname_data << " for data selection.\n"
              << TAG << "Opening " << fname_mc << " for MC selection." << std::endl;

    /*************************************** FLUX *****************************************/
    std::cout << TAG << "Setup Flux Covariance" << std::endl;

    //input File
    TFile* finfluxcov = TFile::Open(parser.flux_cov.fname.c_str(), "READ"); //contains flux systematics info
    std::cout << TAG << "Opening " << parser.flux_cov.fname << " for flux covariance." << std::endl;
    //setup enu bins and covm for flux
    //TH1D* nd_numu_bins_hist = (TH1D*)finfluxcov->Get(parser.flux_cov.binning.c_str());
    //TAxis* nd_numu_bins = nd_numu_bins_hist->GetXaxis();
    TAxis* nd_numu_bins = (TAxis*)finfluxcov->Get("nd5_numode_numu_bins");

    //std::vector<double> enubins;
    //enubins.push_back(nd_numu_bins -> GetBinLowEdge(1));
    //for(int i = 0; i < nd_numu_bins -> GetNbins(); ++i)
    //    enubins.push_back(nd_numu_bins -> GetBinUpEdge(i+1));

    std::cout<< TAG << "Loaded "<<nd_numu_bins -> GetNbins()<<" flux bins."<<std::endl;

    TMatrixDSym *cov_flux_in   = (TMatrixDSym*)finfluxcov->Get("total_flux_cov");
    TMatrixDSym cov_flux(nd_numu_bins->GetNbins());
    std::vector<double> enubins;
    enubins.push_back(nd_numu_bins->GetBinLowEdge(1));
    for(int i=0;i<nd_numu_bins->GetNbins();i++)
    {
      enubins.push_back(nd_numu_bins->GetBinUpEdge(i+1));
      for(int j=0;j<nd_numu_bins->GetNbins();j++)
      {
        cov_flux(i, j) = (*cov_flux_in)(i,j);
      }
    }


    //Cov mat stuff:
    //TMatrixDSym* cov_flux = (TMatrixDSym*)finfluxcov -> Get(parser.flux_cov.matrix.c_str());
    finfluxcov -> Close();

    /*************************************** FLUX END *************************************/
    std::cout << TAG << "Setup Xsec Covariance" << std::endl;
    TFile* file_xsec_cov = TFile::Open(parser.xsec_cov.fname.c_str(), "READ");
    std::cout << TAG << "Opening " << parser.xsec_cov.fname << " for xsec covariance." << std::endl;
    TMatrixDSym* cov_xsec = (TMatrixDSym*)file_xsec_cov -> Get(parser.xsec_cov.matrix.c_str());
    file_xsec_cov -> Close();
    /*(*cov_xsec)(15,15)=100;
    (*cov_xsec)(16,16)=100;
    (*cov_xsec)(17,17)=100;
    (*cov_xsec)(18,18)=100;*/
    /*TMatrixDSym cov_xsec(12);
  // Cov mat for xsec should be an input, not hard coded, WIP!
    cov_xsec(0,0) = 0.4*0.4; //MAQE
    cov_xsec(1,1) = 66.7/(6.*217.)*66.7/(6.*217.); //PF_C
    cov_xsec(2,2) = 0.02*0.02; //CCNuE
    cov_xsec(3,3) = 0.4*0.4; //DISMPiShp
    cov_xsec(4,4) = 0.3*0.3; //NCCoh
    cov_xsec(5,5) = 0.3*0.3; //NCOth
    cov_xsec(6+0,6+0)= 0.17;//FSI
    cov_xsec(6+0,6+1)= -0.002778;
    cov_xsec(6+0,6+2)= 0;
    cov_xsec(6+0,6+3)= 0.02273;
    cov_xsec(6+0,6+4)= 0.005;
    cov_xsec(6+0,6+5)= 0;
    cov_xsec(6+1,6+1)= 0.1142;
    cov_xsec(6+1,6+2)= -0.1667;
    cov_xsec(6+1,6+3)= -0.001263;
    cov_xsec(6+1,6+4)= -0.002083;
    cov_xsec(6+1,6+5)=  -0.09259;
    cov_xsec(6+2,6+2)= 0.25;
    cov_xsec(6+2,6+3)= -5.204e-18;
    cov_xsec(6+2,6+4)= 0;
    cov_xsec(6+2,6+5)= 0.1389;
    cov_xsec(6+3,6+3)= 0.1694;
    cov_xsec(6+3,6+4)= -0.002273;
    cov_xsec(6+3,6+5)= -3.469e-18;
    cov_xsec(6+4,6+4)= 0.3213;
    cov_xsec(6+4,6+5)= 1.735e-18;
    cov_xsec(6+5,6+5)= 0.07716;*/


    TFile* fout = TFile::Open(fname_output.c_str(), "RECREATE");
    std::cout << TAG << "Open output file: " << fname_output << std::endl;

    // Add analysis samples:

    std::vector<AnaSample*> samples;

    for(const auto& opt : parser.samples)
    {
        if(opt.use_sample == true && opt.cut_branch >= 0)
        {
            std::cout << TAG << "Adding new sample to fit.\n"
                      << TAG << "Name: " << opt.name << std::endl
                      << TAG << "CutB: " << opt.cut_branch << std::endl
                      << TAG << "Detector: " << opt.detector << std::endl
                      << TAG << "Use Sample: " << std::boolalpha << opt.use_sample << std::endl;

            auto s = new AnaSample(opt.cut_branch, opt.name, opt.detector, opt.binning, tdata, llh_func);
            s -> SetNorm(potD/potMC);
            //if(opt.cut_branch >= 0)
                samples.push_back(s);
        }
    }

    //read MC events
    AnaTreeMC selTree(fname_mc.c_str(), "selectedEvents", true);
    std::cout << TAG << "Reading and collecting events." << std::endl;
    selTree.GetEvents(samples, parser.signal_definition, false);

    std::cout << TAG << "Getting sample breakdown by reaction." << std::endl;
    for(auto& sample : samples)
        sample -> GetSampleBreakdown(fout, "nominal", topology, false);


    //*************** FITTER SETTINGS **************************
    //In the bit below we choose which params are used in the fit
    //For stats only just use fit params
    //**********************************************************

    //define fit param classes
    std::vector<AnaFitParameters*> fitpara;

    //Fit parameters
    FitParameters sigfitpara("par_fit", false, tkiFit, fixFitParameters,tkiConstraint);
    if(parser.regularise) {
        if (regStrength>0) sigfitpara.SetRegularisation(regStrength, parser.reg_method);
        else sigfitpara.SetRegularisation(parser.reg_strength, parser.reg_method);
    }
    for(const auto& opt : parser.detectors)
    {
        if(opt.use_detector)
            sigfitpara.AddDetector(opt.name, parser.signal_definition);
            //sigfitpara.AddDetector(opt.name, opt.binning);
    }
    sigfitpara.InitEventMap(samples, 0);
    fitpara.push_back(&sigfitpara);

    //Flux parameters
    FluxParameters fluxpara("par_flux");
    if(parser.flux_cov.do_fit)
    {
        fluxpara.SetCovarianceMatrix(cov_flux, parser.flux_cov.decompose);
        fluxpara.SetThrow(parser.flux_cov.do_throw);
        fluxpara.SetInfoFrac(parser.flux_cov.info_frac);
        for(const auto& opt : parser.detectors)
        {
            if(opt.use_detector)
                fluxpara.AddDetector(opt.name, enubins);
        }
        fluxpara.InitEventMap(samples, 0);
        fitpara.push_back(&fluxpara);
    }

    XsecParameters xsecpara("par_xsec");
    if(parser.xsec_cov.do_fit)
    {
        xsecpara.SetCovarianceMatrix(*cov_xsec, parser.xsec_cov.decompose);
        xsecpara.SetThrow(parser.xsec_cov.do_throw);
        for(const auto& opt : parser.detectors)
        {
            if(opt.use_detector)
                xsecpara.AddDetector(opt.name, opt.xsec);
        }
        xsecpara.InitEventMap(samples, 0);
        fitpara.push_back(&xsecpara);
    }

    std::cout << TAG << "Setup Detector Covariance" << std::endl;
    TFile* file_detcov = TFile::Open(parser.det_cov.fname.c_str(), "READ");
    TMatrixDSym* cov_det_in = (TMatrixDSym*)file_detcov -> Get(parser.det_cov.matrix.c_str());
    TMatrixDSym cov_det = *cov_det_in;
    file_detcov -> Close();//TMatrixDSym cov_det;

    DetParameters detpara("par_det");
    if(parser.det_cov.do_fit)
    {
        detpara.SetCovarianceMatrix(cov_det, parser.det_cov.decompose);
        detpara.SetThrow(parser.det_cov.do_throw);
        detpara.SetInfoFrac(parser.det_cov.info_frac);
        for(const auto& opt : parser.detectors)
        {
            if(opt.use_detector) {
                if (useDetBinning) {
                     detpara.AddDetector(opt.name, samples, false, det_binning);
                     detpara.SetUseDetBinning(true);
                }
                else detpara.AddDetector(opt.name, samples, true);
            }
        }
        detpara.InitEventMap(samples, 0);
        fitpara.push_back(&detpara);
    }

    //Instantiate fitter obj
    XsecFitter xsecfit(fout, seed, threads, minosErr, throwPrefit);
    //xsecfit.SetSaveFreq(10000);
    xsecfit.SetMinSettings(parser.min_settings);
    xsecfit.SetPOTRatio(potD/potMC);
    xsecfit.SetTopology(topology);
    xsecfit.SetZeroSyst(parser.zero_syst);

    //init w/ para vector
    xsecfit.InitFitter(fitpara);
    std::cout << TAG << "Fitter initialised." << std::endl;

    /*
    for(int i = 0; i < sigfitpara.GetNpar(); ++i)
    {
        xsecfit.FixParameter("par_fit_ND280_" + std::to_string(i), 1.0);
        xsecfit.FixParameter("par_fit_INGRID_" + std::to_string(i), 1.0);
    }
    */

    bool did_converge = false;
    if(!dry_run)
    {
        did_converge = xsecfit.Fit(samples, parser.fit_type, parser.stat_fluc);
        if(!did_converge)
            std::cout << TAG << "Fit did not coverge." << std::endl;
    }
    fout -> Close();

    std::cout << TAG << "\u3042\u308a\u304c\u3068\u3046\u3054\u3056\u3044\u307e\u3057\u305f\uff01"
              << std::endl;
    return 0;
}
