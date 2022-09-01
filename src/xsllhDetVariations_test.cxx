#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TStyle.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TROOT.h"
#include "TVector.h"
//#include "ROOT/TSeq.hxx"

#include "json.hpp"
using json = nlohmann::json;

#include "BinManager.hh"
#include "ColorOutput.hh"
#include "ProgressBar.hh"

struct FileOptions
{
    std::string fname_input;
    std::string tree_name;
    std::string fname_input_nominal;
    std::string tree_name_nominal;
    std::string detector;
    unsigned int num_samples;
    unsigned int num_toys;
    unsigned int num_syst;
    std::vector<int> cuts;
    std::map<int, std::vector<int>> samples;
    std::vector<BinManager> bin_manager;
};

int main(int argc, char** argv)
{
    const std::string TAG = color::GREEN_STR + "[xsDetVariation]: " + color::RESET_STR;
    const std::string ERR = color::RED_STR + color::BOLD_STR + "[ERROR]: " + color::RESET_STR;

    std::cout << "--------------------------------------------------------\n"
              << TAG << "Welcome to the Super-xsLLh Detector Variation Interface.\n"
              << TAG << "Initializing the variation machinery..." << std::endl;

    ProgressBar pbar(60, "#");
    pbar.SetRainbow();
    pbar.SetPrefix(std::string(TAG + "Reading Events "));

    int sysi=-1;
    std::string outname;

    int nToysOverride = -1;

    bool do_mc_stat = false;

    std::string json_file;
    char option;
    while((option = getopt(argc, argv, "j:h:s:o:t:m")) != -1)
    {
        switch(option)
        {
            case 'j':
                json_file = optarg;
                break;
            case 's':
                sysi = atoi(optarg);
                break;
            case 'o':
                outname = optarg;
                break;
            case 't':
                nToysOverride = atoi(optarg);
                std::cout << "nToysOverride= "<<nToysOverride<<std::endl;
                break;
            case 'm':
                do_mc_stat = true;
                std::cout << "Enabling MC stat errors"<<std::endl;
                break;
            case 'h':
                std::cout << "USAGE: " << argv[0] << "\nOPTIONS:\n"
                          << "-j : JSON input\n";
            default:
                return 0;
        }
    }

    std::fstream f;
    f.open(json_file, std::ios::in);
    std::cout << TAG << "Opening " << json_file << std::endl;
    if(!f.is_open())
    {
        std::cout << ERR << "Unable to open JSON configure file." << std::endl;
        return 1;
    }

    json j;
    f >> j;

    bool do_projection  = j["projection"];
    bool do_single_syst = j["single_syst"];
    bool do_covariance  = j["covariance"];
    bool do_print       = j["pdf_print"];

    unsigned int syst_idx = j["syst_idx"];
    double weight_cut     = j["weight_cut"];

    if (sysi>-1) syst_idx = sysi;

    std::string fname_output  = j["fname_output"];
    std::string variable_plot = j["plot_variable"];
    std::string cov_mat_name  = j["covariance_name"];
    std::string cor_mat_name  = j["correlation_name"];

    if (outname.size()>0) fname_output=outname;

    std::vector<std::string> var_names = j["var_names"].get<std::vector<std::string>>();

    std::vector<FileOptions> v_files;
    for(const auto& file : j["files"])
    {
        if(file["use"])
        {
            FileOptions f;
            f.fname_input = file["fname_input"];
            f.tree_name   = file["tree_name"];
            f.fname_input_nominal = file["fname_input_nominal"];
            f.tree_name_nominal   = file["tree_name_nominal"];
            f.detector    = file["detector"];
            f.num_toys    = file["num_toys"];
            f.num_syst    = file["num_syst"];
            f.num_samples = file["num_samples"];
            f.cuts        = file["cuts"].get<std::vector<int>>();

            std::map<std::string, std::vector<int>> temp_json = file["samples"];
            for(const auto& kv : temp_json)
                f.samples.emplace(std::make_pair(std::stoi(kv.first), kv.second));

            f.bin_manager.resize(f.samples.size());
            std::map<std::string, std::string> temp_bins = file["sample_binning"];
            for(const auto& kv : temp_bins)
                f.bin_manager.at(std::stoi(kv.first)) = std::move(BinManager(kv.second));
                //f.bin_manager.emplace_back(BinManager(kv.second));

            v_files.emplace_back(f);
            std::cout << TAG << "File input : "<<f.fname_input<<std::endl;
        }
    }

    const int nvars = var_names.size();

    std::cout << TAG << "Output ROOT file: " << fname_output << std::endl
              << TAG << "Toy Weight Cut: " << weight_cut << std::endl
              << TAG << "Calculating Covariance: " << std::boolalpha << do_covariance << std::endl;

    std::cout << TAG << "Covariance Variables: ";
    for(const auto& var : var_names)
        std::cout << var << " ";
    std::cout << std::endl;

    int var_plot = -1;
    if(do_projection)
    {
        auto it  = std::find(var_names.begin(), var_names.end(), variable_plot);
        var_plot = std::distance(var_names.begin(), it);
    }

    std::cout << TAG << "Initalizing histograms." << std::endl;
    TH1::SetDefaultSumw2();
    std::vector<std::vector<TH1D>> v_hists;
    std::vector<TH1D> v_avg;
    std::vector<TH1F> v_mc_stat;
    for(const auto& file : v_files)
    {
        for(const auto& kv : file.samples)
        {
            const int sam   = kv.first;
            BinManager bm   = file.bin_manager.at(sam);
            const int nbins = sam==0 ? bm.GetNbins()+1
                                     : bm.GetNbins();
            std::vector<TH1D> v_temp;

            for(unsigned int t = 0; t < file.num_toys; ++t)
            {
                std::stringstream ss;
                ss << file.detector << "_sample" << sam << "_toy" << t;
                if(do_projection)
                {
                    std::vector<double> v_bins = bm.GetBinVector(var_plot);
                    v_temp.emplace_back(
                        TH1D(ss.str().c_str(), ss.str().c_str(), v_bins.size() - 1, &v_bins[0]));
                    if(t == 0)
                    {
                        ss.str("");
                        ss << file.detector << "_sample" << sam << "_avg";
                        v_avg.emplace_back(TH1D(ss.str().c_str(), ss.str().c_str(),
                                                v_bins.size() - 1, &v_bins[0]));
                        ss.str("");
                        ss << "cov_sample" << sam << "_mc_stat";
                        v_mc_stat.emplace_back(
                            TH1F(ss.str().c_str(), ss.str().c_str(), v_bins.size() - 1, &v_bins[0]));

                    }
                }
                else
                {
                    v_temp.emplace_back(TH1D(ss.str().c_str(), ss.str().c_str(), nbins, 0, nbins));
                    if(t == 0)
                    {
                        ss.str("");
                        ss << file.detector << "_sample" << sam << "_avg";
                        v_avg.emplace_back(
                            TH1D(ss.str().c_str(), ss.str().c_str(), nbins, 0, nbins));
                        ss.str("");
                        ss << "cov_sample" << sam << "_mc_stat";
                        v_mc_stat.emplace_back(TH1F(ss.str().c_str(), ss.str().c_str(), nbins, 0, nbins));

                    }
                }
            }
            v_hists.emplace_back(v_temp);
        }
    }

    std::cout << TAG << "Finished initializing histograms" << std::endl
              << TAG << "Reading events from files..." << std::endl;

    unsigned int cov_toys = 0;
    unsigned int offset   = 0;
    int fileCount = 0;
    std::cout << TAG << "v_files.size() = "<<v_files.size() << std::endl;
    for(const auto& file : v_files)
    {   //std::cout << TAG << "Staring file event loop " << file.fname_input << std::endl;
        //if (!file.fname_input.size()>0) continue;
        if (fileCount>=v_files.size()) break;fileCount++;
        std::cout << TAG << "Staring reading file" << file.fname_input << std::endl;
        int accum_level[file.num_toys][file.num_samples];
 //       std::cout << TAG << "checkpoint 1" << std::endl;
        float hist_variables[nvars][file.num_toys];
        float weight_syst_total_noflux[file.num_toys];
        float weight_syst[file.num_toys][file.num_syst];
        float selmu_mom[file.num_toys];
        float selmu_dir[3];
        //float PosPionMom[file.num_toys][1000];
        float PosPionDir[100][3];
        //float TPCProtonMom[file.num_toys][1000];
        float TPCProtonDir[100][3];
//std::cout << TAG << "checkpoint 2" << std::endl;
        int accum_level_nominal[1][file.num_samples];
        float hist_variables_nominal[nvars][1];
        float weight_syst_total_noflux_nominal[10];
        float weight_syst_nominal[file.num_syst];
     //   std::cout << TAG << "checkpoint 4" << std::endl;
        float selmu_mom_nominal[100];
        Float_t selmu_dir_nominal[10];
   //     std::cout << TAG << "checkpoint 5" << std::endl;
        float PosPionMom_nominal[1][500];
        float PosPionDir_nominal[100][3];
 //       std::cout << TAG << "checkpoint 6" << std::endl;
        float TPCProtonMom_nominal[1][500];
        float TPCProtonDir_nominal[100][3];
//std::cout << TAG << "checkpoint 3" << std::endl;

        std::cout << TAG << "Opening file: " << file.fname_input << std::endl
                  << TAG << "Reading tree: " << file.tree_name << std::endl
                  << TAG << "Num Toys: " << file.num_toys << std::endl
                  << TAG << "Num Syst: " << file.num_syst << std::endl;

        std::cout << TAG << "Branch to Sample mapping:" << std::endl;
        for(const auto& kv : file.samples)
        {
            std::cout << TAG << "Sample " << kv.first << ": ";
            for(const auto& b : kv.second)
                std::cout << b << " ";
            std::cout << std::endl;
        }

        //TFile* file_input = TFile::Open(file.fname_input.c_str(), "READ");
        //TTree* tree_event = (TTree*)file_input->Get(file.tree_name.c_str());
//std::cout << TAG << "checkpoint 7" << std::endl;
        TChain* tree_event;
        tree_event = new TChain(file.tree_name.c_str());
        tree_event->Add(file.fname_input.c_str());

        tree_event->SetBranchStatus("*",0); //disable all branches
        tree_event->SetBranchStatus("accum_level",1);
        tree_event->SetBranchStatus("weight_syst",1);
        tree_event->SetBranchStatus("weight_syst_total",1);
        for(unsigned int i = 0; i < nvars; ++i)
            tree_event->SetBranchStatus(var_names[i].c_str(),1);
        tree_event->SetBranchStatus("TPCProtonMom",1);

//std::cout << TAG << "checkpoint 8" << std::endl;
        tree_event->SetBranchAddress("accum_level", &accum_level);
        tree_event->SetBranchAddress("weight_syst", &weight_syst);
        tree_event->SetBranchAddress("weight_syst_total", &weight_syst_total_noflux);
        for(unsigned int i = 0; i < nvars; ++i)
            tree_event->SetBranchAddress(var_names[i].c_str(), &hist_variables[i]);
        //tree_event->SetBranchAddress("selmu_mom", selmu_mom);
        //tree_event->SetBranchAddress("selmu_dir", selmu_dir);
        //std::cout << TAG << "checkpoint 8.5" << std::endl;
        //ULong_t length1 = tree_event->GetBranch("PosPionMom")->GetLeaf("PosPionMom")->GetLenStatic();
        const int nToys = file.num_toys; 
        //const int len1 = length1;
        //float** PosPionMom = new float[nToys][len1];
        //float PosPionMom[nToys][len1];
        //float** PosPionMom = new float*[file.num_toys];
        //for (int i = 0;i<file.num_toys;i++) 
          //  PosPionMom[i] = new float[len1];
        //tree_event->SetBranchAddress("PosPionMom", PosPionMom);
        //tree_event->SetBranchAddress("PosPionDir", PosPionDir);
        //std::cout << TAG << "checkpoint 8.7" << std::endl;
        ULong_t length2 = tree_event->GetBranch("TPCProtonMom")->GetLeaf("TPCProtonMom")->GetLenStatic();
        //std::cout << TAG << "checkpoint 8.8 length2 = " <<length2<< " file.num_toys = "<<file.num_toys<< std::endl;
        const int len2 = length2;
        //float** TPCProtonMom = new float[file.num_toys][len2];
        float TPCProtonMom[file.num_toys][len2];
        //float** TPCProtonMom = new float*[file.num_toys];
        //for (int i = 0;i<file.num_toys;i++) 
          //  TPCProtonMom[i] = new float[len2];
        //std::cout << TAG << "checkpoint 8.9" << std::endl;
        tree_event->SetBranchAddress("TPCProtonMom", &TPCProtonMom);

        //std::cout << TAG << "checkpoint 8.95" << std::endl;
        //tree_event->SetBranchAddress("TPCProtonDir", TPCProtonDir);
//std::cout << TAG << "checkpoint 9" << std::endl;
        TChain* tree_event_nominal;
        tree_event_nominal = new TChain(file.tree_name_nominal.c_str());
        tree_event_nominal->Add(file.fname_input_nominal.c_str());
//std::cout << TAG << "checkpoint 10" << std::endl;
        tree_event_nominal->SetBranchAddress("accum_level", accum_level_nominal);
        tree_event_nominal->SetBranchAddress("weight_syst", weight_syst_nominal);
        tree_event_nominal->SetBranchAddress("weight_syst_total", weight_syst_total_noflux_nominal);
        for(unsigned int i = 0; i < nvars; ++i)
            tree_event_nominal->SetBranchAddress(var_names[i].c_str(), hist_variables_nominal[i]);
        tree_event_nominal->SetBranchAddress("selmu_mom", selmu_mom_nominal);
        tree_event_nominal->SetBranchAddress("selmu_dir", selmu_dir_nominal);
        tree_event_nominal->SetBranchAddress("PosPionMom", PosPionMom_nominal);
        tree_event_nominal->SetBranchAddress("PosPionDir", PosPionDir_nominal);
        tree_event_nominal->SetBranchAddress("TPCProtonMom", TPCProtonMom_nominal);
        tree_event_nominal->SetBranchAddress("TPCProtonDir", TPCProtonDir_nominal);
//std::cout << TAG << "checkpoint 11" << std::endl;

        unsigned int rejected_weights = 0;
        unsigned int total_weights    = 0;
        const unsigned int num_events = tree_event->GetEntries();
        const unsigned int num_events_nominal = tree_event_nominal->GetEntries();

        std::cout << TAG << "Number of events: " << num_events << std::endl;

        int nToysInLoop = file.num_toys;
        if (nToysOverride>-1) nToysInLoop = nToysOverride;

        tree_event->SetMaxVirtualSize(4e9);
        //ROOT::EnableImplicitMT(8);
//#pragma omp parallel for num_threads(8)
        //for (auto i : ROOT::TSeqUL(num_events))
        for(unsigned int i = 0; i < num_events; ++i)
        {   //if (i>1000) continue;
            //std::cout << TAG << "Event: " << i << std::endl;

            tree_event->GetEntry(i);
            /*tree_event->GetBranch("accum_level")->GetEntry(i);
            tree_event->GetBranch("weight_syst")->GetEntry(i);
            tree_event->GetBranch("weight_syst_total")->GetEntry(i);
            for(unsigned int iv = 0; iv < nvars; ++iv)
              tree_event->GetBranch(var_names[iv].c_str())->GetEntry(i);
            tree_event->GetBranch("TPCProtonMom")->GetEntry(i);*/

            //std::cout << TAG << "checkpoint 12" << std::endl;
            if(i % 2000 == 0 || i == (num_events - 1))
                pbar.Print(i, num_events - 1);

            for(unsigned int t = 0; t < nToysInLoop; ++t)
            {   
                //if (nToysOverride>-1 && t>=nToysOverride ) break;
                //std::cout << TAG << "t = " <<t<< std::endl;
                for(const auto& kv : file.samples)
                {
                    unsigned int s = kv.first;
                    //file.bin_manager[s].Print();
                    for(const auto& branch : kv.second)
                    {
                        if(accum_level[t][branch] > file.cuts[branch]-1)
                        {
                            float idx = -1;
                            if(do_projection)
                                idx = hist_variables[var_plot][t];
                            else
                            {
                                //hist_variables[0][t]=0.1;
                                //if (hist_variables[1][t]<-700) hist_variables[1][t]=-699;//overflow for dptt
                                //if (hist_variables[1][t]> 700) hist_variables[1][t]= 699;
                                //if (accum_level[t][branch]==file.cuts[branch]) hist_variables[1][t]=750;//OOPS bin for dptt
                                //if (abs(hist_variables[1][t])>1500) hist_variables[1][t]=1499; //overflow for pN
                                //if (accum_level[t][branch]==file.cuts[branch]) hist_variables[1][t]=1550;//OOPS bin for pN
                                //if (accum_level[t][branch]==file.cuts[branch]) hist_variables[1][t]=185;//OOPS bin for daT
                                //std::vector<double> vars;
                                //std::cout << TAG << "checkpoint 21" << std::endl;
                                //vars.push_back(selmu_mom[t]);vars.push_back(selmu_dir[2]);
                                //std::cout << TAG << "checkpoint 22" << std::endl;
                                //vars.push_back(/*PosPionMom[t][0]*/150.);vars.push_back(PosPionDir[0][2]);
                                //std::cout << TAG << "checkpoint 23" << std::endl;
                                //vars.push_back(TPCProtonMom[t][0]);vars.push_back(TPCProtonDir[0][2]);
                                //std::cout << TAG << "checkpoint 24" << std::endl;
                                //vars.push_back(hist_variables[1][t]);
                                //std::cout<<t<<" "<<selmu_mom[t]<<" "<<selmu_dir[2]<<" "<<PosPionMom[t][0]<<" "<<PosPionDir[0][2]<<" "<<TPCProtonMom[t][0]<<" "<<TPCProtonDir[0][2]<<std::endl;
                                //for(unsigned int v = 0; v < nvars; ++v)
                                //    vars.push_back(hist_variables[v][t]);
                                //idx = file.bin_manager[s].GetBinIndex(vars);
                                if (branch==2 && accum_level[t][branch]==file.cuts[branch]) idx = file.bin_manager[s].GetNbins();
                                else if (accum_level[t][branch] > file.cuts[branch]) {
                                    std::vector<double> vars;
                                    vars.push_back(/*selmu_mom[t]*/250.);vars.push_back(/*selmu_dir[2]*/0.5);
                                    vars.push_back(/*PosPionMom[t][0]*/150.);vars.push_back(/*PosPionDir[0][2]*/0.5);
                                    vars.push_back(TPCProtonMom[t][0]);vars.push_back(/*TPCProtonDir[0][2]*/0.5);
                                    vars.push_back(hist_variables[1][t]);
                                    idx = file.bin_manager[s].GetBinIndex(vars);//
                                }
                            }

                            float weight = do_single_syst ? weight_syst[t][syst_idx]
                                                          : weight_syst_total_noflux[t];
                            if(weight > 0.0 && weight < weight_cut)
                            {
                                v_hists[s + offset][t].Fill(idx, weight);
                                //v_avg[s + offset].Fill(idx, weight / file.num_toys);
                            }
                            else
                                rejected_weights++;
                            total_weights++;
                            break;
                        }
                    }
                }
            }
        }

        std::cout << TAG << "Number of events_nominal: " << num_events_nominal << std::endl;
        for(unsigned int i = 0; i < num_events_nominal; ++i)
        {   //if (i>1000) continue;
            tree_event_nominal->GetEntry(i);
            if(i % 2000 == 0 || i == (num_events_nominal - 1))
                pbar.Print(i, num_events_nominal - 1);

            for(unsigned int t = 0; t < 1; ++t)
            {   
                for(const auto& kv : file.samples)
                {
                    unsigned int s = kv.first;
                    for(const auto& branch : kv.second)
                    {
                        if(accum_level_nominal[t][branch] > file.cuts[branch]-1)
                        {
                            float idx = -1;
                            if(do_projection)
                                idx = hist_variables_nominal[var_plot][t];
                            else
                            {
                                //hist_variables_nominal[0][t]=0.1;
                                //if (hist_variables_nominal[1][t]<-700) hist_variables_nominal[1][t]=-699;//overflow for dptt
                                //if (hist_variables_nominal[1][t]> 700) hist_variables_nominal[1][t]= 699;
                                //if (accum_level_nominal[t][branch]==file.cuts[branch]) hist_variables_nominal[1][t]=750;//OOPS bin for dptt
                                //if (abs(hist_variables_nominal[1][t])>1500) hist_variables_nominal[1][t]=1499; //overflow for pN
                                //if (accum_level_nominal[t][branch]==file.cuts[branch]) hist_variables_nominal[1][t]=1550;//OOPS bin for pN
                                //if (accum_level_nominal[t][branch]==file.cuts[branch]) hist_variables_nominal[1][t]=185;//OOPS bin for daT
                                std::vector<double> vars;
                                //for(unsigned int v = 0; v < nvars; ++v)
                                //    vars.push_back(hist_variables_nominal[v][t]);
                                vars.push_back(selmu_mom_nominal[t]);vars.push_back(selmu_dir_nominal[2]);
                                vars.push_back(PosPionMom_nominal[t][0]);vars.push_back(PosPionDir_nominal[0][2]);
                                vars.push_back(TPCProtonMom_nominal[t][0]);vars.push_back(TPCProtonDir_nominal[0][2]);
                                vars.push_back(hist_variables_nominal[1][t]);
                                //if (branch==2 && accum_level_nominal[t][branch] > file.cuts[branch]) std::cout<<vars[0]<<" "<<vars[1]<<" "<<vars[2]<<" "<<vars[3]<<" "<<vars[4]<<" "<<vars[5]<<" "<<vars[6]<<" "<<weight_syst_total_noflux_nominal[0]<<std::endl;
                                idx = file.bin_manager[s].GetBinIndex(vars);
                                if (branch==2 && accum_level_nominal[t][branch]==file.cuts[branch]) idx = file.bin_manager[s].GetNbins();
                            }

                            float weight = do_single_syst ? weight_syst_nominal[syst_idx]
                                                          : weight_syst_total_noflux_nominal[0];
                            if(weight > 0.0 && weight < weight_cut)
                            {   //std::cout << TAG << "Sth in nominal tree!" << std::endl;
                                //v_hists[s + offset][t].Fill(idx, weight);
                                v_avg[s + offset].Fill(idx, weight);
                            }
                            else
                                rejected_weights++;
                            total_weights++;
                            break;
                        }
                    }
                }
            }
        }



        if(file.num_toys < cov_toys || cov_toys == 0)
            cov_toys = file.num_toys;

        if (nToysOverride>-1) cov_toys=nToysOverride;

        offset += file.samples.size();
        double reject_fraction = (rejected_weights * 1.0) / total_weights;
        std::cout << TAG << "Finished processing events." << std::endl;
        std::cout << TAG << "Total weights: " << total_weights << std::endl;
        std::cout << TAG << "Rejected weights: " << rejected_weights << std::endl;
        std::cout << TAG << "Rejected fraction: " << reject_fraction << std::endl;
    }
    std::cout << TAG << "Starting to calculate covariance matrix." << std::endl;

    unsigned int num_elements = 0;
    TMatrixTSym<double> cov_mat(num_elements);
    TMatrixTSym<double> cor_mat(num_elements);
    std::vector<float> v_mc_error;

    if(do_covariance)
    {
        std::cout << TAG << "Calculating covariance matrix." << std::endl;
        std::vector<std::vector<double>> v_toys;
        std::vector<std::vector<double>> v_nominal;

        for(unsigned int t = 0; t < cov_toys; ++t)
        {
            offset = 0;
            std::vector<double> i_toy;
            std::vector<double> i_nominal;
            for(const auto& file : v_files)
            {
                for(const auto& kv : file.samples)
                {
                    const unsigned int s     = kv.first;
                    const unsigned int nbins = s==0 ? file.bin_manager[s].GetNbins()+1
                                                    : file.bin_manager[s].GetNbins();
                    Double_t* w  = v_avg[s + offset].GetArray();
                    double* w2 = v_avg[s + offset].GetSumw2()->GetArray();
                    for(unsigned int b = 0; b < nbins; ++b) {
                        i_toy.emplace_back(v_hists[s + offset][t].GetBinContent(b + 1));
                        if (t==0) i_nominal.emplace_back(v_avg[s + offset].GetBinContent(b + 1));
                        float rel_error = w2[b+1] / (w[b+1] * w[b+1]);
                        v_mc_error.emplace_back(rel_error);
                    }
                }
                offset += file.samples.size();
            }
            v_toys.emplace_back(i_toy);
            if (t==0) v_nominal.emplace_back(i_nominal);
        }


        std::cout << TAG << "Using " << cov_toys << " toys." << std::endl;
        //num_elements = v_toys.at(0).size();
        std::vector<double> v_mean(num_elements, 0.0);
        cov_mat.ResizeTo(num_elements, num_elements);
        cor_mat.ResizeTo(num_elements, num_elements);
        cov_mat.Zero();
        cor_mat.Zero();

        //for(unsigned int t = 0; t < cov_toys; ++t)
        //{
            for(unsigned int i = 0; i < num_elements; ++i)
                v_mean[i] += v_nominal[0][i] ;
        //}

        for(unsigned int t = 0; t < cov_toys; ++t)
        {
            for(unsigned int i = 0; i < num_elements; ++i)
            {
                for(unsigned int j = 0; j < num_elements; ++j)
                {
                    if(v_mean[i] != 0 && v_mean[j] != 0)
                    {
                        cov_mat(i, j) += (1.0 - v_toys[t][i] / v_mean[i])
                                         * (1.0 - v_toys[t][j] / v_mean[j]) / (1.0 * cov_toys);
                    }
                }
            }
        }

        if(do_mc_stat)
        {
            std::cout << TAG << "Adding MC stat error to covariance." << std::endl;
            for(unsigned int i = 0; i < num_elements; ++i)
            {
                cov_mat(i, i) += v_mc_error[i];
            }
        }


        for(unsigned int i = 0; i < num_elements; ++i)
        {
            if(cov_mat(i, i) <= 0.0)
                cov_mat(i, i) = 1.0;
        }

        for(unsigned int i = 0; i < num_elements; ++i)
        {
            for(unsigned int j = 0; j < num_elements; ++j)
            {
                double bin_i  = cov_mat(i, i);
                double bin_j  = cov_mat(j, j);
                cor_mat(i, j) = cov_mat(i, j) / std::sqrt(bin_i * bin_j);
                if(std::isnan(cor_mat(i, j)))
                    cor_mat(i, j) = 0;
            }
        }
    }

    std::cout << TAG << "Saving to output file." << std::endl;
    TFile* file_output = TFile::Open(fname_output.c_str(), "RECREATE");
    file_output->cd();

    offset = 0;
    gStyle->SetOptStat(0);
    for(const auto& file : v_files)
    {
        for(const auto& kv : file.samples)
        {
            unsigned int s = kv.first;
            std::stringstream ss;
            ss << file.detector << "_sample" << s;
            TCanvas c(ss.str().c_str(), ss.str().c_str(), 1200, 900);
            v_avg[s + offset].Draw("axis");

            for(unsigned int t = 0; t < file.num_toys; ++t)
            {
                v_hists[s + offset][t].SetLineColor(kRed);
                if(do_projection)
                    v_hists[s + offset][t].Scale(1, "width");
                v_hists[s + offset][t].Draw("hist same");
            }

            v_avg[s + offset].SetLineColor(kBlack);
            v_avg[s + offset].SetLineWidth(2);
            if(do_projection)
                v_avg[s + offset].Scale(1, "width");
            v_avg[s + offset].GetYaxis()->SetRangeUser(0, v_avg[s + offset].GetMaximum() * 1.50);
            v_avg[s + offset].Draw("hist same");
            c.Write(ss.str().c_str());

            if(do_print)
                c.Print(std::string(ss.str() + ".pdf").c_str());
        }
        offset += file.samples.size();
    }

    if(do_covariance)
    {
        cov_mat.Write(cov_mat_name.c_str());
        cor_mat.Write(cor_mat_name.c_str());
    }
    if(do_mc_stat)
    {
        int num_elements = v_mc_error.size();
        TVectorT<float> v_mc_root(num_elements, v_mc_error.data());
        v_mc_root.Write("mc_stat_error");
    }


    file_output->Close();

    std::cout << TAG << "Finished." << std::endl;
    return 0;
}
