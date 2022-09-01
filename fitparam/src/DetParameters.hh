#ifndef __DetParameters_hh__
#define __DetParameters_hh__

#include <iomanip>
#include <iostream>
#include <string>

#include "AnaFitParameters.hh"
#include "FitStructs.hh"
#include "BinManager.hh"

class DetParameters : public AnaFitParameters
{
public:
    DetParameters(const std::string& name);
    ~DetParameters();

    void InitParameters();
    void InitEventMap(std::vector<AnaSample*>& sample, int mode);
    int GetBinIndex(const int sam, double D1, double D2) const;
    void ReWeight(AnaEvent* event, const std::string& det, int nsample, int nevent, std::vector<double>& params);
    bool SetBinning(const std::string& file_name, std::vector<xsllh::FitBin>& bins);
    void AddDetector(const std::string& det, std::vector<AnaSample*>& v_sample, bool match_bins, const std::string& detbinning="");
    void SetUseDetBinning (bool usedet) {useDetbinning=usedet;}

private:
    std::map<int, std::vector<xsllh::FitBin>> m_sample_bins;
    std::map<int, int> m_sample_offset;
    std::vector<int> v_samples;
    std::vector<BinManager> m_bin_manager;

    bool useDetbinning = false;

    const std::string TAG = color::GREEN_STR + "[DetParameters]: " + color::RESET_STR;
};

#endif
