#ifndef __FluxParameters_hh__
#define __FluxParameters_hh__

#include "AnaFitParameters.hh"

class FluxParameters : public AnaFitParameters
{
public:
    FluxParameters(const std::string& name);
    ~FluxParameters();

    void InitParameters();
    void InitEventMap(std::vector<AnaSample*>& sample, int mode);
    void ReWeight(AnaEvent* event, const std::string& det, int nsample, int nevent,
                  std::vector<double>& params);
    void AddDetector(const std::string& det, const std::vector<double>& bins);
    int GetBinIndex(const std::string& det, double enu);
    int GetDetectorOffset(const std::string& det) const { return m_det_offset.at(det); };

private:
    std::vector<double> m_enubins;
    std::map<std::string, int> m_det_offset;
    std::map<std::string, std::vector<double>> m_det_bins;
    std::vector<std::string> v_detectors;

    const std::string TAG = color::GREEN_STR + "[FluxParameters]: " + color::RESET_STR;
};

#endif
