#ifndef XSECNORM_H
#define XSECNORM_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <TCollection.h>
#include <TFile.h>
#include <TGraph.h>
#include <TKey.h>

#include "BinManager.hh"
#include "ColorOutput.hh"

class XsecNorm
{
    public:
        XsecNorm(const std::string& dial_name);
        XsecNorm(const std::string& norm_name, const std::vector<int> norm_topologies, const std::vector<std::string> norm_definition, const std::vector<std::vector<double>> norm_def_value, bool norm_fixed);
        XsecNorm(const std::string& dial_name, const std::string& fname_splines,
                 const std::string& fname_binning);

        void SetBinning(const std::string& fname_binning);
        void ReadSplines(const std::string& fname_splines);

        int GetSplineIndex(int topology, int reaction, double q2) const;
        int GetSplineIndex(const std::vector<int>& var, const std::vector<double>& bin, const std::vector<double>& ptheta) const;
        double GetSplineValue(int index, double dial_value) const;
        std::string GetSplineName(int index) const;

        void SetVars(double nominal, double step, double limit_lo, double limit_hi);
        void SetDimensions(int num_top, int num_reac);
        void SetDimensions(const std::vector<int>& dim);
        void Print(bool print_bins = false) const;

        std::string GetName() const { return m_name; }
        double GetNominal() const { return m_nominal; }
        double GetStep() const { return m_step; }
        double GetLimitLow() const { return m_limit_lo; }
        double GetLimitHigh() const { return m_limit_hi; }

        int GetNDef() const {return m_definition.size();};
        std::string GetDef(int i) const {return m_definition.at(i);};
        bool IsApplied(int index, int topo, double val) const;
        bool IsFixed() const {return m_fixed;};

    private:
        int ntop;
        int nreac;
        int nbins;
        std::string m_name;
        double m_nominal;
        double m_step;
        double m_limit_lo;
        double m_limit_hi;
        std::vector<TGraph> v_splines;
        std::vector<int> m_dimensions;
        BinManager bin_manager;

        std::vector<int> m_topologies;
        std::vector<std::string> m_definition;
        std::vector<std::vector<double>> m_def_value;

        bool m_fixed;

        const std::string TAG = color::MAGENTA_STR + "[XsecNorm]: " + color::RESET_STR;
};

#endif
