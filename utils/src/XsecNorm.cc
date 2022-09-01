#include "XsecNorm.hh"

XsecNorm::XsecNorm(const std::string& dial_name)
    : m_name(dial_name)
{
}

XsecNorm::XsecNorm(const std::string& dial_name, const std::vector<int> norm_topologies, const std::vector<std::string> norm_definition, const std::vector<std::vector<double>> norm_def_value, bool norm_fixed)
    : m_name(dial_name), m_topologies(norm_topologies), m_definition(norm_definition), m_def_value(norm_def_value), m_fixed(norm_fixed)
{
}

XsecNorm::XsecNorm(const std::string& dial_name, const std::string& fname_binning,
                   const std::string& fname_splines)
    : m_name(dial_name)
{
    SetBinning(fname_binning);
    ReadSplines(fname_splines);
}

void XsecNorm::SetBinning(const std::string& fname_binning)
{
    bin_manager.SetBinning(fname_binning);
    nbins = bin_manager.GetNbins();
}

void XsecNorm::ReadSplines(const std::string& fname_splines)
{
    v_splines.clear();

    TFile* file_splines = TFile::Open(fname_splines.c_str(), "READ");
    if(file_splines == 0)
    {
        std::cout << "[ERROR]: Failed to open " << fname_splines << std::endl;
        return;
    }

    TIter key_list(file_splines -> GetListOfKeys());
    TKey* key = nullptr;
    while((key = (TKey*)key_list.Next()))
    {
        TGraph* spline = (TGraph*)key -> ReadObj();
        v_splines.emplace_back(*spline);
    }
    file_splines -> Close();
    delete file_splines;

    v_splines.shrink_to_fit();
}

int XsecNorm::GetSplineIndex(int topology, int reaction, double q2) const
{
    const int b = bin_manager.GetBinIndex(std::vector<double>{q2});
    return topology * nreac * nbins + reaction * nbins + b;
}

int XsecNorm::GetSplineIndex(const std::vector<int>& var, const std::vector<double>& bin, const std::vector<double>& ptheta) const
{
    if(var.size() != m_dimensions.size())
        return -1;

//    int idx = bin_manager.GetBinIndex(bin);
//    for(int i = 0; i < var.size(); ++i)
//        idx += var[i] * m_dimensions[i];
    int idxReco = bin_manager.GetBinIndex(std::vector<double>{0.1,bin[0]});
    int idxTrue = bin_manager.GetBinIndex(std::vector<double>{0.1,bin[1]});
    if (idxReco==-1||idxTrue==-1) return -1;
    if (idxReco>=m_dimensions[1]||idxTrue>=m_dimensions[0]) return -1;

    /*if (var[1]!=0) {
      if (idxReco==0) idxReco=2;
      else if (idxReco==1) idxReco=5;
    }*/
    /*std::string histname = Form("RecBin_%d_trueBin_%d_sample_%d_reac_%d",idxReco,idxTrue,var[1],var[0]);

    for (int i=0;i<v_splines.size();i++) {
      if (histname.compare(v_splines.at(i).GetName())==0) return i;
    }
    return -1;*/
    int idx = var[0]+idxTrue*12+idxReco*12*m_dimensions[0]+var[1]*12*m_dimensions[0]*m_dimensions[1]; //12 topo bins, m_dimensions[0] true bins, m_dimensions[1] reco bins, 5+1 samples

    if (m_name.find("PMOM") != std::string::npos) {
      int dialBin = (int)m_name.back()- '0';
      //std::cout<<"m_name = "<<m_name<<" dialBin = "<<dialBin<<std::endl;
      const int nPMombins = 3;
      const double PMombins[nPMombins+1] = {450,700,1000,1200};
      int pMomBin=-1;
      for (int i=0;i<nPMombins;i++) {
          if (ptheta[0]>PMombins[i] && ptheta[0]<PMombins[i+1]) {pMomBin=i;break;}
      }
      if (dialBin!=pMomBin) return -1;
    }
    if (m_name.find("PTHETA") != std::string::npos) {
      int dialBin = (int)m_name.back()- '0';
      const int nPThetabins = 3;
      const double PThetabins[nPThetabins+1] = {0.342,0.6,0.8,1.};
      int pThetaBin=-1;
      for (int i=0;i<nPThetabins;i++) {
          if (ptheta[1]>PThetabins[i] && ptheta[1]<PThetabins[i+1]) {pThetaBin=i;break;}
      }
      if (dialBin!=pThetaBin) return -1;
    }

    return idx;
}

double XsecNorm::GetSplineValue(int index, double dial_value) const
{
    if(index >= 0)
        return v_splines.at(index).Eval(dial_value);
    else
        return 1.0;
}

std::string XsecNorm::GetSplineName(int index) const
{
    return std::string(v_splines.at(index).GetName());
}

void XsecNorm::SetVars(double nominal, double step, double limit_lo, double limit_hi)
{
    m_nominal = nominal;
    m_step = step;
    m_limit_lo = limit_lo;
    m_limit_hi = limit_hi;
}

void XsecNorm::SetDimensions(int num_top, int num_reac)
{
    ntop = num_top;
    nreac = num_reac;
}

void XsecNorm::SetDimensions(const std::vector<int>& dim)
{
    m_dimensions = dim;
}

void XsecNorm::Print(bool print_bins) const
{
    std::cout << TAG << "Name: " <<  m_name << std::endl
              << TAG << "Nominal: " << m_nominal << std::endl
              << TAG << "Step: " << m_step << std::endl
              << TAG << "Limits: [" << m_limit_lo << "," << m_limit_hi << "]" << std::endl
              << TAG << "Splines: " << v_splines.size() << std::endl
              << TAG << "Fixed: " << m_fixed << std::endl;

    std::cout << TAG << "Topologies:";
    for(const auto& topo : m_topologies)
        std::cout << " " << topo;
    std::cout << std::endl;

    std::cout << TAG << "Definition:";
    std::cout << std::endl;
    for(int i=0;i<m_definition.size();i++)
        std::cout << TAG << m_definition.at(i) << ": [" << m_def_value.at(i).at(0) << "," <<m_def_value.at(i).at(1) << "]" << std::endl;

    //if(print_bins)
    //    bin_manager.Print();
}

bool XsecNorm::IsApplied(int index, int topo, double val) const
{
    if (index>=m_definition.size()) {
      std::cout << TAG << "[ERROR]: Out of ranges "<< "Name: " <<  m_name << " size: "<<m_definition.size()<<" index = "<<index<<std::endl;
      return false;
    }
    
    bool trueTopo=false;
    for (int i=0;i<m_topologies.size();i++) {
      if (topo==m_topologies.at(i)) trueTopo=true;
    }
 
    if (trueTopo && val>=m_def_value.at(index).at(0) && val<m_def_value.at(index).at(1)) return true;

    return false;
}
