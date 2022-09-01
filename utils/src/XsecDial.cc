#include "XsecDial.hh"

XsecDial::XsecDial(const std::string& dial_name)
    : m_name(dial_name)
{
}

XsecDial::XsecDial(const std::string& dial_name, const std::string& fname_binning,
                   const std::string& fname_splines)
    : m_name(dial_name)
{
    SetBinning(fname_binning);
    ReadSplines(fname_splines);
}

void XsecDial::SetBinning(const std::string& fname_binning)
{
    bin_manager.SetBinning(fname_binning);
    nbins = bin_manager.GetNbins();
}

void XsecDial::ReadSplines(const std::string& fname_splines)
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

int XsecDial::GetSplineIndex(int topology, int reaction, double q2) const
{
    const int b = bin_manager.GetBinIndex(std::vector<double>{q2});
    return topology * nreac * nbins + reaction * nbins + b;
}

int XsecDial::GetSplineIndex(const std::vector<int>& var, const std::vector<double>& bin, const std::vector<double>& ptheta, int target) const
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

    if (var[1]==5) idxReco=0; //truth tree
    if (var[0]!=2&&var[0]!=3&&var[0]!=11) idxTrue=2; //stupid bug in spline generation

    /*const int nPMombins = 6;
    const double PMombins[nPMombins+1] = {405,575,700,825,950,1075,1320};
    int idxPmom = 0;
    if (var[1]==0) //signal sample
    {
      if (idxReco<m_dimensions[1]-1){
         double pmom=ptheta[1]; //rec pmom
         for (int i=0;i<nPMombins;i++) {
              if (pmom>=PMombins[i] && pmom<PMombins[i+1]) {idxPmom=i;break;}
         }
      }
    }*/
    
    /*if (var[1]!=0) {
      if (idxReco==0) idxReco=2;
      else if (idxReco==1) idxReco=5;
    }*/
    /*std::string histname = Form("RecBin_%d_trueBin_%d_sample_%d_reac_%d",idxReco,idxTrue,var[1],var[0]);

    for (int i=0;i<v_splines.size();i++) {
      if (histname.compare(v_splines.at(i).GetName())==0) return i;
    }
    return -1;*/
    int indexTarget = 0;
    if (target==1) indexTarget=1;
    int idx = indexTarget+var[0]*2+idxTrue*2*12+idxReco*2*12*m_dimensions[0]+var[1]*2*12*m_dimensions[0]*m_dimensions[1]; //2 targets, 12 topo bins, m_dimensions[0] true bins, m_dimensions[1] reco bins, 5+1 samples
    //int idx = idxPmom+indexTarget*nPMombins+var[0]*nPMombins*2+idxTrue*nPMombins*2*12+idxReco*nPMombins*2*12*m_dimensions[0]+var[1]*nPMombins*2*12*m_dimensions[0]*m_dimensions[1]; //6 pbin, 2 targets, 12 topo bins, m_dimensions[0] true bins, m_dimensions[1] reco bins, 5+1 samples

    /*if (m_name.find("PMOM") != std::string::npos) {
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
    }*/

    return idx;
}

double XsecDial::GetSplineValue(int index, double dial_value) const
{
    if(index >= 0)
        return v_splines.at(index).Eval(dial_value);
    else
        return 1.0;
}

std::string XsecDial::GetSplineName(int index) const
{
    return std::string(v_splines.at(index).GetName());
}

void XsecDial::SetVars(double nominal, double step, double limit_lo, double limit_hi)
{
    m_nominal = nominal;
    m_step = step;
    m_limit_lo = limit_lo;
    m_limit_hi = limit_hi;
}

void XsecDial::SetDimensions(int num_top, int num_reac)
{
    ntop = num_top;
    nreac = num_reac;
}

void XsecDial::SetDimensions(const std::vector<int>& dim)
{
    m_dimensions = dim;
}

void XsecDial::Print(bool print_bins) const
{
    std::cout << TAG << "Name: " <<  m_name << std::endl
              << TAG << "Nominal: " << m_nominal << std::endl
              << TAG << "Step: " << m_step << std::endl
              << TAG << "Limits: [" << m_limit_lo << "," << m_limit_hi << "]" << std::endl
              << TAG << "Splines: " << v_splines.size() << std::endl;

    std::cout << TAG << "Dimensions:";
    for(const auto& dim : m_dimensions)
        std::cout << " " << dim;
    std::cout << std::endl;

    if(print_bins)
        bin_manager.Print();
}
