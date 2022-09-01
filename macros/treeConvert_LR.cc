#include</hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/libReadoaAnalysis/ND__NRooTrackerVtx.h>
#include <TTree.h>
#include <TFile.h>

void
treeConvert_LR(TString fname="/scratch/kmtsui/runOutput_run2_4/run8a_aa.root")
{

  Int_t i, j;
  Int_t nevents;

  //gSystem->Load("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/libReadoaAnalysis/libReadoaAnalysis.so");
  
  TTree  *tn;
  ND::NRooTrackerVtx *nvtx;

  TFile f(fname);
  tn = (TTree *)(f.Get("neuttree"));

  nvtx = new ND::NRooTrackerVtx();
  tn->SetBranchAddress("NRooTrackerVtx",&nvtx);

  nevents = tn->GetEntries();nevents=10;

  //exit();

}
  
  
  
