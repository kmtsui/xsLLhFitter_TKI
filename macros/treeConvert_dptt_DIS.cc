/******************************************************

Code to convert a HL2 tree into the format required 
for the fitting code. C

Can't simply read HL2 tree directly since we don't
know what variables will be the tree

Author: Stephen Dolan
Date Created: November 2015

******************************************************/


#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <assert.h>

#include <TCanvas.h>
#include <TH1F.h>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TMath.h>


using namespace std;

bool extraCuts=false;

int treeConvert_dptt_DIS(TString inFileName="/hepstore/kmtsui/T2K/work/highland2/myNumuCCMultiPiAnalysis/v0r0/output/*.root", TString inTreeName="default",  TString inTreeName_T="truth", TString outFileName="../inputs/neutRun8_new_Bin2_50up_Bin5_50down_Bin6_50up_Bin7_50down_Bin8_50down.root",
                Float_t sigWeight=11.53/62.7/*neut Run8 to data POT*/, Int_t EvtStart=0, Int_t EvtEnd=0,
                TString D1NameRec="recDeltapTT_3TPCAllSec", TString D1NameTrue="truthDeltapTT_vtx", TString D1NameTrue_T="truthDeltapTT_vtx",
                TString D2NameRec="selmu_mom", TString D2NameTrue="selmu_truemom", TString D2NameTrue_T="truelepton_mom")
{
  // You need to provide the number of branches in your HL2 tree
  // And the accum_level you want to cut each one at to get your selected events
  // i.e choosing n in accum_level[0][branch_i]>n
  const int nbranches = 10;
  //const int accumToCut[nbranches] = {5,6,7,7,6,4,3,5,7,6};
  const int accumToCut[nbranches] =   {12,12,12,12,12,12,12,9,9,9};
  bool data = false;

  //inFileName="/scratch/kmtsui/prod6D/run*.root";
  //sigWeight=11.53/60.33;
  //outFileName="../inputs/neut_prod6D_dataPOT.root";

  //inFileName="/scratch/kmtsui/runOutput_genie/run*.root";
  //sigWeight=11.53/128.21;
  //outFileName="../inputs/genie_Run2_4_dataPOT.root";

  //inFileName="/hepstore/kmtsui/T2K/work/highland2/myNumuCCMultiPiAnalysis/v0r0/output_data/run*.root";
  //sigWeight=1.;
  //outFileName="../inputs/data_tree.root";
  //data = true;

  //inFileName="/scratch/kmtsui/runOutput_run2_4/run2w*.root";
  //sigWeight=11.53/(12.0);
  //outFileName="../inputs/neut_Run2w_dataPOT.root";

  //inFileName="/scratch/kmtsui/runOutput_run2_4/run8*.root";
  //sigWeight=11.53/(62.7);
  //outFileName="../inputs/neut_Run8_dataPOT.root";

  //inFileName="/scratch/kmtsui/runOutput_run2_4/run8w*.root";
  //sigWeight=11.53/(26.405);
  //outFileName="../inputs/neut_Run8w_dataPOT.root";

  inFileName="/scratch/kmtsui/runOutput_run2_4/run*.root";
  sigWeight=11.53/(62.7+109.7);
  outFileName="../inputs/neut_Run2_8_dataPOT_DIS_mult_0.root";
  //outFileName="../inputs/neut_Run2_8_test_mult_2.root";
  //outFileName="../inputs/test_protonFSI.root";

  //inFileName="/scratch/kmtsui/runOutput_run2_4/run*.root";
  //sigWeight=11.53/(62.7+109.7);
  //outFileName="../inputs/neut_dptt_pmom.root";

  /*TFile *fNuWro = new TFile("/scratch/kmtsui/proton_FSI/nuwro_protonFSI.root");
  TH2D* hPTheta_PreFSI_nuwro = (TH2D*)fNuWro->Get("hPTheta_PreFSI");
  TH2D* hPTheta_PostFSI_nuwro = (TH2D*)fNuWro->Get("hPTheta_PostFSI");
  TFile *fNeut = new TFile("/scratch/kmtsui/proton_FSI/neut_protonFSI.root");
  TH2D* hPTheta_PreFSI_neut = (TH2D*)fNeut->Get("hPTheta_PreFSI");
  TH2D* hPTheta_PostFSI_neut = (TH2D*)fNeut->Get("hPTheta_PostFSI");

  hPTheta_PostFSI_neut->Divide(hPTheta_PreFSI_neut);
  hPTheta_PostFSI_nuwro->Divide(hPTheta_PreFSI_nuwro);
  hPTheta_PostFSI_nuwro->Divide(hPTheta_PostFSI_neut);
  //hPTheta_PostFSI_nuwro->Draw("colz text");
*/

  TFile *fHist2d = new TFile("/scratch/kmtsui/runOutput_run2_4/hist3d.root");
  TH3D* hist3d[4];
  hist3d[0] = (TH3D*)fHist2d->Get("truthDeltapTT_vs_W_Wppi_topology5");
  hist3d[1] = (TH3D*)fHist2d->Get("truthDeltapTT_vs_W_Wppi_topology6");
  hist3d[2] = (TH3D*)fHist2d->Get("truthDeltapTT_vs_W_Wppi_topology7");
  hist3d[3] = (TH3D*)fHist2d->Get("truthDeltapTT_vs_W_Wppi_topology8");

  TAxis *xaxis = hist3d[0]->GetXaxis();
  TAxis *yaxis = hist3d[0]->GetYaxis();
  TAxis *zaxis = hist3d[0]->GetZaxis();

  outFileName="../inputs/neut_Run2_8_dataPOT_DIS_mult_0.root";
  TFile *fHist2d_DIS = new TFile("/scratch/kmtsui/neut_DIS/hist3d_0.root");
  TH3D* hist3d_DIS[4];
  hist3d_DIS[0] = (TH3D*)fHist2d_DIS->Get("truthDeltapTT_vs_W_Wppi_topology5");
  hist3d_DIS[1] = (TH3D*)fHist2d_DIS->Get("truthDeltapTT_vs_W_Wppi_topology6");
  hist3d_DIS[2] = (TH3D*)fHist2d_DIS->Get("truthDeltapTT_vs_W_Wppi_topology7");
  hist3d_DIS[3] = (TH3D*)fHist2d_DIS->Get("truthDeltapTT_vs_W_Wppi_topology8");


  //inFileName="/scratch/kmtsui/runOutput_run4a/*.root";

  TChain* intree;
  intree = new TChain(inTreeName);
  intree->Add(inFileName);
  TChain* intree_T;
  intree_T = new TChain(inTreeName_T);
  intree_T->Add(inFileName);

  //TFile *infile = new TFile(inFileName);
  //TTree *intree = (TTree*)infile->Get(inTreeName);
  //TTree *intree_T = (TTree*)infile->Get(inTreeName_T);

  TFile *outfile = new TFile(outFileName,"recreate");
  TTree *outtree = new TTree("selectedEvents", "selectedEvents");
  TTree *outtree_T = new TTree("trueEvents", "trueEvents");

  // Declaration of leaf types
  Int_t          accum_level[1500][50];
  Int_t          reaction;
  Int_t          cutBranch=-999;
  Int_t          mectopology;
  Int_t          topology_multiproton;
  Int_t          target;
  Float_t        D1true;
  Float_t        D2true;
  Float_t        D1Reco;
  Float_t        D2Reco;
  Float_t        pMomRec;
  Float_t        pMomRecRange;
  Float_t        pThetaRec;
  Float_t        pMomTrue;
  Float_t        pThetaTrue;
  Float_t        muMomRec;
  Float_t        muMomRecRange;
  Float_t        muThetaRec;
  Float_t        muMomTrue;
  Float_t        muThetaTrue;
  Float_t        muCosThetaRec;
  Float_t        muCosThetaTrue;
  Float_t        pCosThetaRec;
  Float_t        pCosThetaTrue;
  Float_t        piMomRec;
  Float_t        piCosThetaRec;
  Float_t        piMomTrue;
  Float_t        piCosThetaTrue;
  Float_t        RecoNuEnergy=0;
  Float_t        TrueNuEnergy=0;
  Float_t        weight;
  Float_t        TPCAllSecDir[100][3];
  Float_t        TPCAllSecTrueDir[100][3];
  Float_t        TrueProtonStartDir[100][3];
  Float_t        selmu_pos[3];
    Float_t        selmu_dir[3];
  Float_t        TPCAllSecMom[100];       
  Float_t        TPCAllSecPosStart[100][3];
  Float_t        TPCAllSecPrPidLik[100];
  Float_t        TPCAllSecPiPidLik[100];
  Float_t        selmu_fgd_V11;
  Float_t        selmu_fgd_V33;
  Float_t        selmu_fgd_V55;
  Float_t        selmu_fgd_V77;
  Float_t        selmu_fgd_VLayer;
  Float_t        recTargetInvMass;
  Float_t        recpTsum;
  Float_t        recpLsum;
  Float_t        recpN;
  Float_t        recdeltaAlphaT;
    Float_t        selmu_mom[1500];
    Float_t        TPCProtonPosStart[100][3];
    Float_t        PosPionPosStart[100][3];
    Float_t        TPCProtonDir[100][3];
    Float_t        PosPionDir[100][3];
    Float_t        TPCProtonMom[1500][100];
    Float_t        PosPionMom[1500][100];
    Int_t          NECalIsoTrack;
    Float_t        ECalIsoTrackEMEnergy;
    Float_t        ECalIsoTrackPIDMipEm;
    Int_t          ECalIsoTrackMostUpStreamLayerHit;
    Int_t          topology_multipi;
    //Int_t          target;
    Float_t recDeltapTT_Toy[1500];
    Float_t recpN_Toy[1500];
    Float_t recdeltaAlphaT_Toy[1500];
    Int_t NPosPion;
    Int_t NTPCproton;
    Int_t NFGDSec;
    Int_t NNegPion;
    Int_t NFGDPi;
    Int_t NME;
    Float_t        truthDeltapTT_vtx;
    Float_t        truthpN;
    Float_t        truthdeltaAlphaT;
    Float_t        selmu_truedir[3];
    Float_t        truepi_mom;
    Float_t        truepi_dir[3];
    Float_t        trueProtonStartDir[100][3];
    Float_t        truelepton_mom;
    Float_t        trueProtonMom[100];


  Int_t          topology;
  Int_t          reaction_T;
  Int_t          mectopology_T;
  Float_t        D1true_T;
  Float_t        D2true_T; 
  Float_t        muMomTrue_T;
  Float_t        pMomTrue_T;
  Float_t        piMomTrue_T;
  Float_t        muCosThetaTrue_T;
  Float_t        pCosThetaTrue_T;
  Float_t        piCosThetaTrue_T;
  Float_t        TrueNuEnergy_T;
  Float_t        weight_T=1.0;
  Int_t          topology_T;
  Float_t        truelepton_dir_T[3];
  Float_t        truepi_dir_T[3];
  Float_t        TrueProtonStartDir_T[100][3];
  Float_t        truelepton_mom_T;
  Float_t        truepi_mom_T;
  Float_t        TrueProtonMom_T[100];
  Int_t          cutBranch_T=-1;

    intree->SetBranchAddress("PosPionPosStart", &PosPionPosStart);
    intree->SetBranchAddress("PosPionDir", &PosPionDir);
    intree->SetBranchAddress("PosPionMom", &PosPionMom);
    intree->SetBranchAddress("TPCProtonPosStart", &TPCProtonPosStart);
    intree->SetBranchAddress("TPCProtonDir", &TPCProtonDir);
    intree->SetBranchAddress("TPCProtonMom", &TPCProtonMom);

    intree->SetBranchAddress("NECalIsoTrack", &NECalIsoTrack);
    intree->SetBranchAddress("ECalIsoTrackEMEnergy", &ECalIsoTrackEMEnergy);
    intree->SetBranchAddress("ECalIsoTrackPIDMipEm", &ECalIsoTrackPIDMipEm);
    intree->SetBranchAddress("ECalIsoTrackMostUpStreamLayerHit", &ECalIsoTrackMostUpStreamLayerHit);

    intree->SetBranchAddress("topology_multipi", &topology_multipi);
    intree->SetBranchAddress("target", &target);
    intree->SetBranchAddress("NPosPion", &NPosPion);
    intree->SetBranchAddress("NTPCproton", &NTPCproton);
    intree->SetBranchAddress("NFGDSec", &NFGDSec);
    intree->SetBranchAddress("NNegPion", &NNegPion);
    intree->SetBranchAddress("NFGDPi", &NFGDPi);
    intree->SetBranchAddress("NME", &NME);

    intree->SetBranchAddress("recDeltapTT_Toy", &recDeltapTT_Toy);
    intree->SetBranchAddress("recpN_Toy", &recpN_Toy);
    intree->SetBranchAddress("recdeltaAlphaT_Toy", &recdeltaAlphaT_Toy);
    intree->SetBranchAddress("truthDeltapTT_vtx", &truthDeltapTT_vtx);
    intree->SetBranchAddress("truthpN", &truthpN);
    intree->SetBranchAddress("truthdeltaAlphaT", &truthdeltaAlphaT);

    intree->SetBranchAddress("selmu_truedir", &selmu_truedir);
    intree->SetBranchAddress("truepi_dir", &truepi_dir);
    intree->SetBranchAddress("trueProtonStartDir", &trueProtonStartDir);
    intree->SetBranchAddress("truelepton_mom", &truelepton_mom);
    intree->SetBranchAddress("truepi_mom", &truepi_mom);
    intree->SetBranchAddress("trueProtonMom", &trueProtonMom);

  intree->SetBranchAddress("accum_level", &accum_level);
  intree->SetBranchAddress("reaction", &reaction);
  intree->SetBranchAddress("mectopology", &mectopology);
  intree->SetBranchAddress("topology_multiproton", &topology_multiproton);
  intree->SetBranchAddress("target", &target);
  //intree->SetBranchAddress(D1NameTrue, &D1true);
  //intree->SetBranchAddress(D2NameTrue, &D2true);
  //intree->SetBranchAddress(D1NameRec, &D1Reco);
  //intree->SetBranchAddress(D2NameRec, &D2Reco);
  //intree->SetBranchAddress("TPCAllSecMom", &pMomRec);
  //intree->SetBranchAddress("TPCAllSec_fgdX", &pMomRecRange);
  //intree->SetBranchAddress("TPCAllSecDir" ,TPCAllSecDir);
  //intree->SetBranchAddress("TPCAllSecTrueMom" ,&pMomTrue);
  //intree->SetBranchAddress("TPCAllSecTrueDir" ,TPCAllSecTrueDir);
  //intree->SetBranchAddress("selmu_mom", &muMomRec);
    intree->SetBranchAddress("selmu_pos", &selmu_pos);
    intree->SetBranchAddress("selmu_dir", &selmu_dir);
    intree->SetBranchAddress("selmu_mom", &selmu_mom);
  intree->SetBranchAddress("selmu_fgd_x", &muMomRecRange);
  intree->SetBranchAddress("selmu_costheta", &muThetaRec);
  //intree->SetBranchAddress("truelepton_mom", &muMomTrue);
  intree->SetBranchAddress("truelepton_costheta", &muCosThetaTrue);
  //intree->SetBranchAddress("nu_trueE", &RecoNuEnergy);
  intree->SetBranchAddress("nu_trueE", &TrueNuEnergy);
  intree->SetBranchAddress("weight_syst_total", &weight);
  intree->SetBranchAddress("TPCAllSecMom", &TPCAllSecMom);
  intree->SetBranchAddress("TPCAllSecPosStart", &TPCAllSecPosStart);
  intree->SetBranchAddress("TPCAllSecPrPidLik", &TPCAllSecPrPidLik);
  intree->SetBranchAddress("TPCAllSecPiPidLik", &TPCAllSecPiPidLik);
  intree->SetBranchAddress("selmu_fgd_V11", &selmu_fgd_V11);
  intree->SetBranchAddress("selmu_fgd_V33", &selmu_fgd_V33);
  intree->SetBranchAddress("selmu_fgd_V55", &selmu_fgd_V55);
  intree->SetBranchAddress("selmu_fgd_V77", &selmu_fgd_V77);
  intree->SetBranchAddress("selmu_fgd_VLayer", &selmu_fgd_VLayer);
  intree->SetBranchAddress("recTargetInvMass", &recTargetInvMass);
  intree->SetBranchAddress("recpTsum", &recpTsum);
  intree->SetBranchAddress("recpLsum", &recpLsum);
  intree->SetBranchAddress("recpN", &recpN);
  intree->SetBranchAddress("recdeltaAlphaT", &recdeltaAlphaT);

  outtree->Branch("reaction", &reaction, "reaction/I");
  outtree->Branch("cutBranch", &cutBranch, "cutBranch/I");
  outtree->Branch("mectopology", &mectopology, "mectopology/I");
  outtree->Branch("topology_multiproton", &topology_multiproton,"topology_multiproton/I");
  outtree->Branch("topology", &topology,"topology/I");
  outtree->Branch("target", &target,"target/I");
  outtree->Branch("D1True", &D1true, ("D1True/F"));
  outtree->Branch("D1Rec", &D1Reco, ("D1Rec/F"));
  outtree->Branch("D2True", &D2true, ("D2True/F"));
  outtree->Branch("D2Rec", &D2Reco, ("D2Rec/F"));
  outtree->Branch("muMomRec", &muMomRec, ("muMomRec/F"));
  outtree->Branch("muMomTrue", &muMomTrue, ("muMomTrue/F"));
  outtree->Branch("muCosThetaRec", &muCosThetaRec, ("muCosThetaRec/F"));
  outtree->Branch("muCosThetaTrue", &muCosThetaTrue, ("muCosThetaTrue/F"));
  outtree->Branch("pMomRec", &pMomRec, ("pMomRec/F"));
  outtree->Branch("pMomTrue", &pMomTrue, ("pMomTrue/F"));
  outtree->Branch("pCosThetaRec", &pCosThetaRec, ("pCosThetaRec/F"));
  outtree->Branch("pCosThetaTrue", &pCosThetaTrue, ("pCosThetaTrue/F"));
  outtree->Branch("piMomRec", &piMomRec, ("piMomRec/F"));
  outtree->Branch("piCosThetaRec", &piCosThetaRec, ("piCosThetaRec/F"));
  outtree->Branch("piMomTrue", &piMomTrue, ("piMomTrue/F"));
  outtree->Branch("piCosThetaTrue", &piCosThetaTrue, ("piCosThetaTrue/F"));
  outtree->Branch("Enureco", &RecoNuEnergy, "Enureco/F");
  outtree->Branch("Enutrue", &TrueNuEnergy, "Enutrue/F");
  outtree->Branch("weight", &weight, "weight/F");

  intree_T->SetBranchAddress("reaction", &reaction_T);
  intree_T->SetBranchAddress("mectopology", &mectopology_T);
  intree_T->SetBranchAddress("topology_multiproton", &topology_multiproton);
  intree_T->SetBranchAddress("topology_multipi", &topology_T);
  intree_T->SetBranchAddress("target", &target);
  intree_T->SetBranchAddress(D1NameTrue_T, &D1true_T);
  intree_T->SetBranchAddress(D2NameTrue_T, &D2true_T);
  //intree_T->SetBranchAddress("TrueProtonMom" ,&pMomTrue_T);
  //intree_T->SetBranchAddress("TrueProtonStartDir" ,TrueProtonStartDir);
  intree_T->SetBranchAddress("truelepton_costheta", &muCosThetaTrue_T);
  intree_T->SetBranchAddress("nu_trueE", &TrueNuEnergy_T);
  intree_T->SetBranchAddress("weight", &weight_T);
  intree_T->SetBranchAddress("truelepton_dir", &truelepton_dir_T);
  intree_T->SetBranchAddress("truepi_dir", &truepi_dir_T);
  intree_T->SetBranchAddress("TrueProtonStartDir", &TrueProtonStartDir_T);
  intree_T->SetBranchAddress("truelepton_mom", &truelepton_mom_T);
  intree_T->SetBranchAddress("truepi_mom", &truepi_mom_T);
  intree_T->SetBranchAddress("TrueProtonMom", &TrueProtonMom_T);

  outtree_T->Branch("reaction", &reaction_T, "reaction/I");
  outtree_T->Branch("mectopology", &mectopology_T, "mectopology/I");
  outtree_T->Branch("topology_multiproton", &topology_multiproton,"topology_multiproton/I");
  outtree_T->Branch("topology", &topology_T,"topology/I");
  outtree_T->Branch("target", &target,"target/I");
  outtree_T->Branch("D1True", &D1true_T, ("D1True/F"));
  outtree_T->Branch("D2True", &D2true_T, ("D2True/F"));
  outtree_T->Branch("muMomTrue", &muMomTrue_T, ("muMomTrue/F"));
  outtree_T->Branch("muCosThetaTrue", &muCosThetaTrue_T, ("muCosThetaTrue/F"));
  outtree_T->Branch("pMomTrue", &pMomTrue_T, ("pMomTrue/F"));
  outtree_T->Branch("pCosThetaTrue", &pCosThetaTrue_T, ("pCosThetaTrue/F"));
  outtree_T->Branch("piMomTrue", &piMomTrue_T, ("piMomTrue/F"));
  outtree_T->Branch("piCosThetaTrue", &piCosThetaTrue_T, ("piCosThetaTrue/F"));
  outtree_T->Branch("Enutrue", &TrueNuEnergy_T, "Enutrue/F");
  outtree_T->Branch("weight", &weight_T, "weight/F");
  outtree_T->Branch("cutBranch", &cutBranch_T, "cutBranch/I");


  Long64_t nentries = intree->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  if(EvtEnd!=0) nentries=EvtEnd;
  int passCount=0;
  std::cout<<"nentries = "<<nentries<<std::endl;
  for (Long64_t jentry=EvtStart; jentry<nentries;jentry++) {
//std::cout<<"jentry = "<<jentry<<std::endl;
      if ( jentry%10000==0 ) { 
	cout << "."  ; 
	cout.flush();
      }
    nb = intree->GetEntry(jentry); nbytes += nb;
    passCount=0;
    RecoNuEnergy=TrueNuEnergy;
    //pCosThetaRec   = TPCAllSecDir[0][2]//TMath::Cos(pThetaRec);
    muMomRec = selmu_mom[0];
    muCosThetaRec  = selmu_dir[2];
    pMomRec = TPCProtonMom[0][0];
    pCosThetaRec = TPCProtonDir[0][2];
    piMomRec = PosPionMom[0][0];
    piCosThetaRec = PosPionDir[0][2];

    muMomTrue = truelepton_mom;
    muCosThetaTrue  = selmu_truedir[2];
    pMomTrue = trueProtonMom[0];
    pCosThetaTrue = trueProtonStartDir[0][2];
    piMomTrue = truepi_mom;
    piCosThetaTrue = truepi_dir[2];
    //muCosThetaTrue = TMath::Cos(muThetaTrue);
    int branches_passed[10]={0};
    //if (recpN<0||recpN>1000||recdeltaAlphaT<0||recdeltaAlphaT>180) continue;
    //for(int i=0; i<nbranches; i++){
      bool passed = 0;
      bool commonVtxCut = (fabs(selmu_pos[0]-TPCProtonPosStart[0][0])<50&&fabs(selmu_pos[0]-PosPionPosStart[0][0])<50&&fabs(TPCProtonPosStart[0][0]-PosPionPosStart[0][0])<50&&fabs(selmu_pos[1]-TPCProtonPosStart[0][1])<50&&fabs(selmu_pos[1]-PosPionPosStart[0][1])<50&&fabs(TPCProtonPosStart[0][1]-PosPionPosStart[0][1])<50&&fabs(selmu_pos[2]-TPCProtonPosStart[0][2])<30&&fabs(selmu_pos[2]-PosPionPosStart[0][2])<30&&fabs(TPCProtonPosStart[0][2]-PosPionPosStart[0][2])<30);
//std::cout<<"commonVtxCut"<<std::endl;
     bool kpsCut = (selmu_dir[2]>0.342&&PosPionDir[0][2]>0.342&&TPCProtonDir[0][2]>0.342&&selmu_mom[0]>225&&selmu_mom[0]<7700&&PosPionMom[0][0]>135&&PosPionMom[0][0]<1320&&TPCProtonMom[0][0]>405&&TPCProtonMom[0][0]<1320);
     //bool kpsCut = (selmu_dir[2]>0.342&&PosPionDir[0][2]>0.342&&TPCProtonDir[0][2]>0.342&&selmu_mom[0]>250&&selmu_mom[0]<7000&&PosPionMom[0][0]>150&&PosPionMom[0][0]<1200&&TPCProtonMom[0][0]>450&&TPCProtonMom[0][0]<1200);
//std::cout<<"kpsCut"<<std::endl;
      bool ecalpi0veto = (!(NECalIsoTrack==1&&ECalIsoTrackEMEnergy>30&&ECalIsoTrackPIDMipEm>0&&ECalIsoTrackMostUpStreamLayerHit<=5));
//std::cout<<"ecalpi0veto"<<std::endl;
      /*if(accum_level[0][1]>6) {
          if (NPosPion==1&&NTPCproton>0&&NFGDSec==0&&commonVtxCut&&ecalpi0veto) {passed=1;cutBranch=0;passCount++;branches_passed[0]++;}
      }
      if(accum_level[0][6]>6) {
          if (NPosPion==1&&NTPCproton>0&&NNegPion+NFGDPi+NME==1&&ecalpi0veto) {passed=1;cutBranch=1;passCount++;branches_passed[1]++;}
          if (NPosPion==1&&NTPCproton>0&&NNegPion+NFGDPi+NME==0&&!ecalpi0veto) {passed=1;cutBranch=2;passCount++;branches_passed[2]++;}
          if (NPosPion>0&&NTPCproton>0&&(NPosPion>1||NNegPion+NFGDPi+NME>0)&&!ecalpi0veto) {passed=1;cutBranch=3;passCount++;branches_passed[3]++;}
          if (NPosPion>0&&NTPCproton>0&&(NPosPion>1||NNegPion+NFGDPi+NME>1)&&ecalpi0veto) {passed=1;cutBranch=4;passCount++;branches_passed[4]++;}
      }*/
      if (accum_level[0][2]>=11) {passed=1;cutBranch=0;passCount++;branches_passed[0]++;}
      if (accum_level[0][11]>8) {passed=1;cutBranch=1;passCount++;branches_passed[1]++;}
      if (accum_level[0][12]>8) {passed=1;cutBranch=2;passCount++;branches_passed[2]++;}
      if (accum_level[0][13]>8) {passed=1;cutBranch=3;passCount++;branches_passed[3]++;}
      if (accum_level[0][14]>8) {passed=1;cutBranch=4;passCount++;branches_passed[4]++;}
//std::cout<<"accum_level"<<std::endl;
      if(!passed) continue;
      D2true=0.1;D2Reco=0.1;
      D1Reco=recDeltapTT_Toy[0];D1true=truthDeltapTT_vtx;
      D2Reco=TPCProtonMom[0][0];

      weight = weight*sigWeight;

    if (D1Reco>700) D1Reco=699;
    if (D1Reco<-700) D1Reco=-699;
    //if (D2Reco>2000) D2Reco=2000;
    //if (!kpsCut && fabs(D1Reco)<700) D1Reco=750;
    if (cutBranch==0 && accum_level[0][2]==11 && fabs(D1Reco)<700) {D1Reco=750;D2Reco=500;}

    if (!data) if (topology_multipi<0) continue;//-1 no truth
    if (topology_multipi<7) topology=topology_multipi; //0 CC0pi, 1 CC1pi0p, 2 CC1pi1p, 3 CC1piNp, 4 CC2pi+, 5 CC1pi+1pi-, 6 CC1pi+Npi0
    else if (topology_multipi==8) topology=7; //7 CC-other-Npi0 
    else if (topology_multipi==10) topology=8; //8 CC-other-0pi0 
    else if (topology_multipi==999) topology=9; //9 backg(NC+antinu+nue)
    else if (topology_multipi==7)   topology=10; //10 OOFV

      double weightMult = 1.;
  
    double Wppi = sqrt((sqrt(938*938+TPCProtonMom[0][0]*TPCProtonMom[0][0])+sqrt(139.6*139.6+PosPionMom[0][0]*PosPionMom[0][0]))*(sqrt(938*938+TPCProtonMom[0][0]*TPCProtonMom[0][0])+sqrt(139.6*139.6+PosPionMom[0][0]*PosPionMom[0][0]))-(TPCProtonMom[0][0]*TPCProtonDir[0][0]+PosPionMom[0][0]*PosPionDir[0][0])*(TPCProtonMom[0][0]*TPCProtonDir[0][0]+PosPionMom[0][0]*PosPionDir[0][0])-(TPCProtonMom[0][0]*TPCProtonDir[0][1]+PosPionMom[0][0]*PosPionDir[0][1])*(TPCProtonMom[0][0]*TPCProtonDir[0][1]+PosPionMom[0][0]*PosPionDir[0][1])-(TPCProtonMom[0][0]*TPCProtonDir[0][2]+PosPionMom[0][0]*PosPionDir[0][2])*(TPCProtonMom[0][0]*TPCProtonDir[0][2]+PosPionMom[0][0]*PosPionDir[0][2]));
    double trueW = sqrt((TrueNuEnergy+938-sqrt(truelepton_mom*truelepton_mom+106*106))*(TrueNuEnergy+938-sqrt(truelepton_mom*truelepton_mom+106*106))-(TrueNuEnergy*TrueNuEnergy+truelepton_mom*truelepton_mom-2*TrueNuEnergy*truelepton_mom*muCosThetaTrue));

      Int_t binx = xaxis->FindBin(Wppi);
      Int_t biny = yaxis->FindBin(trueW);
      Int_t binz = zaxis->FindBin(truthDeltapTT_vtx);
      if (topology>=5&&topology<=8)
        if (hist3d[topology-5]->GetBinContent(binx,biny,binz)>0) 
          weightMult=hist3d_DIS[topology-5]->GetBinContent(binx,biny,binz)/hist3d[topology-5]->GetBinContent(binx,biny,binz);
      weight*=weightMult;

//std::cout<<"topology_multipi"<<std::endl;
    if (topology==2||topology==3) {
      //D2Reco=pMomRec;D2true=trueProtonMom[0];      
      if (!(selmu_truedir[2]>0.342&&truepi_dir[2]>0.342&&trueProtonStartDir[0][2]>0.342&&truelepton_mom>250&&truelepton_mom<7000&&truepi_mom>150&&truepi_mom<1200&&trueProtonMom[0]>450&&trueProtonMom[0]<1200)) topology=11;//OOPS

    }
    if (D1true>700) D1true=699;
    if (D1true<-700) D1true=-699;
//std::cout<<"OOPS"<<std::endl;
    if (topology!=2&&topology!=3&&topology!=11) D1true=0.;
    if (topology==2||topology==3) if (target==1) D2true=1.5;
    if (topology==5) D2true=15;
    if (topology==6) D2true=25;
    if (topology==7) D2true=35;
    if (topology==8) D2true=45;
    if (topology==11) D2true=55;
    //if (topology==11) D2true=1250;

    /*if ((topology==2||topology==3)&&D1true>-100&&D1true<100) weight*=1.5;
    if (topology==5) weight*=0.5;
    if (topology==6) weight*=1.5;
    if (topology==7) weight*=0.5;
    if (topology==8) weight*=1.5;*/

    outtree->Fill();

    if(passCount>1){
      printf("***Warning: More than one cut branch passed***\n");
      for(int j=0;j<10;j++){
        if(branches_passed[j]==1) printf("branch %d passed ...",j);
      }
      printf("\n");
    }
  }

  Long64_t nentries_T = intree_T->GetEntries();
  Long64_t nbytes_T = 0, nb_T = 0;
  std::cout<<"nentries_T = "<<nentries_T<<std::endl;
  for (Long64_t jentry=0; jentry<nentries_T;jentry++) {
    nb_T = intree_T->GetEntry(jentry); nbytes_T += nb_T;
    if (topology_T<7) topology_T=topology_T; //0 CC0pi, 1 CC1pi0p, 2 CC1pi1p, 3 CC1piNp, 4 CC2pi+, 5 CC1pi+1pi-, 6 CC1pi+Npi0
    else if (topology_T==8) topology_T=7; //7 CC-other-Npi0 
    else if (topology_T==10) topology_T=8; //8 CC-other-0pi0 
    else if (topology_T==999) topology_T=9; //9 backg(NC+antinu+nue)
    else if (topology_T==7)   topology_T=10; //10 OOFV
//std::cout<<"topology_multipi"<<std::endl;
    if (topology_T==2||topology_T==3) {
                  if (!(truelepton_dir_T[2]>0.342&&truepi_dir_T[2]>0.342&&TrueProtonStartDir_T[0][2]>0.342&&truelepton_mom_T>250&&truelepton_mom_T<7000&&truepi_mom_T>150&&truepi_mom_T<1200&&TrueProtonMom_T[0]>450&&TrueProtonMom_T[0]<1200)) topology_T=11;//OOPS
    }
    if (D1true_T>700) D1true_T=699;
    if (D1true_T<-700) D1true_T=-699;
    if (topology_T!=2&&topology_T!=3&&topology_T!=11) D1true_T=0.;
    D2true_T=0.1;
    if (topology_T==2||topology_T==3) if (target==1) D2true_T=1.5;
    //D2true_T=TrueProtonMom_T[0];
    if (topology_T==5) D2true_T=15;
    if (topology_T==6) D2true_T=25;
    if (topology_T==7) D2true_T=35;
    if (topology_T==8) D2true_T=45;
    if (topology_T==11) D2true_T=55;

    muMomTrue_T=truelepton_mom_T;
    muCosThetaTrue_T=truelepton_dir_T[2];
    pMomTrue_T=TrueProtonMom_T[0];
    pCosThetaTrue_T=TrueProtonStartDir_T[0][2];
    piMomTrue_T=truepi_mom_T;
    piCosThetaTrue_T=truepi_dir_T[2];

    //if (topology_T==11) D2true_T=1250;
    //pCosThetaTrue_T=TrueProtonStartDir[0][2];
    //weight_T = weight_T*11.53/36.7;//weight_T = weight_T*11.53/62.7;

    /*double weightFSI = 1.;
    int pbin = (TrueProtonMom_T[0]-450.)/75.;
    int thetabin = (TrueProtonStartDir_T[0][2]-0.342)/0.0658;
    if (pbin<hPTheta_PostFSI_nuwro->GetNbinsY()&&thetabin<hPTheta_PostFSI_nuwro->GetNbinsX()) weightFSI=hPTheta_PostFSI_nuwro->GetBinContent(thetabin+1,pbin+1);
    //weight_T = weight_T*weightFSI;*/

    weight_T = weight_T*sigWeight;

    cutBranch_T=-1;
    outtree_T->Fill();
  }
  
  printf("***Output Rec Tree: ***\n");
  outtree->Print();
  printf("***Output True Tree: ***\n");
  outtree_T->Print();
  outfile->Write();

  //delete infile;
  delete outfile;
  return 0;
}
