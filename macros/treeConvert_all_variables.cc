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

bool data = true;

TH2D* h_proton_weight;
TH2D* h_muon_weight;
TH2D* h_pion_weight;
TH1D* h_tki_weight;


enum AnalysisVariable
{
    dptt      = 0,
    pN        = 1,
    dalphaT   = 2,
    W         = 3,
    costhA    = 4,
    phiA      = 5
};

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
  Float_t        weight_syst[16];
  Int_t          NInts;
  Int_t          IntType[100];
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

    Float_t        recW_Toy[1500];
    Float_t        recCosThetaA_Toy[1500];
    Float_t        recPhiA_Toy[1500];
    Float_t        truthW;
    Float_t        truthCosThetaA;
    Float_t        truthPhiA;

  Float_t dptt_rec;
  Float_t dptt_truth;
  Float_t pN_rec;
  Float_t pN_truth;
  Float_t daT_rec;
  Float_t daT_truth;

  Float_t weight_pionSI;

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

void SetInTreeBranches (TChain* intree) {
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

    intree->SetBranchAddress("recW_Toy", &recW_Toy);
    intree->SetBranchAddress("recCosThetaA_Toy", &recCosThetaA_Toy);
    intree->SetBranchAddress("recPhiA_Toy", &recPhiA_Toy);
    intree->SetBranchAddress("truthW", &truthW);
    intree->SetBranchAddress("truthCosThetaA", &truthCosThetaA);
    intree->SetBranchAddress("truthPhiA", &truthPhiA);

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
  intree->SetBranchAddress("weight_syst", &weight_syst);
  //intree->SetBranchAddress("NInts", &NInts);
  //intree->SetBranchAddress("IntType", &IntType);
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
}

void SetOutTreeBranches (TTree *outtree) {
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
  outtree->Branch("weight_pionSI", &weight_pionSI, "weight_pionSI/F");
  outtree->Branch("NTPCproton",&NTPCproton,"NTPCproton/I");
  outtree->Branch("dptt_rec", &dptt_rec, "dptt_rec/F");
  outtree->Branch("dptt_truth", &dptt_truth, "dptt_truth/F");
  outtree->Branch("pN_rec", &pN_rec, "pN_rec/F");
  outtree->Branch("pN_truth", &pN_truth, "pN_truth/F");
  outtree->Branch("daT_rec", &daT_rec, "daT_rec/F");
  outtree->Branch("daT_truth", &daT_truth, "daT_truth/F");
}
  
void  SetInTruthTreeBranches(TChain* intree_T){
  intree_T->SetBranchAddress("reaction", &reaction_T);
  intree_T->SetBranchAddress("mectopology", &mectopology_T);
  intree_T->SetBranchAddress("topology_multiproton", &topology_multiproton);
  intree_T->SetBranchAddress("topology_multipi", &topology_T);
  intree_T->SetBranchAddress("target", &target);
  //intree_T->SetBranchAddress(D1NameTrue_T, &D1true_T);
  //intree_T->SetBranchAddress(D2NameTrue_T, &D2true_T);
  intree_T->SetBranchAddress("truthDeltapTT_vtx",&truthDeltapTT_vtx);
  intree_T->SetBranchAddress("truthpN",&truthpN);
  intree_T->SetBranchAddress("truthdeltaAlphaT",&truthdeltaAlphaT);
  intree_T->SetBranchAddress("truthW", &truthW);
  intree_T->SetBranchAddress("truthCosThetaA", &truthCosThetaA);
  intree_T->SetBranchAddress("truthPhiA", &truthPhiA);
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
}

void SetOutTruthTreeBranches(TTree *outtree_T){
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
  outtree_T->Branch("dptt_truth", &dptt_truth, "dptt_truth/F");
  outtree_T->Branch("pN_truth", &pN_truth, "pN_truth/F");
  outtree_T->Branch("daT_truth", &daT_truth, "daT_truth/F");
}

void SetProtonWeight() {
  TFile* fProtonMom_nom  = new TFile("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/neut/src/output/proton_nominal.root");
  TFile* fProtonMom_fact = new TFile("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/neut/src/output/proton_fact_0.5.root");
  TH2D* hProtonMomTheta_nom  = (TH2D*)fProtonMom_nom->Get("hProtonMomTheta");
  TH2D* hProtonMomTheta_fact = (TH2D*)fProtonMom_fact->Get("hProtonMomTheta");
  hProtonMomTheta_fact->Divide(hProtonMomTheta_nom);
  h_proton_weight = (TH2D*)hProtonMomTheta_fact->Clone();
}

void SetParticleModelWeight(TString fModelPath) {
  TFile* fModel = new TFile(fModelPath);
  TFile* fNom   = new TFile("../inputs/neut_prod6T_dptt.root");

  TTree* tNom = (TTree*)fNom->Get("trueEvents");
  tNom->Draw("pMomTrue:pCosThetaTrue>>truthProton(20,0.342,1,20,450,1200)","topology==2||topology==3");
  TH2D* hNomProton = (TH2D*)gDirectory->Get("truthProton");
  hNomProton->Scale(1./hNomProton->Integral());
  tNom->Draw("piMomTrue:piCosThetaTrue>>truthPion(20,0.342,1,20,150,1200)","topology==2||topology==3");
  TH2D* hNomPion = (TH2D*)gDirectory->Get("truthPion");
  hNomPion->Scale(1./hNomPion->Integral());
  tNom->Draw("muMomTrue:muCosThetaTrue>>truthMuon(20,0.342,1,20,250,7000)","topology==2||topology==3");
  TH2D* hNomMuon = (TH2D*)gDirectory->Get("truthMuon");
  hNomMuon->Scale(1./hNomMuon->Integral());

  TTree* tModel = (TTree*)fModel->Get("TKITree");
  tModel->Draw("pmom:pcostheta>>truthProtonModel(20,0.342,1,20,450,1200)","EvtWght");
  TH2D* hModelProton = (TH2D*)gDirectory->Get("truthProtonModel");
  hModelProton->Scale(1./hModelProton->Integral());
  tModel->Draw("pimom:picostheta>>truthPionModel(20,0.342,1,20,150,1200)","EvtWght");
  TH2D* hModelPion = (TH2D*)gDirectory->Get("truthPionModel");
  hModelPion->Scale(1./hModelPion->Integral());
  tModel->Draw("mumom:mucostheta>>truthMuonModel(20,0.342,1,20,250,7000)","EvtWght");
  TH2D* hModelMuon = (TH2D*)gDirectory->Get("truthMuonModel");
  hModelMuon->Scale(1./hModelMuon->Integral());

  hModelProton->Divide(hNomProton);
  hModelPion->Divide(hNomPion);
  hModelMuon->Divide(hNomMuon);
  h_proton_weight = (TH2D*)hModelProton->Clone();
  h_pion_weight = (TH2D*)hModelPion->Clone();
  h_muon_weight = (TH2D*)hModelMuon->Clone();
}

void SetParticleModelWeight(TString fModelPath, int tkiVar) {
  TFile* fModel = new TFile(fModelPath);
  TFile* fNom;//   = new TFile("../inputs/neut_prod6T_dptt.root");
  TH1D* hTKINom;
  TH1D* hTKIModel;

  switch (tkiVar) {
      case dptt:
        fNom = new TFile("../inputs/neut_prod6T_dptt.root");
        break;
      case pN:
        fNom = new TFile("../inputs/neut_prod6T_pN.root");
        break;
      case dalphaT:
        fNom = new TFile("../inputs/neut_prod6T_dat.root");
        break;
      default:
        printf("***Warning: Not a valid TKI variable***\n");
        break;
  }

  TTree* tNom = (TTree*)fNom->Get("trueEvents");
  tNom->Draw("pMomTrue:pCosThetaTrue>>truthProton(20,0.342,1,20,450,1200)","(topology==2||topology==3)&&target==6");
  TH2D* hNomProton = (TH2D*)gDirectory->Get("truthProton");
  hNomProton->Scale(1./hNomProton->Integral());
  tNom->Draw("piMomTrue:piCosThetaTrue>>truthPion(20,0.342,1,20,150,1200)","(topology==2||topology==3)&&target==6");
  TH2D* hNomPion = (TH2D*)gDirectory->Get("truthPion");
  hNomPion->Scale(1./hNomPion->Integral());
  tNom->Draw("muMomTrue:muCosThetaTrue>>truthMuon(20,0.342,1,20,250,7000)","(topology==2||topology==3)&&target==6");
  TH2D* hNomMuon = (TH2D*)gDirectory->Get("truthMuon");
  hNomMuon->Scale(1./hNomMuon->Integral());

  TTree* tModel = (TTree*)fModel->Get("TKITree");
  tModel->Draw("pmom:pcostheta>>truthProtonModel(20,0.342,1,20,450,1200)");
  TH2D* hModelProton = (TH2D*)gDirectory->Get("truthProtonModel");
  hModelProton->Scale(1./hModelProton->Integral());
  tModel->Draw("pimom:picostheta>>truthPionModel(20,0.342,1,20,150,1200)");
  TH2D* hModelPion = (TH2D*)gDirectory->Get("truthPionModel");
  hModelPion->Scale(1./hModelPion->Integral());
  tModel->Draw("mumom:mucostheta>>truthMuonModel(20,0.342,1,20,250,7000)");
  TH2D* hModelMuon = (TH2D*)gDirectory->Get("truthMuonModel");
  hModelMuon->Scale(1./hModelMuon->Integral());

  switch (tkiVar) {
      case dptt:
        tNom->Draw("D1True>>truthD1(20,-700,700)","(topology==2||topology==3)&&target==6"); 
        hTKINom = (TH1D*)gDirectory->Get("truthD1");
        tModel->Draw("dptt>>D1Model(20,-700,700)");
        hTKIModel = (TH1D*)gDirectory->Get("D1Model");
        hTKINom->Scale(1./hTKINom->Integral());
        hTKIModel->Scale(1./hTKIModel->Integral());
        break;
      case pN:
        tNom->Draw("D1True>>truthD1(20,0,1500)","(topology==2||topology==3)&&target==6"); 
        hTKINom = (TH1D*)gDirectory->Get("truthD1");
        tModel->Draw("pN>>D1Model(20,0,1500)");
        hTKIModel = (TH1D*)gDirectory->Get("D1Model");
        hTKINom->Scale(1./hTKINom->Integral());
        hTKIModel->Scale(1./hTKIModel->Integral());
        break;
      case dalphaT:
        tNom->Draw("D1True>>truthD1(20,0,180)","(topology==2||topology==3)&&target==6"); 
        hTKINom = (TH1D*)gDirectory->Get("truthD1");
        tModel->Draw("dat>>D1Model(20,0,180)");
        hTKIModel = (TH1D*)gDirectory->Get("D1Model");
        hTKINom->Scale(1./hTKINom->Integral());
        hTKIModel->Scale(1./hTKIModel->Integral());
        break;
      default:
        printf("***Warning: Not a valid TKI variable***\n");
        break;
  }

  hModelProton->Divide(hNomProton);
  hModelPion->Divide(hNomPion);
  hModelMuon->Divide(hNomMuon);
  hTKIModel->Divide(hTKINom);
  h_proton_weight = (TH2D*)hModelProton->Clone();
  h_pion_weight = (TH2D*)hModelPion->Clone();
  h_muon_weight = (TH2D*)hModelMuon->Clone();
  h_tki_weight = (TH1D*)hTKIModel->Clone();
}

double GetParticleModelWeight(double pmom, double pcostheta, double pimom, double picostheta, double mumom, double mucostheta){
  double wgt = 1.;

  int ybin = h_proton_weight->GetYaxis()->FindBin(pmom);
  int xbin = h_proton_weight->GetXaxis()->FindBin(pcostheta);
  double weight = h_proton_weight->GetBinContent(xbin,ybin);
  //std::cout<<"weight = "<<weight<<std::endl;
  if (weight>0 && weight<10) wgt*=weight;

  ybin = h_pion_weight->GetYaxis()->FindBin(pimom);
  xbin = h_pion_weight->GetXaxis()->FindBin(picostheta);
  weight = h_pion_weight->GetBinContent(xbin,ybin);
  if (weight>0 && weight<10) wgt*=weight;

  ybin = h_muon_weight->GetYaxis()->FindBin(mumom);
  xbin = h_muon_weight->GetXaxis()->FindBin(mucostheta);
  weight = h_muon_weight->GetBinContent(xbin,ybin);
  if (weight>0 && weight<10) wgt*=weight;

  return wgt;
}

double GetParticleModelWeight(double pmom, double pcostheta, double pimom, double picostheta, double mumom, double mucostheta, double tkiVar){
  double wgt = 1.;

  int ybin = h_proton_weight->GetYaxis()->FindBin(pmom);
  int xbin = h_proton_weight->GetXaxis()->FindBin(pcostheta);
  double weight = h_proton_weight->GetBinContent(xbin,ybin);
  //std::cout<<"weight = "<<weight<<std::endl;
  if (weight>0 && weight<10) wgt*=weight;

  ybin = h_pion_weight->GetYaxis()->FindBin(pimom);
  xbin = h_pion_weight->GetXaxis()->FindBin(picostheta);
  weight = h_pion_weight->GetBinContent(xbin,ybin);
  if (weight>0 && weight<10) wgt*=weight;

  ybin = h_muon_weight->GetYaxis()->FindBin(mumom);
  xbin = h_muon_weight->GetXaxis()->FindBin(mucostheta);
  weight = h_muon_weight->GetBinContent(xbin,ybin);
  if (weight>0 && weight<10) wgt*=weight;

  xbin = h_tki_weight->GetXaxis()->FindBin(tkiVar);
  weight = h_tki_weight->GetBinContent(xbin);
  if (weight>0 && weight<10) wgt*=weight;

  return wgt;
}

std::vector<double> enubins;
std::vector<double> enuWeights;
void SetFluxWeight() {
  TFile* fFluxCov  = new TFile("../inputs/flux_covariance_banff_run1-9d_v2.root");
  TAxis* nd_numu_bins = (TAxis*)fFluxCov->Get("nd5_numode_numu_bins");
  enubins.push_back(nd_numu_bins->GetBinLowEdge(1));
  for(int i=0;i<nd_numu_bins->GetNbins();i++)
    {
      enubins.push_back(nd_numu_bins->GetBinUpEdge(i+1));
      if (i%2==0) enuWeights.push_back(1.5);
      else enuWeights.push_back(1.5);
    }
}

double GetFluxWeight(double enu) {
  enu=enu/1000;//Mev to GeV
  int idx=-1;
  for (int i =0;i<enubins.size()-1;i++) {
    if (enu>enubins[i]&&enu<enubins[i+1]) {
       idx=i; break;
    }
  }
  if (idx>-1) return enuWeights[idx];
  else return 1.0;
}

double GetProtonWeight(double pmom, double pcostheta) {
  int xbin = h_proton_weight->GetXaxis()->FindBin(pmom);
  int ybin = h_proton_weight->GetYaxis()->FindBin(pcostheta);
  double weight = h_proton_weight->GetBinContent(xbin,ybin);
  //std::cout<<"weight = "<<weight<<std::endl;
  if (weight>0 && !std::isnan(weight)) return weight;
  else return 1.;
}

double runweight[7]={0.43/12, 0.36/16.8, 1.59/30.8, 1.70/36.1, 1.79/36.1, 1.58/25.4, 4.15/36.1};
double GetRunWeight(const char* fname){
  std::string filename = fname;
  if (filename.find("run2w")!=std::string::npos) return runweight[0];
  else if (filename.find("run2a")!=std::string::npos) return runweight[1];
  else if (filename.find("run3a")!=std::string::npos) return runweight[2];
  else if (filename.find("run4w")!=std::string::npos) return runweight[3];
  else if (filename.find("run4a")!=std::string::npos) return runweight[4];
  else if (filename.find("run8w")!=std::string::npos) return runweight[5];
  else if (filename.find("run8a")!=std::string::npos) return runweight[6];
  else {
    std::cout<<"wrong run!!!"<<std::endl; return 0.;
  }
}

void RunInTree(TChain* intree, TTree *outtree, int TKIvar, double sigWeight, bool sigOnly = false, bool bkgOnly = false, bool sigSampleOnly = false, bool controlSampleOnly = false) {
  Long64_t nentries = intree->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  int passCount=0;
  std::cout<<"nentries = "<<nentries<<std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    nb = intree->GetEntry(jentry); nbytes += nb;
    if ( jentry%10000==0 ) { 
	    cout << "."  ; 
	    cout.flush();
      //cout<<intree->GetFile()->GetName()<<" weight = "<<GetRunWeight(intree->GetFile()->GetName())<<endl;
    }
    //std::cout<<" jentry = "<<jentry<<std::endl;
    passCount=0;
    RecoNuEnergy=TrueNuEnergy;
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

    dptt_rec = recDeltapTT_Toy[0];
    dptt_truth = truthDeltapTT_vtx;
    pN_rec = recpN_Toy[0];
    pN_truth = truthpN;
    daT_rec = recdeltaAlphaT_Toy[0];
    daT_truth = truthdeltaAlphaT;

    int branches_passed[10]={0};
    bool passed = 0;
    if (accum_level[0][2]>=11) {passed=1;cutBranch=0;passCount++;branches_passed[0]++;}
    if (accum_level[0][11]>8) {passed=1;cutBranch=1;passCount++;branches_passed[1]++;}
    if (accum_level[0][12]>8) {passed=1;cutBranch=2;passCount++;branches_passed[2]++;}
    if (accum_level[0][13]>8) {passed=1;cutBranch=3;passCount++;branches_passed[3]++;}
    if (accum_level[0][14]>8) {passed=1;cutBranch=4;passCount++;branches_passed[4]++;}
    if(!passed) continue;
    D2true=0.1;D2Reco=0.1;
    D2Reco=TPCProtonMom[0][0];

    //weight = weight*sigWeight;
    weight = weight*GetRunWeight(intree->GetFile()->GetName());
    weight_pionSI = weight_syst[13];

    if (!data) if (topology_multipi<0) continue;//-1 no truth
    if (topology_multipi<7) topology=topology_multipi; //0 CC0pi, 1 CC1pi0p, 2 CC1pi1p, 3 CC1piNp, 4 CC2pi+, 5 CC1pi+1pi-, 6 CC1pi+Npi0
    else if (topology_multipi==8) topology=7; //7 CC-other-Npi0 
    else if (topology_multipi==10) topology=8; //8 CC-other-0pi0 
    else if (topology_multipi==999) topology=9; //9 backg(NC+antinu+nue)
    else if (topology_multipi==7)   topology=10; //10 OOFV
    else if (topology_multipi==11)  topology=11; //OOPS
    if (topology==2||topology==3) {   
      if (!(selmu_truedir[2]>0.342&&truepi_dir[2]>0.342&&trueProtonStartDir[0][2]>0.342&&truelepton_mom>250&&truelepton_mom<7000&&truepi_mom>150&&truepi_mom<1200&&trueProtonMom[0]>450&&trueProtonMom[0]<1200)) topology=11;
    }

    if (sigOnly) if (topology!=2&&topology!=3&&topology!=11) continue;
    if (bkgOnly) if (topology==2||topology==3||topology==11) continue;
    if (sigSampleOnly)     if (cutBranch!=0) continue;
    if (controlSampleOnly) if (cutBranch==0) continue;

    switch (TKIvar) {
      case dptt:
        D1Reco=recDeltapTT_Toy[0];D1true=truthDeltapTT_vtx;
        if (D1Reco>700) D1Reco=699;
        if (D1Reco<-700) D1Reco=-699;
        if (cutBranch==0 && accum_level[0][2]==11 && fabs(D1Reco)<700) 
          {
            D1Reco=750;D2Reco=500;
            if (pMomRec>1320&&piMomRec<1320) D2Reco=600;
            else if (pMomRec<1320&&piMomRec>1320) D2Reco=700;
            else if (pMomRec>1320&&piMomRec>1320) D2Reco=800;
          }
        if (D1true>700) D1true=699;
        if (D1true<-700) D1true=-699;
        if (topology!=2&&topology!=3&&topology!=11) D1true=0.1;
        break;
      case pN:
        D1Reco=recpN_Toy[0];D1true=truthpN;
        if (D1Reco>1500) D1Reco=1499; //also act as overflow bins
        if (D1true>1500) D1true=1499;
        if (cutBranch==0 && accum_level[0][2]==11 && D1Reco>0 && D1Reco<1500) {D1Reco=1550;D2Reco=500;}
        if (topology!=2&&topology!=3&&topology!=11) D1true=0.1;
        break;
      case dalphaT:
        D1Reco=recdeltaAlphaT_Toy[0];D1true=truthdeltaAlphaT;
        if (cutBranch==0 && accum_level[0][2]==11) {D1Reco=185;D2Reco=500;}
        if (topology!=2&&topology!=3&&topology!=11) D1true=0.1;
        break;
      case W:
        D1Reco=recW_Toy[0];D1true=truthW;
        if (D1Reco>1800) D1Reco=1799;
        if (D1Reco<1080) D1Reco=1081;
        if (cutBranch==0 && accum_level[0][2]==11) {D1Reco=1850;D2Reco=500;}
        if (D1true>1800) D1true=1799;
        if (D1true<1080) D1true=1081;
        if (topology!=2&&topology!=3&&topology!=11) D1true=1081;
        break;
      case costhA:
        D1Reco=recCosThetaA_Toy[0];D1true=truthCosThetaA;
        if (cutBranch==0 && accum_level[0][2]==11) {D1Reco=1.1;D2Reco=500;}
        if (topology!=2&&topology!=3&&topology!=11) D1true=0.1;
        break;
      case phiA:
        D1Reco=recPhiA_Toy[0];D1true=truthPhiA;
        if (cutBranch==0 && accum_level[0][2]==11) {D1Reco=365;D2Reco=500;}
        if (topology!=2&&topology!=3&&topology!=11) D1true=0.1;
        break;
      default:
        printf("***Warning: Not a valid TKI variable***\n");
        break;
    }

    //if (topology==6||topology==7) weight*=1.5; // +50% pi0
    //if (topology==11) weight*=1.5; // +50% OOPS
    //weight*=GetProtonWeight(trueProtonMom[0],trueProtonStartDir[0][2]);//proton weight
    //weight*=GetFluxWeight(TrueNuEnergy);
    //if (topology==2||topology==3) weight*=GetParticleModelWeight(pMomTrue,pCosThetaTrue,piMomTrue,piCosThetaTrue,muMomTrue,muCosThetaTrue);
    //if ((topology==2||topology==3)&&target==6) weight*=GetParticleModelWeight(pMomTrue,pCosThetaTrue,piMomTrue,piCosThetaTrue,muMomTrue,muCosThetaTrue, D1true);

    if (TKIvar<3) if (topology==2||topology==3) if (target==1) D2true=1.5;//H-flag only for TKI
    if (topology==5) D2true=15;
    if (topology==6) D2true=25;
    if (topology==7) D2true=35;
    if (topology==8) D2true=45;
    if (topology==11) {
      D2true=55;
      if (pMomTrue>1200&&piMomTrue<1200) D2true=65;
      else if (pMomTrue<1200&&piMomTrue>1200) D2true=75;
      else if (pMomTrue>1200&&piMomTrue>1200) D2true=85;
    }

    outtree->Fill();
    if(passCount>1){
      printf("***Warning: More than one cut branch passed***\n");
      for(int j=0;j<10;j++){
        if(branches_passed[j]==1) printf("branch %d passed ...",j);
      }
      printf("\n");
    }
  }
}

void RunTruthTree(TChain* intree_T,TTree *outtree_T,int TKIvar, double sigWeight){
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
    else if (topology_T==11)   topology_T=11; //11 OOPS
    if (topology_T==2||topology_T==3) {
                  if (!(truelepton_dir_T[2]>0.342&&truepi_dir_T[2]>0.342&&TrueProtonStartDir_T[0][2]>0.342&&truelepton_mom_T>250&&truelepton_mom_T<7000&&truepi_mom_T>150&&truepi_mom_T<1200&&TrueProtonMom_T[0]>450&&TrueProtonMom_T[0]<1200)) topology_T=11;//OOPS
    }

    dptt_truth = truthDeltapTT_vtx;
    pN_truth = truthpN;
    daT_truth = truthdeltaAlphaT;
    switch (TKIvar) {
      case dptt:
        D1true_T = truthDeltapTT_vtx;
        if (D1true_T>700) D1true_T=699;
        if (D1true_T<-700) D1true_T=-699;
        if (topology_T!=2&&topology_T!=3&&topology_T!=11) D1true_T=0.1;
        break;
      case pN:
        D1true_T = truthpN;
        if (D1true_T>1500) D1true_T=1499;//overflow
        if (topology_T!=2&&topology_T!=3&&topology_T!=11) D1true_T=0.1; //no truth for bkg
        break;
      case dalphaT:
        D1true_T = truthdeltaAlphaT;
        if (topology_T!=2&&topology_T!=3&&topology_T!=11) D1true_T=0.1;
        break;
      case W:
        D1true_T=truthW;
        if (D1true_T>1800) D1true_T=1799;
        if (D1true_T<1080) D1true_T=1081;
        if (topology_T!=2&&topology_T!=3&&topology_T!=11) D1true_T=1081;
        break;
      case costhA:
        D1true_T=truthCosThetaA;
        if (topology_T!=2&&topology_T!=3&&topology_T!=11) D1true_T=0.1;
        break;
      case phiA:
        D1true_T=truthPhiA;
        if (topology_T!=2&&topology_T!=3&&topology_T!=11) D1true_T=0.1;
        break;
      default:
        printf("***Warning: Not a valid TKI variable***\n");
        break;
    }

    D2true_T=0.1;
    //D2true_T=TrueProtonMom_T[0];
    if (TKIvar<3) if (topology_T==2||topology_T==3) if (target==1) D2true_T=1.5;//H-flag only for TKI
    if (topology_T==5) D2true_T=15;
    if (topology_T==6) D2true_T=25;
    if (topology_T==7) D2true_T=35;
    if (topology_T==8) D2true_T=45;
    if (topology_T==11) {
      D2true_T=55;
      if (TrueProtonMom_T[0]>1200&&truepi_mom_T<1200) D2true_T=65;
      else if (TrueProtonMom_T[0]<1200&&truepi_mom_T>1200) D2true_T=75;
      else if (TrueProtonMom_T[0]>1200&&truepi_mom_T>1200) D2true_T=85;
    }
    //if (topology_T==11) D2true_T=1250;

    muMomTrue_T=truelepton_mom_T;
    muCosThetaTrue_T=truelepton_dir_T[2];
    pMomTrue_T=TrueProtonMom_T[0];
    pCosThetaTrue_T=TrueProtonStartDir_T[0][2];
    piMomTrue_T=truepi_mom_T;
    piCosThetaTrue_T=truepi_dir_T[2];

    //weight_T = weight_T*sigWeight;
    weight_T = weight_T*GetRunWeight(intree_T->GetFile()->GetName());
    //weight_T*=GetProtonWeight(TrueProtonMom_T[0],TrueProtonStartDir_T[0][2]);//proton weight
    //weight_T*=GetFluxWeight(TrueNuEnergy_T);
    //if (topology_T==2||topology_T==3) weight_T*=GetParticleModelWeight(pMomTrue_T,pCosThetaTrue_T,piMomTrue_T,piCosThetaTrue_T,muMomTrue_T,muCosThetaTrue_T);
    //if ((topology_T==2||topology_T==3)&&target==6) weight_T*=GetParticleModelWeight(pMomTrue_T,pCosThetaTrue_T,piMomTrue_T,piCosThetaTrue_T,muMomTrue_T,muCosThetaTrue_T,D1true_T);

    cutBranch_T=-1;
    outtree_T->Fill();
  }

}

int treeConvert_all_variables(TString inFileName="/hepstore/kmtsui/T2K/work/highland2/myNumuCCMultiPiAnalysis/v0r0/output/*.root", TString inTreeName="default",  TString inTreeName_T="truth", TString outFileName="../inputs/neutRun8_new_Bin2_50up_Bin5_50down_Bin6_50up_Bin7_50down_Bin8_50down.root",
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

  double bkgWeight=0;

  //TString inFileName_signal="/scratch/kmtsui/prod6T/single_sys/run*.root";
  //sigWeight=11.53/(131.9);

  //TString inFileName_signal="/scratch/kmtsui/nuwro_prod/nuwro_noFSI.root";
  //sigWeight=11.53/(37.67);

  //inFileName="/scratch/kmtsui/runOutput_run2_4/run*.root";
  //sigWeight=11.53/(62.7+109.7);

  //TString inFileName_signal="/scratch/kmtsui/nuwro_prod/nuwro_sf.root";
  //sigWeight=11.53/(40.46);

  //inFileName="/scratch/kmtsui/prod6D/run*.root";
  //sigWeight=11.53/60.33;

  //TString inFileName_signal="/hepstore/kmtsui/T2K/output/runOutput_run2_4/run*.root";
  //sigWeight=11.53/(62.7+109.7);

  //TString inFileName_bkg="/hepstore/kmtsui/T2K/output/runOutput_run2_4/run*.root";
  //bkgWeight=11.53/(62.7+109.7);

  //TString inFileName_signal="/bundle/data/T2K/users/kmtsui/runOutput_genie/run*.root";
  //sigWeight=11.53/128.21;
  //TString inFileName_bkg="/hepstore/kmtsui/T2K/output/runOutput_genie/run*.root";
  //bkgWeight=11.53/128.21;

  //TString inFileName_bkg="/hepstore/kmtsui/T2K/output/prod6D/run*.root";
  //bkgWeight=11.53/60.33;

  //TString inFileName_bkg="/scratch/kmtsui/prod6T/single_sys/run*.root";
  //bkgWeight=11.53/(131.9);

  //TString inFileName_bkg="/hepstore/kmtsui/T2K/output/nuwro_prod/nuwro_sf.root";
  //bkgWeight=11.53/(40.46);

  //TString inFileName_bkg="/hepstore/kmtsui/T2K/output/nuwro_prod/nuwro_noFSI.root";
  //bkgWeight=11.53/(37.67);

  TString inFileName_signal="/hepstore/kmtsui/T2K/output/prod6T_geombugfix/run*.root";
  //TString inFileName_signal="/scratch/kmtsui/irods/output/run*.root";
  sigWeight=11.60/(193.4);
  //runweight={0.43/12, 0.36/16.8, 1.58/30.8, 1.64/36.1, 1.78/36.1, 1.58/27.2, 4.15/36.1};

  //TString inFileName_signal="/hepstore/kmtsui/T2K/output/runOutput_data_prod6T/run*.root";
  //sigWeight=1.;

  TString inFileName_bkg="/hepstore/kmtsui/T2K/output/prod6T_geombugfix/run*.root";
  bkgWeight=1.;

  //TString inFileName_bkg="/scratch/kmtsui/prod6T/run*.root";
  //bkgWeight=11.53/(195.1);

  outFileName="../inputs/neut_prod6T_dptt_geomfixed.root";
  //outFileName="../inputs/genie_allTKI_geomfixed.root";
  //outFileName="../inputs/data_dptt_prod6T.root";
  //sigWeight*=1.2;
  //int TKIvar = dptt;//0
  //int TKIvar = pN;//1
  //int TKIvar = dalphaT;//2
  int TKIvar = 0;

  TChain* intree;
  intree = new TChain(inTreeName);
  intree->Add(inFileName_signal);
  TChain* intree_T;
  intree_T = new TChain(inTreeName_T);
  intree_T->Add(inFileName_signal);

  TChain* intree_bkg;
  intree_bkg = new TChain(inTreeName);
  intree_bkg->Add(inFileName_bkg);

  TFile *outfile = new TFile(outFileName,"recreate");
  TTree *outtree = new TTree("selectedEvents", "selectedEvents");
  TTree *outtree_T = new TTree("trueEvents", "trueEvents");

  //SetProtonWeight();
  //SetFluxWeight();
  //SetParticleModelWeight("/hepstore/kmtsui/T2K/output/proton_model/TKI_nuisevt_lfgvalenciarpa.root");
  //SetParticleModelWeight("/hepstore/kmtsui/T2K/NuWro/nuwro/output/TKI_out_SF.root", TKIvar);

  SetOutTreeBranches(outtree);

  SetInTreeBranches(intree);
  RunInTree(intree,    outtree,TKIvar,sigWeight,false,false,false,false);

  //SetInTreeBranches(intree_bkg);
  //RunInTree(intree_bkg,outtree,TKIvar,bkgWeight,false,false,false,true);

  SetInTruthTreeBranches(intree_T);
  SetOutTruthTreeBranches(outtree_T);
  RunTruthTree(intree_T,outtree_T,TKIvar,sigWeight);

  outfile->cd();
  printf("***Output Rec Tree: ***\n");
  outtree->Print();
  printf("***Output True Tree: ***\n");
  outtree_T->Print();
  outfile->Write();

  //h_proton_weight->Write();
  //h_pion_weight->Write();
  //h_muon_weight->Write();
  //h_tki_weight->Write();

  //delete infile;
  delete outfile;
  return 0;
}
