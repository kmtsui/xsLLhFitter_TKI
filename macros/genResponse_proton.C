//#include "genResponse_test.h"

void genResponse(string microtreeIn="/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root", string rwpath="", string dialName="CA5", double nom=1., double err=0.12/1.01, int DISMultiPi = -1) 
{ 

  // declare variables for fit parameters, names, and values
  vector<string> param_names;
  vector<double> param_values;
  double bestfit;
  double error;
  //double sigma;


  // define inputs
  string microtree = microtreeIn;
  string rwfile_default = rwpath+"rwfile_default.root";
  string rwfile_truth = rwpath+"rwfile_truth.root";
  string flux = "/data/t2k/dolan/xsToolBasedMEC/CC0PiAnl/fluxFiles/nd5_tuned11bv3.2_11anom_run2_fine.root";
  string fluxname = "enu_nd5_tuned11b_numu";

  //string outname = Form("_allsample.root",sam,reac);
  string outname = "_allsample.root";
  if (DISMultiPi>-1) outname = Form("_%d_allsample.root",DISMultiPi);
  //string outputFileName = "./output/resFunc"+dialName+outname;
  string outputFileName = dialName+outname;
  cout << "Output file name is: " << outputFileName << endl; 
  
  string dialGraphName = dialName+"Graph";
  
  // set up unfolding variables, signal, and cut
  string mc_variable = "trueDpT";
  string data_variable = "Sum$(recDpT[xstool_throw])";
  string signal, cut, bcut; // these will be defined in the main loop
  //string psCut = "&& (Sum$(selp_mom[xstool_throw])>450) && (Sum$(selmu_mom[xstool_throw])>150) && (Sum$(selp_theta[xstool_throw])>-3.15) && (Sum$(selmu_theta[xstool_throw])>1.57)";
  
  // set up recon binning
  const int nRbins = 10;
  //const double Rbins[nRbins+1] = { 0.0, 0.125, 0.2, 0.28, 0.37, 0.47, 0.6, 0.83, 2.5 };
  //const double Rbins[nRbins+1] = { -300, -279, -259, -238, -217, -197, -176, -155, -134, -114, -93, -72, -52, -31, -10, 10, 31, 52, 72, 93 ,114, 134, 155, 176, 197, 217, 238, 259, 279, 300};
  const double Rbins[nRbins+1] = {0,200,400,600,800,1000,1200,1400,1600,1800,2000};

  // set up true binning
  const int nTbins = 10;
  //const double Tbins[nTbins+1] = { 0.0, 0.125, 0.2, 0.28, 0.37, 0.47, 0.6, 0.83, 2.5 };
 // const double Tbins[nTbins+1] = { -300, -279, -259, -238, -217, -197, -176, -155, -134, -114, -93, -72, -52, -31, -10, 10, 31, 52, 72, 93 ,114, 134, 155, 176, 197, 217, 238, 259, 279, 300};
  const double Tbins[nTbins+1] = {0,200,400,600,800,1000,1200,1400,1600,1800,2000};

  // You need to provide the number of branches in your HL2 tree
  // And the accum_level you want to cut each one at to get your selected events
  // i.e choosing n in accum_level[0][branch_i]>n
  const int nbranches = 10;
  const int accumToCut[nbranches] =   {7,8,9,8,7,5,4,7,8,7};

  
  // delcare xsBinning, they will be used always the same in the main loop
  //xsBinning *binningTru = new xsBinning();
  //xsBinning *binningRec = new xsBinning();
  //binningTru->addDimension(mc_variable,nTbins,Tbins);
  //binningRec->addDimension(data_variable,nRbins,Rbins);
  
  // setup on-the-fly reweighting parameter
  param_names.push_back(dialName);
  bestfit = nom;
  error = err;

  //Write TGraph in the output file
  TFile *output = new TFile(outputFileName.c_str(),"RECREATE");   

  // Total number of topologies
  const int nTopologies = 90;// nbranches;
  const int nReactions = 12;
  // Total number of weights
  const int nWeights = 7;
  // reweighted reconstrcted distribution for each weight, true bin, and topology

 for (int sample =0;sample<5;sample++) {

  TH1 *hreco[nWeights][nTbins][nReactions];
  TH1 *htrue[nWeights];
  for (int i=0;i<nWeights;i++) for (int j=0;j<nTbins;j++) for (int k=0;k<nReactions;k++) hreco[i][j][k]=new TH1D("","",nRbins,Rbins);
  for (int i=0;i<nWeights;i++) htrue[i]=new TH1D("","",nTbins,Tbins);

  Double_t dptt[nTopologies], dpttTrue[nTopologies], FluxWeight[nTopologies], pmomTrue[nTopologies], pmomRec[nTopologies];
  Int_t topo_np[nTopologies], topo_npi[nTopologies], topo_target[nTopologies], SelectedSample;

  TChain* sumT;
  sumT = new TChain("sample_sum");
  sumT->Add(microtreeIn.data());

  //TFile *fIn = new TFile(microtreeIn.data(),"READ");
  //TTree *sumT = (TTree*)fIn->Get("sample_sum");
  sumT->SetBranchAddress("SelectedSample", &SelectedSample);
  sumT->SetBranchAddress("dptt", &dptt);
  sumT->SetBranchAddress("dpttTrue", &dpttTrue);
  sumT->SetBranchAddress("topo_np", &topo_np);
  sumT->SetBranchAddress("topo_npi", &topo_npi);
  sumT->SetBranchAddress("topo_target", &topo_target);
  sumT->SetBranchAddress("FluxWeight", &FluxWeight);
  sumT->SetBranchAddress("pmomTrue", &pmomTrue);
  sumT->SetBranchAddress("pmomRec", &pmomRec);
  TClonesArray *arr = new TClonesArray("TGraph");
  sumT->GetBranch(dialGraphName.data())->SetAutoDelete(kFALSE);
  sumT->SetBranchAddress(dialGraphName.data(),&arr);
  //TGraph dialGraph;
  //sumT->SetBranchAddress(dialGraphName.data(),&dialGraph);
 
  Long64_t nentries = sumT->GetEntries();
  std::cout<<"nentries = "<<nentries<<std::endl;
  for (Long64_t ev=0;ev<nentries;ev++) {
    arr->Clear();
    sumT->GetEntry(ev);
    if (SelectedSample!=sample) continue;
    //if (SelectedSample<0||SelectedSample>4) continue;
    //std::cout<<"at entry "<<ev<<std::endl;    
    //std::cout<<"nGraph = "<<arr->GetEntriesFast()<<std::endl;
    int rBin = -1;
    int tBin = -1;
    for (int i=0;i<nRbins;i++) {
      if (pmomRec[SelectedSample]>Rbins[i] && pmomRec[SelectedSample]<Rbins[i+1]) {rBin=i;break;}
    }
    for (int i=0;i<nTbins;i++) {
      if (pmomTrue[SelectedSample]>Tbins[i] && pmomTrue[SelectedSample]<Tbins[i+1]) {tBin=i;break;}
    }
    int evtTopo=topo_npi[SelectedSample];
    if (topo_npi[SelectedSample]<7) evtTopo=topo_npi[SelectedSample]; //0 CC0pi, 1 CC1pi0p, 2 CC1pi1p, 3 CC1piNp, 4 CC2pi+, 5 CC1pi+1pi-, 6 CC1pi+Npi0
    else if (topo_npi[SelectedSample]==8) evtTopo=7; //7 CC-other-Npi0 
    else if (topo_npi[SelectedSample]==10) evtTopo=8; //8 CC-other-0pi0 
    else if (topo_npi[SelectedSample]==999) evtTopo=9; //9 backg(NC+antinu+nue)
    else if (topo_npi[SelectedSample]==7)   evtTopo=10; //10 OOFV
    else if (topo_npi[SelectedSample]==88)   evtTopo=11; //11 OOPS
    if(evtTopo<0||evtTopo>11) {cout << "*** Warning: evtTopology < 0 or > 10, evt Topology is " << evtTopo << endl;continue;}
    if (evtTopo!=2&&evtTopo!=3&&evtTopo!=11) tBin=2;
    Double_t x[10], y[10];
    TGraph *dialGraph=(TGraph*)arr->At(SelectedSample);
    //std::cout<<"dialGraph->GetN()="<<dialGraph->GetN()<<std::endl;
    for(int w = 0; w < dialGraph->GetN(); w++){
      dialGraph->GetPoint(w,x[w],y[w]);  
      //if (dialName.find("FSI") != string::npos) if (evtTopo==2||evtTopo==3||evtTopo==11) y[w]=1.;//For FSI parameters: do not apply weight for signal events
      if (DISMultiPi>-1) { //separte DIS dials for bkg topo
        if (DISMultiPi==0) {if (evtTopo==5||evtTopo==6||evtTopo==7||evtTopo==8) y[w]=1.;} //this DIS dial does not contain the target bkg topo
        else if (DISMultiPi!=evtTopo) y[w]=1.;
        else if (DISMultiPi==evtTopo) y[w]=x[w]+1;
      }
      evtTopo=0;
      hreco[w][0][0]->Fill(pmomRec[SelectedSample],y[w]*FluxWeight[SelectedSample]);
      htrue[w]->Fill(pmomTrue[SelectedSample],y[w]*FluxWeight[SelectedSample]);
    }
    /*if (dialGraph->GetN()>1) {
      //const static int nP = (int)dialGraph->GetN();
      Double_t x[10], y[10];
      for (Int_t i =0;i<dialGraph->GetN();i++) {
         dialGraph->GetPoint(i,x[i],y[i]);
         std::cout<<"Point "<<i<<" x = "<<x[i]<<" y ="<<y[i]<<std::endl; 
      }
    dialGraph->SetMarkerSize(1);
    dialGraph->SetMarkerStyle(8);
    TCanvas *c1 = new TCanvas();
    dialGraph->Draw();
    break;
    }*/
    //dialGraph->Draw();
  } 
  /*TCanvas *c1 = new TCanvas();
  hreco[0][2][2]->Draw();

  TCanvas *c1 = new TCanvas();
  hreco[1][2][2]->Draw();

  TCanvas *c2 = new TCanvas();
  hreco[2][2][2]->Draw();

  TCanvas *c2 = new TCanvas();
  hreco[3][2][2]->Draw();*/
  sumT->Delete();
  //fIn->Close();

  /*for(const int w = 0; w < nWeights; w++){
    for(const int bt = 0; bt < nTbins; bt++){
      int t = topo;
      int r = reac;
      //for(const int t = 0; t < nTopologies; t++){
      
        if(t==0 || t==4 || r==6) continue; // Ignore branches with no proton
        // define xsInput
        xsInputdataHighland2 *inpRwFly = new xsInputdataHighland2();
        inpRwFly->SetupFromHighland(microtree,microtree);
        inpRwFly->SetFlux(flux, fluxname);
        inpRwFly->SetRooTrackerVtxFiles(rwfile_default,rwfile_truth);
        inpRwFly->SetTrueBinning(binningTru);
        inpRwFly->SetReconstructedBinning(binningRec);

        // set cut definition, for each topology
        cut = Form("((Sum$(accum_level[xstool_throw][%d]) > %d) && (mectopology==%d))", t, accumToCut[t], r);

        // set signal definition, for each true bin
        signal = Form("(trueDpT) > %f && (trueDpT) < %f",Tbins[bt],Tbins[bt+1]);
        bcut = Form("((Sum$(trueDpT[xstool_throw]) > %f) && (Sum$(trueDpT[xstool_throw])) < %f)",Tbins[bt],Tbins[bt+1]);

        inpRwFly->SetSignal(signal.c_str());
        inpRwFly->SetCut((cut+" && "+bcut).c_str());

        // define xsReweightParameter
        xsReweightParameter *rwPar = new xsReweightParameter(param_names);
        inpRwFly->AddReweighting(rwPar);


        // set parameter value
        param_values.clear();    
        //param_values.push_back((-3+w)*error/bestfit);
        param_values.push_badialGraphck((bestfit-1)+((-3+w)*error));
        string name = rwPar->SetParameterValue(param_values);

        // reweighting with the given value
        inpRwFly->GenerateThrow(name,0,true);

        // get reconstructed distribution
        hreco[w][bt][t] = inpRwFly->GetReconstructed();

        // delete xsInput
        delete inpRwFly;
	
      //}
    }
  }
  
    
  //  }*/

 
  //Write TGraph in the output file
  //output->mkdir(Form("sample_%d",sample));

  char dir[200];

  TGraph *ReWeight[nRbins][nTbins][nTopologies];
  
  double MA[7];
  double MA_PDD[7]={-1+1,-0.5+1,-0.05+1,0+1,0.05+1,0.5+1,1+1};
  for(int w = 0; w < 7; w++){
      MA[w]=bestfit-error*(3-w);
      if (dialName.compare("PDD_C")==0) MA[w]=MA_PDD[w];
  }  
  //if (dialName.compare("PDD_C")==0) MA={-1,-0.5,-0.05,0,0.05,0.5,1};

  for (int w=0;w<nWeights;w++){
    char nameHisto[256];
    sprintf(nameHisto,"Rec_Weight_%d_sample_%d",w,sample);
    hreco[w][0][0]->SetName(nameHisto);
    hreco[w][0][0]->SetTitle(nameHisto);
    hreco[w][0][0]->SetMarkerStyle(20);
    hreco[w][0][0]->SetMarkerColor(2);
    if (w!=3) hreco[w][0][0]->Divide(hreco[3][0][0]);
    hreco[w][0][0]->Write();
    sprintf(nameHisto,"True_Weight_%d_sample_%d",w,sample);
    htrue[w]->SetName(nameHisto);
    htrue[w]->SetTitle(nameHisto);
    htrue[w]->SetMarkerStyle(20);
    htrue[w]->SetMarkerColor(2);
    if (w!=3) htrue[w]->Divide(htrue[3]);
    htrue[w]->Write();
  }

  /*for(int br = 0; br<nRbins; br++){//reco kinematics bin
    for(int bt = 0; bt < nTbins; bt++){//true kinematics bin
      //if(fabs(br-bt)>20) continue;  //save memory if reco bin and true bin very far away
      //for(int t = 0; t < nTopologies; t++){//topology
      for (int r=0;r<nReactions;r++) {
        int t=sample;
      //int r=reac;

	       //sprintf(dir,"sample_%d",t);
	       //output->cd(dir);
  	     //for(int r=0;r<5;r++){ //reaction
	       char nameHisto[256];
	       sprintf(nameHisto,"RecBin_%d_trueBin_%d_sample_%d_reac_%d",br,bt,t,r);
  	     ReWeight[br][bt][t] = new TGraph(7);
  	     ReWeight[br][bt][t]->SetName(nameHisto);
  	     ReWeight[br][bt][t]->SetTitle(nameHisto);
  	     ReWeight[br][bt][t]->SetMarkerStyle(20);
  	     ReWeight[br][bt][t]->SetMarkerColor(2);
	       for(int w=0;w<nWeights;w++){
	         if(hreco[3][bt][r]->GetBinContent(br+1) !=0 ){
	           ReWeight[br][bt][t]->SetPoint(w,MA[w],hreco[w][bt][r]->GetBinContent(br+1)/hreco[3][bt][r]->GetBinContent(br+1));
                   //if (abs(hreco[w][bt][r]->GetBinContent(br+1)/hreco[3][bt][r]->GetBinContent(br+1)-1)>0.01) std::cout<<"br "<<br<<" bt "<<bt<<" r "<<r<<std::endl;
	          }
	          else{
	            ReWeight[br][bt][t]->SetPoint(w,MA[w],1);
	          }
	          ReWeight[br][bt][t]->GetYaxis()->SetTitle("weight");	   
	       }
               //output->cd();
	       ReWeight[br][bt][t]->Write();
	     //}
       }
     }
   }*/
 } 
   output->Close();
  
}


void genResponse_proton(){
  /*genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","MAQE",1.21,0.4/1.2*1.21);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","PF_C",475./2.,66.7/(6.));
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","CA5",1.01,0.12/1.01);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","MARES",0.95,0.15/0.95);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","BgSclRes",1.3,0.2/1.3);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","CCNuE",1,0.02/1.);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","DISMPiShp",1,0.4/1.);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","NCCoh",1,0.3/1.);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","NCOth",1,0.3/1.);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","FSI_PI_ABS",1.1,0.267*1.1);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","FSI_CEX_LO",1,0.267);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","FSI_INEL_LO",1,0.267);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","FSI_PI_PROD",1.,0.267*1.);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","FSI_CEX_HI",1.8,0.267*1.8);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","FSI_INEL_HI",1.8,0.267*1.8);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","MEC_C",1,1./3.);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","PDD_C",1,1.);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","CCCoh",1,0.3);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","DISMPiShp",1,0.4/1.,0);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","DISMPiShp",1,0.4/1.,5);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","DISMPiShp",1,0.4/1.,6);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","DISMPiShp",1,0.4/1.,7);
  genResponse("/hepstore/kmtsui/T2K/work/GlobalAnalysisTools/T2KReWeight/app/run*.root","","DISMPiShp",1,0.4/1.,8);*/
  genResponse("/scratch/kmtsui/prod6D/reweight/run*.root","","Nucleon_MFP",1.0,0.2);
}

