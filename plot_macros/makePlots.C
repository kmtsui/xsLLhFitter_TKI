#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>

void makePlots(){
  //**** Set Style for Plots ****
  TStyle *t2kstyle = new TStyle("T2K","T2K approved plots style");
  SetStyleVariables(t2kstyle);
  gROOT->SetStyle("T2K");

  TFile* fIn = new TFile("xsec_pN_neut.root");
  TMatrixDSym* xsec_cor = (TMatrixDSym*)fIn -> Get("xsec_cor");
  TMatrixDSym* postfit_cor = (TMatrixDSym*)fIn -> Get("postfit_cor");
  TMatrixDSym* postfit_cov = (TMatrixDSym*)fIn -> Get("postfit_cov");
  TVectorD* postfit_param = (TVectorD*)fIn -> Get("postfit_param");
  TH1D* sel_best_fit = (TH1D*)fIn -> Get("sel_best_fit");
  TH1D* h_Diag_Err = new TH1D("diagonal_Error","Fractional error of diagonal elements",postfit_param->GetNoElements(),0,postfit_param->GetNoElements());
  for (int i=0;i<postfit_param->GetNoElements();i++) {
    h_Diag_Err->SetBinContent(i+1,sqrt((*postfit_cov)[i][i])/(*postfit_param)[i]);
  }
  TH1D* hist_par_xsec_error_final = (TH1D*)fIn -> Get("hist_par_xsec_error_final");
  TH1D* hist_par_xsec_error_prior = (TH1D*)fIn -> Get("hist_par_xsec_error_prior");
  hist_par_xsec_error_final->Divide(hist_par_xsec_error_prior);
  hist_par_xsec_error_final->SetTitle("Final error/Prior error");

  //const int nDpttbins = 5;
  //const double Dpttbins[nDpttbins+1] = {-700,-300,-100,100,300,700};
  const int nBins = 4;
  const double Bins[nBins+1] = {0,150,300,600,1500};

  TH1D* h_unfolded_xsec = new TH1D("","",nBins,Bins);
  for (int i=0;i<nBins;i++){
    h_unfolded_xsec->SetBinContent(i+1,sel_best_fit->GetBinContent(i+1));
    h_unfolded_xsec->SetBinError(i+1,sel_best_fit->GetBinError(i+1));
  }
  

  TCanvas *c1 = new TCanvas();
  xsec_cor->Draw("colz");
  TH2D *h_xsec_cor = (TH2D*) c1->GetPrimitive("TMatrixDBase");
  h_xsec_cor->GetZaxis()->SetRangeUser(-1,1);

  TCanvas *c2 = new TCanvas();
  postfit_cor->Draw("colz");
  TH2D *h_postfit_cor = (TH2D*) c2->GetPrimitive("TMatrixDBase");
  h_postfit_cor->GetZaxis()->SetRangeUser(-1,1);

  TCanvas *c3 = new TCanvas();
  h_Diag_Err->Draw();

  TCanvas *c4 = new TCanvas();
  //h_unfolded_xsec->GetXaxis()->SetTitle("#deltap_{TT} (MeV)");
  //h_unfolded_xsec->GetYaxis()->SetTitle("#frac{d#sigma}{d#deltap_{TT}}(Nucleon^{-1}cm^{2}MeV^{-1})");
  h_unfolded_xsec->GetXaxis()->SetTitle("p_{N} (MeV)");
  h_unfolded_xsec->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{N}}(Nucleon^{-1}cm^{2}MeV^{-1})");
  h_unfolded_xsec->Draw("E1");

  TCanvas *c5 = new TCanvas();
  hist_par_xsec_error_final->Draw();
  TPaveText *t = new TPaveText(0,0.9,0.3,1.0);
  t->AddText("Final error/prior error");
  t->Draw();
}

//#############################################################################
void SetStyleVariables(TStyle *t2kStyle){

  // use plain black on white colors
  t2kStyle->SetFrameBorderMode(0);
  t2kStyle->SetCanvasBorderMode(0);
  t2kStyle->SetPadBorderMode(0);
  t2kStyle->SetPadColor(0);
  t2kStyle->SetCanvasColor(0);
  t2kStyle->SetStatColor(0);
  //t2kStyle->SetFillColor(0);
  t2kStyle->SetLegendBorderSize(1);

  // set the paper & margin sizes
  t2kStyle->SetPaperSize(20,26);
  t2kStyle->SetPadTopMargin(0.1);
  t2kStyle->SetPadRightMargin(0.15); //0.05 
  t2kStyle->SetPadBottomMargin(0.15);
  t2kStyle->SetPadLeftMargin(0.15);

  // use large Times-Roman fonts
  t2kStyle->SetTextFont(132);
  t2kStyle->SetTextSize(0.08);
  t2kStyle->SetLabelFont(132,"x");
  t2kStyle->SetLabelFont(132,"y");
  t2kStyle->SetLabelFont(132,"z");
  t2kStyle->SetLabelSize(0.05,"x");
  t2kStyle->SetTitleSize(0.06,"x");
  t2kStyle->SetLabelSize(0.05,"y");
  t2kStyle->SetTitleSize(0.06,"y");
  t2kStyle->SetLabelSize(0.05,"z");
  t2kStyle->SetTitleSize(0.06,"z");
  t2kStyle->SetLabelFont(132,"t");
  t2kStyle->SetTitleFont(132,"x");
  t2kStyle->SetTitleFont(132,"y");
  t2kStyle->SetTitleFont(132,"z");
  t2kStyle->SetTitleFont(132,"t");
  t2kStyle->SetTitleFillColor(0);
  t2kStyle->SetTitleX(0.25);
  t2kStyle->SetTitleFontSize(0.08);
  t2kStyle->SetTitleFont(132,"pad");

  //t2kStyle->SetPadGridX(true);
  //t2kStyle->SetPadGridY(true);

  // use bold lines and markers
  //  t2kStyle->SetMarkerStyle(20);
  t2kStyle->SetHistLineWidth(1.85);
  t2kStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  //  t2kStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  t2kStyle->SetOptTitle(0);
  t2kStyle->SetOptStat(0);
  t2kStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  //t2kStyle->SetPadTickX(1);
  t2kStyle->SetPadTickY(1);

  t2kStyle->SetPalette(1,0);  // use the nice red->blue palette
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                   NCont);
  t2kStyle->SetNumberContours(NCont);

  // End of definition of t2kStyle
}
