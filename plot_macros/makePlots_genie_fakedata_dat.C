#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>

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

//#############################################################################

void makePlots(){
  //**** Set Style for Plots ****
  TStyle *t2kstyle = new TStyle("T2K","T2K approved plots style");
  SetStyleVariables(t2kstyle);
  gROOT->SetStyle("T2K");

  TFile* fIn = new TFile("xsec_neutMC_genieData_dat.root");
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
  //hist_par_xsec_error_final->Divide(hist_par_xsec_error_prior);
  //hist_par_xsec_error_final->SetTitle("Final error/Prior error");

  TH1D* hist_prefit_par_all = (TH1D*)fIn -> Get("hist_prefit_par_all");
  TMatrixDSym* prefit_cov = (TMatrixDSym*)fIn -> Get("prefit_cov");

  static const int Npar = postfit_param->GetNoElements();
  TBox *box_err[Npar];

  const int nBins = 3;
  const double Bins[nBins+1] = {0,60,120,180};

  TH1D* chi2_tot_periter = (TH1D*)fIn -> Get("chi2_tot_periter");
  double chi2_prefit = chi2_tot_periter->GetBinContent(1);
  double chi2_postfit = chi2_tot_periter->GetBinContent(chi2_tot_periter->GetNbinsX()-1);

  TCanvas *c1 = new TCanvas();
  TH1D* hist_par_all = new TH1D("","",postfit_param->GetNoElements(),0,postfit_param->GetNoElements());
  hist_par_all->Draw();
  TLine *line = new TLine(0,1,Npar,1);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  //line->Draw("same");
  for (int i=0;i<postfit_param->GetNoElements();i++) {
    hist_par_all->SetBinContent(i+1,(*postfit_param)[i]/hist_prefit_par_all->GetBinContent(i+1));
    hist_par_all->SetBinError(i+1,sqrt((*postfit_cov)[i][i])/hist_prefit_par_all->GetBinContent(i+1));
    if (i<nBins+1) continue;
    box_err[i] = new TBox(hist_par_all->GetBinLowEdge(i+1),1-sqrt((*prefit_cov)[i][i])/hist_prefit_par_all->GetBinContent(i+1),
                          hist_par_all->GetBinLowEdge(i+2),1+sqrt((*prefit_cov)[i][i])/hist_prefit_par_all->GetBinContent(i+1));
    box_err[i]->SetFillColor(856);
    box_err[i]->Draw("same");
  }
  gStyle->SetErrorX(0);
  hist_par_all->SetMarkerStyle(kFullCircle);
  hist_par_all->Draw("E1 X0 same");
  line->Draw("same");

  TCanvas *c2 = new TCanvas();
  postfit_cor->Draw("colz");
  TH2D *h_postfit_cor = (TH2D*) c2->GetPrimitive("TMatrixDBase");
  h_postfit_cor->GetZaxis()->SetRangeUser(-1,1);

  TCanvas *c3 = new TCanvas();
  TH1D* evhist_prefit_all = (TH1D*)fIn -> Get("evhist_prefit_all");
  TH1D* evhist_postfit_all = (TH1D*)fIn -> Get("evhist_postfit_all");
  TH1D* evhist_data_all = (TH1D*)fIn -> Get("evhist_data_all");
  evhist_prefit_all->SetLineColor(kRed);
  evhist_postfit_all->SetLineColor(kBlue);
  evhist_data_all->SetLineColor(kBlack);
  evhist_data_all->SetMarkerStyle(kFullCircle);
  evhist_postfit_all->Draw();
  evhist_prefit_all->Draw("same");
  evhist_data_all->Draw("E1 X0 same");
  auto legend = new TLegend(0.15,0.72,0.75,0.87);
  legend->AddEntry(evhist_data_all,"Data","lep");
  legend->AddEntry(evhist_prefit_all,Form("Pre-fit, #chi^{2}=%4.1f",chi2_prefit),"l");
  legend->AddEntry(evhist_postfit_all,Form("Post-fit, #chi^{2}=%4.1f",chi2_postfit),"l");
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->Draw();

  TCanvas *c4 = new TCanvas();
  TH1D* sel_best_fit_results = (TH1D*)fIn -> Get("sel_best_fit");
  TFile* fIn_neut = new TFile("xsec_neut_dat.root");
  TH1D* sel_best_fit_neut = (TH1D*)fIn_neut -> Get("sel_best_fit");
  TFile* fIn_genie = new TFile("xsec_genie_dat.root");
  TH1D* sel_best_fit_genie = (TH1D*)fIn_genie -> Get("sel_best_fit");
  TH1D* h_unfolded_xsec = new TH1D("","",nBins,Bins);
  TH1D* h_xsec_neut = new TH1D("","",nBins,Bins);
  TH1D* h_xsec_genie = new TH1D("","",nBins,Bins);
  for (int i=0;i<nBins;i++){
    h_unfolded_xsec->SetBinContent(i+1,sel_best_fit_results->GetBinContent(i+1));
    h_unfolded_xsec->SetBinError(i+1,sel_best_fit_results->GetBinError(i+1));
    h_xsec_neut->SetBinContent(i+1,sel_best_fit_neut->GetBinContent(i+1));
    h_xsec_genie->SetBinContent(i+1,sel_best_fit_genie->GetBinContent(i+1));
  }
  h_xsec_neut->SetLineColor(kBlue);
  h_xsec_genie->SetLineColor(kGreen);
  h_unfolded_xsec->SetLineColor(kRed);
  h_unfolded_xsec->SetMarkerStyle(kFullCircle);
  h_unfolded_xsec->GetXaxis()->SetTitle("#delta#alpha_{T} (degree)");
  h_unfolded_xsec->GetYaxis()->SetTitle("#frac{d#sigma}{d#delta#alpha_{T}}(Nucleon^{-1}cm^{2})");
  h_unfolded_xsec->GetYaxis()->SetRangeUser(0,90e-42);
  h_unfolded_xsec->Draw("E0 X0)");
  h_xsec_neut->Draw("same");
  h_xsec_genie->Draw("same");
  auto legend2 = new TLegend(0.15,0.78,0.72,0.89);
  legend2->AddEntry(h_xsec_neut,"Input Monte Carlo (NEUT)","l");
  legend2->AddEntry(h_xsec_genie,"Fake Data (GENIE)","l");
  legend2->AddEntry(h_unfolded_xsec,"Result","lep");
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->Draw();
}

void makePlots_genie_fakedata_dat(){
  makePlots();
}
