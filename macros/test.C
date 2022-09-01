void test(){
  double VarBins[100];
  const int nDpttbins = 6;
  double Dpttbins[nDpttbins+1] = {-700,-300,-100,100,300,700,800};
  for (int j=0;j<nDpttbins+1;j++) VarBins[j] = Dpttbins[j];
  TH1D* hTest = new TH1D("","",nDpttbins,VarBins);
  hTest->Draw();
}
