double myfunc(double* x, double* par)
{
  double xx = x[0];

  double x0 = par[4];

  double A = par[0];
  double C = par[1];
  double D = par[2];
  double E = par[3];
  double B = D*E*exp(E*x0) - 2.*C*x0;
  double F = A + B*x0 + C*x0*x0 - D*exp(E*x0);

  if(xx < x0) return par[5]*(A + B*xx + C*xx*xx);
  else        return par[5]*(D*exp(E*xx)+F);
}



void drawDCRVsVov()
{
  TFile* inFile_DCRParams = TFile::Open("../data/DCRParams_new.root","READ");
  TGraph* g_DCR_par0 = (TGraph*)( inFile_DCRParams->Get("g_par0"));
  TGraph* g_DCR_par1 = (TGraph*)( inFile_DCRParams->Get("g_par1"));
  TGraph* g_DCR_par2 = (TGraph*)( inFile_DCRParams->Get("g_par2"));
  TGraph* g_DCR_par3 = (TGraph*)( inFile_DCRParams->Get("g_par3"));
  TGraph* g_DCR_par4 = (TGraph*)( inFile_DCRParams->Get("g_par4"));

  float fluence = 1E13;
  float DCRRef_par0 = g_DCR_par0->Eval(fluence);
  float DCRRef_par1 = g_DCR_par1->Eval(fluence);
  float DCRRef_par2 = g_DCR_par2->Eval(fluence);
  float DCRRef_par3 = g_DCR_par3->Eval(fluence);
  float DCRRef_par4 = g_DCR_par4->Eval(fluence);
  TF1* f_DCRRef_vs_Vov = new TF1("f_DCRRef_vs_Vov",myfunc,0.,5.,6);
  f_DCRRef_vs_Vov -> SetParameters(DCRRef_par0,DCRRef_par1,DCRRef_par2,DCRRef_par3,DCRRef_par4,1.);
  //f_DCRRef_vs_Vov -> SetParameter(5,15000./f_DCRRef_vs_Vov->Eval(4));
  
  f_DCRRef_vs_Vov -> Draw();


  fluence = 2E14;
  DCRRef_par0 = g_DCR_par0->Eval(fluence);
  DCRRef_par1 = g_DCR_par1->Eval(fluence);
  DCRRef_par2 = g_DCR_par2->Eval(fluence);
  DCRRef_par3 = g_DCR_par3->Eval(fluence);
  DCRRef_par4 = g_DCR_par4->Eval(fluence);
  TF1* f_DCRRef_vs_Vov_2 = new TF1("f_DCRRef_vs_Vov_2",myfunc,0.,5.,6);
  f_DCRRef_vs_Vov_2 -> SetParameters(DCRRef_par0,DCRRef_par1,DCRRef_par2,DCRRef_par3,DCRRef_par4,1.);
  f_DCRRef_vs_Vov_2 -> SetParameter(5,f_DCRRef_vs_Vov->Eval(1)/f_DCRRef_vs_Vov_2->Eval(1));
  f_DCRRef_vs_Vov_2 -> SetLineColor(kBlue);
  
  f_DCRRef_vs_Vov_2 -> Draw("same");
}

