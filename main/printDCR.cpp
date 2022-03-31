#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/Functions.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"


int main(int argc, char** argv)
{
  float fluence = atof(argv[1]);
  float Vov = atof(argv[2]);
  float DCRScale = atof(argv[3]);
  float VovScale = atof(argv[4]);
  
  TFile* inFile_DCRParams = TFile::Open("data/DCRParams_new_Dic2021.root","READ");
  int nParDCR = 8;
  std::map<int,TGraph*> g_DCR_pars;
  for(int iPar = 0; iPar < nParDCR; ++iPar)
  {
    g_DCR_pars[iPar] = (TGraph*)( inFile_DCRParams->Get(Form("g_par%d",iPar)) );
  }
  
  
  TF1* f_DCRRef_vs_Vov = new TF1("f_DCRRef_vs_Vov",myfunc_DCR,0.,7.,nParDCR);
  float* DCRRef_pars = new float[nParDCR];
  for(int iPar = 0; iPar < nParDCR; ++iPar)
  {    
    DCRRef_pars[iPar] = g_DCR_pars[iPar]->Eval(fluence);
    f_DCRRef_vs_Vov -> SetParameter(iPar,DCRRef_pars[iPar]);
  }

  std::cout << "DCR[" << Vov << "] = " << DCRScale * f_DCRRef_vs_Vov->Eval(Vov)/f_DCRRef_vs_Vov->Eval(VovScale) << std::endl;
  
  
  return 0;
}
