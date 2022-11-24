#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/Functions.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TBox.h"
#include "TLatex.h"

#define timeStep 60. // in minutes



int main(int argc, char** argv)
{
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);

  double T_ref = 21.;
  double T_ann = opts.GetOpt<double>("Input.Tann");
  
  double totFluence = opts.GetOpt<double>("Input.totFluence");
  
  
  TGraph g_instLumi_vs_time, g_alpha_vs_time, g_alphaNorm_vs_time, g_dcr_vs_time;
  TGraph g_fluence_vs_time;
  TGraph g_intLumi_vs_time;
  TGraph g_temp_vs_time;
  
  int year, month, irr;
  long double lumi;
  long double alpha = 0., norm = 0.;
  long double alpha0[20], tau[20];
  for(int i = 0; i < 20; ++i)
  {
    tau[i] = pow(10,i-5);
    alpha0[i] = (alpha_p0+alpha_p1*(i-5)+alpha_p2*pow(i-5,2)+alpha_p3*pow(i-5,3)) / 20.;
  }
  
  long double vett[20];
  for(int j = 0; j < 20; ++j)
    vett[j] = 0;
  
  
  
  //-------------------------------
  // define irradiation temperature
  double T_irr = 20.;
  
  
  
  //--------------------------------------------------------------------------------------
  // further update luminosity and temperature  within a month (with timestep granularity)
  for(int jmin = timeStep; jmin < minInMonth*12; jmin+=timeStep) // loop over 60*24*30*12=518400 minutes in a year
  {
    float effLumi = 0.;
    float effTemp = 21.;
    float theta = 1.;
    
    //--- irradiation (20 hours at 20° C)
    if( jmin < 1200 )
    {
      effLumi = 1.;
      effTemp = 20.;
      theta = timeScale(effTemp, T_ref);
    }

    //--- cool down time (2 weeks at -30° C)
    else if( jmin < 21360 )
    {
      effLumi = 0.;
      effTemp = -30.;
      theta = timeScale(effTemp, T_ref);
    }
    
    //--- 40 min at 20° C for safety inspection
    else if( jmin < 21400 ) 
    {
      effLumi = 0.;
      effTemp = 20.;
      theta = timeScale(effTemp, T_ref);
    }

    //--- shipping to CERN (1 day at 0° C)
    else if( jmin < 22840 )
    {
      effLumi = 0.;
      effTemp = 0.;
      theta = timeScale(effTemp, T_ref);
    }
    
    //--- storage at CERN (1 week at -30° C)
    else if( jmin < 32920 )
    {
      effLumi = 0.;
      effTemp = -30.;
      theta = timeScale(effTemp, T_ref);
    }

    //--- annealing 40 min 70° C
    else if( jmin < 32960 )
    {
      effLumi = 0.;
      effTemp = 70.;
      theta = timeScale(effTemp, T_ref);
    }

    //--- measurement
    else if( jmin < 35120 )
    {
      effLumi = 0.;
      effTemp = 20.;
      theta = timeScale(effTemp, T_ref);
    }

    //--- annealing 3 days at 90° C
    else if( jmin < 39440 )
    {
      effLumi = 0.;
      effTemp = 90.;
      theta = timeScale(effTemp, T_ref);
    }
    
    //--- measurement
    else if( jmin < 41600 )
    {
      effLumi = 0.;
      effTemp = 20.;
      theta = timeScale(effTemp, T_ref);
    }

    //--- annealing 3 days at 110° C
    else if( jmin < 45920 )
    {
      effLumi = 0.;
      effTemp = 110.;
      theta = timeScale(effTemp, T_ref);
    }
    
    //--- measurement
    else if( jmin < 48080 )
    {
      effLumi = 0.;
      effTemp = 20.;
      theta = timeScale(effTemp, T_ref);
    }
    
    //--- annealing
    else
    {
      effLumi = 0.;
      effTemp = T_ann;
      theta = timeScale(effTemp, T_ref);      
    }
    

    g_instLumi_vs_time.SetPoint(g_instLumi_vs_time.GetN(), jmin/minInDay, effLumi);
    g_temp_vs_time.SetPoint(g_temp_vs_time.GetN(),jmin/minInDay,effTemp);
    
    norm += (effLumi*timeStep);
    
    for(int i = 0; i < 20; ++i)
    {
      vett[i] = alpha0[i]*tau[i]/(theta)*effLumi*(1-expl(-theta*timeStep/tau[i])) + (vett[i])*expl(-theta*timeStep/tau[i]);
      alpha += vett[i];
    }
    //std::cout << "t_month: " << t_month << " (" << (jmin+t_month)/minInYear << ") " << "   effTemp: " << effTemp << "   theta: " << theta << "   alpha: " << alpha << std::endl;
    
    g_alpha_vs_time.SetPoint(g_alpha_vs_time.GetN(),     jmin/minInDay, alpha);
    g_fluence_vs_time.SetPoint(g_fluence_vs_time.GetN(), jmin/minInDay,  norm);
    g_intLumi_vs_time.SetPoint(g_intLumi_vs_time.GetN(), jmin/minInDay,  norm);
    
    alpha = 0.;
  } // loop over 60*24*30*12 minutes in a year
  
  
  
  //---------------------------
  // normalize to total fluence
  for(int point = 0; point < g_alpha_vs_time.GetN(); ++point)
  {
    double x,y;
    g_alpha_vs_time.GetPoint(point,x,y);
    g_alphaNorm_vs_time.SetPoint(point,x,y/norm);
    
    g_fluence_vs_time.GetPoint(point,x,y);
    g_fluence_vs_time.SetPoint(point,x,y/norm*totFluence);
    
    g_intLumi_vs_time.GetPoint(point,x,y);
    g_intLumi_vs_time.SetPoint(point,x,y/norm*3000);
  }
  
  
  
  //------------
  // save graphs
  TFile* outFile = TFile::Open(Form("plots/computeAlpha_Tann_%.0f.root",T_ann),"RECREATE");
  outFile -> cd();
  
  g_alphaNorm_vs_time.GetXaxis()->SetTitle("days");
  g_alphaNorm_vs_time.GetYaxis()->SetTitle("#alpha [10^{-17} A/cm]");
  g_alphaNorm_vs_time.SetMarkerStyle(20);
  g_alphaNorm_vs_time.SetMarkerSize(0.5);
  g_alphaNorm_vs_time.Write("g_alphaNorm_vs_time");

  g_fluence_vs_time.GetXaxis()->SetTitle("days");
  g_fluence_vs_time.GetYaxis()->SetTitle("fluence [cm^{-2}]");
  g_fluence_vs_time.SetMarkerStyle(20);
  g_fluence_vs_time.SetMarkerSize(0.5);
  g_fluence_vs_time.Write("g_fluence_vs_time");
  
  g_intLumi_vs_time.GetXaxis()->SetTitle("days");
  g_intLumi_vs_time.GetYaxis()->SetTitle("int. luminosity [fb^{-1}]");
  g_intLumi_vs_time.SetMarkerStyle(20);
  g_intLumi_vs_time.SetMarkerSize(0.5);
  g_intLumi_vs_time.Write("g_intLumi_vs_time");
  
  g_instLumi_vs_time.GetXaxis()->SetTitle("days");
  g_instLumi_vs_time.GetYaxis()->SetTitle("inst. luminosity [a.u.]");
  g_instLumi_vs_time.SetMarkerStyle(20);
  g_instLumi_vs_time.SetMarkerSize(0.5);
  g_instLumi_vs_time.Write("g_instLumi_vs_time");
  
  g_temp_vs_time.GetXaxis()->SetTitle("days");
  g_temp_vs_time.GetYaxis()->SetTitle("temperature [#circ C]");
  g_temp_vs_time.SetMarkerStyle(20);
  g_temp_vs_time.SetMarkerSize(0.5);
  g_temp_vs_time.Write("g_temp_vs_time");
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  return 0;
}
