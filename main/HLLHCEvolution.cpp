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
  
  std::string SiPMType = opts.GetOpt<std::string>("Input.SiPMType");
    
  double T_ref = 21.;
  double T_ann = opts.GetOpt<double>("Input.Tann");
  double T_ann2 = opts.GetOpt<double>("Input.Tann2");
  double T_CO2 = opts.GetOpt<double>("Input.TCO2");
  
  int useTECs = opts.GetOpt<int>("Input.useTECs");
  double TECsDeltaT = opts.GetOpt<double>("Input.TECsDeltaT");
  
  int interfillAnnealing = opts.GetOpt<int>("Input.interfillAnnealing");
  double interfillTemp = opts.GetOpt<double>("Input.interfillTemp");
  
  double maxPower = opts.GetOpt<double>("Input.maxPower");
  
  double noiseTerm = opts.GetOpt<double>("Input.noiseTerm");
  double LO = opts.GetOpt<double>("Input.LO");
  double taud = opts.GetOpt<double>("Input.taud");
  double dropPDE = opts.GetOpt<double>("Input.dropPDE");
  double dropGain = opts.GetOpt<double>("Input.dropGain");

  double totFluence = opts.GetOpt<double>("Input.totFluence");

  double DCRScale = opts.GetOpt<double>("Scalings.DCRScale");
  double LOScale = opts.GetOpt<double>("Scalings.LOScale");
  double gainScale = opts.GetOpt<double>("Scalings.gainScale");
  double PDEScale = opts.GetOpt<double>("Scalings.PDEScale");
  
  double Ncells = opts.GetOpt<double>("SiPM.Ncells");
  double rechargeTime = opts.GetOpt<double>("SiPM.rechargeTime");

  std::string outputLabel = opts.GetOpt<std::string>("Input.outputLabel");
  
  
  TGraph g_instLumi_vs_time, g_alpha_vs_time, g_alphaNorm_vs_time, g_instLumidark_vs_time, g_dcr_vs_time;
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
    //if( i != 7 ) alpha0[i] = 0.;
  }
  
  long double vett[20];
  for(int j = 0; j < 20; ++j)
    vett[j] = 0;
  
  
  
  //------------------------------------------------------
  // define operating temperature (i.e. during collisions)
  double T_op = T_CO2;
  if( useTECs ) T_op -= TECsDeltaT;
  
  
  
  //----------------------
  // parse HL-LHC schedule
  std::string HLLHCScheduleLabel = opts.GetOpt<std::string>("Input.HLLHCSchedule");
  std::ifstream HLLHCSchedule;
  HLLHCSchedule.open(Form("data/%s.txt",HLLHCScheduleLabel.c_str()), std::ios::in);
  
  std::string line;
  while(true)
  {
    getline(HLLHCSchedule,line);
    if( !HLLHCSchedule.good() ) break;
    
    if( line[0] == '#' ) continue;
    
    std::stringstream ss(line);
    ss >> year >> month >> irr >> lumi;
    double t_month = (year*12+month) * minInMonth; // in minutes
    
    
    //---------------------------------------------------------------------------
    // define luminosity and temperature for each period (with month granularity)
    double theta = 1.;
    double temp;
    if( irr == 0 ) // technical stops - assume temperature T_ann
    {
      lumi = 0.;
      theta = timeScale(T_ann, T_ref);
      temp = T_ann;
    }
    
    else if( irr == 4 ) // technical stops - assume annealing temperature T_ann2
    {
      lumi = 0.;
      theta = timeScale(T_ann2, T_ref);
      temp = T_ann2;
    }
    
    else if( irr == 1 ) // beam commissioning - temperature is T_op
    {
      lumi = lumi*0.1;
      theta = timeScale(T_op, T_ref);
      temp = T_op;
    }
    
    else if( irr == 3 ) // Pb-Pb collisions - temperature is T_op
    {
      lumi = lumi*0.001;
      theta=timeScale(T_op, T_ref);
      temp = T_op;
    }
    
    else if( irr == 2 || irr == 9 || irr == 91 ) // pp collisions - temperature is T_op
    {
      theta = timeScale(T_op, T_ref);
      temp = T_op;
    }

    else if( irr == 10 || irr == 11 || irr == 12 || irr == 13 || irr == 14 ) // fake
    {
      lumi = 0.;
      theta = timeScale(20., T_ref);
      temp = 20.;
    }
    
    
    //----------------------------------------------------------------------------------
    // further update luminosity and temperature  within a month (with timestep granularity)
    for(int jmin = timeStep; jmin < minInMonth; jmin+=timeStep) // loop over 60*24*30=43200 minutes in a month 
    {
      float effLumi = lumi;
      float effTemp = temp;

      //--- if TECs, assume 6 periods of 12 hours / month (10% of a month) at higher temperature for annealing
      
      if( interfillAnnealing && irr == 2 )
      {
        if(
          ( (jmin >= 6600  && jmin <  7320) ||
            (jmin >= 13920 && jmin < 14640) ||
            (jmin >= 21240 && jmin < 21960) ||
            (jmin >= 28560 && jmin < 29280) ||
            (jmin >= 35880 && jmin < 36600) ) &&
          useTECs
          )
        {
          effLumi = 0;
          theta = timeScale(interfillTemp, T_ref);
          effTemp = interfillTemp;
        }
        else
        {
          effLumi = lumi;
          theta = timeScale(T_op, T_ref);
          effTemp = T_op;
        }
      }

      //--- assume 4 days for annealing at T_ann2 during the data taking, include interfill annealings
      
      if( irr == 9 || irr == 91 )
      {
        // 26 days with normal data taking
        if( !interfillAnnealing && jmin < 37440 )
        {
          effLumi = lumi;
          theta = timeScale(T_op, T_ref);
          effTemp = T_op;
        }
        
        if( interfillAnnealing && jmin < 37440 )
        {
          if(
            ( (jmin >= 6600  && jmin <  7320) ||
              (jmin >= 13920 && jmin < 14640) ||
              (jmin >= 21240 && jmin < 21960) ||
              (jmin >= 28560 && jmin < 29280) ||
              (jmin >= 35880 && jmin < 36600) ) &&
            useTECs
            )
          {
            effLumi = 0;
            theta = timeScale(interfillTemp, T_ref);
            effTemp = interfillTemp;
          }
          else
          {
            effLumi = lumi;
            theta = timeScale(T_op, T_ref);
            effTemp = T_op;
          }          
        }
        
        // 4 days for annealing
        if( jmin >= 37440 )
        {
          effLumi = 0.;
          if( useTECs )
          {
            if( irr == 9 )
            {
              theta = timeScale(interfillTemp, T_ref);
              effTemp = interfillTemp;
            }
            if( irr == 91 )
            {
              theta = timeScale(T_ann2, T_ref);
              effTemp = T_ann2;
            }
          }
          else
          {
            if( irr == 9 )
            {
              theta = timeScale(T_op, T_ref);
              effTemp = T_op;
            }
            if( irr == 91 )
            {
              theta = timeScale(0., T_ref);
              effTemp = 0.;
            }
          }
        }
      }

      //--- fake, emulate fast annealing in the lab
      
      if( irr == 11 ) // 40 min at 70° C
      {
        if( jmin > 43160 )
        {
          effLumi = 0.;
          theta = timeScale(70., T_ref);
          effTemp = 70.;
        }
        else
        {
          effLumi = 0.;
          theta = timeScale(20., T_ref);
          effTemp = 20.;
        }
      }

      if( irr == 12 ) // 3.5 days at 110° C
      {
        if( jmin > 38160 )
        {
          effLumi = 0.;
          theta = timeScale(110., T_ref);
          effTemp = 110.;
        }
        else
        {
          effLumi = 0.;
          theta = timeScale(20., T_ref);
          effTemp = 20.;
        }
      }
      
      if( irr == 13 ) // 4 days at 90° C
      {
        if( jmin > 37440 )
        {
          effLumi = 0.;
          theta = timeScale(90., T_ref);
          effTemp = 90.;
        }
        else
        {
          effLumi = 0.;
          theta = timeScale(20., T_ref);
          effTemp = 20.;
        }
      }

      if( irr == 14 ) // 5 days at 90° C
      {
        if( jmin > 36000 )
        {
          effLumi = 0.;
          theta = timeScale(90., T_ref);
          effTemp = 90.;
        }
        else
        {
          effLumi = 0.;
          theta = timeScale(20., T_ref);
          effTemp = 20.;
        }
      }

      //---
      
      g_instLumi_vs_time.SetPoint(g_instLumi_vs_time.GetN(), (jmin+t_month)/minInYear, effLumi);
      g_temp_vs_time.SetPoint(g_temp_vs_time.GetN(),(jmin+t_month)/minInYear,effTemp);
      
      norm += (effLumi*timeStep);
      
      for(int i = 0; i < 20; ++i)
      {
        vett[i] = alpha0[i]*tau[i]/(theta)*effLumi*(1-expl(-theta*timeStep/tau[i])) + (vett[i])*expl(-theta*timeStep/tau[i]);
        alpha += vett[i];
      }

      //std::cout << "t_month: " << t_month << " (" << (jmin+t_month)/minInYear << ") " << "   effTemp: " << effTemp << "   theta: " << theta << "   alpha: " << alpha << std::endl;

      g_alpha_vs_time.SetPoint(g_alpha_vs_time.GetN(),     (jmin+t_month)/minInYear, alpha);
      g_fluence_vs_time.SetPoint(g_fluence_vs_time.GetN(), (jmin+t_month)/minInYear,  norm);
      g_intLumi_vs_time.SetPoint(g_intLumi_vs_time.GetN(), (jmin+t_month)/minInYear,  norm);
      
      alpha = 0.;
      
    } // loop over 60*24*30 minutes in a month 
    
  } // loop over months as in the HL-LHC schedule file
  
  
  
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
  
  
  
  
  //------------------------------
  // compute optimum working point
  std::string TECsLabel = useTECs ? "_TECs" : "";
  std::string interfillLabel = interfillAnnealing ? "_interfillAnnealing" : "";
  TFile* outFile = TFile::Open(Form("plots/outFile__LO%d_tau%.1f__%s_dropPDE%.2f_dropGain%.2f%s__noise%.0f__Top_%d_Tann1_%d_Tann2_%d%s__%s__maxPower%.0f__%s.root",int(LO),taud,SiPMType.c_str(),dropPDE,dropGain,TECsLabel.c_str(),noiseTerm,int(T_op),int(T_ann),int(T_ann2),interfillLabel.c_str(),HLLHCScheduleLabel.c_str(),maxPower,outputLabel.c_str()),"RECREATE");
  
  TGraph g_tResBest_vs_time;
  TGraph g_tResBest_stoch_vs_time;
  TGraph g_tResBest_noise_vs_time;
  TGraph g_tResBest_DCR_vs_time;
  TGraph g_VbiasBest_vs_time;
  TGraph g_VovBest_vs_time;
  TGraph g_DCRBest_vs_time;
  TGraph g_PDEBest_vs_time;
  TGraph g_nPEBest_vs_time;
  TGraph g_occupancyBest_vs_time;
  TGraph g_gainBest_vs_time;
  TGraph g_SoNBest_vs_time;
  TGraph g_powerBest_vs_time;
  TGraph g_staticPowerBest_vs_time;
  TGraph g_dynamicPowerBest_vs_time;
  TGraph g_TECsPowerBest_vs_time;
  TGraph g_currentBest_vs_time;
  TGraph g_staticCurrentBest_vs_time;
  TGraph g_dynamicCurrentBest_vs_time;
  TGraph g_tempBest_vs_time;
  
  TGraph g_tResBest_vs_intLumi;
  TGraph g_tResBest_stoch_vs_intLumi;
  TGraph g_tResBest_noise_vs_intLumi;
  TGraph g_tResBest_DCR_vs_intLumi;
  TGraph g_VbiasBest_vs_intLumi;
  TGraph g_VovBest_vs_intLumi;
  TGraph g_DCRBest_vs_intLumi;
  TGraph g_PDEBest_vs_intLumi;
  TGraph g_nPEBest_vs_intLumi;
  TGraph g_occupancyBest_vs_intLumi;
  TGraph g_gainBest_vs_intLumi;
  TGraph g_SoNBest_vs_intLumi;
  TGraph g_powerBest_vs_intLumi;
  TGraph g_staticPowerBest_vs_intLumi;
  TGraph g_dynamicPowerBest_vs_intLumi;
  TGraph g_TECsPowerBest_vs_intLumi;
  TGraph g_currentBest_vs_intLumi;
  TGraph g_staticCurrentBest_vs_intLumi;
  TGraph g_dynamicCurrentBest_vs_intLumi;
  TGraph g_tempBest_vs_intLumi;
  
  TGraph g_tRes_vs_Vov_EoO;
  TGraph g_tRes_stoch_vs_Vov_EoO;
  TGraph g_tRes_noise_vs_Vov_EoO;
  TGraph g_tRes_DCR_vs_Vov_EoO;
  TGraph g_DCR_vs_Vov_EoO;
  TGraph g_PDE_vs_Vov_EoO;
  TGraph g_nPE_vs_Vov_EoO;
  TGraph g_occupancy_vs_Vov_EoO;
  TGraph g_gain_vs_Vov_EoO;
  TGraph g_power_vs_Vov_EoO;
  TGraph g_dynamicPower_vs_Vov_EoO;
  TGraph g_staticPower_vs_Vov_EoO;
  TGraph g_TECsPower_vs_Vov_EoO;
  
  TF1* f_PDE_default = new TF1("f_PDE_default","[0]*(1-exp(-1.*[1]*x))",0.,10.);
  if( SiPMType == "HPK" ) f_PDE_default -> SetParameters(1.0228*0.384,0.583);
  if( SiPMType == "FBK" ) f_PDE_default -> SetParameters(0.8847*0.466,0.314);
  TF1* f_PDE = new TF1("f_PDE","[0]*(1-exp(-1.*[1]*x))",0.,10.);
  if( SiPMType == "HPK" ) f_PDE -> SetParameters(PDEScale*1.0228*0.384,PDEScale*0.583);
  if( SiPMType == "FBK" ) f_PDE -> SetParameters(PDEScale*0.8847*0.466,PDEScale*0.314);
  //if( SiPMType == "HPK" ) f_PDE -> SetParameters(PDEScale*1.0228*0.384,0.583);
  //if( SiPMType == "FBK" ) f_PDE -> SetParameters(PDEScale*0.8847*0.466,0.314);
  
  TFile* sipmParams;
  if( SiPMType == "HPK" ) sipmParams = TFile::Open("data/sipm_spec_input_HDR2-015-v2-1e13.root","READ");
  if( SiPMType == "FBK" ) sipmParams = TFile::Open("data/sipm_spec_input_FBK-W7C-1e13.root","READ");
  TF1* f_gain = (TF1*)( sipmParams->Get("fGain_vs_OV") );
  TF1* f_ENF = (TF1*)( sipmParams->Get("fENF_vs_OV") );
  TF1* turnOn_PDE = new TF1("turnOn_PDE",Form("1.-%.2f/2E14*x",dropPDE),0.,2E14);
  TF1* turnOn_gain = new TF1("turnOn_gain",Form("1-%.2f/2E14*x",dropGain),0.,2E14);
  
  TFile* inFile_DCRParams = TFile::Open("data/DCRParams_new_Dic2021.root","READ");
  int nParDCR = 8;
  std::map<int,TGraph*> g_DCR_pars;
  for(int iPar = 0; iPar < nParDCR; ++iPar)
  {
    g_DCR_pars[iPar] = (TGraph*)( inFile_DCRParams->Get(Form("g_par%d",iPar)) );
  }
  
  TFile* inFile_TECsPower = TFile::Open("data/TECsPower.root");
  TF1* f_TECsPower_noSiPMLoad = (TF1*)( inFile_TECsPower->Get("f_noSiPMLoad") );

  std::string slewRateLabel = opts.GetOpt<std::string>("SiPM.slewRateLabel");
  TFile* inFile_slewRate_vs_amp = TFile::Open(slewRateLabel.c_str(),"READ");
  TGraph* g_slewRate_vs_amp = (TGraph*)( inFile_slewRate_vs_amp -> Get("g_SR") );
  TF1* f_slewRate_vs_amp = new TF1("f_slewRate_vs_amp",myfunc_amp,0.,10.,4);
  f_slewRate_vs_amp -> SetParameters(5.32470e-01,0.,2.92152e+01,7.79368e+00);
    
  
  float tResAvg = 0.;
  int nTResAvg = 0;
  for(int point = 0; point < g_fluence_vs_time.GetN(); ++point)
  {
    double time, fluence, alpha, instLumi, intLumi;
    g_instLumi_vs_time.GetPoint(point,time,instLumi);
    g_fluence_vs_time.GetPoint(point,time,fluence);
    g_intLumi_vs_time.GetPoint(point,time,intLumi);
    g_alphaNorm_vs_time.GetPoint(point,time,alpha);
    
    std::cout << "time: " << std::fixed << std::setprecision(3) << std::setw(5) << time
              << " y   fluence: " << std::fixed << std::scientific << std::setprecision(3) << std::setw(5) << fluence << "\r" << std::flush;
    if( instLumi == 0. ) continue;
    
    bool isEoO = false;
    if( time > 11.42 && g_tRes_vs_Vov_EoO.GetN() == 0 )
      isEoO = true;
    
    
    //------------------
    // evaluate SiPM Vbr
    float intercept = 0.;
    float slope = 0.;

    if( SiPMType == "HPK" )
    { 
      intercept = 37.95 + 1.19E-14*fluence;   // at 0° C, from Carlos data
      slope = 37.5 + 1.1E-14*fluence;         //  mV/° C, from Carlos data
    }
    if( SiPMType == "FBK" )
    {
      intercept = 32.04 + 0.28E-14*fluence;   // at 0° C, from SiPM specs document
      slope = 31.2 + 1.1E-14*fluence;         //  mV/° C, from SiPM specs document
    }
    
    float Vbr = intercept + 0.001*slope*(T_op-0.);
    
    TF1* f_DCRRef_vs_Vov = new TF1("f_DCRRef_vs_Vov",myfunc_DCR,0.,7.,nParDCR);
    float* DCRRef_pars = new float[nParDCR];
    for(int iPar = 0; iPar < nParDCR; ++iPar)
    {    
      DCRRef_pars[iPar] = g_DCR_pars[iPar]->Eval(fluence);
      f_DCRRef_vs_Vov -> SetParameter(iPar,DCRRef_pars[iPar]);
    }
    //float DCRRef = f_DCRRef_vs_Vov -> Eval(1.); // DCR at 1 V and - 30° C as per Carlos data
    float DCRRef = f_DCRRef_vs_Vov -> Eval(0.8); // DCR at 0.8 V and - 40° C as per TB data
    
    float tResBest = 999999.;
    float tResBest_stoch = 999999.;
    float tResBest_noise = 999999.;
    float tResBest_DCR = 999999.;
    float VbiasBest = -1.;
    float VovBest = -1.;
    float DCRBest = -1.;
    float PDEBest = 0.;
    float nPEBest = -1.;
    float occupancyBest = 0.;
    float gainBest = -1.;
    float SoNBest = -1.;
    float powerBest = 0.;
    float dynamicPowerBest = 0.;
    float staticPowerBest = 0.;
    float TECsPowerBest = 0.;
    float currentBest = 0.;
    float dynamicCurrentBest = 0.;
    float staticCurrentBest = 0.;
    float tempBest = 999.;
    
    //for(float effTECsDeltaT = 0; effTECsDeltaT < 20; effTECsDeltaT+=1)
    float effTECsDeltaT = TECsDeltaT;
    for(float Vov = 0.2; Vov < 5.; Vov += 0.01)
    {
      T_op = T_CO2 - effTECsDeltaT;
      
      //-----------------------------------------------
      // evaluate effective DCR for a given OV and T_op
      float DCR = ( 16300./8. * alpha * 1.E-17 * totFluence ) * // from HPK2E14 used at TB (after 3.5d at 110° C, alpha is 0.89, assume 8% gain loss (i.e. 29. GHz at Vov = 1.5 V, DCR(1.5)/DCR(0.8)~8)
                  (f_DCRRef_vs_Vov->Eval(Vov)/DCRRef) *         // morphing vs. OV using CERN Oct. TB data
                  ( DCRScale ) *
                  f_PDE->Eval(Vov)/f_PDE_default->Eval(Vov);
        
      if( SiPMType == "FBK" ) DCR *= 1.10;
      
      float B = -0.00416498*alpha + 0.0798623;   // DCR scaling with temperature, including dependence of scaling factor on alpha
      //DCR = DCR * exp(B*(T_op-(-35.)));
      DCR = DCR * exp(B*(T_op-(-40.)));
      

      // cell occupancy due to DCR -- assuming here 2 tau_R
      float occupancy = 2. * rechargeTime * DCR / Ncells;
      
      
      //-----------------------------------------
      // evaluate total power dissipated per SiPM
      float staticCurrent = DCR*1E09 * f_ENF->Eval(Vov) * gainScale * f_gain->Eval(Vov)*(turnOn_gain->Eval(fluence)) * 1.602E-19;
      float staticPower = staticCurrent * (Vbr + Vov) * 1000.;   // in mW
      if( staticPower > maxPower ) continue;
      
      float dynamicCurrent = (2.3*1E06 * instLumi * 5/7.5) * (4.2 * LO * LOScale * (1-occupancy) * f_PDE->Eval(Vov)*(turnOn_PDE->Eval(fluence))/f_PDE_default->Eval(3.5)) * f_ENF->Eval(Vov) * gainScale * f_gain->Eval(Vov)*(turnOn_gain->Eval(fluence)) * 1.602E-19;   // 2.3 MHz equivalent MIP rate at 200 PU.
      float dynamicPower = dynamicCurrent * (Vbr + Vov) * 1000.;   // in mW
      if( (staticPower+dynamicPower) > maxPower ) continue;

      float TECsPower = 0.; // in mW per channel
      if( useTECs )
      {
        f_TECsPower_noSiPMLoad -> SetParameter(0.,(4./420.)*16.*(staticPower+dynamicPower));
        TECsPower = f_TECsPower_noSiPMLoad -> Eval(-1.*effTECsDeltaT) / 16.;
      }
      if( (staticPower+dynamicPower+TECsPower) > maxPower ) continue;
      
      
      //-------------------------
      // evaluate time resolution
      float nPE = 4.2 * LO * LOScale * (1-occupancy) * f_PDE->Eval(Vov)*(turnOn_PDE->Eval(fluence))/f_PDE_default->Eval(3.5);
      
      float sigma_stoch = 27. * sqrt(7000./(nPE*38.5/taud));
      //float sigma_noise = sqrt( pow(noiseTerm/1.2/g_slewRate_vs_amp->Eval(nPE*gainScale*f_gain->Eval(Vov)*(turnOn_gain->Eval(fluence))/f_gain->Eval(3.5)/9500.),2) + pow(16.7,2) )/sqrt(2);
      float sigma_noise = sqrt( pow(noiseTerm/1.2/g_slewRate_vs_amp->Eval(nPE*gainScale*f_gain->Eval(Vov)*(turnOn_gain->Eval(fluence))/f_gain->Eval(3.5)/9500.),2) + pow(16.7,2) )/sqrt(2);
      float sigma_DCR   = 40. * 6000./(nPE*38.5/taud) * pow(DCR/30.,0.41);
      float sigma_clock = 15.;
      
      float tResCurr = sqrt( pow(sigma_stoch,2) + pow(sigma_noise,2) + pow(sigma_DCR,2) + pow(sigma_clock,2) );
      
      if( tResCurr < tResBest )
      {
        tResBest = tResCurr;
        tResBest_stoch = sigma_stoch;
        tResBest_noise = sigma_noise;
        tResBest_DCR = sigma_DCR;
        VbiasBest = Vov + Vbr;
        VovBest = Vov;
        DCRBest = DCR;
        PDEBest = f_PDE->Eval(Vov)*(turnOn_PDE->Eval(fluence));
        nPEBest = nPE;
        occupancyBest = occupancy;
        gainBest = gainScale*f_gain->Eval(Vov)*(turnOn_gain->Eval(fluence));
        SoNBest = f_PDE->Eval(Vov)*(turnOn_PDE->Eval(fluence))/sqrt(DCR);
        powerBest = staticPower+dynamicPower+TECsPower;
        dynamicPowerBest = dynamicPower;
        staticPowerBest = staticPower;
        TECsPowerBest = TECsPower;
        currentBest = (staticPower+dynamicPower) / (Vbr + Vov);
        dynamicCurrentBest = (dynamicPower) / (Vbr + Vov);
        staticCurrentBest = (staticPower) / (Vbr + Vov);
        tempBest = T_op;
      }

      if( isEoO )
      {
        g_tRes_vs_Vov_EoO.SetPoint(g_tRes_vs_Vov_EoO.GetN(),Vov,tResCurr);
        g_tRes_stoch_vs_Vov_EoO.SetPoint(g_tRes_stoch_vs_Vov_EoO.GetN(),Vov,sigma_stoch);
        g_tRes_noise_vs_Vov_EoO.SetPoint(g_tRes_noise_vs_Vov_EoO.GetN(),Vov,sigma_noise);
        g_tRes_DCR_vs_Vov_EoO.SetPoint(g_tRes_DCR_vs_Vov_EoO.GetN(),Vov,sigma_DCR);
        g_DCR_vs_Vov_EoO.SetPoint(g_DCR_vs_Vov_EoO.GetN(),Vov,DCR);
        g_PDE_vs_Vov_EoO.SetPoint(g_PDE_vs_Vov_EoO.GetN(),Vov,f_PDE->Eval(Vov)*(turnOn_PDE->Eval(fluence)));
        g_nPE_vs_Vov_EoO.SetPoint(g_nPE_vs_Vov_EoO.GetN(),Vov,nPE);
        g_occupancy_vs_Vov_EoO.SetPoint(g_occupancy_vs_Vov_EoO.GetN(),Vov,occupancy);
        g_gain_vs_Vov_EoO.SetPoint(g_gain_vs_Vov_EoO.GetN(),Vov,gainScale*f_gain->Eval(Vov));
        g_power_vs_Vov_EoO.SetPoint(g_power_vs_Vov_EoO.GetN(),Vov,staticPower+dynamicPower+TECsPower);
        g_dynamicPower_vs_Vov_EoO.SetPoint(g_dynamicPower_vs_Vov_EoO.GetN(),Vov,dynamicPower);
        g_staticPower_vs_Vov_EoO.SetPoint(g_staticPower_vs_Vov_EoO.GetN(),Vov,staticPower);
        g_TECsPower_vs_Vov_EoO.SetPoint(g_TECsPower_vs_Vov_EoO.GetN(),Vov,TECsPower);
      }
    }
    
    //std::cout << "Vov: " << VovBest << "   nPEBest: " << nPEBest << "   gainBest: " << gainBest << "   DCRBest: " << DCRBest << "   tResBest: " << tResBest << std::endl;
    g_tResBest_vs_time.SetPoint(g_tResBest_vs_time.GetN(),time,tResBest);
    g_tResBest_stoch_vs_time.SetPoint(g_tResBest_stoch_vs_time.GetN(),time,tResBest_stoch);
    g_tResBest_noise_vs_time.SetPoint(g_tResBest_noise_vs_time.GetN(),time,tResBest_noise);
    g_tResBest_DCR_vs_time.SetPoint(g_tResBest_DCR_vs_time.GetN(),time,tResBest_DCR);
    g_VbiasBest_vs_time.SetPoint(g_VbiasBest_vs_time.GetN(),time,VbiasBest);
    g_VovBest_vs_time.SetPoint(g_VovBest_vs_time.GetN(),time,VovBest);
    g_DCRBest_vs_time.SetPoint(g_DCRBest_vs_time.GetN(),time,DCRBest);
    g_PDEBest_vs_time.SetPoint(g_PDEBest_vs_time.GetN(),time,PDEBest);
    g_nPEBest_vs_time.SetPoint(g_nPEBest_vs_time.GetN(),time,nPEBest);
    g_occupancyBest_vs_time.SetPoint(g_occupancyBest_vs_time.GetN(),time,occupancyBest);
    g_gainBest_vs_time.SetPoint(g_gainBest_vs_time.GetN(),time,gainBest);
    g_SoNBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,SoNBest);
    g_powerBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,powerBest);
    g_dynamicPowerBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,dynamicPowerBest);
    g_staticPowerBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,staticPowerBest);
    g_TECsPowerBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,TECsPowerBest);
    g_currentBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,currentBest);
    g_dynamicCurrentBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,dynamicCurrentBest);
    g_staticCurrentBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,staticCurrentBest);
    g_tempBest_vs_time.SetPoint(g_tempBest_vs_time.GetN(),time,tempBest);
    
    g_tResBest_vs_intLumi.SetPoint(g_tResBest_vs_intLumi.GetN(),intLumi,tResBest);
    g_tResBest_stoch_vs_intLumi.SetPoint(g_tResBest_vs_intLumi.GetN(),intLumi,tResBest_stoch);
    g_tResBest_noise_vs_intLumi.SetPoint(g_tResBest_vs_intLumi.GetN(),intLumi,tResBest_noise);
    g_tResBest_DCR_vs_intLumi.SetPoint(g_tResBest_vs_intLumi.GetN(),intLumi,tResBest_DCR);
    g_VbiasBest_vs_intLumi.SetPoint(g_VbiasBest_vs_intLumi.GetN(),intLumi,VbiasBest);
    g_VovBest_vs_intLumi.SetPoint(g_VovBest_vs_intLumi.GetN(),intLumi,VovBest);
    g_DCRBest_vs_intLumi.SetPoint(g_DCRBest_vs_intLumi.GetN(),intLumi,DCRBest);
    g_PDEBest_vs_intLumi.SetPoint(g_PDEBest_vs_intLumi.GetN(),intLumi,PDEBest);
    g_nPEBest_vs_intLumi.SetPoint(g_nPEBest_vs_intLumi.GetN(),intLumi,nPEBest);
    g_occupancyBest_vs_intLumi.SetPoint(g_occupancyBest_vs_intLumi.GetN(),intLumi,occupancyBest);
    g_gainBest_vs_intLumi.SetPoint(g_gainBest_vs_intLumi.GetN(),intLumi,gainBest);
    g_SoNBest_vs_intLumi.SetPoint(g_SoNBest_vs_intLumi.GetN(),intLumi,SoNBest);
    g_powerBest_vs_intLumi.SetPoint(g_powerBest_vs_intLumi.GetN(),intLumi,powerBest);
    g_dynamicPowerBest_vs_intLumi.SetPoint(g_dynamicPowerBest_vs_intLumi.GetN(),intLumi,dynamicPowerBest);
    g_staticPowerBest_vs_intLumi.SetPoint(g_staticPowerBest_vs_intLumi.GetN(),intLumi,staticPowerBest);
    g_TECsPowerBest_vs_intLumi.SetPoint(g_TECsPowerBest_vs_intLumi.GetN(),intLumi,TECsPowerBest);
    g_currentBest_vs_intLumi.SetPoint(g_currentBest_vs_intLumi.GetN(),intLumi,currentBest);
    g_dynamicCurrentBest_vs_intLumi.SetPoint(g_dynamicCurrentBest_vs_intLumi.GetN(),intLumi,dynamicCurrentBest);
    g_staticCurrentBest_vs_intLumi.SetPoint(g_staticCurrentBest_vs_intLumi.GetN(),intLumi,staticCurrentBest);
    g_tempBest_vs_intLumi.SetPoint(g_tempBest_vs_intLumi.GetN(),intLumi,tempBest);
    
    tResAvg += tResBest;
    ++nTResAvg;
  }
  
  tResAvg /= nTResAvg;
  
  


  //------------
  // save graphs
  outFile -> cd();
  
  g_alphaNorm_vs_time.GetXaxis()->SetTitle("years from 2027");
  g_alphaNorm_vs_time.GetYaxis()->SetTitle("#alpha [10^{-17} A/cm]");
  g_alphaNorm_vs_time.SetMarkerStyle(20);
  g_alphaNorm_vs_time.SetMarkerSize(0.5);
  g_alphaNorm_vs_time.Write("g_alphaNorm_vs_time");

  g_fluence_vs_time.GetXaxis()->SetTitle("years from 2027");
  g_fluence_vs_time.GetYaxis()->SetTitle("fluence [cm^{-2}]");
  g_fluence_vs_time.SetMarkerStyle(20);
  g_fluence_vs_time.SetMarkerSize(0.5);
  g_fluence_vs_time.Write("g_fluence_vs_time");
  
  g_intLumi_vs_time.GetXaxis()->SetTitle("years from 2027");
  g_intLumi_vs_time.GetYaxis()->SetTitle("int. luminosity [fb^{-1}]");
  g_intLumi_vs_time.SetMarkerStyle(20);
  g_intLumi_vs_time.SetMarkerSize(0.5);
  g_intLumi_vs_time.Write("g_intLumi_vs_time");
  
  g_instLumi_vs_time.GetXaxis()->SetTitle("years from 2027");
  g_instLumi_vs_time.GetYaxis()->SetTitle("inst. luminosity [a.u.]");
  g_instLumi_vs_time.SetMarkerStyle(20);
  g_instLumi_vs_time.SetMarkerSize(0.5);
  g_instLumi_vs_time.Write("g_instLumi_vs_time");
  
  g_temp_vs_time.GetXaxis()->SetTitle("years from 2027");
  g_temp_vs_time.GetYaxis()->SetTitle("temperature [#circ C]");
  g_temp_vs_time.SetMarkerStyle(20);
  g_temp_vs_time.SetMarkerSize(0.5);
  g_temp_vs_time.Write("g_temp_vs_time");
  
  g_tResBest_vs_time.SetTitle(";years from 2027;#sigma_{t}^{best} [ps]");
  g_tResBest_vs_time.Write("g_tResBest_vs_time");
  g_tResBest_stoch_vs_time.SetTitle(";years from 2027;#sigma_{t}^{best} [ps]");
  g_tResBest_stoch_vs_time.Write("g_tResBest_stoch_vs_time");
  g_tResBest_noise_vs_time.SetTitle(";years from 2027;#sigma_{t}^{best} [ps]");
  g_tResBest_noise_vs_time.Write("g_tResBest_noise_vs_time");
  g_tResBest_DCR_vs_time.SetTitle(";years from 2027;#sigma_{t}^{best} [ps]");
  g_tResBest_DCR_vs_time.Write("g_tResBest_DCR_vs_time");
  g_VbiasBest_vs_time.SetTitle(";years from 2027;V_{bias} at best #sigma_{t} [V]");
  g_VbiasBest_vs_time.Write("g_VbiasBest_vs_time");
  g_VovBest_vs_time.SetTitle(";years from 2027;V_{OV} at best #sigma_{t} [V]");
  g_VovBest_vs_time.Write("g_VovBest_vs_time");
  g_DCRBest_vs_time.SetTitle(";years from 2027;DCR at best #sigma_{t} [GHz]");
  g_DCRBest_vs_time.Write("g_DCRBest_vs_time");
  g_PDEBest_vs_time.SetTitle(";years from 2027;PDE at best #sigma_{t}");
  g_PDEBest_vs_time.Write("g_PDEBest_vs_time");
  g_nPEBest_vs_time.SetTitle(";years from 2027;N_{p.e.} at best #sigma_{t}");
  g_nPEBest_vs_time.Write("g_nPEBest_vs_time");
  g_occupancyBest_vs_time.SetTitle(";years from 2027;SiPM occupancy at best #sigma_{t}");
  g_occupancyBest_vs_time.Write("g_occupancyBest_vs_time");
  g_gainBest_vs_time.SetTitle(";years from 2027;SiPM gain at best #sigma_{t}");
  g_gainBest_vs_time.Write("g_gainBest_vs_time");
  g_SoNBest_vs_time.SetTitle(";years from 2027;S/N at best #sigma_{t}");
  g_SoNBest_vs_time.Write("g_SoNBest_vs_time");
  g_powerBest_vs_time.SetTitle(";years from 2027;total power per ch. at best #sigma_{t} [mW]");
  g_powerBest_vs_time.Write("g_powerBest_vs_time");
  g_dynamicPowerBest_vs_time.SetTitle(";years from 2027;dynamic power per ch. at best #sigma_{t} [mW]");
  g_dynamicPowerBest_vs_time.Write("g_dynamicPowerBest_vs_time");
  g_staticPowerBest_vs_time.SetTitle(";years from 2027;static power per ch. at best #sigma_{t} [mW]");
  g_staticPowerBest_vs_time.Write("g_staticPowerBest_vs_time");
  g_TECsPowerBest_vs_time.SetTitle(";years from 2027;TECs power per ch. at best #sigma_{t} [mW]");
  g_TECsPowerBest_vs_time.Write("g_TECsPowerBest_vs_time");
  g_currentBest_vs_time.SetTitle(";years from 2027;total current per ch. at best #sigma_{t} [mA]");
  g_currentBest_vs_time.Write("g_currentBest_vs_time");
  g_dynamicCurrentBest_vs_time.SetTitle(";years from 2027;dynamic current per ch. at best #sigma_{t} [mA]");
  g_dynamicCurrentBest_vs_time.Write("g_dynamicCurrentBest_vs_time");
  g_staticCurrentBest_vs_time.SetTitle(";years from 2027;static current per ch. at best #sigma_{t} [mA]");
  g_staticCurrentBest_vs_time.Write("g_staticCurrentBest_vs_time");
  g_tempBest_vs_time.SetTitle(";years from 2027;temperature at best #sigma_{t} [mA]");
  g_tempBest_vs_time.Write("g_tempBest_vs_time");
  
  
  g_tResBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];#sigma_{t}^{best} [ps]");
  g_tResBest_vs_intLumi.Write("g_tResBest_vs_intLumi");
  g_tResBest_stoch_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];#sigma_{t}^{best} [ps]");
  g_tResBest_stoch_vs_intLumi.Write("g_tResBest_stoch_vs_intLumi");
  g_tResBest_noise_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];#sigma_{t}^{best} [ps]");
  g_tResBest_noise_vs_intLumi.Write("g_tResBest_noise_vs_intLumi");
  g_tResBest_DCR_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];#sigma_{t}^{best} [ps]");
  g_tResBest_DCR_vs_intLumi.Write("g_tResBest_DCR_vs_intLumi");
  g_VbiasBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];V_{OV} at best #sigma_{t} [V]");
  g_VbiasBest_vs_intLumi.Write("g_VbiasBest_vs_intLumi");
  g_VovBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];V_{OV} at best #sigma_{t} [V]");
  g_VovBest_vs_intLumi.Write("g_VovBest_vs_intLumi");
  g_DCRBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];DCR at best #sigma_{t} [GHz]");
  g_DCRBest_vs_intLumi.Write("g_DCRBest_vs_intLumi");
  g_PDEBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];PDE at best #sigma_{t}");
  g_PDEBest_vs_intLumi.Write("g_PDEBest_vs_intLumi");
  g_nPEBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];N_{p.e.} at best #sigma_{t}");
  g_nPEBest_vs_intLumi.Write("g_nPEBest_vs_intLumi");
  g_occupancyBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];SiPM occupancy at best #sigma_{t}");
  g_occupancyBest_vs_intLumi.Write("g_occupancyBest_vs_intLumi");
  g_gainBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];SiPM gain at best #sigma_{t}");
  g_gainBest_vs_intLumi.Write("g_gainBest_vs_intLumi");
  g_SoNBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];S/N at best #sigma_{t}");
  g_SoNBest_vs_intLumi.Write("g_SoNBest_vs_intLumi");
  g_powerBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];total power per ch. at best #sigma_{t} [mW]");
  g_powerBest_vs_intLumi.Write("g_powerBest_vs_intLumi");
  g_dynamicPowerBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];dynamic power per ch. at best #sigma_{t} [mW]");
  g_dynamicPowerBest_vs_intLumi.Write("g_dynamicPowerBest_vs_intLumi");
  g_staticPowerBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];static power per ch. at best #sigma_{t} [mW]");
  g_staticPowerBest_vs_intLumi.Write("g_staticPowerBest_vs_intLumi");
  g_TECsPowerBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];TECs power per ch. at best #sigma_{t} [mW]");
  g_TECsPowerBest_vs_intLumi.Write("g_TECsPowerBest_vs_intLumi");
  g_currentBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];total current per ch. at best #sigma_{t} [mA]");
  g_currentBest_vs_intLumi.Write("g_currentBest_vs_intLumi");
  g_dynamicCurrentBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];dynamic current per ch. at best #sigma_{t} [mA]");
  g_dynamicCurrentBest_vs_intLumi.Write("g_dynamicCurrentBest_vs_intLumi");
  g_staticCurrentBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];static current per ch. at best #sigma_{t} [mA]");
  g_staticCurrentBest_vs_intLumi.Write("g_staticCurrentBest_vs_intLumi");
  g_tempBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];temp at best #sigma_{t} [mA]");
  g_tempBest_vs_intLumi.Write("g_tempBest_vs_intLumi");
  
  g_tRes_vs_Vov_EoO.Write("g_tRes_vs_Vov_EoO");
  g_tRes_stoch_vs_Vov_EoO.Write("g_tRes_stoch_vs_Vov_EoO");
  g_tRes_noise_vs_Vov_EoO.Write("g_tRes_noise_vs_Vov_EoO");
  g_tRes_DCR_vs_Vov_EoO.Write("g_tRes_DCR_vs_Vov_EoO");
  g_DCR_vs_Vov_EoO.Write("g_DCR_vs_Vov_EoO");
  g_PDE_vs_Vov_EoO.Write("g_PDE_vs_Vov_EoO");
  g_nPE_vs_Vov_EoO.Write("g_nPE_vs_Vov_EoO");
  g_occupancy_vs_Vov_EoO.Write("g_occupancy_vs_Vov_EoO");
  g_gain_vs_Vov_EoO.Write("g_gain_vs_Vov_EoO");
  g_power_vs_Vov_EoO.Write("g_power_vs_Vov_EoO");
  g_dynamicPower_vs_Vov_EoO.Write("g_dynamicPower_vs_Vov_EoO");
  g_staticPower_vs_Vov_EoO.Write("g_staticPower_vs_Vov_EoO");
  g_TECsPower_vs_Vov_EoO.Write("g_TECsPower_vs_Vov_EoO");
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  return 0;
}
