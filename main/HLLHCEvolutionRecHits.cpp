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
#include "TTree.h"

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
    
    
    //--- define luminosity and temperature for each period (with month granularity)
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
    
    
    //--- further update luminosity and temperature  within a month (with timestep granularity)
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
  
  
  
  
  //--------------------------
  // define various parameters
  std::string TECsLabel = useTECs ? "_TECs" : "";
  std::string interfillLabel = interfillAnnealing ? "_interfillAnnealing" : "";
  TFile* outFile = TFile::Open(Form("plots/outFileRecHits__LO%d_tau%.1f__%s_dropPDE%.2f_dropGain%.2f%s__noise%.0f__Top_%d_Tann1_%d_Tann2_%d%s__%s__maxPower%.0f__%s.root",int(LO),taud,SiPMType.c_str(),dropPDE,dropGain,TECsLabel.c_str(),noiseTerm,int(T_op),int(T_ann),int(T_ann2),interfillLabel.c_str(),HLLHCScheduleLabel.c_str(),maxPower,outputLabel.c_str()),"RECREATE");
  
  TF1* f_PDE_default = new TF1("f_PDE_default","[0]*(1-exp(-1.*[1]*x))",0.,10.);
  if( SiPMType == "HPK" ) f_PDE_default -> SetParameters(1.0228*0.384,0.583);
  if( SiPMType == "FBK" ) f_PDE_default -> SetParameters(0.8847*0.466,0.314);
  TF1* f_PDE = new TF1("f_PDE","[0]*(1-exp(-1.*[1]*x))",0.,10.);
  if( SiPMType == "HPK" ) f_PDE -> SetParameters(PDEScale*1.0228*0.384,PDEScale*0.583);
  if( SiPMType == "FBK" ) f_PDE -> SetParameters(PDEScale*0.8847*0.466,PDEScale*0.314);
  
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
  
  std::vector<std::string> slewRateLabel = opts.GetOpt<std::vector<std::string> >("SiPM.slewRateLabel");
  std::map<int,TGraph*> g_slewRate_vs_amp;
  for(int RUIt = 1; RUIt <= 6; ++RUIt)
  {
    TFile* inFile_slewRate_vs_amp = TFile::Open((slewRateLabel.at(RUIt-1)).c_str(),"READ");
    g_slewRate_vs_amp[RUIt] = (TGraph*)( inFile_slewRate_vs_amp -> Get("g_SR") );
    inFile_slewRate_vs_amp -> Close();
  }
  
  std::vector<int> types = opts.GetOpt<std::vector<int> >("RUTypes.types");
  std::vector<float> loscales = opts.GetOpt<std::vector<float> >("RUTypes.LOScales");
  std::map<int,int> RUTypes;
  std::map<int,float> LOScales;
  for(int RUIt = 1; RUIt <= 6; ++RUIt)
  {
    RUTypes[RUIt] = types.at(RUIt-1);
    LOScales[RUIt] = loscales.at(RUIt-1);
  }

  TF1* f_fluence_vs_eta = new TF1("f_fluence_vs_eta","pol2",0.,1.5);
  f_fluence_vs_eta -> SetParameters(1.66071E14,1.08123e+11,1.0406e+13);
  
  

  
  //------------
  // define time
  double v_time, v_fluence, v_alpha, v_instLumi, v_intLumi;
  for(int point = 0; point < g_fluence_vs_time.GetN(); ++point)
  {
    g_instLumi_vs_time.GetPoint(point,v_time,v_instLumi);
    g_fluence_vs_time.GetPoint(point,v_time,v_fluence);
    g_intLumi_vs_time.GetPoint(point,v_time,v_intLumi);
    g_alphaNorm_vs_time.GetPoint(point,v_time,v_alpha);
    
    if( v_instLumi == 0. ) continue;
    
    if( v_time > 11.42 ) break;
  }    
  std::cout << "time: " << std::fixed << std::setprecision(3) << std::setw(5) << v_time << " y   "
            << "fluence: " << std::fixed << std::scientific << std::setprecision(3) << std::setw(5) << v_fluence
            << std::endl;
  

  //---------------------------
  // get the tracks and recHits
  TFile* inFile_recHits = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/data/TDR/CMSSW_10_4_0_mtd5/ntuple_MinBias_14TeV_barphiflat.root","READ");
  TTree* tree = (TTree*)( inFile_recHits->Get("FTLDumpHits/DumpHits") );
  int nEntries = tree->GetEntries();
  
  tree -> SetBranchStatus("*",0);
  
  std::vector<float>* tracks_pt  = new std::vector<float>;
  std::vector<float>* tracks_eta = new std::vector<float>;
  std::vector<float>* tracks_phi = new std::vector<float>;
  std::vector<float>* tracks_mcMatch_genPt = new std::vector<float>;
  std::vector<float>* tracks_mcMatch_DR = new std::vector<float>;
  tree -> SetBranchStatus("track_pt",           1); tree -> SetBranchAddress("track_pt",           &tracks_pt);
  tree -> SetBranchStatus("track_eta",          1); tree -> SetBranchAddress("track_eta",          &tracks_eta);
  tree -> SetBranchStatus("track_phi",          1); tree -> SetBranchAddress("track_phi",          &tracks_phi);
  tree -> SetBranchStatus("track_mcMatch_genPt",1); tree -> SetBranchAddress("track_mcMatch_genPt",&tracks_mcMatch_genPt);
  tree -> SetBranchStatus("track_mcMatch_DR",   1); tree -> SetBranchAddress("track_mcMatch_DR",   &tracks_mcMatch_DR);
  
  std::vector<std::vector<int> >*   matchedRecHits_det = new std::vector<std::vector<int> >;
  std::vector<std::vector<float> >* matchedRecHits_energy = new std::vector<std::vector<float> >;
  tree -> SetBranchStatus("matchedRecHits_det",1);        tree -> SetBranchAddress("matchedRecHits_det",        &matchedRecHits_det);
  tree -> SetBranchStatus("matchedRecHits_energy",1);     tree -> SetBranchAddress("matchedRecHits_energy",     &matchedRecHits_energy);

  
  //------------------
  // evaluate SiPM Vbr
  float intercept = 0.;
  float slope = 0.;
  
  if( SiPMType == "HPK" )
  { 
    intercept = 37.95 + 1.19E-14*v_fluence;   // at 0° C, from Carlos data
    slope = 37.5 + 1.1E-14*v_fluence;         //  mV/° C, from Carlos data
  }
  if( SiPMType == "FBK" )
  {
    intercept = 32.04 + 0.28E-14*v_fluence;   // at 0° C, from SiPM specs document
    slope = 31.2 + 1.1E-14*v_fluence;         //  mV/° C, from SiPM specs document
  }
  
  TF1* f_DCRRef_vs_Vov = new TF1("f_DCRRef_vs_Vov",myfunc_DCR,0.,7.,nParDCR);
  float* DCRRef_pars = new float[nParDCR];
  for(int iPar = 0; iPar < nParDCR; ++iPar)
  {    
    DCRRef_pars[iPar] = g_DCR_pars[iPar]->Eval(v_fluence);
    f_DCRRef_vs_Vov -> SetParameter(iPar,DCRRef_pars[iPar]);
  }
  //float DCRRef = f_DCRRef_vs_Vov -> Eval(1.); // DCR at 1 V and - 30° C as per Carlos data
  float DCRRef = f_DCRRef_vs_Vov -> Eval(0.8); // DCR at 0.8 V and - 40° C as per TB data


  //----------------------
  // optimize Vov and temp
  std::map<int,float> tResBest_vs_RU;
  std::map<int,std::pair<float,float> > condBest_vs_RU;
  for(int RUIt = 1; RUIt <= 6; ++RUIt)
  {
    tResBest_vs_RU[RUIt] = 999999.;
    condBest_vs_RU[RUIt] = std::make_pair<float,float>(-999.,-999.);
  }
  
  
  //float effTECsDeltaT = TECsDeltaT;
  for(float effTECsDeltaT = 5; effTECsDeltaT <= 15; effTECsDeltaT+=0.5)
  {
    std::cout << ">>> effTECsDeltaT: " << effTECsDeltaT << std::endl;
    
    for(float Vov = 0.2; Vov < 5.; Vov += 0.01)
    {
      std::cout << ">>>>>> Vov: " << Vov << "\r" << std::flush;
      
      T_op = T_CO2 - effTECsDeltaT;
      if( T_op < -50. ) continue;
      
      
      //-----------------------------------------------
      // evaluate effective DCR for a given OV and T_op
      float DCR = ( 16300./8. * v_alpha * 1.E-17 * totFluence ) * // from HPK2E14 used at TB (after 3.5d at 110° C, alpha is 0.89, assume 8% gain loss (i.e. 29. GHz at Vov = 1.5 V, DCR(1.5)/DCR(0.8)~8)
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
      float Vbr = intercept + 0.001*slope*(T_op-0.);
      
      float staticCurrent = DCR*1E09 * f_ENF->Eval(Vov) * gainScale * f_gain->Eval(Vov)*(turnOn_gain->Eval(v_fluence)) * 1.602E-19;
      float staticPower = staticCurrent * (Vbr + Vov) * 1000.;   // in mW
      if( staticPower > maxPower ) continue;
      if( staticCurrent*1000. > 50./16. ) continue; // current in mA, 50 mA is the ALDO limit
      if( staticCurrent > 2./48./16. ) continue; // current in A, 2A is the PS limit, as per Krzysztof's plots
      
      float dynamicCurrent = (2.3*1E06 * v_instLumi * 5/7.5) * (4.2 * LO * LOScale * (1-occupancy) * f_PDE->Eval(Vov)*(turnOn_PDE->Eval(v_fluence))/f_PDE_default->Eval(3.5)) * f_ENF->Eval(Vov) * gainScale * f_gain->Eval(Vov)*(turnOn_gain->Eval(v_fluence)) * 1.602E-19;   // 2.3 MHz equivalent MIP rate at 200 PU.
      float dynamicPower = dynamicCurrent * (Vbr + Vov) * 1000.;   // in mW
      if( (staticPower+dynamicPower) > maxPower ) continue;
      if( (staticCurrent+dynamicCurrent)*1000. > 50./16. ) continue; // current in mA, 50 mA is the ALDO limit
      if( (staticCurrent+dynamicCurrent) > 2./48./16. ) continue; // current in A, 2A is the PS limit, as per Krzysztof's plots
      
      float TECsPower = 0.; // in mW per channel
      if( useTECs )
      {
        f_TECsPower_noSiPMLoad -> SetParameter(0.,(4./420.)*16.*(staticPower+dynamicPower));
        TECsPower = f_TECsPower_noSiPMLoad -> Eval(-1.*effTECsDeltaT) / 16.;
      }
      if( TECsPower > 125. ) continue; // limit to 2W/array as per Krzysztof's plots, with some margin
      
      if( (staticPower+dynamicPower+TECsPower) > maxPower ) continue;
      
      
      //-------------------------
      // evaluate time resolution
      TProfile* p_tRes_vs_RU = new TProfile("p_tRes_vs_RU","",6,0.5,6.5);
      
      for(int entry = 0; entry < 10000; ++entry)
      {
        //if( entry%1000 == 0 ) std::cout << ">>> reading entry " << entry << " / " << nEntries << "\r" << std::endl;
        tree -> GetEntry(entry);
        
        for(unsigned int trackIt = 0; trackIt < tracks_pt->size(); ++trackIt)
        {
          //std::cout << ">>>>>> trackIt " << trackIt << " / " << tracks_pt->size() << std::endl;
          
          float pt = tracks_pt->at(trackIt);
          float feta = fabs(tracks_eta->at(trackIt));
          float genPt = tracks_mcMatch_genPt->at(trackIt);
          float DR = tracks_mcMatch_DR->at(trackIt);
          
          if( DR > 0.01 ) continue;
          if( fabs(pt/genPt-1.) > 0.05 ) continue; 
          
          int RUId = -1.;
          float DCRCorr = 1.;
          float LOCorr = 1.;
          if( feta > 0.000 && feta < 0.347 ) { RUId = 1; };
          if( feta > 0.347 && feta < 0.661 ) { RUId = 2; };
          if( feta > 0.661 && feta < 0.927 ) { RUId = 3; };
          if( feta > 0.927 && feta < 1.148 ) { RUId = 4; };
          if( feta > 1.148 && feta < 1.338 ) { RUId = 5; };
          if( feta > 1.338 && feta < 1.498 ) { RUId = 6; };
          
          if( RUId <= 0. ) continue;
          
          if( RUTypes[RUId] == 1 ) { DCRCorr *= 1.25; LOCorr *= 1.05; }
          if( RUTypes[RUId] == 3 ) { DCRCorr *= 0.8;  LOCorr *= 0.95; }
          LOCorr *= LOScales[RUId];
          
          DCRCorr *= f_fluence_vs_eta -> Eval(feta) / 2.E14;
          
          std::vector<float> Npes;
          std::vector<float> recHitEs;
          float totEnergy = 0;
          for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(trackIt)).size(); ++recHitIt)
          {
            if( (matchedRecHits_det->at(trackIt)).at(recHitIt) != 1 ) continue;
            float recHitE = (matchedRecHits_energy->at(trackIt)).at(recHitIt);
            if( recHitE*LOScales[RUId] < 1. ) continue;
            
            Npes.push_back( recHitE * LO * LOCorr * LOScale * (1-occupancy) * f_PDE->Eval(Vov)*(turnOn_PDE->Eval(v_fluence))/f_PDE_default->Eval(3.5) );
            recHitEs.push_back( recHitE );
            
            totEnergy += recHitE;
          }
          if( totEnergy < 1. ) continue;
          //std::cout << ">>>>>> trackIt " << trackIt << " / " << tracks_pt->size() << "   feta: " << feta << "   RUId: " << RUId << "   totEnergy: " << totEnergy << std::endl;
          
          float sigma_weighted = 0.;
          int jj = 0;
          for(auto nPE : Npes)
          {
            float sigma_stoch = 27. * sqrt(7000./(nPE*38.5/taud));
            float sigma_noise = sqrt( pow(noiseTerm/1.2/g_slewRate_vs_amp[RUId]->Eval(nPE*gainScale*f_gain->Eval(Vov)*(turnOn_gain->Eval(v_fluence))/f_gain->Eval(3.5)/9500.),2) + pow(16.7,2) )/sqrt(2);
            float sigma_DCR   = 40. * 6000./(nPE*38.5/taud) * pow(DCR*DCRCorr/30.,0.41);
            float sigma_clock = 15.;
            float sigma_tot = sqrt( pow(sigma_stoch,2) + pow(sigma_noise,2) + pow(sigma_DCR,2) + pow(sigma_clock,2) );
            sigma_weighted += 1. / pow(sigma_tot,2.);
            
            //std::cout << ">>>>>>>>> recHit E: " << recHitEs.at(jj) << "   Npe: " << nPE << "   sigma_tot: " << sigma_tot << std::endl;
            ++jj;
          }
          sigma_weighted = sqrt( 1./sigma_weighted );
          p_tRes_vs_RU -> Fill(RUId,sigma_weighted);
          //std::cout << ">>>>>> RUId" << RUId << "   sigma: " << sigma_weighted << std::endl;
          
        } // loop over tracks
        
      } // loop over events
      
      for(int RUIt = 1; RUIt <= 6; ++RUIt)
      {
        if( p_tRes_vs_RU->GetBinContent(RUIt) < tResBest_vs_RU[RUIt] )
        {
          tResBest_vs_RU[RUIt] = p_tRes_vs_RU->GetBinContent(RUIt);
          condBest_vs_RU[RUIt].first = effTECsDeltaT;
          condBest_vs_RU[RUIt].second = Vov;
        }
      }
      delete p_tRes_vs_RU;
      
    } // loop over Vov
    std::cout << std::endl;
    
  } // loop over TECs deltaT

  for(int RUIt = 1; RUIt <= 6; ++RUIt)
    std::cout << "RU" << RUIt << "   TECsDeltaT: " << condBest_vs_RU[RUIt].first << "   Vov: " << condBest_vs_RU[RUIt].second << "   tResBest: " << tResBest_vs_RU[RUIt] << std::endl;
  
  

  
  //------------
  // final plots
  outFile -> cd();

  TProfile* p_tResBest_vs_eta = new TProfile("p_tResBest_vs_eta","",48,0.,1.52);
  TProfile* p_tResBest_stoch_vs_eta = new TProfile("p_tResBest_stoch_vs_eta","",48,0.,1.52);
  TProfile* p_tResBest_noise_vs_eta = new TProfile("p_tResBest_noise_vs_eta","",48,0.,1.52);
  TProfile* p_tResBest_DCR_vs_eta = new TProfile("p_tResBest_DCR_vs_eta","",48,0.,1.52);

  TProfile* p_VbiasBest_vs_RU = new TProfile("p_VbiasBest_vs_RU","",6,0.5,6.5);
  TProfile* p_VovBest_vs_RU = new TProfile("p_VovBest_vs_RU","",6,0.5,6.5);
  TProfile* p_DCRBest_vs_RU = new TProfile("p_DCRBest_vs_RU","",6,0.5,6.5);
  TProfile* p_PDEBest_vs_RU = new TProfile("p_PDEBest_vs_RU","",6,0.5,6.5);
  TProfile* p_nPEBest_vs_RU = new TProfile("p_nPEBest_vs_RU","",6,0.5,6.5);
  TProfile* p_occupancyBest_vs_RU = new TProfile("p_occupancyBest_vs_RU","",6,0.5,6.5);
  TProfile* p_gainBest_vs_RU = new TProfile("p_gainBest_vs_RU","",6,0.5,6.5);
  TProfile* p_powerBest_vs_RU = new TProfile("p_powerBest_vs_RU","",6,0.5,6.5);
  TProfile* p_dynamicPowerBest_vs_RU = new TProfile("p_dynamicPowerBest_vs_RU","",6,0.5,6.5);
  TProfile* p_staticPowerBest_vs_RU = new TProfile("p_staticPowerBest_vs_RU","",6,0.5,6.5);
  TProfile* p_TECsPowerBest_vs_RU = new TProfile("p_TECsPowerBest_vs_RU","",6,0.5,6.5);
  TProfile* p_currentBest_vs_RU = new TProfile("p_currentBest_vs_RU","",6,0.5,6.5);
  TProfile* p_staticCurrentBest_vs_RU = new TProfile("p_staticCurrentBest_vs_RU","",6,0.5,6.5);
  TProfile* p_dynamicCurrentBest_vs_RU = new TProfile("p_dynamicCurrentBest_vs_RU","",6,0.5,6.5);
  TProfile* p_tempBest_vs_RU = new TProfile("p_tempBest_vs_RU","",6,0.5,6.5);  
  
  for(int RUIt = 1; RUIt <= 6; ++RUIt)
  {
    float effTECsDeltaT = condBest_vs_RU[RUIt].first;
    float Vov = condBest_vs_RU[RUIt].second;
    T_op = T_CO2 - effTECsDeltaT;
    
    
    //--- evaluate effective DCR for a given OV and T_op
    float DCR = ( 16300./8. * v_alpha * 1.E-17 * totFluence ) * // from HPK2E14 used at TB (after 3.5d at 110° C, alpha is 0.89, assume 8% gain loss (i.e. 29. GHz at Vov = 1.5 V, DCR(1.5)/DCR(0.8)~8)
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
    float Vbr = intercept + 0.001*slope*(T_op-0.);
    
    float staticCurrent = DCR*1E09 * f_ENF->Eval(Vov) * gainScale * f_gain->Eval(Vov)*(turnOn_gain->Eval(v_fluence)) * 1.602E-19;
    float staticPower = staticCurrent * (Vbr + Vov) * 1000.;   // in mW
    if( staticPower > maxPower ) continue;
    if( staticCurrent*1000. > 50./16. ) continue; // current in mA, 50 mA is the ALDO limit
    if( staticCurrent > 2./48./16. ) continue; // current in A, 2A is the PS limit, as per Krzysztof's plots
    
    float dynamicCurrent = (2.3*1E06 * v_instLumi * 5/7.5) * (4.2 * LO * LOScale * (1-occupancy) * f_PDE->Eval(Vov)*(turnOn_PDE->Eval(v_fluence))/f_PDE_default->Eval(3.5)) * f_ENF->Eval(Vov) * gainScale * f_gain->Eval(Vov)*(turnOn_gain->Eval(v_fluence)) * 1.602E-19;   // 2.3 MHz equivalent MIP rate at 200 PU.
    float dynamicPower = dynamicCurrent * (Vbr + Vov) * 1000.;   // in mW
    if( (staticPower+dynamicPower) > maxPower ) continue;
    if( (staticCurrent+dynamicCurrent)*1000. > 50./16. ) continue; // current in mA, 50 mA is the ALDO limit
    if( (staticCurrent+dynamicCurrent) > 2./48./16. ) continue; // current in A, 2A is the PS limit, as per Krzysztof's plots
    
    float TECsPower = 0.; // in mW per channel
    if( useTECs )
    {
      f_TECsPower_noSiPMLoad -> SetParameter(0.,(4./420.)*16.*(staticPower+dynamicPower));
      TECsPower = f_TECsPower_noSiPMLoad -> Eval(-1.*effTECsDeltaT) / 16.;
    }
    if( TECsPower > 125. ) continue; // limit to 2W/array as per Krzysztof's plots, with some margin
    
    if( (staticPower+dynamicPower+TECsPower) > maxPower ) continue;
    
    
    //--- evaluate time resolution
    for(int entry = 0; entry < 10000; ++entry)
    {
      //if( entry%1000 == 0 ) std::cout << ">>> reading entry " << entry << " / " << nEntries << "\r" << std::endl;
      tree -> GetEntry(entry);
      
      for(unsigned int trackIt = 0; trackIt < tracks_pt->size(); ++trackIt)
      {
        //std::cout << ">>>>>> trackIt " << trackIt << " / " << tracks_pt->size() << std::endl;
        
        float pt = tracks_pt->at(trackIt);
        float feta = fabs(tracks_eta->at(trackIt));
        float genPt = tracks_mcMatch_genPt->at(trackIt);
        float DR = tracks_mcMatch_DR->at(trackIt);
        
        if( DR > 0.01 ) continue;
        if( fabs(pt/genPt-1.) > 0.05 ) continue; 
        
        int RUId = -1.;
        float DCRCorr = 1.;
        float LOCorr = 1.;
        if( feta > 0.000 && feta < 0.347 ) { RUId = 1; };
        if( feta > 0.347 && feta < 0.661 ) { RUId = 2; };
        if( feta > 0.661 && feta < 0.927 ) { RUId = 3; };
        if( feta > 0.927 && feta < 1.148 ) { RUId = 4; };
        if( feta > 1.148 && feta < 1.338 ) { RUId = 5; };
        if( feta > 1.338 && feta < 1.498 ) { RUId = 6; };
        
        if( RUId != RUIt ) continue;
        
        if( RUTypes[RUId] == 1 ) { DCRCorr *= 1.25; LOCorr *= 1.05; }
        if( RUTypes[RUId] == 3 ) { DCRCorr *= 0.8;  LOCorr *= 0.95; }
        LOCorr *= LOScales[RUId];
        
        DCRCorr *= f_fluence_vs_eta -> Eval(feta) / 2.E14;
        
        std::vector<float> Npes;
        std::vector<float> recHitEs;
        float totEnergy = 0;
        for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(trackIt)).size(); ++recHitIt)
        {
          if( (matchedRecHits_det->at(trackIt)).at(recHitIt) != 1 ) continue;
          float recHitE = (matchedRecHits_energy->at(trackIt)).at(recHitIt);
          if( recHitE*LOScales[RUId] < 1. ) continue;
          
          Npes.push_back( recHitE * LO * LOCorr * LOScale * (1-occupancy) * f_PDE->Eval(Vov)*(turnOn_PDE->Eval(v_fluence))/f_PDE_default->Eval(3.5) );
          recHitEs.push_back( recHitE );
          
          totEnergy += recHitE;
        }
        if( totEnergy < 1. ) continue;
        //std::cout << ">>>>>> trackIt " << trackIt << " / " << tracks_pt->size() << "   feta: " << feta << "   RUId: " << RUId << "   totEnergy: " << totEnergy << std::endl;
        
        float sigma_weighted = 0.;
        float sigma_stoch_weighted = 0.;
        float sigma_noise_weighted = 0.;
        float sigma_DCR_weighted = 0.;
        float nPETot = 0;
        for(auto nPE : Npes)
        {
          nPETot += nPE;
          float sigma_stoch = 27. * sqrt(7000./(nPE*38.5/taud));
          float sigma_noise = sqrt( pow(noiseTerm/1.2/g_slewRate_vs_amp[RUId]->Eval(nPE*gainScale*f_gain->Eval(Vov)*(turnOn_gain->Eval(v_fluence))/f_gain->Eval(3.5)/9500.),2) + pow(16.7,2) )/sqrt(2);
          float sigma_DCR   = 40. * 6000./(nPE*38.5/taud) * pow(DCR*DCRCorr/30.,0.41);
          float sigma_clock = 15.;
          float sigma_tot = sqrt( pow(sigma_stoch,2) + pow(sigma_noise,2) + pow(sigma_DCR,2) + pow(sigma_clock,2) );
          sigma_weighted += 1. / pow(sigma_tot,2.);
          sigma_stoch_weighted += 1. / pow(sigma_stoch,2.);
          sigma_noise_weighted += 1. / pow(sigma_noise,2.);
          sigma_DCR_weighted += 1. / pow(sigma_DCR,2.);
        }
        sigma_weighted = sqrt( 1./sigma_weighted );
        sigma_stoch_weighted = sqrt( 1./sigma_stoch_weighted );
        sigma_noise_weighted = sqrt( 1./sigma_noise_weighted );
        sigma_DCR_weighted = sqrt( 1./sigma_DCR_weighted );
        
        p_tResBest_vs_eta -> Fill(feta,sigma_weighted);
        p_tResBest_stoch_vs_eta -> Fill(feta,sigma_stoch_weighted);
        p_tResBest_noise_vs_eta -> Fill(feta,sigma_noise_weighted);
        p_tResBest_DCR_vs_eta -> Fill(feta,sigma_DCR_weighted);
        
        p_VbiasBest_vs_RU -> Fill(RUId,Vbr+Vov);
        p_VovBest_vs_RU -> Fill(RUId,Vov);
        p_DCRBest_vs_RU -> Fill(RUId,DCR*DCRCorr);
        p_PDEBest_vs_RU -> Fill(RUId,f_PDE->Eval(Vov)*(turnOn_PDE->Eval(v_fluence)));
        p_nPEBest_vs_RU -> Fill(RUId,nPETot);
        p_occupancyBest_vs_RU -> Fill(RUId,occupancy);
        p_gainBest_vs_RU -> Fill(RUId,gainScale*f_gain->Eval(Vov)*(turnOn_gain->Eval(v_fluence)));
        p_powerBest_vs_RU -> Fill(RUId,staticPower+dynamicPower+TECsPower);
        p_staticPowerBest_vs_RU -> Fill(RUId,staticPower);
        p_dynamicPowerBest_vs_RU -> Fill(RUId,dynamicPower);
        p_TECsPowerBest_vs_RU -> Fill(RUId,TECsPower);
        p_currentBest_vs_RU -> Fill(RUId,(staticPower+dynamicPower)/(Vbr+Vov));
        p_staticCurrentBest_vs_RU -> Fill(RUId,(staticPower)/(Vbr+Vov));
        p_dynamicCurrentBest_vs_RU -> Fill(RUId,(dynamicPower)/(Vbr+Vov));
        p_tempBest_vs_RU -> Fill(RUId,T_op);
        
        //std::cout << ">>>>>> RUId" << RUId << "   sigma: " << sigma_weighted << std::endl;
        
      } // loop over tracks
      
    } // loop over events
    
  } // loop over RUs
  
  
  
  
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
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  return 0;
}
