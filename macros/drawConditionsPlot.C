
std::pair<std::vector<TGraph*>,float> get_tRes_graphs(TGraph* g_tRes, TGraph* g_temp)
{
  std::vector<TGraph*> vec;

  float tRes_ave = 0.;
  int tRes_n = 0;
  
  bool isTechStop = true;
  for(int point = 0; point < g_temp->GetN(); ++point)
  {
    double time, temp;
    g_temp -> GetPoint(point,time,temp);
    
    if( temp < -20. && isTechStop == true )
    {
      isTechStop = false;
      vec.push_back(new TGraph());
      
      TGraph* g = vec.at(vec.size()-1);
      float val = g_tRes->Eval(time);
      g -> SetPoint(g->GetN(),time,val);

      tRes_ave += val;
      ++tRes_n;
    }

    if( temp < -20. && isTechStop == false )
    {
      TGraph* g = vec.at(vec.size()-1);
      float val = g_tRes->Eval(time);
      g -> SetPoint(g->GetN(),time,val);

      tRes_ave += val;
      ++tRes_n;
    }

    if( temp > -20. && isTechStop == false )
    {
      isTechStop = true;
      continue;
      
    }
    
    if( isTechStop == true )
    {
      continue;
    }
    
  }

  std::pair<std::vector<TGraph*>,float> res(vec,tRes_ave/tRes_n);
  return res;
}




void makePlot(const std::vector<std::string>& inFileNames, const std::vector<int>& colors, const std::vector<std::string>& labels,
              const std::string& graphName, const std::string& label,
              const float& xMin, const float& yMin, const float& xMax, const float& yMax,
              const std::string& xAxisTitle, const std::string& yAxisTitle)
{
  TCanvas* c = new TCanvas(Form("c_%s",label.c_str()),Form("c_%s",label.c_str()),1400,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(xMin,yMin,xMax,yMax) );
  hPad -> SetTitle(Form(";%s;%s",xAxisTitle.c_str(),yAxisTitle.c_str()));
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.90,0.99-0.05*inFileNames.size(),0.99,0.99);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.03);
  
  int it = 0;
  for(auto inFileName : inFileNames)
  {
    TFile* inFile = TFile::Open(inFileName.c_str());
    
    TGraph* g = (TGraph*)( inFile->Get(graphName.c_str()) );
    g -> SetLineColor(colors.at(it));
    g -> SetLineWidth(3);
    g -> Draw("L,same");
    
    if( graphName == "g_powerBest_vs_time" )
    {
      TGraph* g1 = (TGraph*)( inFile->Get("g_dynamicPowerBest_vs_time") );
      g1 -> SetLineColor(colors.at(it)-4);
      g1 -> SetLineStyle(7);
      g1 -> SetLineWidth(3);
      g1 -> Draw("L,same");

      TGraph* g2 = (TGraph*)( inFile->Get("g_staticPowerBest_vs_time") );
      g2 -> SetLineColor(colors.at(it)-7);
      g2 -> SetLineStyle(10);
      g2 -> SetLineWidth(3);
      g2 -> Draw("L,same");
      
      TGraph* g3 = (TGraph*)( inFile->Get("g_TECsPowerBest_vs_time") );
      g3 -> SetLineColor(kGreen+1);
      g3 -> SetLineStyle(1);
      g3 -> SetLineWidth(3);
      g3 -> Draw("L,same");
    }

    if( graphName == "g_powerBest_vs_intLumi" )
    {
      TGraph* g1 = (TGraph*)( inFile->Get("g_dynamicPowerBest_vs_intLumi") );
      g1 -> SetLineColor(colors.at(it)-4);
      g1 -> SetLineStyle(7);
      g1 -> SetLineWidth(3);
      g1 -> Draw("L,same");

      TGraph* g2 = (TGraph*)( inFile->Get("g_staticPowerBest_vs_intLumi") );
      g2 -> SetLineColor(colors.at(it)-7);
      g2 -> SetLineStyle(10);
      g2 -> SetLineWidth(3);
      g2 -> Draw("L,same");
      
      TGraph* g3 = (TGraph*)( inFile->Get("g_TECsPowerBest_vs_intLumi") );
      g3 -> SetLineColor(kGreen+1);
      g3 -> SetLineStyle(1);
      g3 -> SetLineWidth(3);
      g3 -> Draw("L,same");      
    }
    
    
    if( graphName == "g_currentBest_vs_time" )
    {
      TGraph* g1 = (TGraph*)( inFile->Get("g_dynamicCurrentBest_vs_time") );
      g1 -> SetLineColor(colors.at(it)-4);
      g1 -> SetLineStyle(7);
      g1 -> SetLineWidth(3);
      g1 -> Draw("L,same");

      TGraph* g2 = (TGraph*)( inFile->Get("g_staticCurrentBest_vs_time") );
      g2 -> SetLineColor(colors.at(it)-7);
      g2 -> SetLineStyle(10);
      g2 -> SetLineWidth(3);
      g2 -> Draw("L,same");
    }

    if( graphName == "g_currentBest_vs_intLumi" )
    {
      TGraph* g1 = (TGraph*)( inFile->Get("g_dynamicCurrentBest_vs_intLumi") );
      g1 -> SetLineColor(colors.at(it)-4);
      g1 -> SetLineStyle(7);
      g1 -> SetLineWidth(3);
      g1 -> Draw("L,same");

      TGraph* g2 = (TGraph*)( inFile->Get("g_staticCurrentBest_vs_intLumi") );
      g2 -> SetLineColor(colors.at(it)-7);
      g2 -> SetLineStyle(10);
      g2 -> SetLineWidth(3);
      g2 -> Draw("L,same");
    }
    
    legend -> AddEntry(g,(labels.at(it)).c_str(),"L");
    ++it;
  }
  
  legend -> Draw("same");

  c -> Print(Form("HLLHC_HPK_%s.png",label.c_str()));
}



void drawConditionsPlot(const int& TECs)
{
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleOffset(0.8, "Y");
  
  std::vector<std::string> inFileNames;
  std::vector<int> colors;
  std::vector<std::string> labels;

  std::string TECsLabel;
  
  if( TECs == 1 )
  {
    inFileNames.push_back("../plots/outFile__LO1600__HPK_dropPDE0.15_dropGain0.30_TECs__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower80.root");
    colors.push_back(kRed);
    labels.push_back("HPK");
    
    //inFileNames.push_back("plots/outFile_FBK_Top-45_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower65.root");
    //colors.push_back(kBlue);
    //labels.push_back("FBK");

    TECsLabel = "withTECs_";
  }
  else
  {
    inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-35_Tann15_HLLHCSchedule8_maxPower65.root");
    colors.push_back(kRed);
    labels.push_back("HPK");
    
    inFileNames.push_back("plots/outFile_FBK_Top-35_Tann15_HLLHCSchedule8_maxPower65.root");
    colors.push_back(kBlue);
    labels.push_back("FBK");    
    
    TECsLabel = "noTECs_";
  }
  

  
  makePlot(inFileNames, colors, labels,
           "g_nPEBest_vs_time", TECsLabel+"nPE_vs_time",
           0., 0., 13., 12000,
           "years from 2027", "N_{p.e.}");
  makePlot(inFileNames, colors, labels,
           "g_nPEBest_vs_intLumi", TECsLabel+"nPE_vs_intLumi",
           0., 0., 3000., 12000,
           "integrated luminosity [fb^{-1}]", "N_{p.e.}");
  
  makePlot(inFileNames, colors, labels,
           "g_PDEBest_vs_time", TECsLabel+"PDE_vs_time",
           0., 0., 13., 0.5,
           "years from 2027", "PDE");
  makePlot(inFileNames, colors, labels,
           "g_PDEBest_vs_intLumi", TECsLabel+"PDE_vs_intLumi",
           0., 0., 3000., 0.5,
           "integrated luminosity [fb^{-1}]", "PDE");
  
  makePlot(inFileNames, colors, labels,
           "g_VovBest_vs_time", TECsLabel+"Vov_vs_time",
           0., 0., 13., 6.,
           "years from 2027", "V_{OV} [V]");
  makePlot(inFileNames, colors, labels,
           "g_VovBest_vs_intLumi", TECsLabel+"Vov_vs_intLumi",
           0., 0., 3000., 6.,
           "integrated luminosity [fb^{-1}]", "V_{OV} [V]");
  
  makePlot(inFileNames, colors, labels,
           "g_VbiasBest_vs_time", TECsLabel+"Vbias_vs_time",
           0., 30., 13., 45.,
           "years from 2027", "V_{bias} [V]");
  makePlot(inFileNames, colors, labels,
           "g_VbiasBest_vs_intLumi", TECsLabel+"Vbias_vs_intLumi",
           0., 30., 3000., 45.,
           "integrated luminosity [fb^{-1}]", "V_{bias} [V]");
  
  makePlot(inFileNames, colors, labels,
           "g_DCRBest_vs_time", TECsLabel+"DCR_vs_time",
           0., 0., 13., 100.,
           "years from 2027", "DCR [GHz]");
  makePlot(inFileNames, colors, labels,
           "g_DCRBest_vs_intLumi", TECsLabel+"DCR_vs_intLumi",
           0., 0., 3000., 100.,
           "integrated luminosity [fb^{-1}]", "DCR [GHz]");
  
  makePlot(inFileNames, colors, labels,
           "g_gainBest_vs_time", TECsLabel+"gain_vs_time",
           0., 0., 13., 8E05,
           "years from 2027", "gain");
  makePlot(inFileNames, colors, labels,
           "g_gainBest_vs_intLumi", TECsLabel+"gain_vs_intLumi",
           0., 0., 3000., 8E05,
           "integrated luminosity [fb^{-1}]", "gain");
  
  makePlot(inFileNames, colors, labels,
           "g_powerBest_vs_time", TECsLabel+"power_vs_time",
           0., 0., 13., 80,
           "years from 2027", "total power / ch. [mW]");
  makePlot(inFileNames, colors, labels,
           "g_powerBest_vs_intLumi", TECsLabel+"power_vs_intLumi",
           0., 0., 3000., 80,
           "integrated luminosity [fb^{-1}]", "total power / ch. [mW]");
  
  makePlot(inFileNames, colors, labels,
           "g_currentBest_vs_time", TECsLabel+"current_vs_time",
           0., 0., 13., 2.5,
           "years from 2027", "total current / ch. [mA]");
  makePlot(inFileNames, colors, labels,
           "g_currentBest_vs_intLumi", TECsLabel+"current_vs_intLumi",
           0., 0., 3000., 2.5
           ,
           "integrated luminosity [fb^{-1}]", "total current / ch. [mA]");
}

