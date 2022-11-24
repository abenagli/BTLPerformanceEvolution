
std::pair<std::vector<TGraphAsymmErrors*>,float> get_tRes_graphs(TGraph* g_tRes, TGraph* g_temp)
{
  std::vector<TGraphAsymmErrors*> vec;

  float tRes_ave = 0.;
  int tRes_n = 0;

  float tRes_TDR_ave = 0.;
  int tRes_TDR_n = 0;
  
  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  //f_TDR -> Draw("same");
  
  bool isTechStop = true;
  for(int point = 0; point < g_temp->GetN(); ++point)
  {
    double time, temp;
    g_temp -> GetPoint(point,time,temp);
    
    if( temp < -20. && isTechStop == true )
    {
      isTechStop = false;
      vec.push_back(new TGraphAsymmErrors());
      
      TGraph* g = vec.at(vec.size()-1);
      float val = g_tRes->Eval(time);
      g -> SetPoint(g->GetN(),time,val);

      tRes_ave += val;
      ++tRes_n;
      
      tRes_TDR_ave += f_TDR->Eval(time) ;
      ++tRes_TDR_n;
    }
    
    if( temp < -20. && isTechStop == false )
    {
      TGraph* g = vec.at(vec.size()-1);
      float val = g_tRes->Eval(time);
      g -> SetPoint(g->GetN(),time,val);

      tRes_ave += val;
      ++tRes_n;

      tRes_TDR_ave += f_TDR->Eval(time) ;
      ++tRes_TDR_n;
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

  std::cout << "TDR: " << 1.*tRes_TDR_ave/tRes_TDR_n << std::endl;
  
  std::pair<std::vector<TGraphAsymmErrors*>,float> res(vec,tRes_ave/tRes_n);
  return res;
}



void drawFinalPlot()
{
  gStyle -> SetPadLeftMargin(0.15);
  gStyle -> SetPadRightMargin(0.05);
  gStyle -> SetTitleOffset(1.20, "Y");
  
  std::vector<std::pair<int,int> >vals;
  std::vector<int> TECFlags;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  vals.push_back(std::make_pair<int,int>(-35,15));  TECFlags.push_back(0); colors.push_back(kBlack);   labels.push_back("T_{op} = -35#circ C   T_{ann} = 15#circ C   power budget: 65 mW (no TECs)");
  vals.push_back(std::make_pair<int,int>(-35,40));  TECFlags.push_back(0); colors.push_back(kGray);    labels.push_back("T_{op} = -35#circ C   T_{ann} = 40#circ C   power budget: 65 mW (no TECs)");
  vals.push_back(std::make_pair<int,int>(-45,15));  TECFlags.push_back(0); colors.push_back(99);       labels.push_back("T_{op} = -45#circ C   T_{ann} = 15#circ C   power budget: 65 mW (no TECs)");
  vals.push_back(std::make_pair<int,int>(-45,40));  TECFlags.push_back(0); colors.push_back(95);       labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   power budget: 65 mW (no TECs)");
  vals.push_back(std::make_pair<int,int>(-45,70));  TECFlags.push_back(0); colors.push_back(90);       labels.push_back("T_{op} = -45#circ C   T_{ann} = 70#circ C   power budget: 65 mW (no TECs)");
  vals.push_back(std::make_pair<int,int>(-43,15));  TECFlags.push_back(1); colors.push_back(70);       labels.push_back("T_{op} = -43#circ C   T_{ann} = 15#circ C   power budget: 65 mW (w/ TECs)");
  vals.push_back(std::make_pair<int,int>(-43,40));  TECFlags.push_back(1); colors.push_back(65);       labels.push_back("T_{op} = -43#circ C   T_{ann} = 40#circ C   power budget: 65 mW (w/ TECs)");
  vals.push_back(std::make_pair<int,int>(-43,70));  TECFlags.push_back(1); colors.push_back(60);       labels.push_back("T_{op} = -43#circ C   T_{ann} = 70#circ C   power budget: 65 mW (w/ TECs)");
  
  
  TCanvas* c = new TCanvas("c","",1200,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,20.,10.,140.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.18,0.64,0.80,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  int it = 0;
  for(auto val : vals)
  {
    TFile* inFile;
    if( TECFlags.at(it) == 0 ) inFile = TFile::Open(Form("plots/outFile_Top%d_Tann%d_HLLHCSchedule3_PDELoss15perc_maxPower65.root",val.first,val.second));
    else                       inFile = TFile::Open(Form("plots/outFile_Top%d_Tann%d_TECs_HLLHCSchedule3_PDELoss15perc_maxPower65.root",val.first,val.second));
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(3);
      graph -> Draw("L,same");

      /* if( TECFlags.at(it) == 0 && val.first < -40. ) graph -> SetLineStyle(3); */
    }
    
    legend -> AddEntry(graphs[0],(labels.at(it)+Form("   #LT #sigma_{t} #GT = %.0f ps",res.second)).c_str(),"L");
    ++it;
  }

  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  f_TDR -> Draw("same");

  legend -> AddEntry(f_TDR,"TDR","L");
  
  legend -> Draw("same");

  c -> Print("HLLHC_tRes.png");
}



void drawFinalPlot_TECs_vsPower()
{
  std::vector<int>vals;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  vals.push_back(30);  colors.push_back(83);   labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   power budget:  30 mW (w/ TECs)");
  vals.push_back(50);  colors.push_back(87);   labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   power budget:  50 mW (w/ TECs)");
  vals.push_back(65);  colors.push_back(91);   labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   power budget:  65 mW (w/ TECs)");
  vals.push_back(80);  colors.push_back(95);   labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   power budget:  80 mW (w/ TECs)");
  vals.push_back(100); colors.push_back(kRed); labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   power budget: 100 mW (w/ TECs)");
  
  
  TCanvas* c = new TCanvas("c","",1200,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,20.,10.,140.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.18,0.64,0.80,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  int it = 0;
  for(auto val : vals)
  {
    TFile* inFile = TFile::Open(Form("plots/outFile_Top-45_Tann40_TECs_HLLHCSchedule3_PDELoss15perc_maxPower%d.root",val));
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(3);
      graph -> Draw("L,same");
      
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(3);
      graph -> Draw("L,same");
    }
    
    legend -> AddEntry(graphs[0],(labels.at(it)+Form("   #LT #sigma_{t} #GT = %.0f ps",res.second)).c_str(),"L");
    ++it;
  }

  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  f_TDR -> Draw("same");

  legend -> AddEntry(f_TDR,"TDR","L");
  
  legend -> Draw("same");

  c -> Print("HLLHC_tRes_TECs_vsPower.png");
}



void drawFinalPlot_noTECs_vsPower()
{
  std::vector<int>vals;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  vals.push_back(30);  colors.push_back(83);   labels.push_back("T_{op} = -35#circ C   T_{ann} = 40#circ C   power budget:  30 mW (no TECs)");
  vals.push_back(50);  colors.push_back(87);   labels.push_back("T_{op} = -35#circ C   T_{ann} = 40#circ C   power budget:  50 mW (no TECs)");
  vals.push_back(65);  colors.push_back(91);   labels.push_back("T_{op} = -35#circ C   T_{ann} = 40#circ C   power budget:  65 mW (no TECs)");
  vals.push_back(80);  colors.push_back(95);   labels.push_back("T_{op} = -35#circ C   T_{ann} = 40#circ C   power budget:  80 mW (no TECs)");
  vals.push_back(100); colors.push_back(kRed); labels.push_back("T_{op} = -35#circ C   T_{ann} = 40#circ C   power budget: 100 mW (no TECs)");
  
  
  TCanvas* c = new TCanvas("c","",1200,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,20.,10.,140.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.18,0.64,0.80,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  int it = 0;
  for(auto val : vals)
  {
    TFile* inFile = TFile::Open(Form("plots/outFile_Top-35_Tann40_HLLHCSchedule3_PDELoss15perc_maxPower%d.root",val));
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(3);
      graph -> Draw("L,same");
      
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(3);
      graph -> Draw("L,same");
    }
    
    legend -> AddEntry(graphs[0],(labels.at(it)+Form("   #LT #sigma_{t} #GT = %.0f ps",res.second)).c_str(),"L");
    ++it;
  }

  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  f_TDR -> Draw("same");

  legend -> AddEntry(f_TDR,"TDR","L");
  
  legend -> Draw("same");

  c -> Print("HLLHC_tRes_noTECs_vsPower.png");
}



void drawFinalPlot_vsTECsDeltaT()
{
  std::vector<int>vals;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  vals.push_back(2);  colors.push_back(kRed); labels.push_back("T_{op} = -37#circ C   T_{ann} = 40#circ C   power budget: 65 mW (w/ TECs)");
  vals.push_back(5);  colors.push_back(92);   labels.push_back("T_{op} = -40#circ C   T_{ann} = 40#circ C   power budget: 65 mW (w/ TECs)");
  vals.push_back(8);  colors.push_back(84);   labels.push_back("T_{op} = -43#circ C   T_{ann} = 40#circ C   power budget: 65 mW (w/ TECs)");
  vals.push_back(10); colors.push_back(76);   labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   power budget: 65 mW (w/ TECs)");
  vals.push_back(12); colors.push_back(68);   labels.push_back("T_{op} = -47#circ C   T_{ann} = 40#circ C   power budget: 65 mW (w/ TECs)");
  vals.push_back(15); colors.push_back(60);   labels.push_back("T_{op} = -50#circ C   T_{ann} = 40#circ C   power budget: 65 mW (w/ TECs)");
  
  
  TCanvas* c = new TCanvas("c","",1200,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,20.,10.,140.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.18,0.64,0.80,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  int it = 0;
  for(auto val : vals)
  {
    TFile* inFile = TFile::Open(Form("plots/outFile_Top%d_Tann40_TECs_HLLHCSchedule3_PDELoss15perc_maxPower65.root",-35-val));
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(3);
      graph -> Draw("L,same");
      
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(3);
      graph -> Draw("L,same");
    }
    
    legend -> AddEntry(graphs[0],(labels.at(it)+Form("   #LT #sigma_{t} #GT = %.0f ps",res.second)).c_str(),"L");
    ++it;
  }

  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  f_TDR -> Draw("same");

  legend -> AddEntry(f_TDR,"TDR","L");
  
  legend -> Draw("same");

  c -> Print("HLLHC_tRes_vsTECsDeltaT.png");
}






void drawFinalPlot_vsInfraYearAnnealing()
{
  std::vector<int>vals;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  vals.push_back(3);  colors.push_back(kBlack); labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C");
  vals.push_back(4);  colors.push_back(95);     labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   2 weeks / year at  0#circ C");
  vals.push_back(5);  colors.push_back(85);     labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   2 weeks / year at 15#circ C");
  vals.push_back(6);  colors.push_back(75);     labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   2 weeks / year at 40#circ C");
  vals.push_back(7);  colors.push_back(65);     labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   2 weeks / year at 70#circ C");
  vals.push_back(41); colors.push_back(95);     labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   4 weeks / year at  0#circ C");
  vals.push_back(51); colors.push_back(85);     labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   4 weeks / year at 15#circ C");
  vals.push_back(61); colors.push_back(75);     labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   4 weeks / year at 40#circ C");
  vals.push_back(71); colors.push_back(65);     labels.push_back("T_{op} = -45#circ C   T_{ann} = 40#circ C   4 weeks / year at 70#circ C");

  
  
  TCanvas* c = new TCanvas("c","",1200,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,20.,10.,140.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.18,0.64,0.80,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  int it = 0;
  for(auto val : vals)
  {
    TFile* inFile = TFile::Open(Form("plots/outFile_Top-45_Tann40_TECs_HLLHCSchedule%d_PDELoss15perc_maxPower65.root",val));
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      if( val > 10 ) graph -> SetLineStyle(7);
      graph -> SetLineWidth(3);
      graph -> Draw("L,same");
    }
    
    legend -> AddEntry(graphs[0],(labels.at(it)+Form("   #LT #sigma_{t} #GT = %.0f ps",res.second)).c_str(),"L");
    ++it;
  }

  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  f_TDR -> Draw("same");

  legend -> AddEntry(f_TDR,"TDR","L");
  
  legend -> Draw("same");

  c -> Print("HLLHC_tRes_vsInfraYearAnnealing.png");
}






void drawFinalPlot_withManyInfraYearAnnealings()
{
  std::vector<std::pair<int,int> > vals;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  vals.push_back(std::make_pair<int,int>(0,0));  colors.push_back(kBlack);    labels.push_back("T^{op.} = -35#circ C   T^{ann.} = 15#circ C   #rightarrow  ");
  vals.push_back(std::make_pair<int,int>(1,0));  colors.push_back(kRed);      labels.push_back("T^{op.} = -45#circ C   T^{ann.} = 40#circ C   4 days/6 weeks at 0#circ C   #rightarrow  ");
  vals.push_back(std::make_pair<int,int>(1,1));  colors.push_back(kBlue);     labels.push_back("T^{op.} = -45#circ C   T^{ann.} = 40#circ C   4 days/6 weeks at 0#circ C,  inter-fills at 0#circ C  #rightarrow  ");
  vals.push_back(std::make_pair<int,int>(1,2));  colors.push_back(kGreen+1);  labels.push_back("T^{op.} = -45#circ C   T^{ann.} = 40#circ C   4 days/6 weeks at 40#circ C, inter-fills at 0#circ C, 80 mW/ch, PDE 7.5% #rightarrow  ");
  
  
  
  TCanvas* c = new TCanvas("c","",1200,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,20.,10.,140.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.18,0.64,0.85,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.02);
  
  int it = 0;
  for(auto val : vals)
  {
    std::string TECLabel = "";
    if( val.first == 1 ) TECLabel = "TECs_";
    
    TFile* inFile;
    if( val.first == 0) inFile = TFile::Open(Form("plots/outFile_Top-35_Tann15_HLLHCSchedule8_PDELoss15perc_maxPower65.root"));
    else
    {
      if(val.second == 0 ) inFile = TFile::Open(Form("plots/outFile_Top-45_Tann40_%sHLLHCSchedule8_PDELoss15perc_maxPower65.root",TECLabel.c_str()));
      if(val.second == 1 ) inFile = TFile::Open(Form("plots/outFile_Top-45_Tann40_%sinterfillAnnealing_HLLHCSchedule8_PDELoss15perc_maxPower65.root",TECLabel.c_str()));
      if(val.second == 2 ) inFile = TFile::Open(Form("plots/outFile_Top-45_Tann40_%sinterfillAnnealing_HLLHCSchedule8_PDELoss7p5perc_maxPower80.root",TECLabel.c_str()));
    }
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(3);
      graph -> Draw("L,same");
    }
    
    legend -> AddEntry(graphs[0],(labels.at(it)+Form("#LT #sigma_{t} #GT = %.0f ps",res.second)).c_str(),"L");
    ++it;
  }

  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  f_TDR -> Draw("same");

  legend -> AddEntry(f_TDR,"TDR","L");
  
  legend -> Draw("same");

  c -> Print("HLLHC_tRes_optimisticTECScenario.png");
}



void drawFinalPlot_HPK_vs_FBK()
{
  std::vector<std::pair<int,int> > vals;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  vals.push_back(std::make_pair<int,int>(1,1));  colors.push_back(kRed);  labels.push_back("T^{op.} = -45#circ C   T^{ann.} = 40#circ C   4 days/6 weeks at 40#circ C, inter-fills at 0#circ C, 60 mW/ch, HPK #rightarrow  ");
  vals.push_back(std::make_pair<int,int>(1,2));  colors.push_back(kBlue); labels.push_back("T^{op.} = -45#circ C   T^{ann.} = 40#circ C   4 days/6 weeks at 40#circ C, inter-fills at 0#circ C, 60 mW/ch, FBK #rightarrow  ");
  
  
  
  TCanvas* c = new TCanvas("c","",1200,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,20.,10.,140.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.18,0.64,0.85,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.02);
  
  int it = 0;
  for(auto val : vals)
  {
    std::string TECLabel = "";
    if( val.first == 1 ) TECLabel = "TECs_";
    
    TFile* inFile;
    if( val.first == 0 ) inFile = TFile::Open(Form("plots/outFile_Top-35_Tann15_HLLHCSchedule8_PDELoss15perc_maxPower65.root"));
    else
    {
      if(val.second == 1 ) inFile = TFile::Open(Form("plots/outFile_HPK_Top-45_Tann40_%sinterfillAnnealing_HLLHCSchedule8_PDELoss7p5perc_maxPower80.root",TECLabel.c_str()));
      if(val.second == 2 ) inFile = TFile::Open(Form("plots/outFile_FBK_Top-45_Tann40_%sinterfillAnnealing_HLLHCSchedule8_PDELoss7p5perc_maxPower80.root",TECLabel.c_str()));
    }
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(3);
      graph -> Draw("L,same");
    }
    
    legend -> AddEntry(graphs[0],(labels.at(it)+Form("#LT #sigma_{t} #GT = %.0f ps",res.second)).c_str(),"L");
    ++it;
  }

  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  f_TDR -> Draw("same");

  legend -> AddEntry(f_TDR,"TDR","L");
  
  legend -> Draw("same");

  c -> Print("HLLHC_tRes_HPK_FBK.png");
}






void drawFinalPlot_final_HPK()
{
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleOffset(0.8, "Y");
  
  std::vector<std::string> inFileNames;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-35_Tann15_HLLHCSchedule8_maxPower65.root");
  colors.push_back(kBlack);
  labels.push_back("no TECs: T^{op.} = -35#circ C T^{ann.} = 15#circ C,   65 mW/ch   #rightarrow  ");

  //inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-35_Tann15_HLLHCSchedule81_maxPower65.root");
  //colors.push_back(kGray+2);
  //labels.push_back("no TECs: T^{op.} = -35#circ C T^{ann.} = 15#circ C,   4 days/6 weeks at  0#circ C,   65 mW/ch   #rightarrow  ");
  
  //inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-35_Tann15_HLLHCSchedule81_maxPower80.root");
  //colors.push_back(kGray+1);
  //labels.push_back("no TECs: T^{op.} = -35#circ C T^{ann.} = 15#circ C,   4 days/6 weeks at  0#circ C,   80 mW/ch   #rightarrow  ");
  
  //inFileNames.push_back("plots/outFile_HPK_PDELoss07p5_Top-35_Tann15_HLLHCSchedule81_maxPower80.root");
  //colors.push_back(kGray);
  //labels.push_back("no TECs: T^{op.} = -35#circ C T^{ann.} = 15#circ C,   4 days/6 weeks at  0#circ C,   80 mW/ch,   PDE loss -50%   #rightarrow  ");
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-45_Tann40_TECs_HLLHCSchedule8_maxPower65.root");
  colors.push_back(kRed);
  labels.push_back("w/ TECs: T^{op.} = -45#circ C T^{ann.} = 40#circ C,   4 days/6 weeks at  0#circ C,   65 mW/ch   #rightarrow  ");
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-45_Tann40_TECs_HLLHCSchedule81_maxPower65.root");
  colors.push_back(kOrange);
  labels.push_back("w/ TECs: T^{op.} = -45#circ C T^{ann.} = 40#circ C,   4 days/6 weeks at 40#circ C,   65 mW/ch   #rightarrow  ");
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-45_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower65.root");
  colors.push_back(kAzure+10);
  labels.push_back("w/ TECs: T^{op.} = -45#circ C T^{ann.} = 40#circ C,   4 days/6 weeks at 40#circ C,   interfills at 0#circ C,   65 mW/ch   #rightarrow  ");
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-45_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower80.root");
  colors.push_back(kBlue);
  labels.push_back("w/ TECs: T^{op.} = -45#circ C T^{ann.} = 40#circ C,   4 days/6 weeks at 40#circ C,   interfills at 0#circ C,   80 mW/ch   #rightarrow  ");
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss07p5_Top-45_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower80.root");
  colors.push_back(kGreen+1);
  labels.push_back("w/ TECs: T^{op.} = -45#circ C T^{ann.} = 40#circ C,   4 days/6 weeks at 40#circ C,   interfills at 0#circ C,   80 mW/ch,   PDE loss -50%   #rightarrow  ");
  
  
  
  TCanvas* c = new TCanvas("c","",1400,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,20.,10.,130.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.15,0.94-0.04*inFileNames.size(),0.50,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.02);
  
  int it = 0;
  for(auto inFileName : inFileNames)
  {
    TFile* inFile = TFile::Open(inFileName.c_str());
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(3);
      graph -> Draw("L,same");
    }
    
    legend -> AddEntry(graphs[0],(labels.at(it)+Form("#LT #sigma_{t} #GT = %.0f ps",res.second)).c_str(),"L");
    ++it;
  }

  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  f_TDR -> Draw("same");

  legend -> AddEntry(f_TDR,"TDR","L");
  
  legend -> Draw("same");

  c -> Print("HLLHC_final_HPK.png");
}






void drawFinalPlot_final_HPK_tommaso()
{
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleOffset(0.8, "Y");
  
  std::vector<std::string> inFileNames;
  std::vector<std::string> inFileNames_up;
  std::vector<std::string> inFileNames_down;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-35_Tann15_HLLHCSchedule8_maxPower65.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kRed);
  labels.push_back("Post-TDR w/o TECs");
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-45_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower80.root");
  inFileNames_up.push_back("plots/outFile_HPK_PDELoss15p0_Top-45_Tann40_TECs_HLLHCSchedule81_maxPower65.root");
  inFileNames_down.push_back("plots/outFile_HPK_PDELoss07p5_Top-45_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower80.root");
  colors.push_back(kBlue);
  labels.push_back("Post-TDR w/ TECs [for a range of operation scenarios]");

  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-45_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower80_reducedLO.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kRed+2);
  labels.push_back("Post-TDR w/ TECs and -30% LO");
  
  
  
  TCanvas* c = new TCanvas("c","",1200,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,20.,10.,120.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.15,0.94-0.10*inFileNames.size(),0.85,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.04);
  
  int it = 0;
  for(auto inFileName : inFileNames)
  {
    std::cout << labels.at(it) << std::endl;
    
    TFile* inFile = TFile::Open(inFileName.c_str());
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    
    TGraph* g_up;
    TGraph* g_down;
    if( inFileNames_up.at(it) != "NULL" )
    {
      TFile* inFile_up = TFile::Open(inFileNames_up.at(it).c_str());
      g_up = (TGraph*)( inFile_up->Get("g_tResBest_vs_time") );
    }
    if( inFileNames_down.at(it) != "NULL" )
    {
      TFile* inFile_down = TFile::Open(inFileNames_down.at(it).c_str());
      g_down = (TGraph*)( inFile_down->Get("g_tResBest_vs_time") );
    }
    
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(4);
      
      if( inFileNames_up.at(it) != "NULL" )
      {
        for(int point = 0; point < graph->GetN(); ++point)
        {
          double x = graph -> GetPointX(point);
          double y = graph -> GetPointY(point);
          double y_up = g_up -> Eval(x);
          graph -> SetPointEYhigh(point,y_up-y);
        }
      }
      if( inFileNames_down.at(it) != "NULL" )
      {
        for(int point = 0; point < graph->GetN(); ++point)
        {
          double x = graph -> GetPointX(point);
          double y = graph -> GetPointY(point);
          double y_down = g_down -> Eval(x);
          graph -> SetPointEYlow(point,y-y_down);
        }
      }

      graph -> SetFillColor(colors.at(it));
      graph -> SetFillStyle(3001);
      graph -> Draw("L3,same");
    }
    
    legend -> AddEntry(graphs[0],(labels.at(it)).c_str(),"L");
    ++it;
  }

  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  f_TDR -> Draw("same");

  legend -> AddEntry(f_TDR,"TDR goal","L");
  
  legend -> Draw("same");

  c -> Print("HLLHC_final_HPK_tommaso.png");
  c -> Print("HLLHC_final_HPK_tommaso.C");
}



void drawFinalPlot_final_HPK_thermalInterfaces()
{
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleOffset(0.8, "Y");
  
  std::vector<std::string> inFileNames;
  std::vector<std::string> inFileNames_up;
  std::vector<std::string> inFileNames_down;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-35_Tann15_HLLHCSchedule8_maxPower65.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kRed);
  labels.push_back("Post-TDR w/o TECs");
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-45_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower80.root");
  inFileNames_up.push_back("plots/outFile_HPK_PDELoss15p0_Top-45_Tann40_TECs_HLLHCSchedule81_maxPower65.root");
  inFileNames_down.push_back("plots/outFile_HPK_PDELoss07p5_Top-45_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower80.root");
  colors.push_back(kBlue);
  labels.push_back("Post-TDR w/ TECs [for a range of operation scenarios]");
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-42_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower80.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kRed+1);
  labels.push_back("Post-TDR w/ TECs and #DeltaT = 3#circ C lost in thermal int.");
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-39_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower80.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kRed+2);
  labels.push_back("Post-TDR w/ TECs and #DeltaT = 6#circ C lost in thermal int.");
  
  inFileNames.push_back("plots/outFile_HPK_PDELoss15p0_Top-36_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower80.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kRed+3);
  labels.push_back("Post-TDR w/ TECs and #DeltaT = 9#circ C lost in thermal int.");
  
  
  
  
  TCanvas* c = new TCanvas("c","",1200,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,20.,10.,120.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.15,0.94-0.05*inFileNames.size(),0.85,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.04);
  
  int it = 0;
  for(auto inFileName : inFileNames)
  {
    std::cout << labels.at(it) << std::endl;
    
    TFile* inFile = TFile::Open(inFileName.c_str());
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    
    TGraph* g_up;
    TGraph* g_down;
    if( inFileNames_up.at(it) != "NULL" )
    {
      TFile* inFile_up = TFile::Open(inFileNames_up.at(it).c_str());
      g_up = (TGraph*)( inFile_up->Get("g_tResBest_vs_time") );
    }
    if( inFileNames_down.at(it) != "NULL" )
    {
      TFile* inFile_down = TFile::Open(inFileNames_down.at(it).c_str());
      g_down = (TGraph*)( inFile_down->Get("g_tResBest_vs_time") );
    }
    
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(4);
      
      if( inFileNames_up.at(it) != "NULL" )
      {
        for(int point = 0; point < graph->GetN(); ++point)
        {
          double x = graph -> GetPointX(point);
          double y = graph -> GetPointY(point);
          double y_up = g_up -> Eval(x);
          graph -> SetPointEYhigh(point,y_up-y);
        }
      }
      if( inFileNames_down.at(it) != "NULL" )
      {
        for(int point = 0; point < graph->GetN(); ++point)
        {
          double x = graph -> GetPointX(point);
          double y = graph -> GetPointY(point);
          double y_down = g_down -> Eval(x);
          graph -> SetPointEYlow(point,y-y_down);
        }
      }

      graph -> SetFillColor(colors.at(it));
      graph -> SetFillStyle(3001);
      graph -> Draw("L3,same");
    }
    
    legend -> AddEntry(graphs[0],(labels.at(it)).c_str(),"L");
    ++it;
  }

  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  f_TDR -> Draw("same");

  legend -> AddEntry(f_TDR,"TDR goal","L");
  
  legend -> Draw("same");

  c -> Print("HLLHC_final_HPK_thermalInterfaces.png");
}






void drawFinalPlot_final_FBK()
{
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleOffset(0.8, "Y");
  
  std::vector<std::string> inFileNames;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  inFileNames.push_back("plots/outFile_FBK_Top-35_Tann15_HLLHCSchedule8_maxPower65.root");
  colors.push_back(kBlack);
  labels.push_back("no TECs: T^{op.} = -35#circ C T^{ann.} = 15#circ C,   65 mW/ch   #rightarrow  ");
  
  //inFileNames.push_back("plots/outFile_FBK_Top-35_Tann15_HLLHCSchedule81_maxPower65.root");
  //colors.push_back(kGray+2);
  //labels.push_back("no TECs: T^{op.} = -35#circ C T^{ann.} = 15#circ C,   4 days/6 weeks at  0#circ C,   65 mW/ch   #rightarrow  ");

  //inFileNames.push_back("plots/outFile_FBK_Top-35_Tann15_HLLHCSchedule81_maxPower80.root");
  //colors.push_back(kGray+1);
  //labels.push_back("no TECs: T^{op.} = -35#circ C T^{ann.} = 15#circ C,   4 days/6 weeks at  0#circ C,   80 mW/ch   #rightarrow  ");
  
  inFileNames.push_back("plots/outFile_FBK_Top-45_Tann40_TECs_HLLHCSchedule8_maxPower65.root");
  colors.push_back(kRed);
  labels.push_back("w/ TECs: T^{op.} = -45#circ C T^{ann.} = 40#circ C,   4 days/6 weeks at  0#circ C,   65 mW/ch   #rightarrow  ");
  
  inFileNames.push_back("plots/outFile_FBK_Top-45_Tann40_TECs_HLLHCSchedule81_maxPower65.root");
  colors.push_back(kOrange);
  labels.push_back("w/ TECs: T^{op.} = -45#circ C T^{ann.} = 40#circ C,   4 days/6 weeks at 40#circ C,   65 mW/ch   #rightarrow  ");
  
  inFileNames.push_back("plots/outFile_FBK_Top-45_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower65.root");
  colors.push_back(kAzure+10);
  labels.push_back("w/ TECs: T^{op.} = -45#circ C T^{ann.} = 40#circ C,   4 days/6 weeks at 40#circ C,   interfills at 0#circ C,   65 mW/ch   #rightarrow  ");
  
  inFileNames.push_back("plots/outFile_FBK_Top-45_Tann40_TECs_interfillAnnealing_HLLHCSchedule81_maxPower80.root");
  colors.push_back(kBlue);
  labels.push_back("w/ TECs: T^{op.} = -45#circ C T^{ann.} = 40#circ C,   4 days/6 weeks at 40#circ C,   interfills at 0#circ C,   80 mW/ch   #rightarrow  ");
  
  
  
  TCanvas* c = new TCanvas("c","",1400,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,20.,10.,130.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.15,0.94-0.04*inFileNames.size(),0.50,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.02);
  
  int it = 0;
  for(auto inFileName : inFileNames)
  {
    TFile* inFile = TFile::Open(inFileName.c_str());
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(3);
      graph -> Draw("L,same");
    }
    
    legend -> AddEntry(graphs[0],(labels.at(it)+Form("#LT #sigma_{t} #GT = %.0f ps",res.second)).c_str(),"L");
    ++it;
  }

  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  f_TDR -> Draw("same");

  legend -> AddEntry(f_TDR,"TDR","L");
  
  legend -> Draw("same");

  c -> Print("HLLHC_final_FBK.png");
}






void drawAlpha()
{
  gStyle -> SetPadLeftMargin(0.15);
  gStyle -> SetPadRightMargin(0.05);
  gStyle -> SetTitleOffset(1.20, "Y");
  
  std::vector<std::pair<int,int> >vals;
  std::vector<int> TECFlags;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  vals.push_back(std::make_pair<int,int>(-35,-35)); TECFlags.push_back(0); colors.push_back(kBlack); labels.push_back("T_{op} = -35#circ C   T_{ann} = -35#circ C");
  vals.push_back(std::make_pair<int,int>(-35,15));  TECFlags.push_back(0); colors.push_back(57);     labels.push_back("T_{op} = -35#circ C   T_{ann} = +15#circ C");
  vals.push_back(std::make_pair<int,int>(-35,30));  TECFlags.push_back(0); colors.push_back(67);     labels.push_back("T_{op} = -35#circ C   T_{ann} = +30#circ C");
  vals.push_back(std::make_pair<int,int>(-35,40));  TECFlags.push_back(0); colors.push_back(77);     labels.push_back("T_{op} = -35#circ C   T_{ann} = +40#circ C");
  vals.push_back(std::make_pair<int,int>(-35,50));  TECFlags.push_back(0); colors.push_back(87);     labels.push_back("T_{op} = -35#circ C   T_{ann} = +50#circ C");
  vals.push_back(std::make_pair<int,int>(-35,70));  TECFlags.push_back(0); colors.push_back(97);     labels.push_back("T_{op} = -35#circ C   T_{ann} = +70#circ C");
  
  vals.push_back(std::make_pair<int,int>(-45,15));  TECFlags.push_back(1); colors.push_back(57);     labels.push_back("T_{op} = -45#circ C   T_{ann} = +15#circ C");
  vals.push_back(std::make_pair<int,int>(-45,40));  TECFlags.push_back(1); colors.push_back(77);     labels.push_back("T_{op} = -45#circ C   T_{ann} = +40#circ C");
  vals.push_back(std::make_pair<int,int>(-45,70));  TECFlags.push_back(1); colors.push_back(97);     labels.push_back("T_{op} = -45#circ C   T_{ann} = +70#circ C");
  
  
  TCanvas* c = new TCanvas("c","",1500,700);
  c -> Divide(2,1);
  
  c -> cd(1);
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,10.) );
  hPad -> SetTitle(";years from 2027;#alpha [10^{-17} A/cm]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.18,0.54,0.50,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);

  std::map<int,TGraph*> g;
  int it = 0;
  for(auto val : vals)
  {
    TFile* inFile;
    if( TECFlags.at(it) == 0 ) inFile = TFile::Open(Form("plots/outFile_Top%d_Tann%d_HLLHCSchedule3_PDELoss15perc_maxPower65.root",val.first,val.second));
    else                       inFile = TFile::Open(Form("plots/outFile_Top%d_Tann%d_TECs_HLLHCSchedule3_PDELoss15perc_maxPower65.root",val.first,val.second));

    g[it] = (TGraph*)( inFile->Get("g_alphaNorm_vs_time") );
    
    g[it] -> SetLineColor(colors.at(it));
    g[it] -> SetLineWidth(3);
    g[it] -> Draw("L,same");

    if( TECFlags.at(it) == 1 ) g[it] -> SetLineStyle(2);

    legend -> AddEntry(g[it],(labels.at(it)).c_str(),"L");
    
    ++it;
  }

  legend -> Draw("same");
  
  
  c -> cd(2);
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad2 = (TH1F*)( gPad->DrawFrame(0.,0.,10.,2.) );
  hPad2 -> SetTitle(";years from 2027;#alpha / #alpha_{-35#circ C,40#circ C}");
  hPad2 -> Draw();
  
  std::map<int,TGraph*> g_ratio;
  
  it = 0;
  for(auto val : vals)
  {
    TFile* inFile;
    if( TECFlags.at(it) == 0 ) inFile = TFile::Open(Form("plots/outFile_Top%d_Tann%d_HLLHCSchedule3_PDELoss15perc_maxPower65.root",val.first,val.second));
    else                       inFile = TFile::Open(Form("plots/outFile_Top%d_Tann%d_TECs_HLLHCSchedule3_PDELoss15perc_maxPower65.root",val.first,val.second));

    g_ratio[it] = new TGraph();

    for(int point = 0; point < g[it]->GetN(); ++point)
    {
      double x,y;
      g[it] -> GetPoint(point,x,y);
      double xref,yref;
      g[3] -> GetPoint(point,xref,yref);

      if( yref == 0 ) continue;
      
      g_ratio[it] -> SetPoint(g_ratio[it]->GetN(),xref,y/yref);
    }

    
    g_ratio[it] -> SetLineColor(colors.at(it));
    g_ratio[it] -> SetLineWidth(3);
    g_ratio[it] -> Draw("L,same");

    if( TECFlags.at(it) == 1 ) g_ratio[it] -> SetLineStyle(2);

    ++it;
  }

  c -> Print("alphaPlot.png");
}




void drawFinalPlot_final_HPK_dec21()
{
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleOffset(0.8, "Y");
  
  std::vector<std::string> inFileNames;
  std::vector<std::string> inFileNames_up;
  std::vector<std::string> inFileNames_down;
  std::vector<int> colors;
  std::vector<std::string> labels;
  
  inFileNames.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.00_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower80__type2.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(1);
  labels.push_back("RU 2-3");
  
  inFileNames.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.00_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower80__type1.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kRed);
  labels.push_back("RU 0-1");
  
  inFileNames.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.00_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower80__type3.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kBlue);
  labels.push_back("RU 4-5");

  inFileNames.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.00_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower80__type2mod.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kGray);
  labels.push_back("RU 2-3 (type1 thickness)");
  
  inFileNames.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.00_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower80__type3mod.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kBlue-7);
  labels.push_back("RU 4-5 (type1 thickness)");
  
  /*
  inFileNames.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.30_TECs__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower80.root");
  //inFileNames_up.push_back("../plots/outFile__LO1600_tau38.5__HPK_dropPDE0.15_dropGain0.30_TECs__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower50.root");
  //inFileNames_down.push_back("../plots/outFile__LO1600_tau38.5__HPK_dropPDE0.07_dropGain0.30_TECs__Top_-45_Tann1_15_Tann2_50_interfillAnnealing__HLLHCSchedule81__maxPower80.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(50+7*1);
  labels.push_back("T_{annealing} = 40#circ C");
  
  inFileNames.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.30_TECs__Top_-45_Tann1_15_Tann2_50_interfillAnnealing__HLLHCSchedule81__maxPower80.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(50+7*2);
  labels.push_back("T_{annealing} = 50#circ C");
  
  inFileNames.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.30_TECs__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower80.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(50+7*3);
  labels.push_back("T_{annealing} = 60#circ C");
  
  inFileNames.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.30_TECs__Top_-45_Tann1_15_Tann2_70_interfillAnnealing__HLLHCSchedule81__maxPower80.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(50+7*4);
  labels.push_back("T_{annealing} = 70#circ C");
  
  inFileNames.push_back("../plots/outFile__LO1600_tau38.5__HPK_dropPDE0.15_dropGain0.30_TECs__Top_-45_Tann1_15_Tann2_80_interfillAnnealing__HLLHCSchedule81__maxPower80.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(50+7*5);
  labels.push_back("T_{annealing} = 80#circ C");
  
  inFileNames.push_back("../plots/outFile__LO1600_tau38.5__HPK_dropPDE0.15_dropGain0.30_TECs__Top_-45_Tann1_15_Tann2_90_interfillAnnealing__HLLHCSchedule81__maxPower80.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(50+7*6);
  labels.push_back("T_{annealing} = 90#circ C");
  
  inFileNames.push_back("../plots/outFile__LO1600_tau38.5__HPK_dropPDE0.15_dropGain0.30_TECs__Top_-45_Tann1_15_Tann2_100_interfillAnnealing__HLLHCSchedule81__maxPower80.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(50+7*7);
  labels.push_back("T_{annealing} = 100#circ C");
  */
  
  TCanvas* c = new TCanvas("c","",1200,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,20.,13.,120.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.15,0.94-0.04*(inFileNames.size()+4),0.55,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  int it = 0;
  for(auto inFileName : inFileNames)
  {
    std::cout << labels.at(it) << std::endl;
    
    TFile* inFile = TFile::Open(inFileName.c_str());
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    
    TGraph* g_up;
    TGraph* g_down;
    if( inFileNames_up.at(it) != "NULL" )
    {
      TFile* inFile_up = TFile::Open(inFileNames_up.at(it).c_str());
      g_up = (TGraph*)( inFile_up->Get("g_tResBest_vs_time") );
    }
    if( inFileNames_down.at(it) != "NULL" )
    {
      TFile* inFile_down = TFile::Open(inFileNames_down.at(it).c_str());
      g_down = (TGraph*)( inFile_down->Get("g_tResBest_vs_time") );
    }
    
    int it2 = 0;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(4);
      if( it2 == 0 )
        legend -> AddEntry(graph,(labels.at(it)).c_str(),"L");
      ++it2;
      
      if( inFileNames_up.at(it) != "NULL" )
      {
        for(int point = 0; point < graph->GetN(); ++point)
        {
          double x = graph -> GetPointX(point);
          double y = graph -> GetPointY(point);
          double y_up = g_up -> Eval(x);
          graph -> SetPointEYhigh(point,y_up-y);
        }
      }
      if( inFileNames_down.at(it) != "NULL" )
      {
        for(int point = 0; point < graph->GetN(); ++point)
        {
          double x = graph -> GetPointX(point);
          double y = graph -> GetPointY(point);
          double y_down = g_down -> Eval(x);
          graph -> SetPointEYlow(point,y-y_down);
        }
      }

      graph -> SetFillColor(colors.at(it));
      graph -> SetFillStyle(3001);
      graph -> Draw("L3,same");
    }
    
    /*
    g = (TGraph*)( inFile->Get("g_tResBest_stoch_vs_time") );
    res = get_tRes_graphs(g,g_temp);
    graphs = res.first;
    int it = 0;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(kGreen);
      graph -> SetLineWidth(4);
      graph -> Draw("L,same");
v      if( it == 0 )
        legend -> AddEntry(graph,"stochastics","L");
      ++it;
    }
    
    
    g = (TGraph*)( inFile->Get("g_tResBest_noise_vs_time") );
    res = get_tRes_graphs(g,g_temp);
    graphs = res.first;
    it = 0;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(kBlue);
      graph -> SetLineWidth(4);
      graph -> Draw("L,same");
      if( it == 0 )
        legend -> AddEntry(graph,"noise","L");
      ++it;
    }
    

    g = (TGraph*)( inFile->Get("g_tResBest_DCR_vs_time") );
    res = get_tRes_graphs(g,g_temp);
    graphs = res.first;
    it = 0;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(kRed);
      graph -> SetLineWidth(4);
      graph -> Draw("L,same");
      if( it == 0 )
        legend -> AddEntry(graph,"DCR","L");
      ++it;
    }
    */
    
    ++it;
  }

  TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.);
  f_TDR -> SetLineColor(kBlack);
  f_TDR -> SetLineWidth(3);
  f_TDR -> SetLineStyle(7);
  f_TDR -> Draw("same");

  legend -> AddEntry(f_TDR,"TDR goal","L");
  
  legend -> Draw("same");

  c -> Print("HLLHC_final_HPK_tommaso.png");
  c -> Print("HLLHC_final_HPK_tommaso.C");
}












void drawFinalPlot_final_HPK_apr22()
{
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleOffset(0.8, "Y");
  
  std::vector<std::string> inFileNames;
  std::vector<std::string> inFileNames_up;
  std::vector<std::string> inFileNames_down;
  std::vector<int> colors;
  std::vector<std::string> labels;

  /*
  inFileNames.push_back("../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower72__15um_final.root");
  inFileNames_up.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower65__15um_final.root");
  inFileNames_down.push_back("../plots/outFile__LO1300_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower80__15um_final.root");
  colors.push_back(51+0);
  labels.push_back("15 #mum cell size - MS SiPMs");
  
  inFileNames.push_back("../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower72__20um_final.root");
  inFileNames_up.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower65__20um_final.root");
  inFileNames_down.push_back("../plots/outFile__LO1300_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower80__20um_final.root");
  colors.push_back(51+15);
  labels.push_back("20 #mum cell size");
  
  inFileNames.push_back("../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower72__25um_final.root");
  inFileNames_up.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower65__25um_final.root");
  inFileNames_down.push_back("../plots/outFile__LO1300_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower80__25um_final.root");
  colors.push_back(51+15*2);
  labels.push_back("25 #mum cell size");

  inFileNames.push_back("../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower72__30um_final.root");
  inFileNames_up.push_back("../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower65__30um_final.root");
  inFileNames_down.push_back("../plots/outFile__LO1300_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower80__30um_final.root");
  colors.push_back(51+15*3);
  labels.push_back("30 #mum cell size");
  */

  /*
  inFileNames.push_back("../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower72__25um_final.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kBlue+1);
  labels.push_back("25 #mum cell size - 40#circ C annealing");

  inFileNames.push_back("../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um_final.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kRed+1);
  labels.push_back("25 #mum cell size - 60#circ C annealing");
  */
  
  inFileNames.push_back("../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower72__20um_final.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kBlue);
  labels.push_back("20 #mum cell size - Cg = 100 pF");

  inFileNames.push_back("../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower72__20um_Cg20_final.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kRed);
  labels.push_back("20 #mum cell size - Cg = 20 pF");

  /*
  inFileNames.push_back("../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise420__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower72__25um_final.root");
  inFileNames_up.push_back("NULL");
  inFileNames_down.push_back("NULL");
  colors.push_back(kRed+1);
  labels.push_back("25 #mum cell size - bad TOFHIR2B");
  */
  
  TCanvas* c = new TCanvas("c","",1200,700);
  c -> cd();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,30.,13.,120.) );
  hPad -> SetTitle(";years from 2027;#sigma_{t} [ps]");
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.15,0.94-0.04*(inFileNames.size()+4),0.45,0.94);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.05);
  
  int it = 0;
  for(auto inFileName : inFileNames)
  {
    std::cout << labels.at(it) << std::endl;
    
    TFile* inFile = TFile::Open(inFileName.c_str());
    
    TGraph* g = (TGraph*)( inFile->Get("g_tResBest_vs_time") );
    TGraph* g_temp = (TGraph*)( inFile->Get("g_temp_vs_time") );
    
    float tResAve = 0.;
    int N_tResAve = 0;
    std::pair<std::vector<TGraphAsymmErrors*>,float> res = get_tRes_graphs(g,g_temp);
    std::vector<TGraphAsymmErrors*> graphs = res.first;
    
    TGraph* g_up;
    TGraph* g_down;
    if( inFileNames_up.at(it) != "NULL" )
    {
      TFile* inFile_up = TFile::Open(inFileNames_up.at(it).c_str());
      g_up = (TGraph*)( inFile_up->Get("g_tResBest_vs_time") );
    }
    if( inFileNames_down.at(it) != "NULL" )
    {
      TFile* inFile_down = TFile::Open(inFileNames_down.at(it).c_str());
      g_down = (TGraph*)( inFile_down->Get("g_tResBest_vs_time") );
    }
    
    int it2 = 0;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(colors.at(it));
      graph -> SetLineWidth(5);
      if( it2 == 0 )
        legend -> AddEntry(graph,(labels.at(it)).c_str(),"L");
      ++it2;
      
      if( inFileNames_up.at(it) != "NULL" )
      {
        for(int point = 0; point < graph->GetN(); ++point)
        {
          double x = graph -> GetPointX(point);
          double y = graph -> GetPointY(point);
          double y_up = g_up -> Eval(x);
          graph -> SetPointEYhigh(point,y_up-y);
        }
      }
      if( inFileNames_down.at(it) != "NULL" )
      {
        for(int point = 0; point < graph->GetN(); ++point)
        {
          double x = graph -> GetPointX(point);
          double y = graph -> GetPointY(point);
          double y_down = g_down -> Eval(x);
          graph -> SetPointEYlow(point,y-y_down);
        }
      }

      graph -> SetFillColor(colors.at(it));
      graph -> SetFillStyle(3001);
      graph -> Draw("L3,same");
    }
    
    /*
    g = (TGraph*)( inFile->Get("g_tResBest_stoch_vs_time") );
    res = get_tRes_graphs(g,g_temp);
    graphs = res.first;
    int it = 0;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(kGreen);
      graph -> SetLineWidth(4);
      graph -> Draw("L,same");
v      if( it == 0 )
        legend -> AddEntry(graph,"stochastics","L");
      ++it;
    }
    
    
    g = (TGraph*)( inFile->Get("g_tResBest_noise_vs_time") );
    res = get_tRes_graphs(g,g_temp);
    graphs = res.first;
    it = 0;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(kBlue);
      graph -> SetLineWidth(4);
      graph -> Draw("L,same");
      if( it == 0 )
        legend -> AddEntry(graph,"noise","L");
      ++it;
    }
    

    g = (TGraph*)( inFile->Get("g_tResBest_DCR_vs_time") );
    res = get_tRes_graphs(g,g_temp);
    graphs = res.first;
    it = 0;
    for(auto graph : graphs)
    {
      graph -> SetLineColor(kRed);
      graph -> SetLineWidth(4);
      graph -> Draw("L,same");
      if( it == 0 )
        legend -> AddEntry(graph,"DCR","L");
      ++it;
    }
    */
    
    ++it;
  }
  
  /* TF1* f_TDR = new TF1("f_TDR","30+28./10.*x",0.,20.); */
  /* f_TDR -> SetLineColor(kBlack); */
  /* f_TDR -> SetLineWidth(3); */
  /* f_TDR -> SetLineStyle(7); */
  /* f_TDR -> Draw("same"); */

  /* legend -> AddEntry(f_TDR,"TDR goal","L"); */
  
  legend -> Draw("same");

  c -> Print("HLLHC_final_HPK_tommaso.png");
  c -> Print("HLLHC_final_HPK_tommaso.C");
}
