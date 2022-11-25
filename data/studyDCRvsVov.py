#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time

import ROOT

#set the tdr style
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetTitleOffset(1.05,'Y')
ROOT.gStyle.SetLabelSize(0.04)
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
#ROOT.gROOT.SetBatch(True)
ROOT.gROOT.SetBatch(False)


def myfunc_DCR(x, par):
    xx = x[0]
    x0 = par[4]
    A = par[0]
    C = par[1]
    D = par[2]
    E = par[3]
    B = D*E*math.exp(E*x0) - 2.*C*x0
    F = A + B*x0 + C*x0*x0 - D*math.exp(E*x0)
    if xx < x0 :
        return par[5] * ( A + B*xx + C*xx*xx )
    else :
        return par[5] * ( D*math.exp(E*xx) + F )
    
    

c = ROOT.TCanvas('c1','c1',1600,1000)
c.Divide(2,2)



fluence = 1E13
Vbr = 37.62
list_Vov = [1.5, 1.65, 2, 2.5, 2.8, 3.0, 3.2]

g_I_Vset_1E13 = ROOT.TGraphErrors()
g_I_Veff_1E13 = ROOT.TGraphErrors()

rootFile = ROOT.TFile('logs_3.00/logIV_ASIC2_ALDOA_ch3_time_2021-10-16_10:34:48.root','READ')
graph = rootFile.Get('g_IV')
for point in range(0, graph.GetN()):
    g_I_Vset_1E13.SetPoint(g_I_Vset_1E13.GetN(),graph.GetPointX(point),graph.GetPointY(point))

for Vov in list_Vov:
    logFiles = glob.glob('logs_3.00/logStressTest_ch3_bv%.2f*.root'%(Vbr+Vov))
    #print logFiles
    for logFile in logFiles:
        rootFile = ROOT.TFile(logFile,'READ') 
        graph = rootFile.Get('g_I')
        fitFunc = ROOT.TF1('fitFunc','pol0',250,300)
        graph.Fit(fitFunc,'QNRS')
        g_I_Vset_1E13.SetPoint(g_I_Vset_1E13.GetN(),Vbr+Vov,fitFunc.GetParameter(0))

for point in range(0, g_I_Vset_1E13.GetN()):
    g_I_Veff_1E13.SetPoint(point,g_I_Vset_1E13.GetPointX(point)-Vbr-g_I_Vset_1E13.GetPointY(point)*25/1E06,g_I_Vset_1E13.GetPointY(point))


norm = g_I_Veff_1E13.Eval(0.8)
for point in range(0,g_I_Veff_1E13.GetN()):
    g_I_Veff_1E13.SetPoint(point,g_I_Veff_1E13.GetPointX(point),g_I_Veff_1E13.GetPointY(point)/norm)
    g_I_Veff_1E13.SetPointError(point,0,g_I_Veff_1E13.GetPointY(point)*0.01)

fitFunc_preBd_1E13 = ROOT.TF1('fitFunc_preBd_1E13','expo',-1.5,10.)
fitFunc_preBd_1E13.SetLineColor(ROOT.kGreen)
fitFunc_preBd_1E13.SetLineStyle(2)
fitFunc_preBd_1E13.SetLineWidth(1)
g_I_Veff_1E13.Fit(fitFunc_preBd_1E13,'NRS','',-1.5,-0.5)

fitFunc_bd_1E13 = ROOT.TF1('fitFunc_bd_1E13','expo',-1.5,10)
fitFunc_bd_1E13.SetLineColor(ROOT.kMagenta)
fitFunc_bd_1E13.SetLineStyle(2)
fitFunc_bd_1E13.SetLineWidth(1)
g_I_Veff_1E13.Fit(fitFunc_bd_1E13,'QNS','',-0.2,0.2)

fitFunc_postBd_1E13 = ROOT.TF1('fitFunc_postBd_1E13','expo',-1.5,10)
fitFunc_postBd_1E13.SetLineColor(ROOT.kCyan)
fitFunc_postBd_1E13.SetLineStyle(2)
fitFunc_postBd_1E13.SetLineWidth(1)
g_I_Veff_1E13.Fit(fitFunc_postBd_1E13,'QNS','',1.,2.)

fitFunc_all_1E13 = ROOT.TF1('fitFunc_all_1E13','expo(2)+expo(4)*0.5*(1.+TMath::Erf(-(x-[0])/[1]))+0.5*(1.+TMath::Erf((x-[0])/[1]))*expo(6)',-1.5,10.)
fitFunc_all_1E13.SetParameters(
    0.3,0.2,
    fitFunc_preBd_1E13.GetParameter(0),fitFunc_preBd_1E13.GetParameter(1),
    fitFunc_bd_1E13.GetParameter(0),fitFunc_bd_1E13.GetParameter(1),
    fitFunc_postBd_1E13.GetParameter(0),fitFunc_postBd_1E13.GetParameter(1)
)
fitFunc_all_1E13.FixParameter(7,fitFunc_postBd_1E13.GetParameter(1))

g_I_Veff_1E13.Fit(fitFunc_all_1E13,'QNS','',-1.5,2.)


c.cd(1)
hPad1 = ROOT.gPad.DrawFrame(-1.5,0.,2.5,10.)
hPad1.SetTitle(";V_{ov}^{eff} [V]; current / current(0.8 V)");
hPad1.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();
g_I_Veff_1E13.SetMarkerColor(ROOT.kBlack)
g_I_Veff_1E13.SetMarkerSize(1.)
g_I_Veff_1E13.SetMarkerStyle(20)
g_I_Veff_1E13.Draw('PL,same')
fitFunc_preBd_1E13.Draw('same')
fitFunc_bd_1E13.Draw('same')
fitFunc_postBd_1E13.Draw('same')
fitFunc_all_1E13.Draw('same')

lat_1E13 = ROOT.TLatex()
lat_1E13.SetNDC()
lat_1E13.SetTextSize( 0.05 )
lat_1E13.DrawLatex( 0.40, 0.91, 'HPK 1E13')

c.cd(3)
hPad3 = ROOT.gPad.DrawFrame(-1.5,0.0001,2.5,100.)
hPad3.SetTitle(";V_{ov}^{eff} [V]; current / current(0.8 V)");
hPad3.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();
ROOT.gPad.SetLogy()
g_I_Veff_1E13.SetMarkerColor(ROOT.kBlack)
g_I_Veff_1E13.SetMarkerSize(1.)
g_I_Veff_1E13.SetMarkerStyle(20)
g_I_Veff_1E13.Draw('PL,same')
fitFunc_preBd_1E13.Draw('same')
fitFunc_bd_1E13.Draw('same')
fitFunc_postBd_1E13.Draw('same')
fitFunc_all_1E13.Draw('same')




fluence = 2E14
Vbr = 38.83
list_Vov = [1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.5, 3.]

g_I_Vset_2E14 = ROOT.TGraphErrors()
g_I_Veff_2E14 = ROOT.TGraphErrors()

rootFile = ROOT.TFile('logs_26.00/logIV_ASIC2_ALDOA_ch3_time_2021-10-26_20:00:32.root','READ')
graph = rootFile.Get('g_IV')
for point in range(0, graph.GetN()):
    g_I_Vset_2E14.SetPoint(g_I_Vset_2E14.GetN(),graph.GetPointX(point),graph.GetPointY(point))

for Vov in list_Vov:
    logFiles = glob.glob('logs_26.00/logStressTest_ch3_bv%.2f_HPK_2E14_T-40C*.root'%(Vbr+Vov))
    #print logFiles
    for logFile in logFiles:
        rootFile = ROOT.TFile(logFile,'READ') 
        graph = rootFile.Get('g_I')
        fitFunc = ROOT.TF1('fitFunc','pol0',250,300)
        graph.Fit(fitFunc,'QNRS')
        g_I_Vset_2E14.SetPoint(g_I_Vset_2E14.GetN(),Vbr+Vov,fitFunc.GetParameter(0))

for point in range(0, g_I_Vset_2E14.GetN()):
    g_I_Veff_2E14.SetPoint(point,g_I_Vset_2E14.GetPointX(point)-Vbr-g_I_Vset_2E14.GetPointY(point)*25/1E06,g_I_Vset_2E14.GetPointY(point))


norm = g_I_Veff_2E14.Eval(0.8)
for point in range(0,g_I_Veff_2E14.GetN()):
    g_I_Veff_2E14.SetPoint(point,g_I_Veff_2E14.GetPointX(point),g_I_Veff_2E14.GetPointY(point)/norm)
    g_I_Veff_2E14.SetPointError(point,0,0.01*g_I_Veff_2E14.GetPointY(point))
    
fitFunc_preBd_2E14 = ROOT.TF1('fitFunc_preBd_2E14','expo',-1.5,10.)
fitFunc_preBd_2E14.SetLineColor(ROOT.kGreen)
fitFunc_preBd_2E14.SetLineStyle(2)
fitFunc_preBd_2E14.SetLineWidth(1)
g_I_Veff_2E14.Fit(fitFunc_preBd_2E14,'NRS','',-1.5,-0.5)

fitFunc_bd_2E14 = ROOT.TF1('fitFunc_bd_2E14','expo',-1.5,10)
fitFunc_bd_2E14.SetLineColor(ROOT.kMagenta)
fitFunc_bd_2E14.SetLineStyle(2)
fitFunc_bd_2E14.SetLineWidth(1)
g_I_Veff_2E14.Fit(fitFunc_bd_2E14,'QNS','',-0.2,0.2)

fitFunc_postBd_2E14 = ROOT.TF1('fitFunc_postBd_2E14','expo',-1.5,10)
fitFunc_postBd_2E14.SetLineColor(ROOT.kCyan)
fitFunc_postBd_2E14.SetLineStyle(2)
fitFunc_postBd_2E14.SetLineWidth(1)
g_I_Veff_2E14.Fit(fitFunc_postBd_2E14,'QNS','',1.,2.)

fitFunc_all_2E14 = ROOT.TF1('fitFunc_all_2E14','expo(2)+expo(4)*0.5*(1.+TMath::Erf(-(x-[0])/[1]))+0.5*(1.+TMath::Erf((x-[0])/[1]))*expo(6)',-1.5,10.)
fitFunc_all_2E14.SetParameters(
    0.3,0.2,
    fitFunc_preBd_2E14.GetParameter(0),fitFunc_preBd_2E14.GetParameter(1),
    fitFunc_bd_2E14.GetParameter(0),fitFunc_bd_2E14.GetParameter(1),
    fitFunc_postBd_2E14.GetParameter(0),fitFunc_postBd_2E14.GetParameter(1)
)
g_I_Veff_2E14.Fit(fitFunc_all_2E14,'QNS','',-1.5,2.)
        

c.cd(2)
hPad2 = ROOT.gPad.DrawFrame(-1.5,0.,2.5,30.)
hPad2.SetTitle(";V_{ov}^{eff} [V]; current / current(0.8 V)");
hPad2.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();
g_I_Veff_2E14.SetMarkerColor(ROOT.kBlack)
g_I_Veff_2E14.SetMarkerSize(1.)
g_I_Veff_2E14.SetMarkerStyle(20)
g_I_Veff_2E14.Draw('PL,same')
fitFunc_preBd_2E14.Draw('same')
fitFunc_bd_2E14.Draw('same')
fitFunc_postBd_2E14.Draw('same')
fitFunc_all_2E14.Draw('same')

lat_1E13 = ROOT.TLatex()
lat_1E13.SetNDC()
lat_1E13.SetTextSize( 0.05 )
lat_1E13.DrawLatex( 0.40, 0.91, 'HPK 2E14')

c.cd(4)
hPad4 = ROOT.gPad.DrawFrame(-1.5,0.0001,2.5,100.)
hPad4.SetTitle(";V_{ov}^{eff} [V]; current / current(0.8 V)");
hPad4.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();
ROOT.gPad.SetLogy()
g_I_Veff_2E14.SetMarkerColor(ROOT.kBlack)
g_I_Veff_2E14.SetMarkerSize(1.)
g_I_Veff_2E14.SetMarkerStyle(20)
g_I_Veff_2E14.Draw('PL,same')
fitFunc_preBd_2E14.Draw('same')
fitFunc_bd_2E14.Draw('same')
fitFunc_postBd_2E14.Draw('same')
fitFunc_all_2E14.Draw('same')

c.Print('plot.png')



c2 = ROOT.TCanvas('c2','c2')
hPad22 = ROOT.gPad.DrawFrame(-1.5,0.,2.5,8.)
hPad22.SetTitle(";V_{ov}^{eff} [V];ratio 2E14 / 1E13");
hPad22.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();

ratio = ROOT.TGraph()
nPoints = 100
for point in range(0, nPoints):
    xx = -1. + (2. - (-1.))/nPoints*point
    ratio.SetPoint(point,xx,fitFunc_all_2E14.Eval(xx)/fitFunc_all_1E13.Eval(xx))

ratio.Draw('L,same')
c2.Print('plot2.png')



g_pars = {}
for iPar in range(0, 8):
    g_pars[iPar] = ROOT.TGraphErrors()
    g_pars[iPar].SetPoint(0,1E13,fitFunc_all_1E13.GetParameter(iPar))
    g_pars[iPar].SetPointError(0,0.,fitFunc_all_1E13.GetParError(iPar))
    g_pars[iPar].SetPoint(1,2E14,fitFunc_all_2E14.GetParameter(iPar))
    g_pars[iPar].SetPointError(1,0.,fitFunc_all_2E14.GetParError(iPar))
    c3 = ROOT.TCanvas('c3','c3')
    g_pars[iPar].Draw('APL')
    c3.Print('plot_iPar%d.png'%iPar)
    c3.Delete()


fluences = [0., 1E13, 5E13, 1E14, 2E14, 3E14]
fitFunc_fluences = {}

for fluence in fluences:
    fitFunc_fluences[fluence] = ROOT.TF1('fitFunc_%.2e'%fluence,'expo(2)+expo(4)*0.5*(1.+TMath::Erf(-(x-[0])/[1]))+0.5*(1.+TMath::Erf((x-[0])/[1]))*expo(6)',-1.5,10.)    
    for iPar in range(0, 8):
        fitFunc_fluences[fluence].SetParameter(iPar,g_pars[iPar].Eval(fluence))

c4 = ROOT.TCanvas('c4','c4',1600,1000)
c4.Divide(2,1)
c4.cd(1)
hPad42 = ROOT.gPad.DrawFrame(-1.5,0.,2.5,30.)
hPad42.SetTitle(";V_{ov}^{eff} [V]; current / current(0.8 V)");
hPad42.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();
fluenceIt = 0
for fluence in fluences:
    print(fluence)
    fitFunc_fluences[fluence].SetLineColor(51+fluenceIt*8)
    fitFunc_fluences[fluence].SetLineWidth(5)
    if abs(1.-fluence/2E14) < 0.1 or abs(1.-fluence/1E13) < 0.1 :
        fitFunc_fluences[fluence].SetLineStyle(7)
    fitFunc_fluences[fluence].Draw('same')
    fluenceIt += 1
c4.cd(2)
hPad4 = ROOT.gPad.DrawFrame(-1.5,0.0001,2.5,100.)
hPad4.SetTitle(";V_{ov}^{eff} [V]; current / current(0.8 V)");
hPad4.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();
ROOT.gPad.SetLogy()
for fluence in fluences:
    fitFunc_fluences[fluence].Draw('same')
c4.Print('plot_fluences.png')


outfile = ROOT.TFile('DCRParams_new_Dic2021','RECREATE')
for iPar in range(0, 8):
    g_pars[iPar].Write('g_par%d'%iPar)
