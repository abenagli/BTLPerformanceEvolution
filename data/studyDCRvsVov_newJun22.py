#! /usr/bin/env python3
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
    
    


fluences = ['1E13', '1E14', '2E14']
annealingFactors = {
    '1E13' : 3.6*1.4,   # 3.6 is the alpha factor after 40 mins at 70C, 1.38 is to account for measurements at -40C, while others are at -45C
    '1E14' : 0.89,
    '2E14' : 0.89
}

g_I_Veff = {}
g_I_Veff_input = {}

for fluence in fluences:
    
    g_I_Veff_input[fluence] = ROOT.TGraphErrors()
    g_I_Veff[fluence] = ROOT.TGraphErrors()
    
    inFile = open('HPK_%s_Jun22.txt'%fluence,'r')
    for line in inFile:
        readings = line.rstrip().split()
        g_I_Veff_input[fluence].SetPoint(g_I_Veff_input[fluence].GetN(),float(readings[0]),float(readings[1]))
        g_I_Veff[fluence].SetPoint(g_I_Veff[fluence].GetN(),float(readings[0]),float(readings[1]))
    
    
    norm = g_I_Veff[fluence].Eval(0.4)
    for point in range(0,g_I_Veff[fluence].GetN()):
        g_I_Veff[fluence].SetPoint(point,g_I_Veff[fluence].GetPointX(point),g_I_Veff[fluence].GetPointY(point)/norm)
        g_I_Veff[fluence].SetPointError(point,0,g_I_Veff[fluence].GetPointY(point)*0.01)
    
    fitFunc_preBd = ROOT.TF1('fitFunc_preBd_%s'%fluence,'expo',-1.5,10.)
    fitFunc_preBd.SetLineColor(ROOT.kGreen)
    fitFunc_preBd.SetLineStyle(2)
    fitFunc_preBd.SetLineWidth(1)
    g_I_Veff[fluence].Fit(fitFunc_preBd,'NRS','',-1.5,-0.5)
    
    fitFunc_bd = ROOT.TF1('fitFunc_bd_%s'%fluence,'expo',-1.5,10)
    fitFunc_bd.SetLineColor(ROOT.kMagenta)
    fitFunc_bd.SetLineStyle(2)
    fitFunc_bd.SetLineWidth(1)
    g_I_Veff[fluence].Fit(fitFunc_bd,'QNS','',-0.2,0.2)
    
    fitFunc_postBd = ROOT.TF1('fitFunc_postBd_%s'%fluence,'expo',-1.5,10)
    fitFunc_postBd.SetLineColor(ROOT.kCyan)
    fitFunc_postBd.SetLineStyle(2)
    fitFunc_postBd.SetLineWidth(1)
    g_I_Veff[fluence].Fit(fitFunc_postBd,'QNS','',1.,2.)
    
    fitFunc_all = ROOT.TF1('fitFunc_all_%s'%fluence,'expo(2)+expo(4)*0.5*(1.+TMath::Erf(-(x-[0])/[1]))+0.5*(1.+TMath::Erf((x-[0])/[1]))*expo(6)',-1.5,10.)
    fitFunc_all.SetParameters(
        0.3,0.2,
        fitFunc_preBd.GetParameter(0),fitFunc_preBd.GetParameter(1),
        fitFunc_bd.GetParameter(0),fitFunc_bd.GetParameter(1),
        fitFunc_postBd.GetParameter(0),fitFunc_postBd.GetParameter(1)
    )
    fitFunc_all.FixParameter(7,fitFunc_postBd.GetParameter(1))
    
    g_I_Veff[fluence].Fit(fitFunc_all,'QNS','',-1.5,2.)
    

    c = ROOT.TCanvas('c1_%s'%fluence,'c1_%s'%fluence,1600,700)
    c.Divide(2,1)
    
    c.cd(1)
    hPad1 = ROOT.gPad.DrawFrame(-0.1,0.,2.5,200.)
    hPad1.SetTitle(";V_{ov}^{eff} [V]; current / current(0.8 V)");
    hPad1.Draw();
    ROOT.gPad.SetGridx();
    ROOT.gPad.SetGridy();
    g_I_Veff[fluence].SetMarkerColor(ROOT.kBlack)
    g_I_Veff[fluence].SetMarkerSize(1.)
    g_I_Veff[fluence].SetMarkerStyle(20)
    g_I_Veff[fluence].Draw('PL,same')
    fitFunc_preBd.Draw('same')
    fitFunc_bd.Draw('same')
    fitFunc_postBd.Draw('same')
    fitFunc_all.Draw('same')
    
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize( 0.05 )
    lat.DrawLatex( 0.40, 0.91, 'HPK %s'%fluence)
    
    c.cd(2)
    hPad3 = ROOT.gPad.DrawFrame(-0.1,0.001,2.5,1000.)
    hPad3.SetTitle(";V_{ov}^{eff} [V]; current / current(0.8 V)");
    hPad3.Draw();
    ROOT.gPad.SetGridx();
    ROOT.gPad.SetGridy();
    ROOT.gPad.SetLogy()
    g_I_Veff[fluence].SetMarkerColor(ROOT.kBlack)
    g_I_Veff[fluence].SetMarkerSize(1.)
    g_I_Veff[fluence].SetMarkerStyle(20)
    g_I_Veff[fluence].Draw('PL,same')
    fitFunc_preBd.Draw('same')
    fitFunc_bd.Draw('same')
    fitFunc_postBd.Draw('same')
    fitFunc_all.Draw('same')
    
    c.Print('c_IV_%s.png'%fluence)



c2 = ROOT.TCanvas('c2','c2',1600,700)
c2.Divide(2,1)

c2.cd(1)
hPad21 = ROOT.gPad.DrawFrame(0.,0.,2,7.)
hPad21.SetTitle(";V_{ov}^{eff} [V];ratio fluence / 1E14");
hPad21.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();
c2.cd(2)
hPad22 = ROOT.gPad.DrawFrame(0.,0.01,2,10.)
hPad22.SetTitle(";V_{ov}^{eff} [V];ratio fluence / 1E14");
hPad22.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();
ROOT.gPad.SetLogy()

ratios = {}
it = 1
leg = ROOT.TLegend(0.15, 0.88-0.05*len(fluences), 0.40, 0.88)
leg.SetFillColor(0)
leg.SetFillStyle(1000)
leg.SetTextFont(42)
leg.SetTextSize(0.04)
for fluence in fluences:
    ratios[fluence] = ROOT.TGraph()
    nPoints = 100
    for point in range(0, nPoints):
        xx = 0. + (2. - (0.))/nPoints*point
        ratios[fluence].SetPoint(point,xx,g_I_Veff_input[fluence].Eval(xx)/annealingFactors[fluence]/g_I_Veff_input['1E14'].Eval(xx)*annealingFactors['1E14'])
    ratios[fluence].SetLineColor(it)
    ratios[fluence].SetMarkerColor(it)
    leg.AddEntry(ratios[fluence],fluence,'L')
    c2.cd(1)
    ratios[fluence].Draw('L,same')
    c2.cd(2)
    ratios[fluence].Draw('L,same')
    it += 1
c2.cd(1)
leg.Draw('same')
c2.cd(2)
leg.Draw('same')
c2.Print('plot2.png')



c3 = ROOT.TCanvas('c3','c3',800,700)
hPad3 = ROOT.gPad.DrawFrame(0.,0.,2,2.)
hPad3.SetTitle(";V_{ov}^{eff} [V];ratio");
hPad3.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();

ratios2 = {}
it = 1
leg = ROOT.TLegend(0.15, 0.88-0.05*len(fluences), 0.40, 0.88)
leg.SetFillColor(0)
leg.SetFillStyle(1000)
leg.SetTextFont(42)
leg.SetTextSize(0.04)
for fluence in fluences:
    ratios2[fluence] = ROOT.TGraph()
    nPoints = 100
    for point in range(0, nPoints):
        xx = 0. + (2. - (0.))/nPoints*point
        ratios2[fluence].SetPoint(point,xx,g_I_Veff_input[fluence].Eval(xx)/annealingFactors[fluence]/float(fluence)/g_I_Veff_input['2E14'].Eval(xx)*annealingFactors['2E14']*float('2E14'))
    ratios2[fluence].SetLineColor(it)
    ratios2[fluence].SetMarkerColor(it)
    leg.AddEntry(ratios[fluence],fluence,'L')
    ratios2[fluence].Draw('L,same')
    it += 1
leg.Draw('same')
c3.Print('plot3.png')





outfile = ROOT.TFile('DCRParams_new_Jun2022','RECREATE')
for iPar in range(0, 8):
    g_pars[iPar].Write('g_par%d'%iPar)
