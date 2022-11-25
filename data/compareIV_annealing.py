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
    
    


labels = [ 'FBK_2E14_3daysat110C', 'FBK_2E14_3daysat90C', 'FBK_2E14_40minat70C', 'FBK_2E14_noAnn' ]

scalings = {
    'FBK_2E14_3daysat110C' : 1.,
    'FBK_2E14_3daysat90C'  : 1.45,
    'FBK_2E14_40minat70C'  : 3.2,
    'FBK_2E14_noAnn'       : 6.5
    }

g_I_Veff = {}
g_I_Veff_scaled = {}

c = ROOT.TCanvas('c1','c1',1600,700)
c.Divide(2,1)
c.cd(1)
hPad1 = ROOT.gPad.DrawFrame(-0.1,0.,2.5,3.)
hPad1.SetTitle(";V_{ov}^{eff} [V]; current [mA]");
hPad1.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();
c.cd(2)
hPad2 = ROOT.gPad.DrawFrame(-0.1,0.,2.5,3.)
hPad2.SetTitle(";V_{ov}^{eff} [V]; current scaled [mA]");
hPad2.Draw();
hPad2.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();

latexs = {}
it = 1
for label in labels:
    c.cd(1)
    g_I_Veff[label] = ROOT.TGraphErrors('/Users/abenagli/BTLPerformanceEvolution/data/'+label+'.txt')
    g_I_Veff[label].SetMarkerColor(51+10*it)
    g_I_Veff[label].SetLineColor(51+10*it)
    g_I_Veff[label].SetMarkerSize(1.)
    g_I_Veff[label].SetMarkerStyle(20)
    g_I_Veff[label].Draw('PL,same')
    c.cd(2)
    g_I_Veff_scaled[label] = ROOT.TGraphErrors()
    for point in range(g_I_Veff[label].GetN()):
        x = g_I_Veff[label].GetPointX(point)
        y = g_I_Veff[label].GetPointY(point)
        g_I_Veff_scaled[label].SetPoint(point,x,y/scalings[label])
    g_I_Veff_scaled[label].SetMarkerColor(51+10*it)
    g_I_Veff_scaled[label].SetLineColor(51+10*it)
    g_I_Veff_scaled[label].SetMarkerSize(1.)
    g_I_Veff_scaled[label].SetMarkerStyle(20)
    g_I_Veff_scaled[label].Draw('PL,same')
    latexs[label] = ROOT.TLatex()
    latexs[label].SetNDC()
    latexs[label].SetTextSize( 0.04 )
    latexs[label].SetTextFont(42)
    latexs[label].SetTextColor(51+10*it)
    latexs[label].DrawLatex( 0.20, 0.90-0.05*it, '%s - scaling: %.1f'%(label,scalings[label]))
    it += 1
c.Print('c_compareIV_annealing.png')
