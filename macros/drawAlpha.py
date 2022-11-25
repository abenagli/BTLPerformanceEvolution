#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time
import numpy as np
from collections import OrderedDict

import ROOT
import tdrstyle

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetTitleOffset(1.05,'Y')
ROOT.gStyle.SetLabelSize(0.04)
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.SetBatch(False)



def GetTimeForAlpha(graph,alpha,tMin):
    for point in range(0,graph.GetN()):
        if graph.GetPointX(point) < tMin:
            continue
        if graph.GetPointY(point) < alpha:
            return graph.GetPointX(point)

    return -999.



c = ROOT.TCanvas('c_tRes_vs_tau','c_tRes_vs_tau',2400,1400)
c.Divide(2,1)
c.cd(1)
hPad1 = ROOT.gPad.DrawFrame(0.,0.,365.,7.)
hPad1.SetTitle(";days;#alpha [10^{-17} A/cm]")
hPad1.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
c.cd(2)
hPad2 = ROOT.gPad.DrawFrame(1.,0.,365.,7.)
hPad2.SetTitle(";days;#alpha [10^{-17} A/cm]")
hPad2.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogx()

temps = [70, 80, 90, 100, 110, 120, 130]

leg = ROOT.TLegend(0.15, 0.94-0.035*len(temps), 0.95, 0.94)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.03)

it = 0
for temp in temps:
    print temp
    infile = ROOT.TFile('../plots/computeAlpha_Tann_%.0f.root'%temp)
    graph = infile.Get('g_alphaNorm_vs_time')
    graph.SetLineColor(99-8*it)
    graph.SetLineWidth(2)
    
    tReq = GetTimeForAlpha(graph,0.67,30.)
    if tReq > 0.:
        leg.AddEntry(graph,'T_{ann.} = %.0f, time for #alpha = 0.67: %.1f days'%(temp,tReq-30.),'L')
    else:
        leg.AddEntry(graph,'T_{ann.} = %.0f, time for #alpha = 0.67: > 1 year'%(temp),'L')
    c.cd(1)
    graph.Draw('L,same')
    c.cd(2)
    graph.Draw('L,same')
    it += 1

c.cd(1)
leg.Draw('same')
c.cd(2)
leg.Draw('same')
c.Print('c_alpha.png')

