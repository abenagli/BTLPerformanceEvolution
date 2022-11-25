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




infilenames = OrderedDict()

infilenames[('1nominal',1)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU1__newHPKPars__final.root'
infilenames[('1nominal',2)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU2__newHPKPars__final.root'
infilenames[('1nominal',3)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU3__newHPKPars__final.root'
infilenames[('1nominal',4)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU4__newHPKPars__final.root'
infilenames[('1nominal',5)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU5__newHPKPars__final.root'
infilenames[('1nominal',6)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU6__newHPKPars__final.root'

infilenames[('2thick',1)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU1__newHPKPars__final.root'
infilenames[('2thick',2)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU2__newHPKPars__final.root'
infilenames[('2thick',3)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU3t1__newHPKPars__final.root'
infilenames[('2thick',4)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU4__newHPKPars__final.root'
infilenames[('2thick',5)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU5t2__newHPKPars__final.root'
infilenames[('2thick',6)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU6t2__newHPKPars__final.root'

infilenames[('3thickest',1)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU1__newHPKPars__final.root'
infilenames[('3thickest',2)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU2__newHPKPars__final.root'
infilenames[('3thickest',3)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU3t1__newHPKPars__final.root'
infilenames[('3thickest',4)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU4t1__newHPKPars__final.root'
infilenames[('3thickest',5)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU5t1__newHPKPars__final.root'
infilenames[('3thickest',6)] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower72__25um__RU6t1__newHPKPars__final.root'

keys = infilenames.keys()
labels = [ it[0] for it in keys]
labels = list(set(labels))
labels = sorted(labels)
print(labels)

g_tResBest_BoO_vs_RU = {}
g_tResBest_EoO_vs_RU = {}
for label in labels:
    g_tResBest_BoO_vs_RU[label] = ROOT.TGraph()
    g_tResBest_EoO_vs_RU[label] = ROOT.TGraph()

infiles = {}
for key in keys:
    label = key[0]
    infiles[key] = ROOT.TFile(infilenames[key])
    
    graph = infiles[key].Get("g_tResBest_vs_time")
    for point in range(graph.GetN()):
        x = graph.GetPointX(point)
        y = graph.GetPointY(point)

        if x >= 11.5:
            g_tResBest_EoO_vs_RU[label].SetPoint(g_tResBest_EoO_vs_RU[label].GetN(),key[1],y)
            break
    for point in range(graph.GetN()):
        x = graph.GetPointX(point)
        y = graph.GetPointY(point)
        
        if x >= 0.75:
            g_tResBest_BoO_vs_RU[label].SetPoint(g_tResBest_BoO_vs_RU[label].GetN(),key[1],y)
            break        


    
c = ROOT.TCanvas('c','c',1200,700)
hPad = ROOT.gPad.DrawFrame(0.,0.,7.,80.)
hPad.SetTitle(";#RU;#sigma_{t} [ps]")
hPad.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()

it = 1
for label in labels:
    g_tResBest_EoO_vs_RU[label].SetMarkerSize(1.2)
    g_tResBest_EoO_vs_RU[label].SetMarkerColor(it)
    g_tResBest_EoO_vs_RU[label].SetLineColor(it)
    g_tResBest_EoO_vs_RU[label].Draw('PL,same')

    g_tResBest_BoO_vs_RU[label].SetMarkerSize(1.2)
    g_tResBest_BoO_vs_RU[label].SetMarkerColor(it)
    g_tResBest_BoO_vs_RU[label].SetMarkerStyle(24)
    g_tResBest_BoO_vs_RU[label].SetLineColor(it)
    g_tResBest_BoO_vs_RU[label].SetLineStyle(2)
    g_tResBest_BoO_vs_RU[label].Draw('PL,same')

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    latex.SetTextColor(it)
    if label == '1nominal':
        latex.DrawLatex(0.30,0.40-0.05*it,'RU1:T1   RU2: T1 /  RU3: T2   RU4: T2 /  RU5: T3   RU6: T3')
    if label == '2thick':
        latex.DrawLatex(0.30,0.40-0.05*it,'RU1:T1   RU2: T1 /  RU3: T1   RU4: T2 /  RU5: T2   RU6: T2')
    if label == '3thickest':
        latex.DrawLatex(0.30,0.40-0.05*it,'RU1:T1   RU2: T1 /  RU3: T1   RU4: T1 /  RU5: T1   RU6: T1')
    latex.Draw('same')
    it += 1

c.Print('c_tRes_vs_RU.png')
