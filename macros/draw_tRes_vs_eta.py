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

infilenames['1nominal'] = '../plots/outFileRecHits__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower75__25um__RUs__newHPKPars__final.root'
#infilenames['1nominalAnn'] = '../plots/outFileRecHits__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower75__25um__RUs__newHPKPars__final.root'
infilenames['2thicker'] = '../plots/outFileRecHits__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower75__25um__RUs__thicker__newHPKPars__final.root'
infilenames['3thickest'] = '../plots/outFileRecHits__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower75__25um__RUs__thickest__newHPKPars__final.root'

#infilenames['1-75 mW'] = '../plots/outFileRecHits__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower75__25um__RUs__newHPKPars__final.root'
#infilenames['2-50 mW'] = '../plots/outFileRecHits__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower50__25um__RUs__newHPKPars__final.root'
#infilenames['3-100 mW'] = '../plots/outFileRecHits__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower100__25um__RUs__newHPKPars__final.root'

keys = infilenames.keys()
labels = [ it for it in keys]
labels = list(set(labels))
labels = sorted(labels)
print(labels)

p_tResBest_vs_eta = {}
infiles = {}
for key in keys:
    label = key
    infiles[key] = ROOT.TFile(infilenames[key])
    p_tResBest_vs_eta[label] = infiles[key].Get("p_tResBest_vs_eta")


c = ROOT.TCanvas('c','c',1100,700)
hPad = ROOT.gPad.DrawFrame(0.,28.,1.52,82.)
hPad.SetTitle(";#eta;#sigma_{t} [ps]")
hPad.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()

etas = [ 0.3520, 0.6645, 0.9310, 1.1550, 1.3440 ]
it = 0
for label in labels:
    p_tResBest_vs_eta[label].SetLineColor(1+it)
    p_tResBest_vs_eta[label].SetLineWidth(3)
    p_tResBest_vs_eta[label].Draw('hist,same')

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    latex.SetTextColor(1+it)
    if label == '1nominal':
        latex.DrawLatex(0.30,0.40-0.05*it,'RUs: 2#timesT1, 2#timesT2, 2#timesT3, annealing 40#circ C')
    if label == '1nominalAnn':
        latex.DrawLatex(0.30,0.40-0.05*it,'T1-T1-T2-T2-T3-T3, annealing 60#circ C')
    if label == '2thicker':
        latex.DrawLatex(0.30,0.40-0.05*it,'RUs: 3#timesT1, 3#timesT2, annealing 40#circ C')
    if label == '3thickest':
        latex.DrawLatex(0.30,0.40-0.05*it,'RUs: 6#timesT1, annealing 40#circ C')        
    if label == '1-75 mW':
        latex.DrawLatex(0.30,0.40-0.05*it,'T1-T1-T2-T2-T3-T3, annealing 60#circ C, 75 mW/ch.')
    if label == '2-50 mW':
        latex.DrawLatex(0.30,0.40-0.05*it,'T1-T1-T2-T2-T3-T3, annealing 60#circ C, 50 mW/ch.')
    if label == '3-100 mW':
        latex.DrawLatex(0.30,0.40-0.05*it,'T1-T1-T2-T2-T3-T3, annealing 60#circ C, 100 mW/ch.')        
    latex.Draw('same')
    it += 1

lines = {}
for eta in etas:
    lines[eta] = ROOT.TLine(eta,28.,eta,82)
    lines[eta].SetLineStyle(7)
    lines[eta].Draw('same')

c.Print('c_tRes_vs_eta.png')
