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
ROOT.gStyle.SetOptStat(11111)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetTitleOffset(1.05,'Y')
ROOT.gStyle.SetLabelSize(0.04)
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
ROOT.gROOT.SetBatch(False)


#11.5 0.75

infile = ROOT.TFile('../data/DQM_V0001_R000000001__RelValMinBias_14TeV__CMSSW_12_6_0_pre4-125X_mcRun4_realistic_v2_2026D88noPU-v1__DQMIO.root')

c1 = ROOT.TCanvas('c1','c1')
c1.cd()
histo = infile.Get('DQMData/Run 1/MTD/Run summary/BTL/LocalReco/BtlCluEnergy')
histo.Draw()
latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextFont(42)
latex.SetTextSize(0.05)
latex.DrawLatex(0.40,0.70,'#LT E #GT = %.1f MeV'%histo.GetMean())
latex.Draw('same')
ROOT.gPad.Update()

c2 = ROOT.TCanvas('c2','c2')
c2.cd()
prof = infile.Get('DQMData/Run 1/MTD/Run summary/BTL/LocalReco/BtlCluEnergyVsEta')
prof.GetYaxis().SetTitle('E_{RECO} [MeV]')
prof.Draw()

etaRanges = [0., 0.35, 0.67, 0.94, 1.15, 1.35, 1.6]
funcs = {}
latexs = {}
for it in range(len(etaRanges)-1):
    eta1 = etaRanges[it]
    eta2 = etaRanges[it+1]
    funcs[it+1] = ROOT.TF1('func_RU%d'%(it+1),'pol0',eta1,eta2)
    prof.Fit(funcs[it+1],'QNRS')
    funcs[it+1].SetLineColor(it+1)
    funcs[it+1].Draw('same')
    latexs[it+1] = ROOT.TLatex()
    latexs[it+1].SetNDC()
    latexs[it+1].SetTextFont(42)
    latexs[it+1].SetTextSize(0.05)
    latexs[it+1].SetTextColor(it+1)
    latexs[it+1].DrawLatex(0.40,0.90-0.05*it,'RU%d: %.2f'%(it+1,funcs[it+1].GetParameter(0)/histo.GetMean()))
    latexs[it+1].Draw('same')
ROOT.gPad.Update()

input('ok?')
