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



def SetGraphStyle(graph, label):
    if label == '15 #mum cell size':
        graph.SetMarkerStyle(20)
    if label == '20 #mum cell size':
        graph.SetMarkerStyle(21)        
    if label == '25 #mum cell size':
        graph.SetMarkerStyle(22)
    if label == '30 #mum cell size':
        graph.SetMarkerStyle(23)

infilenames = OrderedDict()
##infilenames['15 #mum cell size'] = '../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.00_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower65_15um_new5.root'
##infilenames['20 #mum cell size'] = '../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.00_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower65_20um_new5.root'
##infilenames['25 #mum cell size'] = '../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.00_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower65_25um_new5.root'
##infilenames['30 #mum cell size'] = '../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.15_dropGain0.00_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower65_30um_new5.root'
#infilenames['15 #mum cell size'] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower72__15um_final.root'
#infilenames['20 #mum cell size'] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower72__20um_final.root'
infilenames['25 #mum cell size'] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower72__25um_final.root'
#infilenames['30 #mum cell size'] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower72__30um_final.root'

labels = infilenames.keys()

infiles = {}
for label in labels:
    infiles[label] = ROOT.TFile(infilenames[label])



c = ROOT.TCanvas('c_tRes_vs_intLumi','c_tRes_vs_intLumi',1200,700)
#hPad = ROOT.gPad.DrawFrame(0.,0.,12.,90.)
#hPad.SetTitle(";years from 2027;#sigma_{t}^{best} [ps]")
hPad = ROOT.gPad.DrawFrame(0.,0.,3000.,90.)
hPad.SetTitle(";integrated luminosity [fb^{-1}];#sigma_{t}^{best} [ps]")
hPad.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()

leg = ROOT.TLegend(0.17, 0.84-0.035*len(labels)*4, 0.50, 0.84)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.04)

it = 0
for label in labels:
    graph = infiles[label].Get('g_tResBest_vs_intLumi')
    graph_stoch = infiles[label].Get('g_tResBest_stoch_vs_intLumi')
    graph_noise = infiles[label].Get('g_tResBest_noise_vs_intLumi')
    graph_DCR   = infiles[label].Get('g_tResBest_DCR_vs_intLumi')
    print(label)
    c.cd()
    
    #graph.SetMarkerColor(51+it*15)
    #graph.SetLineColor(51+it*15)
    graph.SetMarkerColor(ROOT.kGreen+1)
    graph.SetLineColor(ROOT.kGreen+1)
    graph.SetLineWidth(3)
    #SetGraphStyle(graph,label)
    graph.Draw('PL,same')
    leg.AddEntry(graph,'%s'%'total intLumi resolution','PL')
    
    graph_stoch.SetMarkerColor(ROOT.kYellow+1)
    graph_stoch.SetLineColor(ROOT.kYellow+1)
    graph_stoch.SetLineWidth(3)
    #SetGraphStyle(graph_stoch,label)
    graph_stoch.Draw('PL,same')
    leg.AddEntry(graph_stoch,'%s'%'photostatistics','PL')

    graph_DCR.SetMarkerColor(ROOT.kRed+1)
    graph_DCR.SetLineColor(ROOT.kRed+1)
    graph_DCR.SetLineWidth(3)
    #SetGraphStyle(graph_DCR,label)
    graph_DCR.Draw('PL,same')
    leg.AddEntry(graph_DCR,'%s'%'DCR','PL')

    graph_noise.SetMarkerColor(ROOT.kAzure+1)
    graph_noise.SetLineColor(ROOT.kAzure+1)
    graph_noise.SetLineWidth(3)
    #SetGraphStyle(graph_noise,label)
    graph_noise.Draw('PL,same')
    leg.AddEntry(graph_noise,'%s'%'electronics noise','PL')
    
    it += 1

leg.Draw('same')
c.Print('c_tRes_vs_intLumi_model.png')



c = ROOT.TCanvas('c_power_vs_time','c_power_vs_time',800,700)
hPad = ROOT.gPad.DrawFrame(0.,0.,12.,90.)
hPad.SetTitle(";years from 2027;total power per ch. [mW]")
hPad.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()

leg = ROOT.TLegend(0.15, 0.34-0.035*len(labels), 0.95, 0.34)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.03)

it = 0
for label in labels:
    graph = infiles[label].Get('g_powerBest_vs_time')
    c.cd()
    
    graph.SetMarkerColor(51+it*15)
    graph.SetLineColor(51+it*15)
    SetGraphStyle(graph,label)
    graph.Draw('PL,same')

    leg.AddEntry(graph,'%s'%label,'PL')
    
    it += 1

leg.Draw('same')
c.Print('c_power_vs_time_model.png')



c = ROOT.TCanvas('c_SiPMPower_vs_time','c_SiPMPower_vs_time',800,700)
hPad = ROOT.gPad.DrawFrame(0.,0.,12.,90.)
hPad.SetTitle(";years from 2027;SiPM power per ch. [mW]")
hPad.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()

leg = ROOT.TLegend(0.15, 0.34-0.035*len(labels), 0.95, 0.34)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.03)

graphs = {}
it = 0
for label in labels:
    graph1 = infiles[label].Get('g_dynamicPowerBest_vs_time')
    graph2 = infiles[label].Get('g_staticPowerBest_vs_time')
    c.cd()
    
    graphs[label] = ROOT.TGraph()
    for point in range(0,graph1.GetN()):
        graphs[label].SetPoint(point,graph1.GetPointX(point),graph1.GetPointY(point)+graph2.GetPointY(point))
    
    graphs[label].SetMarkerColor(51+it*15)
    graphs[label].SetLineColor(51+it*15)
    SetGraphStyle(graphs[label],label)
    graphs[label].Draw('PL,same')

    leg.AddEntry(graphs[label],'%s'%label,'PL')
    
    it += 1

leg.Draw('same')
c.Print('c_SiPMPower_vs_time_model.png')
