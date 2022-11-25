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
infilenames['20 #mum cell size'] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__20um_final.root'
infilenames['30 #mum cell size'] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__30um_final.root'
infilenames['15 #mum cell size'] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__15um_final.root'
infilenames['25 #mum cell size'] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__25um_final.root'

infilenames_up = OrderedDict()
infilenames_up['15 #mum cell size'] = '../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__15um_final.root'
infilenames_up['20 #mum cell size'] = '../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__20um_final.root'
infilenames_up['25 #mum cell size'] = '../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__25um_final.root'
infilenames_up['30 #mum cell size'] = '../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__30um_final.root'

infilenames_down = OrderedDict()
infilenames_down['15 #mum cell size'] = '../plots/outFile__LO1300_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower800__15um_final.root'
infilenames_down['20 #mum cell size'] = '../plots/outFile__LO1300_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower800__20um_final.root'
infilenames_down['25 #mum cell size'] = '../plots/outFile__LO1300_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower800__25um_final.root'
infilenames_down['30 #mum cell size'] = '../plots/outFile__LO1300_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower800__30um_final.root'

colors = {}
colors['15 #mum cell size'] = 51
colors['20 #mum cell size'] = 51+15
colors['25 #mum cell size'] = 51+15+15
colors['30 #mum cell size'] = 51+15+15+15

labels = infilenames.keys()

infiles = {}
infiles_up = {}
infiles_down = {}
for label in labels:
    infiles[label] = ROOT.TFile(infilenames[label])
    infiles_up[label] = ROOT.TFile(infilenames_up[label])
    infiles_down[label] = ROOT.TFile(infilenames_down[label])        




Vov_limit = {}
for label in labels:
    
    c_it = ROOT.TCanvas('c_power_vs_Vov_%s'%label,'c_power_vs_Vov_%s'%label,2400,1400)
    hPad_it = ROOT.gPad.DrawFrame(0.,0.,2.,100.)
    hPad_it.SetTitle(";V_{OV} [V];power/ch [mW]")
    hPad_it.Draw()
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    
    graph = infiles[label].Get('g_power_vs_Vov_EoO')
    graph.SetLineColor(ROOT.kBlack)
    graph.SetLineWidth(3)
    graph.Draw('L,same')
    
    for point in range(0,graph.GetN()):
        if graph.GetPointY(point) > 80.:
            Vov_limit[label] = graph.GetPointX(point)
            break
    
    graph_TECs = infiles[label].Get('g_TECsPower_vs_Vov_EoO')
    graph_TECs.SetLineColor(ROOT.kGreen+1)
    graph_TECs.SetLineWidth(3)
    graph_TECs.Draw('L,same')
    
    graph_dynamic = infiles[label].Get('g_dynamicPower_vs_Vov_EoO')
    graph_dynamic.SetLineColor(ROOT.kBlue-4)
    graph_dynamic.SetLineWidth(3)
    graph_dynamic.Draw('L,same')

    graph_static = infiles[label].Get('g_staticPower_vs_Vov_EoO')
    graph_static.SetLineColor(ROOT.kOrange-3)
    graph_static.SetLineWidth(3)
    graph_static.Draw('PL,same')

    graph_SiPM = ROOT.TGraph()
    for point in range(0,graph_static.GetN()):
        graph_SiPM.SetPoint(point,graph_static.GetPointX(point),graph_static.GetPointY(point)+graph_dynamic.GetPointY(point))
    graph_SiPM.SetLineColor(ROOT.kRed+1)
    graph_SiPM.SetLineWidth(3)
    graph_SiPM.Draw('L,same')
    
    leg_it = ROOT.TLegend(0.15, 0.94-0.035*(len(labels)+1), 0.95, 0.94)
    leg_it.SetFillStyle(0)
    leg_it.SetTextFont(42)
    leg_it.SetTextSize(0.03)
    leg_it.AddEntry(graph,'total','L')
    leg_it.AddEntry(graph_TECs,'TECs','L')
    leg_it.AddEntry(graph_dynamic,'dynamic','L')
    leg_it.AddEntry(graph_static,'static','L')
    leg_it.AddEntry(graph_SiPM,'SiPM total','L')
    leg_it.Draw('same')
    
    latex_it = ROOT.TLatex()
    latex_it.SetNDC()
    latex_it.SetTextFont(42)
    latex_it.SetTextSize(0.04)
    latex_it.DrawLatex(0.16,0.96,label)
    
    c_it.Print('c_power_vs_Vov_model_%s.png'%label)



c = ROOT.TCanvas('c_tRes_vs_Vov','c_tRes_vs_Vov',2000,1400)
hPad = ROOT.gPad.DrawFrame(0.2,0.,1.8,150.)
hPad.SetTitle(";V_{OV} [V];#sigma_{t} [ps]")
hPad.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()

leg = ROOT.TLegend(0.15, 0.34-0.035*len(labels), 0.95, 0.34)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.04)

graphsErr = {}
lines = {}
it = 0
for label in labels:
    graph = infiles[label].Get('g_tRes_vs_Vov_EoO')
    graph_up = infiles_up[label].Get('g_tRes_vs_Vov_EoO')
    graph_down = infiles_down[label].Get('g_tRes_vs_Vov_EoO')
    graphsErr[label] = ROOT.TGraphErrors()
    #pointOpt = -1
    #tResOpt = 9999.
    #for point in range(0,graph.GetN()):
    #    if graph.GetPointY(point) < tResOpt:
    #        tResOpt = graph.GetPointY(point)
    #        pointOpt = point
    for point in range(0,graph.GetN()):
        x = graph.GetPointX(point)
        y_up = graph_up.GetPointY(point)
        y_down = graph_down.GetPointY(point)
        #y_up = graph_up.GetPointY(pointOpt)
        #y_down = graph_down.GetPointY(pointOpt)
        graphsErr[label].SetPoint(point,x,0.5*(y_up+y_down))
        graphsErr[label].SetPointError(point,0,0.5*(y_up-y_down))
        graph.SetPoint(point,x,0.5*(y_up+y_down))
    lines[label] = ROOT.TLine(Vov_limit[label],0.,Vov_limit[label],graph.Eval(Vov_limit[label]))
    lines[label].SetLineStyle(10)
    lines[label].SetLineWidth(5)
    
    c_it = ROOT.TCanvas('c_tRes_vs_Vov_%s'%label,'c_tRes_vs_Vov_%s'%label,2400,1400)
    hPad_it = ROOT.gPad.DrawFrame(0.2,0.,1.8,150.)
    hPad_it.SetTitle(";V_{OV} [V];#sigma_{t} [ps]")
    hPad_it.Draw()
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    
    graph.SetLineColor(ROOT.kBlack)
    graph.SetLineWidth(6)
    graph.Draw('L,same')
    
    graph_stoch = infiles[label].Get('g_tRes_stoch_vs_Vov_EoO')
    graph_stoch.SetLineColor(ROOT.kGreen+1)
    graph_stoch.SetLineStyle(7)
    graph_stoch.SetLineWidth(6)
    graph_stoch.Draw('L,same')
    
    graph_noise = infiles[label].Get('g_tRes_noise_vs_Vov_EoO')
    graph_noise.SetLineColor(ROOT.kBlue-4)
    graph_noise.SetLineStyle(7)
    graph_noise.SetLineWidth(6)
    graph_noise.Draw('L,same')
    
    graph_DCR = infiles[label].Get('g_tRes_DCR_vs_Vov_EoO')
    graph_DCR.SetLineColor(ROOT.kOrange-3)
    graph_DCR.SetLineStyle(7)
    graph_DCR.SetLineWidth(6)
    graph_DCR.Draw('L,same')
    
    leg_it = ROOT.TLegend(0.60, 0.94-0.045*len(labels), 0.80, 0.94)
    leg_it.SetFillStyle(1001)
    leg_it.SetTextFont(42)
    leg_it.SetTextSize(0.05)
    leg_it.AddEntry(graph,'total','L')
    leg_it.AddEntry(graph_stoch,'stochastic','L')
    leg_it.AddEntry(graph_noise,'elec. noise','L')
    leg_it.AddEntry(graph_DCR,'DCR','L')
    leg_it.Draw('same')
    
    lines[label].SetLineColor(ROOT.kBlack)
    lines[label].Draw('same')
    
    latex_it = ROOT.TLatex()
    latex_it.SetNDC()
    latex_it.SetTextFont(42)
    latex_it.SetTextSize(0.05)
    latex_it.DrawLatex(0.16,0.96,label)
    
    latex_line = ROOT.TLatex()
    latex_line.SetTextFont(42)
    latex_line.SetTextSize(0.04)
    latex_line.DrawLatex(Vov_limit[label]+0.025,10,'80 mW/ch limit')
    
    c_it.Print('c_tRes_vs_Vov_model_%s.png'%label)
    
    c.cd()
    
    graphsErr[label].SetFillStyle(3002)
    graphsErr[label].SetFillColor(colors[label])
    graphsErr[label].SetMarkerColor(colors[label])
    graphsErr[label].SetMarkerStyle(20)
    graph.SetMarkerColor(colors[label])
    graph.SetLineColor(colors[label])
    SetGraphStyle(graph,label)
    graph.Draw('PL,same')
    if label == '25 #mum cell size':
        graphsErr[label].Draw('E3,same')
    lines[label].SetLineColor(colors[label])
    lines[label].Draw('same')

    #if label == '15 #mum cell size':
    #    leg.AddEntry(graphsErr[label],'%s - MS SiPMs'%label,'PF')
    #else:
    #    leg.AddEntry(graphsErr[label],'%s'%label,'PF')
        
    it += 1

leg.AddEntry(graphsErr['15 #mum cell size'],'%s - MS SiPMs'%'15 #mum cell size','PF')
leg.AddEntry(graphsErr['20 #mum cell size'],'%s'%'20 #mum cell size','PF')
leg.AddEntry(graphsErr['25 #mum cell size'],'%s'%'25 #mum cell size','PF')
leg.AddEntry(graphsErr['30 #mum cell size'],'%s'%'30 #mum cell size','PF')

leg.Draw('same')

infile_TB = ROOT.TFile('f_tRes_vs_Vov.root')
graph_TB = infile_TB.Get('graph_TB')
graph_TB.SetMarkerStyle(20)
graph_TB.SetMarkerSize(1.5)
graph_TB.Draw('P,same')
graph_model = infile_TB.Get('graph_tRes_tot')
graph_model.SetLineWidth(3)
graph_model.Draw('L,same')

c.cd()


c.Print('c_tRes_vs_Vov_model.png')
