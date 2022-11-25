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



def myfunc_amp(x, par):
    xx = x[0]
    x0 = par[0]
    
    if xx < x0:
        return (par[2]*xx + par[1])
    else:
        return (par[3]*math.log(xx) + par[2]*x0 - par[3]*log(x0))


def myfunc_DCR(x, par):
    xx = x[0]
    
    turnOff = 0.5 * ( 1. + ROOT.TMath.Erf(-1.*(xx-par[0])/par[1]) )
    turnOn  = 0.5 * ( 1. + ROOT.TMath.Erf(+1.*(xx-par[0])/par[1]) )
    return math.exp(par[2]+par[3]*xx) + turnOff*math.exp(par[4]+par[5]*xx) + turnOn*math.exp(par[6]+par[7]*xx)


def PDE(ov, sipm):
    if (sipm=='HPK'):
        return 1.0228*0.384 * ( 1. - math.exp(-1.*0.583*ov) ) # 1.0228 factor to account for LYSO emission spectrum
    if (sipm=='FBK'):
        return 0.8847*0.466 * ( 1. - math.exp(-1.*0.314*ov) ) # 0.8847 factor to account for LYSO emission spectrum


def Gain(ov, sipm):
    if (sipm=='HPK'):
        return (97602.*(ov+0.377962))

    

def SetGraphStyle(graph, label):
    if label == '15 #mum cell size':
        graph.SetMarkerStyle(24)
    if label == '20 #mum cell size':
        graph.SetMarkerStyle(24)        
    if label == '25 #mum cell size':
        graph.SetMarkerStyle(24)
    if label == '30 #mum cell size':
        graph.SetMarkerStyle(24)

infilenames = OrderedDict()
infilenames['15 #mum cell size'] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__15um_final.root'
#infilenames['20 #mum cell size'] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__20um_final.root'
#infilenames['25 #mum cell size'] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__25um_final.root'
#infilenames['30 #mum cell size'] = '../plots/outFile__LO1265_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__30um_final.root'

infilenames_up = OrderedDict()
infilenames_up['15 #mum cell size'] = '../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__15um_final_PDEslow.root'
#infilenames_up['20 #mum cell size'] = '../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__20um_final_PDEslow.root'
#infilenames_up['25 #mum cell size'] = '../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__25um_final_PDEslow.root'
#infilenames_up['30 #mum cell size'] = '../plots/outFile__LO1200_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise335__Top_-45_Tann1_15_Tann2_40_interfillAnnealing__HLLHCSchedule81__maxPower800__30um_final_PDEslow.root'

infilenames_down = OrderedDict()
infilenames_down['15 #mum cell size'] = '../plots/outFile__LO1300_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower800__15um_final.root'
#infilenames_down['20 #mum cell size'] = '../plots/outFile__LO1300_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower800__20um_final.root'
#infilenames_down['25 #mum cell size'] = '../plots/outFile__LO1300_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower800__25um_final.root'
#infilenames_down['30 #mum cell size'] = '../plots/outFile__LO1300_tau38.5__HPK_dropPDE0.22_dropGain0.08_TECs__noise294__Top_-45_Tann1_15_Tann2_60_interfillAnnealing__HLLHCSchedule81__maxPower800__30um_final.root'

labels = infilenames.keys()

infiles = {}
infiles_up = {}
infiles_down = {}
for label in labels:
    infiles[label] = ROOT.TFile(infilenames[label])
    infiles_up[label] = ROOT.TFile(infilenames_up[label])
    infiles_down[label] = ROOT.TFile(infilenames_down[label])        



SiPMType = 'HPK'
LO = 1150 #pe/MeV @ 3.5 V
taud = 38.5 # ns
noiseTerm = 420 # uA/ns
fluence = 2.E14
DCRNorm = 29
PDEDrop = 0.22

f_slewRate_vs_amp = ROOT.TF1("f_slewRate_vs_amp",myfunc_amp,0.,10.,4);
f_slewRate_vs_amp.SetParameters(5.32470e-01,0.,2.92152e+01,7.79368e+00);

nParDCR = 8
inFile_DCRParams = ROOT.TFile.Open("../data/DCRParams_new_Dic2021.root","READ");
g_DCR_pars = {}
for iPar in range(0,nParDCR):
    g_DCR_pars[iPar] = inFile_DCRParams.Get('g_par%d'%iPar)
    
f_DCRRef_vs_Vov = ROOT.TF1("f_DCRRef_vs_Vov",myfunc_DCR,0.,7.,nParDCR)
DCRRef_pars = {}
for iPar in range(0,nParDCR):
    print(iPar,g_DCR_pars[iPar].Eval(fluence))
    DCRRef_pars[iPar] = g_DCR_pars[iPar].Eval(fluence)
    f_DCRRef_vs_Vov.SetParameter(iPar,DCRRef_pars[iPar])
DCRRef = f_DCRRef_vs_Vov.Eval(1.5);


graph_DCR = ROOT.TGraphErrors()
graph_tRes_tot = ROOT.TGraphErrors()
graph_tRes_stoch = ROOT.TGraphErrors()
graph_tRes_noise = ROOT.TGraphErrors()
graph_tRes_DCR = ROOT.TGraphErrors()

for Vov in np.arange(0.5, 2.5, 0.01):
    
    nPE = 4.2 * LO * PDE(Vov,SiPMType)*(1.-PDEDrop)/PDE(3.5,'HPK')
    
    DCR = DCRNorm*f_DCRRef_vs_Vov.Eval(Vov)/DCRRef;
                  
    sigma_stoch = 27. * math.sqrt(7000./(nPE*38.5/taud))
    sigma_noise = math.sqrt( pow(noiseTerm/1.2/f_slewRate_vs_amp.Eval(nPE*Gain(Vov,SiPMType)/Gain(3.5,SiPMType)/9500.),2) + pow(16.7,2) )/math.sqrt(2)
    sigma_DCR   = 40. * 6000./(nPE*38.5/taud) * pow(DCR/30.,0.5)
    sigma_clock = 15.
    
    tResCurr = math.sqrt( math.pow(sigma_stoch,2) + math.pow(sigma_noise,2) + math.pow(sigma_DCR,2) + math.pow(sigma_clock,2) );
    
    graph_DCR.SetPoint(graph_DCR.GetN(),Vov,DCR)
    graph_tRes_tot.SetPoint(graph_tRes_tot.GetN(),Vov,tResCurr)
    graph_tRes_stoch.SetPoint(graph_tRes_stoch.GetN(),Vov,sigma_stoch)
    graph_tRes_noise.SetPoint(graph_tRes_noise.GetN(),Vov,sigma_noise)
    graph_tRes_DCR.SetPoint(graph_tRes_DCR.GetN(),Vov,sigma_DCR)


Vov_limit = {}
for label in labels:
    
    c_it = ROOT.TCanvas('c_power_vs_Vov_%s'%label,'c_power_vs_Vov_%s'%label,2400,1400)
    hPad_it = ROOT.gPad.DrawFrame(0.8,0.,2.,100.)
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
    leg_it.SetTextSize(0.04)
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
hPad = ROOT.gPad.DrawFrame(0.8,0.,1.8,150.)
hPad.SetTitle(";V_{OV} [V];#sigma_{t} [ps]")
hPad.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()

leg = ROOT.TLegend(0.18, 0.34-0.035*(len(labels)+4), 0.38, 0.34)
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
    for point in range(0,graph_up.GetN()):
        x = graph_up.GetPointX(point)
        y_up = graph_up.GetPointY(point)
        y_down = graph_down.GetPointY(point)
        graphsErr[label].SetPoint(point,x,0.5*(y_up+y_down))
        graphsErr[label].SetPointError(point,0,0.5*(y_up-y_down))
        #graph.SetPoint(point,x,0.5*(y_up+y_down))
    lines[label] = ROOT.TLine(Vov_limit[label],0.,Vov_limit[label],graph.Eval(Vov_limit[label]))
    lines[label].SetLineStyle(10)
    lines[label].SetLineWidth(3)
    
    c_it = ROOT.TCanvas('c_tRes_vs_Vov_%s'%label,'c_tRes_vs_Vov_%s'%label,2400,1400)
    hPad_it = ROOT.gPad.DrawFrame(0.8,0.,2.,150.)
    hPad_it.SetTitle(";V_{OV} [V];#sigma_{t} [ps]")
    hPad_it.Draw()
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    
    graph.SetLineColor(ROOT.kBlack)
    graph.SetLineWidth(3)
    graph.Draw('L,same')
    
    graph_stoch = infiles[label].Get('g_tRes_stoch_vs_Vov_EoO')
    graph_stoch.SetLineColor(ROOT.kGreen+1)
    graph_stoch.SetLineStyle(7)
    graph_stoch.SetLineWidth(3)
    graph_stoch.Draw('L,same')
    
    graph_noise = infiles[label].Get('g_tRes_noise_vs_Vov_EoO')
    graph_noise.SetLineColor(ROOT.kBlue-4)
    graph_noise.SetLineStyle(7)
    graph_noise.SetLineWidth(3)
    graph_noise.Draw('L,same')
    
    graph_DCR = infiles[label].Get('g_tRes_DCR_vs_Vov_EoO')
    graph_DCR.SetLineColor(ROOT.kOrange-3)
    graph_DCR.SetLineStyle(7)
    graph_DCR.SetLineWidth(3)
    graph_DCR.Draw('L,same')

    leg_it = ROOT.TLegend(0.15, 0.34-0.035*len(labels), 0.95, 0.34)
    leg_it.SetFillStyle(0)
    leg_it.SetTextFont(42)
    leg_it.SetTextSize(0.03)
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
    latex_it.SetTextSize(0.04)
    latex_it.DrawLatex(0.16,0.96,label)
    
    latex_line = ROOT.TLatex()
    latex_line.SetTextFont(42)
    latex_line.SetTextSize(0.04)
    latex_line.DrawLatex(Vov_limit[label]+0.025,10,'80 mW/ch limit')
    
    c_it.Print('c_tRes_vs_Vov_model_%s.png'%label)
    
    c.cd()
    
    graphsErr[label].SetFillStyle(3002)
    graphsErr[label].SetMarkerStyle(24)
    graphsErr[label].SetMarkerColor(ROOT.kGray)
    graphsErr[label].SetFillColor(ROOT.kGray)
    graph.SetMarkerColor(ROOT.kBlack)
    graph.SetMarkerStyle(24)
    graph.SetLineColor(ROOT.kGray)
    SetGraphStyle(graph,label)
    graph.Draw('PL,same')
    graphsErr[label].Draw('E3,same')
    lines[label].SetLineColor(ROOT.kGray)
    #lines[label].Draw('same')
    
    leg.AddEntry(graphsErr[label],'BTL EoO, TOFHIR2B 30% noise reduction, 15% LO gain','PF')
    
    it += 1


c.cd()

graph_tRes_tot.SetLineColor(ROOT.kBlack)
graph_tRes_tot.SetLineWidth(3)
graph_tRes_tot.SetLineStyle(1)
graph_tRes_tot.Draw('L,same')
graph_tRes_stoch.SetLineColor(ROOT.kGreen)
graph_tRes_stoch.SetLineWidth(3)
graph_tRes_stoch.SetLineStyle(7)
graph_tRes_stoch.Draw('L,same')
graph_tRes_noise.SetLineColor(ROOT.kBlue)
graph_tRes_noise.SetLineWidth(3)
graph_tRes_noise.SetLineStyle(7)
graph_tRes_noise.Draw('L,same')
graph_tRes_DCR.SetLineColor(ROOT.kOrange)
graph_tRes_DCR.SetLineWidth(3)
graph_tRes_DCR.SetLineStyle(7)
graph_tRes_DCR.Draw('L,same')

infile_TB = ROOT.TFile('plots_timeResolution_2E14_mipPeak.root')
graph_TB = infile_TB.Get('g_Data_vs_Vov_HPK_2E14_T-40C_bar10')
graph_TB.SetMarkerStyle(20)
graph_TB.SetMarkerSize(2)
graph_TB.Draw('P,same')

graph_TB_noise = infile_TB.Get('g_Noise_vs_Vov_HPK_2E14_T-40C_bar10')
graph_TB_noise.SetMarkerStyle(20)
graph_TB_noise.SetMarkerColor(ROOT.kBlue)
graph_TB_noise.SetMarkerSize(2)
graph_TB_noise.Draw('P,same')

graph_TB_stoch = infile_TB.Get('g_Stoch_vs_Vov_HPK_2E14_T-40C_bar10')
graph_TB_stoch.SetMarkerStyle(20)
graph_TB_stoch.SetMarkerColor(ROOT.kGreen)
graph_TB_stoch.SetMarkerSize(2)
graph_TB_stoch.Draw('P,same')

graph_TB_DCR = infile_TB.Get('g_DCR_vs_Vov_HPK_2E14_T-40C_bar10')
graph_TB_DCR.SetMarkerStyle(20)
graph_TB_DCR.SetMarkerColor(ROOT.kOrange)
graph_TB_DCR.SetMarkerSize(2)
graph_TB_DCR.Draw('P,same')

leg.AddEntry(graph_TB,"TB data (irradiated modules)")
leg.AddEntry(graph_TB_noise,"TB data - elec. noise")
leg.AddEntry(graph_TB_DCR,"TB data  - DCR")
leg.AddEntry(graph_TB_stoch,"TB data - stochastic")
leg.Draw('same')

c.Print('c_money_plot.png')
