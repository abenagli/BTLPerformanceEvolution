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
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.SetBatch(False)




SiPMType = 'HPK'
LO = 1150*1.15 #pe/MeV @ 3.5 V
taud = 38.5 # ns
noiseTerm = 420 # uA/ns
fluence = 2.E14
DCRNorm = 17.5
PDEDrop = 0.10


'''
SiPMType = 'HPK'
LO = 1265 #pe/MeV @ 3.5 V
taud = 38.5 # ns
noiseTerm = 294 # uA/ns
fluence = 2.E14
DCRNorm = 29*1.5
PDEDrop = 0.22
'''


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
    print iPar,g_DCR_pars[iPar].Eval(fluence)
    DCRRef_pars[iPar] = g_DCR_pars[iPar].Eval(fluence)
    f_DCRRef_vs_Vov.SetParameter(iPar,DCRRef_pars[iPar])
DCRRef = f_DCRRef_vs_Vov.Eval(1.53);


graph_DCR = ROOT.TGraphErrors()
graph_tRes_tot = ROOT.TGraphErrors()
graph_tRes_stoch = ROOT.TGraphErrors()
graph_tRes_noise = ROOT.TGraphErrors()
graph_tRes_DCR = ROOT.TGraphErrors()

#infile_DCR_TB = ROOT.TFile('logIV_ASIC0_ALDOA_ch5_time_2022-06-12_00:00:01.root')
#Vbr = 38.39
infile_DCR_TB = ROOT.TFile('logIV_ASIC0_ALDOA_ch5_time_2022-06-17_12:31:37.root')
Vbr = 37.58
graph_IV_TB = infile_DCR_TB.Get('g_IV')
graph_DCR_TB = ROOT.TGraph()
for point in range(graph_IV_TB.GetN()):
    x = graph_IV_TB.GetPointX(point)
    y = graph_IV_TB.GetPointY(point)
    gain = (((x-Vbr)-y*1E-06*25)+0.25)*(12.7+3.2)/1.602/0.0001
    print ((x-Vbr),gain,y)
    graph_DCR_TB.SetPoint(point,(x-Vbr)-y*1E-06*25,y*1E-06/16./gain/1.602E-19/1E09)
for Vov in np.arange(0.5, 2.5, 0.01):
    
    nPE = 4.2 * LO * PDE(Vov,SiPMType)*(1.-PDEDrop)/PDE(3.5,'HPK')
    
    #DCR = DCRNorm*f_DCRRef_vs_Vov.Eval(Vov)/DCRRef;
    DCR = graph_DCR_TB.Eval(Vov)
                  
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



c = ROOT.TCanvas('c_tRes_vs_tau','c_tRes_vs_tau',1200,700)
c.Divide(2,1)
c.cd(1)
hPad1 = ROOT.gPad.DrawFrame(0.,0.,2.6,100)
hPad1.SetTitle(";V_{OV} [V];DCR [GHz]")
hPad1.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
graph_DCR.Draw('PL,same')
graph_DCR_TB.SetMarkerColor(ROOT.kRed)
graph_DCR_TB.Draw('PL,same')
c.cd(2)
hPad2 = ROOT.gPad.DrawFrame(1.,20.,2.5,180)
hPad2.SetTitle(";V_{OV} [V];#sigma_{t}^{bar} [ps]")
hPad2.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
graph_tRes_tot.SetLineColor(ROOT.kBlack)
graph_tRes_tot.SetLineWidth(2)
graph_tRes_tot.Draw('L,same')
graph_tRes_stoch.SetLineColor(ROOT.kGreen)
graph_tRes_stoch.SetLineWidth(2)
graph_tRes_stoch.Draw('L,same')
graph_tRes_noise.SetLineColor(ROOT.kBlue)
graph_tRes_noise.SetLineWidth(2)
graph_tRes_noise.Draw('L,same')
graph_tRes_DCR.SetLineColor(ROOT.kOrange)
graph_tRes_DCR.SetLineWidth(2)
graph_tRes_DCR.Draw('L,same')

#infile_TB = ROOT.TFile('plots_timeResolution_2E14_mipPeak.root')
infile_TB = ROOT.TFile('HPK_2E14_LYSO796_T-40C_summary.root')
#graph_TB = infile_TB.Get('g_Data_vs_Vov_HPK_2E14_T-40C_bar10')
graph_TB = infile_TB.Get('g_deltaT_energyRatioCorr_vs_vov_bar07_th07_enBin01')
graph_TB.SetMarkerStyle(20)
graph_TB.SetMarkerColor(ROOT.kBlack)
graph_TB.SetMarkerSize(1)
graph_TB.Draw('P,same')

'''graph_TB_noise = infile_TB.Get('g_Noise_vs_Vov_HPK_2E14_T-40C_bar10')
graph_TB_noise.SetMarkerStyle(20)
graph_TB_noise.SetMarkerColor(ROOT.kBlue)
graph_TB_noise.SetMarkerSize(1)
graph_TB_noise.Draw('P,same')

graph_TB_stoch = infile_TB.Get('g_Stoch_vs_Vov_HPK_2E14_T-40C_bar10')
graph_TB_stoch.SetMarkerStyle(20)
graph_TB_stoch.SetMarkerColor(ROOT.kGreen)
graph_TB_stoch.SetMarkerSize(1)
graph_TB_stoch.Draw('P,same')

graph_TB_DCR = infile_TB.Get('g_DCR_vs_Vov_HPK_2E14_T-40C_bar10')
graph_TB_DCR.SetMarkerStyle(20)
graph_TB_DCR.SetMarkerColor(ROOT.kOrange)
graph_TB_DCR.SetMarkerSize(1)
graph_TB_DCR.Draw('P,same')
'''
c.Print('c_tRes_vs_Vov.png')

outfile = ROOT.TFile('f_tRes_vs_Vov.root','RECREATE')
graph_TB.Write('graph_TB')
#graph_TB_noise.Write('graph_TB_noise')
#graph_TB_stoch.Write('graph_TB_stoch')
#graph_TB_DCR.Write('graph_TB_DCR')
graph_tRes_tot.Write('graph_tRes_tot')
graph_tRes_noise.Write('graph_tRes_noise')
graph_tRes_DCR.Write('graph_tRes_DCR')
graph_tRes_stoch.Write('graph_tRes_stoch')
outfile.Close()
