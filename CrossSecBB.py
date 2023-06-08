#!/usr/bin/env python
# -*-coding:Latin-1 -*
from math import *

import ROOT
import ROOT as root
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import atlasplots as aplt
import pandas as pd
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TH1
from array import array
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TLatex, TChain

def histtoarray(hist):
    vector = []
    for i in range(0, hist.GetNbinsX()):
        vector.append(hist.GetBinContent(i+1))

    array = np.array(vector)
    return array
def rebinTheMultiJetBackground(nominal, histogram_MJ):

    nominal.SetBinContent(1, 0)
    nominal.SetBinContent(2, histogram_MJ.GetBinContent(26) + histogram_MJ.GetBinContent(27) + histogram_MJ.GetBinContent(28) + histogram_MJ.GetBinContent(29) + histogram_MJ.GetBinContent(30))
    nominal.SetBinContent(3, histogram_MJ.GetBinContent(31) + histogram_MJ.GetBinContent(32) + histogram_MJ.GetBinContent(33) + histogram_MJ.GetBinContent(34) + histogram_MJ.GetBinContent(35))
    nominal.SetBinContent(4, histogram_MJ.GetBinContent(36) + histogram_MJ.GetBinContent(37) + histogram_MJ.GetBinContent(38) + histogram_MJ.GetBinContent(39) + histogram_MJ.GetBinContent(40))
    nominal.SetBinContent(5, histogram_MJ.GetBinContent(41) + histogram_MJ.GetBinContent(42) + histogram_MJ.GetBinContent(43) + histogram_MJ.GetBinContent(44) + histogram_MJ.GetBinContent(45))
    nominal.SetBinContent(6, histogram_MJ.GetBinContent(46) + histogram_MJ.GetBinContent(47) + histogram_MJ.GetBinContent(48) + histogram_MJ.GetBinContent(49) + histogram_MJ.GetBinContent(50))

    sum1 = 0
    sum2 = 0
    sum3 = 0
    for i in range(51, 61):     sum1 = sum1 + histogram_MJ.GetBinContent(i)
    for i in range(61, 81):     sum2 = sum2 + histogram_MJ.GetBinContent(i)
    for i in range(81, 101):    sum3 = sum3 + histogram_MJ.GetBinContent(i)

    nominal.SetBinContent(7, sum1)
    nominal.SetBinContent(8, sum2)
    nominal.SetBinContent(9, sum3)

    nominal.SetBinError(1, 0)
    nominal.SetBinError(2, sqrt(pow(histogram_MJ.GetBinError(26), 2) + pow(histogram_MJ.GetBinError(27), 2) + pow(histogram_MJ.GetBinError(28), 2) + pow(histogram_MJ.GetBinError(29), 2) + pow(histogram_MJ.GetBinError(30), 2)))
    nominal.SetBinError(3, sqrt(pow(histogram_MJ.GetBinError(31), 2) + pow(histogram_MJ.GetBinError(32), 2) + pow(histogram_MJ.GetBinError(33), 2) + pow(histogram_MJ.GetBinError(34), 2) + pow(histogram_MJ.GetBinError(35), 2)))
    nominal.SetBinError(4, sqrt(pow(histogram_MJ.GetBinError(36), 2) + pow(histogram_MJ.GetBinError(37), 2) + pow(histogram_MJ.GetBinError(38), 2) + pow(histogram_MJ.GetBinError(39), 2) + pow(histogram_MJ.GetBinError(40), 2)))
    nominal.SetBinError(5, sqrt(pow(histogram_MJ.GetBinError(41), 2) + pow(histogram_MJ.GetBinError(42), 2) + pow(histogram_MJ.GetBinError(43), 2) + pow(histogram_MJ.GetBinError(44), 2) + pow(histogram_MJ.GetBinError(45), 2)))
    nominal.SetBinError(6, sqrt(pow(histogram_MJ.GetBinError(46), 2) + pow(histogram_MJ.GetBinError(47), 2) + pow(histogram_MJ.GetBinError(48), 2) + pow(histogram_MJ.GetBinError(49), 2) + pow(histogram_MJ.GetBinError(50), 2)))

    errsum1 = 0
    errsum2 = 0
    errsum3 = 0
    for i in range(51, 61):     errsum1 = errsum1 + pow(histogram_MJ.GetBinError(i), 2)
    for i in range(61, 81):     errsum2 = errsum2 + pow(histogram_MJ.GetBinError(i), 2)
    for i in range(81, 101):    errsum3 = errsum3 + pow(histogram_MJ.GetBinError(i), 2)

    nominal.SetBinError(7, sqrt(errsum1))
    nominal.SetBinError(8, sqrt(errsum2))
    nominal.SetBinError(9, sqrt(errsum3))

    return nominal
def CrossSection(  data, Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, lepton, var, multi, Energy, IdUncertai, IsoUncertai, TrigUncertai, recoSFUncertai, LepCalibration, RecoilUncert, NormalisationSyst, MultiJetCalibrea, ttvaUncertai, MuoCalibration):

    HistData   = data.Get(channel+"Selection/Lepton_pT_Reco_cut7")
    HistReco   = Signal.Get(channel+"Selection/Lepton_pT_Reco_cut7")
    HistTruth  = Signal.Get("TruthSelection/Lepton_pT_Truth_cut4")

    if ("enu" in channel)  : lepton = "el" 
    if ("munu" in channel) : lepton = "mu"

    HistBkg1   = Background_Top.Get(channel+"Selection/Lepton_pT_Reco_cut7")
    HistBkg2   = Background_diboson.Get(channel+"Selection/Lepton_pT_Reco_cut7")
    HistBkg3   = Background_W.Get(channel+"Selection/Lepton_pT_Reco_cut7")
    HistBkg4   = Background_Z.Get(channel+"Selection/Lepton_pT_Reco_cut7")
    HistBkg5   = Background_MiltiJet.Get("hist/"+lepton+"Pt")

    # sum of all the backgrounds:
    BkgTotal   = HistBkg1.Integral() + HistBkg2.Integral() + HistBkg3.Integral() + HistBkg4.Integral() + HistBkg5.Integral()

    # calculate the cross sections:
    CrossSecti = ((HistData.Integral() - BkgTotal)/lum)*(HistTruth.Integral()/HistReco.Integral())

    # the statistical uncertainties:
    StatisticalUncert = CrossSecti*sqrt(HistData.Integral()- BkgTotal) / (HistData.Integral() - BkgTotal)

    # the SF and claibration uncertainties
    MuonCalibrationUn   = 0
    ElecCalibrationUn   = 0
    SFIdUncertainties   = 0
    SFTTvaUncertainties = 0

    if ("enu"  in channel) : 
        for i in range(0, len (IdUncertai)):
            SFIdUncertainties      = SFIdUncertainties   + pow(abs(CrossSecti * IdUncertai[i]      / HistReco.Integral()),2)
        for i in range(0, len (LepCalibration)):
            ElecCalibrationUn      = ElecCalibrationUn   + pow(abs(CrossSecti * LepCalibration[i]  / HistReco.Integral()),2)

    if ("munu" in channel) : 
        for i in range(0, len (ttvaUncertai)):
            SFTTvaUncertainties    = SFTTvaUncertainties + pow(abs(CrossSecti * ttvaUncertai[i]    / HistReco.Integral()),2)
        for i in range(0, len (MuoCalibration)):
            MuonCalibrationUn      = MuonCalibrationUn   + pow(abs(CrossSecti * MuoCalibration[i]  / HistReco.Integral()),2)

    MuonCalibrationUn   = sqrt(MuonCalibrationUn)
    ElecCalibrationUn   = sqrt(ElecCalibrationUn)
    SFIdUncertainties   = sqrt(SFIdUncertainties)
    SFTTvaUncertainties = sqrt(SFTTvaUncertainties)

    # define the iso, reco and trig uncertainties
    SFIsoUncertainties   = 0
    SFrecoUncertainties  = 0
    SFTrigUncertainties  = 0

    for i in range(0, len(IsoUncertai)):
        SFIsoUncertainties   = SFIsoUncertainties  + pow(abs(CrossSecti * IsoUncertai[i]    / HistReco.Integral()),2)

    for i in range(0, len(recoSFUncertai)):
        SFrecoUncertainties  = SFrecoUncertainties + pow(abs(CrossSecti * recoSFUncertai[i] / HistReco.Integral()),2)

    for i in range(0, len(TrigUncertai)):
        SFTrigUncertainties  = SFTrigUncertainties + pow(abs(CrossSecti * TrigUncertai[i]   / HistReco.Integral()),2)

    SFIsoUncertainties  = sqrt(SFIsoUncertainties)
    SFrecoUncertainties = sqrt(SFrecoUncertainties)
    SFTrigUncertainties = sqrt(SFTrigUncertainties)

    # the Recoil uncertainties:
    RecoilCalibrationUn  = 0
    for i in range(0, len(RecoilUncert)):
        RecoilCalibrationUn = RecoilCalibrationUn + pow( abs(CrossSecti * RecoilUncert[i] / HistReco.Integral()), 2 )

    RecoilCalibrationUn  = sqrt(RecoilCalibrationUn)

    # the normalisation uncertainties:
    NormSyst    = 0
    for i in range(0, len(NormalisationSyst)):
        NormSyst = NormSyst + pow(CrossSecti * NormalisationSyst[i] / HistReco.Integral(), 2)
    NormSyst = sqrt(NormSyst)

    # the MJ uncertainties:
    MJUncertainties = CrossSecti * MultiJetCalibrea / HistReco.Integral()

    # Luminosity Uncert:
    if(Energy == "5TeV"):
        ErrForLumi = 0.016
    if(Energy == "13TeV"):
        ErrForLumi = 0.015
    LuminUncertain =  CrossSecti * ErrForLumi 


    # the total uncertainty:
    if ("enu"  in channel) : 
        TotalUncertainty = sqrt(pow(SFIdUncertainties, 2)   + pow(SFIsoUncertainties, 2) + pow(SFrecoUncertainties, 2) + pow(SFTrigUncertainties, 2) + pow(ElecCalibrationUn, 2) + pow(RecoilCalibrationUn, 2) + pow(NormSyst, 2) + pow(MJUncertainties, 2))
    if ("munu"  in channel) : 
        TotalUncertainty = sqrt(pow(SFTTvaUncertainties, 2) + pow(SFIsoUncertainties, 2) + pow(SFrecoUncertainties, 2) + pow(SFTrigUncertainties, 2) + pow(MuonCalibrationUn, 2) + pow(RecoilCalibrationUn, 2) + pow(NormSyst, 2) + pow(MJUncertainties, 2))


    print(" The cross section                   : ",                      CrossSecti)
    print(" The statistical uncert              : ",                      StatisticalUncert)
    if ("enu"  in channel) : 
        print(" The ID  SF uncert               : ",                      SFIdUncertainties)
    if ("munu" in channel) : 
        print(" The TTva  SF uncert             : ",                      SFTTvaUncertainties)
    print(" The Iso SF uncert                   : ",                      SFIsoUncertainties)
    print(" The reco SF uncert                  : ",                      SFrecoUncertainties)
    print(" The Trig SF uncert                  : ",                      SFTrigUncertainties)
    if ("enu"  in channel) : 
        print(" The Lepton Calibration uncert   : ",                      ElecCalibrationUn)
    if ("munu"  in channel) : 
        print(" The Lepton Calibration uncert   : ",                      MuonCalibrationUn)
    print(" The Recoil Calibration uncert       : ",                      RecoilCalibrationUn)
    print(" The Normalisation uncert            : ",                      NormSyst)
    print(" The MJ uncert                       : ",                      MJUncertainties)
    print(" The total uncertainties             : ",                      TotalUncertainty)
    print(" The Luminsity uncert                : ",                      LuminUncertain)



    latexFile = open("LatexTables/Fiducial_Section_"+channel+"_"+Energy+".tex","w+")
    latexFile.write("\\documentclass[12pt]{article} \n")
    latexFile.write("\\usepackage{amsmath}\n")
    latexFile.write("\\usepackage{graphicx}\n")
    latexFile.write("\\usepackage{hyperref}\n")
    latexFile.write("\\usepackage{hyperref}\n")
    latexFile.write("\\usepackage[latin1]{inputenc}\n")
    latexFile.write("\\begin{document}\n")

    latexFile.write("\\begin{table}[ht]\n")
    latexFile.write("\\begin{tabular}{c|c|}\n")
    latexFile.write("\\cline{2-2}\n")
    latexFile.write("                                                                   &    %s  \\\ \\hline \\hline \n"%(Indice))
    #latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ (Unfolding)}         &    %5.3f   $\\pm$ %5.2f(Stat) $\\pm$ %5.2f(Syst) $\\pm$ %5.2f(Unf) $\\pm$ %5.2f(Lum)     \\\ \\hline \n" %( Xs,    Stat,    Syst,     Bias, 1.6))
    latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ }        &    %5.3f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Lum)     \\\ \\hline \n" %( CrossSecti, StatisticalUncert, TotalUncertainty,  LuminUncertain))
    latexFile.write("\\end{tabular}\n")
    latexFile.write("\\end{table}\n")
    latexFile.write("\\end{document}\n")
    latexFile.close()

    latexFile = open("LatexTables/BreakDown_Section_"+channel+"_"+Energy+".tex","w+")
    latexFile.write("\\documentclass[12pt]{article} \n")
    latexFile.write("\\usepackage{amsmath}\n")
    latexFile.write("\\usepackage{graphicx}\n")
    latexFile.write("\\usepackage{hyperref}\n")
    latexFile.write("\\usepackage{hyperref}\n")
    latexFile.write("\\usepackage[latin1]{inputenc}\n")
    latexFile.write("\\begin{document}\n")

    latexFile.write("\\begin{table}[ht]\n")
    latexFile.write("\\begin{tabular}{c|c|}\n")
    latexFile.write("\\cline{2-2}\n")
    latexFile.write("                                                                   &    %s  \\\ \\hline \\hline \n"%(Indice))
    #latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ (Unfolding)}         &    %5.3f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Unf) $\\pm$ %5.3f(Lum)     \\\ \\hline \n" %( Xs,    Stat,    Syst,     Bias, 1.6))
    #latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ (Bin-by-Bin)}        &    %5.3f   $\\pm$ %5.3f(Stat) $\\pm$ %5.3f(Syst) $\\pm$ %5.3f(Unf) $\\pm$ %5.3f(Lum)     \\\ \\hline \n" %( CrossSecti, CrossSecti, CrossSecti,  CrossSecti, 1.6))
    latexFile.write("\\multicolumn{1}{|l|}{$\\sigma_{fid}$ }                            &    %5.2f  \\\ \\hline \\hline \n" %( CrossSecti       ))
    latexFile.write("\\multicolumn{1}{|l|}{Data statistical}                &    %5.2f  \\\ \\hline \\hline \n" %( StatisticalUncert))
    if ("enu"  in channel) : 
        latexFile.write("\\multicolumn{1}{|l|}{ID SF}                       &    %5.2f  \\\ \\hline \n" %( SFIdUncertainties))
    if ("munu" in channel) : 
        latexFile.write("\\multicolumn{1}{|l|}{TTva SF}                     &    %5.2f  \\\ \\hline \n" %( SFTTvaUncertainties))
    latexFile.write("\\multicolumn{1}{|l|}{Iso  SF}                         &    %5.2f  \\\ \\hline \n" %( SFIsoUncertainties))
    latexFile.write("\\multicolumn{1}{|l|}{reco SF}                         &    %5.2f  \\\ \\hline \n" %( SFrecoUncertainties))
    latexFile.write("\\multicolumn{1}{|l|}{Trig SF}                         &    %5.2f  \\\ \\hline \\hline \n" %( SFTrigUncertainties))
    if ("enu"  in channel) : 
        latexFile.write("\\multicolumn{1}{|l|}{Lepton Calib.}          &    %5.2f  \\\ \\hline \n" %( ElecCalibrationUn))
    if ("munu" in channel) : 
        latexFile.write("\\multicolumn{1}{|l|}{Lepton Calib.}          &    %5.2f  \\\ \\hline \n" %( MuonCalibrationUn))
    latexFile.write("\\multicolumn{1}{|l|}{Recoil Calib.}              &    %5.2f  \\\ \\hline \\hline \n" %( RecoilCalibrationUn))
    latexFile.write("\\multicolumn{1}{|l|}{MJy}                              &    %5.2f  \\\ \\hline \n" %( MJUncertainties))
    latexFile.write("\\multicolumn{1}{|l|}{Normalisation}                   &    %5.2f  \\\ \\hline \\hline \\hline \n" %( NormSyst))
    latexFile.write("\\multicolumn{1}{|l|}{Total Syst.}                &    %5.2f  \\\ \\hline \n" %( TotalUncertainty))
    latexFile.write("\\end{tabular}\n")
    latexFile.write("\\end{table}\n")
    latexFile.write("\\end{document}\n")
    latexFile.close()
	



def CrossSectionNormalisa(channel, NormalisationPath):

    lepton = "ee"
    if "munu" in channel: lepton = "mumu"
    NameVector = ["ttbar", "singletop", "Wtop", "Wminusenu", "Wplusenu", "Wminusmunu", "Wplusmunu", "Wplustaunu", "Wminustaunu", "Z"+lepton, "Ztautau", "ZZ", "WZ", "WW", "multijet"]

    Variation   = []
    DireName    = channel + "Selection_NormVariations"
    files       = os.listdir(NormalisationPath)
    for i in files:
        rootFile    = ROOT.TFile.Open(NormalisationPath+"/"+i)
        HistNominal = rootFile.Get(DireName+"/Lepton_pT_Reco_cut7_Nom")
        for j in NameVector:
            HistVaried   = rootFile.Get(DireName+"/Lepton_pT_Reco_cut7_Xsec_up_"+j)
            if HistVaried == None: continue
            else: break
        Variation.append(HistVaried.Integral()-HistNominal.Integral())
    
    return Variation
def CrossSectionSF(Signal_SF, channel):

    # HistVariation:
    IsoVersion1 = []
    DireName = channel + "Selection_WeightVariations"
    directory = Signal_SF.GetDirectory(DireName)

    for key in directory.GetListOfKeys():
        hist = key.ReadObj()
        if (hist.ClassName() == 'TH1D' and (hist.GetName()).find("Lepton_pT_Reco_cut7") != -1 and (hist.GetName()).find("up") == -1):
            IsoVersion1.append( hist.Integral())

    return IsoVersion1
def CrossSectionMJ(Background_MiltiJet, channel, Energy, DirectoName):

    Dir_Nom         = "hist"
    Dir_uTSliceErr  = "uTSliceErr_uphist"
    Dir_ExtrapErr   = "ExtrapErr_uphist"
    Dir_shape       = "shape_uphist"

    if "enu" in channel:
        Hist_Nom        = Background_MiltiJet.Get(  Dir_Nom         +   "/elPt" )
        Hist_uTSliceErr = Background_MiltiJet.Get(  Dir_uTSliceErr  +   "/elPt" )
        Hist_ExtrapErr  = Background_MiltiJet.Get(  Dir_ExtrapErr   +   "/elPt" )
        Hist_shape      = Background_MiltiJet.Get(  Dir_shape       +   "/elPt" )
    if "munu" in channel:
        Hist_Nom        = Background_MiltiJet.Get(  Dir_Nom         +   "/muPt" )
        Hist_uTSliceErr = Background_MiltiJet.Get(  Dir_uTSliceErr  +   "/muPt" )
        Hist_ExtrapErr  = Background_MiltiJet.Get(  Dir_ExtrapErr   +   "/muPt" )
        Hist_shape      = Background_MiltiJet.Get(  Dir_shape       +   "/muPt" )

    Error1 = sqrt(Hist_Nom.Integral())
    Error2 = Hist_uTSliceErr.Integral() - Hist_Nom.Integral()
    Error3 = Hist_ExtrapErr.Integral()  - Hist_Nom.Integral()
    Error4 = Hist_shape.Integral()      - Hist_Nom.Integral()

    MJUncer = sqrt( pow(Error1, 2) +  pow(Error2, 2) +  pow(Error3, 2) +  pow(Error4, 2))
    return MJUncer
def CrossSectionMuCalibration(Signal, channel, Energy, DirectoName):

    rec_hist  = Signal.Get(DirectoName+"Selection/Lepton_pT_Reco_cut7")

    CalibSystFiles = []
    CalibSystVariation = []

    CalibSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_'+str(Energy)+'/MuCalibVar/merge/MC_varMUON_ID_1down.root'))
    CalibSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_'+str(Energy)+'/MuCalibVar/merge/MC_varMUON_MS_1down.root'))
    CalibSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_'+str(Energy)+'/MuCalibVar/merge/MC_varMUON_SCALE_1down.root'))
    CalibSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_'+str(Energy)+'/MuCalibVar/merge/MC_varSagittaBiasOffsetDown.root'))
    for i in range(1, 25):
        CalibSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_'+str(Energy)+'/MuCalibVar/merge/mc16_'+str(Energy)+'.varSagittaBiasstatDownbin'+str(i)+'.root'))

    for i in range(0, len(CalibSystFiles)):
        CalibSystVariation.append( (CalibSystFiles[i].Get(DirectoName+"Selection/Lepton_pT_Reco_cut7")).Integral() - rec_hist.Integral() )
    
    return CalibSystVariation
def CrossSectionCalibration(Signal, channel, Energy, DirectoName):

    rec_hist  = Signal.Get(DirectoName+"Selection/Lepton_pT_Reco_cut7")

    CalibSystFiles = []
    CalibSystVariation = []

    for i in range(1, 25):
        CalibSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_'+str(Energy)+'/ElCalibVar/merge/mc16_'+str(Energy)+'.varscaleDownbin'+str(i)+'.root'))
    for i in range(1, 25):
        CalibSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel +'_MC_' +str(Energy)+'/ElCalibVar/merge/mc16_'+str(Energy)+'.varcDownbin'+str(i)+'.root'))

    for i in range(0, len(CalibSystFiles)):
        CalibSystVariation.append( (CalibSystFiles[i].Get(DirectoName+"Selection/Lepton_pT_Reco_cut7")).Integral() - rec_hist.Integral() )
    
    return CalibSystVariation
def CrossSectionRecoil(Signal, channel, Energy, DirectoName):

    rec_hist  = Signal.Get(DirectoName+"Selection/Lepton_pT_Reco_cut7")

    RecoilSystFiles = []
    RecoilSystVariation = []

    if(Energy == "5TeV"):

        RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel+'_MC_5TeV/RecoilVar/merge/mc16_5TeV.varRESPONSE_EXTSYS_DOWNbin1.root'))
        RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel+'_MC_5TeV/RecoilVar/merge/mc16_5TeV.varRESOLUTION_EXTSYS_DOWNbin1.root'))
        RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel+'_MC_5TeV/RecoilVar/merge/mc16_5TeV.varRESPONSE_SYS_DOWNbin1.root'))

        for i in range(1, 20):  # 1, 20
            RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_5TeV/RecoilVar/merge/mc16_5TeV.varRESPONSE_STAT0_DOWNbin'+str(i)+'.root'))
            RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_5TeV/RecoilVar/merge/mc16_5TeV.varRESPONSE_STAT1_DOWNbin'+str(i)+'.root'))  # error

        for i in range(1, 13):  # 1, 13
            RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel +'_MC_5TeV/RecoilVar/merge/mc16_5TeV.varRESOLUTION_STAT0_DOWNbin'+str(i)+'.root'))
            RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel +'_MC_5TeV/RecoilVar/merge/mc16_5TeV.varRESOLUTION_STAT1_DOWNbin'+str(i)+'.root'))

    if(Energy == "13TeV"):

        RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_13TeV/RecoilVar/merge/mc16_13TeV.varRESPONSE_EXTSYS_DOWNbin1.root'))
        RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_13TeV/RecoilVar/merge/mc16_13TeV.varRESOLUTION_EXTSYS_DOWNbin1.root'))
        RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_13TeV/RecoilVar/merge/mc16_13TeV.varRESPONSE_SYS_DOWNbin1.root'))

        for i in range(1, 16):
            RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_13TeV/RecoilVar/merge/mc16_13TeV.varRESPONSE_STAT0_DOWNbin'+str(i)+'.root'))
            RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_13TeV/RecoilVar/merge/mc16_13TeV.varRESPONSE_STAT1_DOWNbin'+str(i)+'.root'))

        for i in range(1, 15):
            RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_13TeV/RecoilVar/merge/mc16_13TeV.varRESOLUTION_STAT0_DOWNbin'+str(i)+'.root'))
            RecoilSystFiles.append(ROOT.TFile('/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w' + channel + '_MC_13TeV/RecoilVar/merge/mc16_13TeV.varRESOLUTION_STAT1_DOWNbin'+str(i)+'.root'))



    for i in range(0, len(RecoilSystFiles)):
        RecoilSystVariation.append( abs((RecoilSystFiles[i].Get(DirectoName+"Selection/Lepton_pT_Reco_cut7")).Integral() - rec_hist.Integral() ))
    
    return RecoilSystVariation
    
	# define signal, data and background
    data_hist = data.Get(channel+"Selection/Lepton_"+var+"_Reco_cut7")
    rec_hist  = Signal.Get(channel+"Selection/Lepton_"+var+"_Reco_cut7")
    tru_hist  = Signal.Get("TruthSelection/Lepton_pT_Truth_cut4")
    bkg_hist1 = Background_Top.Get(channel+"Selection/Lepton_"+var+"_Reco_cut7")
    bkg_hist2 = Background_diboson.Get(channel+"Selection/Lepton_"+var+"_Reco_cut7")
    bkg_hist3 = Background_W.Get(channel+"Selection/Lepton_"+var+"_Reco_cut7")
    bkg_hist4 = Background_Z.Get(channel+"Selection/Lepton_"+var+"_Reco_cut7")

    if "enu" in channel:
        bkg_MJ = (Background_MiltiJet.Get("hist/elPt")).Clone("bkg_MJ")
    if "munu" in channel:
        bkg_MJ = (Background_MiltiJet.Get("hist/muPt")).Clone("bkg_MJ")

    # sum backgroun
    bkg_sum = bkg_hist1.Integral()+bkg_hist2.Integral()+bkg_hist3.Integral()+bkg_hist4.Integral()+bkg_MJ.Integral()

    # cross section
    if Energy == "5TeV":  lum = 256.827
    if Energy == "13TeV": lum = 335.180
    Norma_uncert    = 0
    Cross_sections  = (tru_hist.Integral()/rec_hist.Integral())*(data_hist.Integral() - bkg_sum)/lum
    stat_uncert     = abs(  Cross_sections  * ( sqrt(data_hist.Integral() - bkg_sum)  / (data_hist.Integral() - bkg_sum ) ))

    sumNorm = 0
    for i in range(0, len(NormalisationSyst)):
        Norma_uncert    = Cross_sections * NormalisationSyst[i] / (data_hist.Integral() - bkg_sum)
        sumNorm         = sumNorm + pow(Norma_uncert,2)
    sumNorm             = sqrt(sumNorm)

    Mj_uncert       = abs(  Cross_sections  * ( abs(MultiJetCalibrea)                 / (data_hist.Integral() - bkg_sum ) ))
    Id_uncert       = abs(  Cross_sections  * ( IdUncertai                            / rec_hist.Integral() ) )
    Iso_uncert      = abs(  Cross_sections  * ( IsoUncertai                           / rec_hist.Integral() ) )
    ttva_uncert     = abs(  Cross_sections  * ( ttvaUncertai                          / rec_hist.Integral() ) )
    Trig_uncert     = abs(  Cross_sections  * ( TrigUncertai                          / rec_hist.Integral() ) )
    Reco_uncert     = abs(  Cross_sections  * ( recoSFUncertai                        / rec_hist.Integral() ) )
    Calib_uncert    = abs(  Cross_sections  * ( LepCalibration                        / rec_hist.Integral() ) )
    MuCalib_uncert  = abs(  Cross_sections  * ( MuoCalibration                        / rec_hist.Integral() ) )
    Recoil_uncert   = abs(  Cross_sections  * ( RecoilUncert                          / rec_hist.Integral() ) )

    total           = sqrt(    pow(Norma_uncert,2) + pow(Mj_uncert,2)    + pow(Id_uncert, 2)   + pow(Recoil_uncert, 2)
                           +   pow(Iso_uncert,2)   + pow(Trig_uncert, 2) + pow(Reco_uncert, 2) + pow(Calib_uncert, 2))
 

    print( "Data   : ", data_hist.Integral() )
    print( "Reco   : ", rec_hist.Integral()  )
    print( "Truth  : ", tru_hist.Integral()  ) 
    print( "Bkgr   : ", bkg_sum              )
    print( "Xs     : ", Cross_sections       )
    print( "StatE  : ", stat_uncert          )
    print( "IdSF   : ", Id_uncert            )
    print( "IsoSF  : ", Iso_uncert           )
    print( "ttvaSF : ", ttva_uncert          )
    print( "TriSF  : ", Trig_uncert          )
    print( "RecSF  : ", Reco_uncert          )
    print( "Normal : ", sumNorm              )
    print( "MJ     : ", Mj_uncert            )
    print( "Calib  : ", Calib_uncert         )
    print( "MuCal  : ", MuCalib_uncert       )
    print( "Recoil : ", Recoil_uncert        )
    print( "Total  : ", total                )   

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define the cross sections input:
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if sys.argv[2] == "minusenu"    : MJv =  "em"
if sys.argv[2] == "minusmunu"   : MJv =  "mum"
if sys.argv[2] == "plusenu"     : MJv =  "ep"
if sys.argv[2] == "plusmunu"    : MJv =  "mup"

channel              = "W"+sys.argv[2]
Energy               = sys.argv[3]
if Energy == "5TeV":  lum = 256.827
if Energy == "13TeV": lum = 335.180
Indice               = sys.argv[4]
Indic                = "W^{-}#rightarrow #mu^{-}#nu"

data                = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[2]+"_DATA_"+Energy+"/Nominal/DATA.root")
Signal              = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[2]+"_MC_"+Energy+"/Nominal/MC.root")
Background_Top      = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[2]+"_MC_"+Energy+"/Nominal/merge/Nominal_Top.root")
Background_diboson  = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[2]+"_MC_"+Energy+"/Nominal/merge/Nominal_dilepton.root")
Background_W        = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[2]+"_MC_"+Energy+"/Nominal/merge/Nominal_W.root")
Background_Z        = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[2]+"_MC_"+Energy+"/Nominal/merge/Nominal_Z.root")
#Background_Normal   = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[2]+"_MC_"+Energy+"/Normalization/merge/Normalization.root")
Background_MiltiJet = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/MJbkg/"+Energy+"/W"+MJv+".root")

#define the normalisation systematics
NormalisationPath   = "/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w"+sys.argv[2]+"_MC_"+Energy+"/Normalization/"
NormalisationSyst   = CrossSectionNormalisa(channel, NormalisationPath)

# define the different sources of uncert:
if "enu"  in channel: Lepton ="Elec"
if "munu" in channel: Lepton ="Muon"

Signal_Iso    = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[2] + "_MC_" + sys.argv[3]+"/"+Lepton+"SF_iso/MC.root")
Signal_Trig   = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[2] + "_MC_" + sys.argv[3]+"/"+Lepton+"SF_trig/MC.root")
Signal_recoSF = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[2] + "_MC_" + sys.argv[3]+"/"+Lepton+"SF_reco/MC.root")

# define the SF systematics:
IsoUncertai       = CrossSectionSF(Signal_Iso,   channel)
TrigUncertai      = CrossSectionSF(Signal_Trig,   channel)
recoSFUncertai    = CrossSectionSF(Signal_recoSF, channel)

# define the Recoil systematic:
RecoilCalibrea    = CrossSectionRecoil(Signal, sys.argv[2], Energy, channel)

# define the Multijet systematic:
MultiJetCalibrea  = CrossSectionMJ(Background_MiltiJet, sys.argv[2], Energy, channel)

if "enu" in channel:
    Signal_Id       = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[2] + "_MC_" + sys.argv[3]+"/ElecSF_id/MC.root")
    IdUncertai      = CrossSectionSF(Signal_Id,   channel)
    LepCalibration  = CrossSectionCalibration(Signal, sys.argv[2], Energy, channel)
    ttvaUncertai    = 0
    MuoCalibration  = 0
if "munu" in channel:
    Signal_ttva     = ROOT.TFile.Open("/sps/atlas/h/hatmani/DATA/DATA_Xs/WCrossSections_w" + sys.argv[2] + "_MC_" + sys.argv[3]+"/MuonSF_ttva/MC.root")
    ttvaUncertai    = CrossSectionSF(Signal_ttva,   channel)
    MuoCalibration  = CrossSectionMuCalibration(Signal, sys.argv[2], Energy, channel)
    IdUncertai      = 0
    LepCalibration  = 0

print(" \n *********************************** ", channel, " ******************", Energy, " *********************************** \n")
CrossSection( data, Signal, Background_Top, Background_diboson, Background_W, Background_Z, Background_MiltiJet, Indice, channel, "mu", "pT", "Pt", Energy, IdUncertai, IsoUncertai, TrigUncertai, recoSFUncertai, LepCalibration, RecoilCalibrea, NormalisationSyst, MultiJetCalibrea, ttvaUncertai, MuoCalibration)

print("end")
os._exit(0)
