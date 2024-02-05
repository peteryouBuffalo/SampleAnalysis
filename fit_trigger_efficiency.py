## == IMPORT STATEMENTS =======================================================

import ROOT
import sys,os
import copy
import math
import ConfigParser
import array
from math import *

from my_funcs import makePlot

## == Useful Methods ==========================================================

##############################################
## Scale a histogram to the proper luminosity
##############################################
def scaleToLumi1(fName, xSec, lumi):
  f = ROOT.TFile.Open(fName, 'read')
  hTmp = f.Get('Nevt')
  nP = hTmp.GetBinContent(3)
  nN = hTmp.GetBinContent(1)
  return lumi*xSec/(nP-nN)

####################################################
## Get the histograms for a given sample & variable
####################################################
def getHist(pN, sample_name, fH, lS, printSamples=True):
  hOut = {}
  
  ## Go through each year we're interested in.
  for y in years:
    ## Get the first sample
    if printSamples: print sample_name[0], pN, y
    hOut[y] = fH[sample_name[0]][y][0].Get(pN).Clone()
    if sample_name[0] not in ['JetHT', 'SingleMuon']:
      hOut[y].Scale(lS[sample_name[0]][y][0])
    
    ## Add the other samples 
    for iS in range(len(sample_name)):
      for fi in range(len(fH[sample_name[iS]][y])):
        if iS == 0 and fi == 0: continue
        h = fH[sample_name[iS]][y][fi].Get(pN).Clone()
        if sample_name[iS] not in ['JetHT', 'SingleMuon']:
          h.Scale(lS[sample_name[iS]][y][fi])
        hOut[y].Add(h)
    
  if printSamples: 
    print "==============================================="
    
  return hOut

###############################################################################
## MAIN CODE
###############################################################################

###############################
## These you can edit / change
###############################

## Years to run over
years = ['16', '17', '18']
#years = ['16', '18']
#years = ['17']

## regions
regions = [ "_ideal", "_3B", "_2b2c" ]
region_name = {
  "_ideal": "Tagging: 4 b-jets",
  "_3B": "Tagging: at least 3 b-jets",
  "_2b2c": "Tagging: 2 b-jets, 2 c-jets"
}

## Do we want to use a log scale?
useLogY = False

## Input & Output
resultpath = '../condor_results/
