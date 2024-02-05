import os,sys
import time
import ROOT

#sNames = ['SingleElectron_DATA_2016']
sNames = ['DY_madgraphMLM_MC_2016','WW_MC_2016','WZ_MC_2016','ZZ_MC_2016']

folder = '/uscms/home/duong/Scratch/Output_ZplusB_withGenDr/condor_run_NONE/'
eosFolder = '/eos/uscms/store/user/duong/Output_ZplusB_withGenDr/NONE/'

compileCmd = 'make FORMAT=DATA_2016 FORMATPROCESSING=NORMAL INPUT=TCHAIN PDFSCALESYST=NONE'
runOpt = ' -data 1 -year 2016 -syst NONE -centralGenWeight 0'

for sName in sNames:
  continue
  work_dir = folder+'/'+sName+'_test/'
  os.chdir(work_dir)
  os.system(compileCmd)
  os.system('ls sampleList_*.txt >| tmp.txt')
  time.sleep(1)
  lines_tmp = open('tmp.txt', 'r').readlines()
  nJob = len(lines_tmp)
  print nJob
  
  for i in range(nJob):
    f = ROOT.TFile.Open(eosFolder+'/'+sName+'/output_'+str(i)+'.root','READ')
    badFile = False
    try:
      h = f.Get('Nevt')
      #print h.GetEntries()
      if h.GetEntries() < 1:
        badFile = True
    except:
      badFile = True
    if badFile:
      try:
        print 'Job ', i, ' was failed. Rerun now'
        cmd = './main -filelist sampleList_'+str(i)+'.txt -out output_'+str(i)+'.root ' + runOpt
        print cmd
        os.system(cmd)
        cmd = 'xrdcp output_'+str(i)+'.root root://cmseos.fnal.gov/'+eosFolder+'/'+sName+'/'
        print cmd
        os.system(cmd)
      except:
        print 'Crash in reprocessing ', i, ' .Skip'


  #nameTmp = 'outFileListTmp.txt'
  #  cmd_tmp = 'xrdfs root://cmseos.fnal.gov/ ls -u ' + dir_final_rootFile + '/ | grep ".root" >| ' + nameTmp
   # os.system(cmd_tmp)

