import os,sys

#work_dir = 'condor_run_NONE/SingleMuon_DATA_2017_test'
work_dir = 'condor_run_NONE/DY_madgraphMLM_MC_2016_test'
failJobs = ['9']

os.chdir(work_dir)
for failJob in failJobs:
  #create workdir
  work_dir_1 = 'Rerun_'+failJob
  os.system('rm -r ' + work_dir_1)
  os.system('mkdir ' + work_dir_1)
  os.system('cp input.tar ' + work_dir_1)
  #write executive file
  #exe_NONE_SingleMuon_DATA_2017.sh'
  lines = open('exe_NONE_'+work_dir.split('/')[1].replace('_test','')+'.sh').readlines()
  exeFileName = 'exe_'+failJob+'.sh'
  exeFile = open(exeFileName,'w')
  for l in lines:
    if 'tar -xvf input.tar' in l or 'make FORMAT' in l or './main -filelist' in l or 'xrdcp *.root' in l:
      if './main -filelist' in l:
        l = l.replace('$1',failJob)
      exeFile.write(l)
  exeFile.close()
  os.system('mv '+exeFileName+' '+work_dir_1)
  #sampleList_118.txt
  os.system('cp sampleList_'+failJob+'.txt '+work_dir_1)
  os.chdir(work_dir_1)
  os.system('chmod +x '+exeFileName)
  os.system('./'+exeFileName+'>&log.txt&')
  os.chdir('../')
