#!/usr/bin/env python 
# catGetDatasetInfo v7-4-4 # to make dataset lists
# sed -i 's/^\/store/root:\/\/cms-xrdr.sdfarm.kr:1094\/\/xrd\/store/g' *

import os, json, array, sys
import numpy as np
from math import ceil       

analysis = sys.argv[1]

mcFiles_h2mumu = ['ttH',
                  "WWTo2L2Nu", "WZTo3LNu_amcatnlo", "WZTo2LQQ", "ZZTo2L2Nu", "ZZTo2L2Q",
                  "TTZToLLNuNu", "ZZTo4L_powheg", "ttWToLNu",
                  "WWW", "WWZ", "WZZ", "ZZZ",
                  "ZZ", "WZ", "WW"
]

mcFiles_topmass = [
                   'TT_powheg',
                #   'SingleTop_tW', 'SingleTbar_tW',
                #   'DYJets','DYJets_10to50',
                #   'WJets', 'ZZ', 'WW', 'WZ',
                #   'TT_powheg_mtop1665', 'TT_powheg_mtop1695', 'TT_powheg_mtop1715',
                #   'TT_powheg_mtop1735', 'TT_powheg_mtop1755', 'TT_powheg_mtop1785',
                  ]

mcFiles_vts = ["TT_powheg", 'DYJets', 'DYJets_MG_10to50']

dataFiles = ['SingleMuon_Run2016', 'SingleEG_Run2016',
             'DoubleMuon_Run2016', 'DoubleEG_Run2016',
             'MuonEG_Run2016', 'MuonEG_Run2016',
             'DoubleMuon_Run2016', 'DoubleEG_Run2016']
dataFiles = [data+period for period in ["B","Bv1","C","D","E","F","G","H"] for data in dataFiles]

if   analysis == 'h2mumu' : RunFiles = mcFiles_h2mumu  + dataFiles; analyser = "h2muAnalyser";
elif analysis == 'topMass': RunFiles = mcFiles_topmass + dataFiles; analyser = "topAnalyser";
elif analysis == 'vts'    : RunFiles = mcFiles_vts     + dataFiles; analyser = "vtsAnalyser";
elif analysis == 'cutbased': RunFiles = mcFiles_topmass; analyser = "cutbased";
else: print "put right name of analysis (h2mumu/topMass/hadron/vts/cutbased)"
#RunFiles = ['WW'] # for test
RunFiles = ['tsw'] # for test

maxFiles = 300
SetDir = "test"
datadir = '{}/src/nano/nanoAOD/data/dataset/dataset_'.format(os.environ['CMSSW_BASE'])

for datasetName in RunFiles:
    fileList = datadir + datasetName + '.txt'
    jobName = analysis+'_'+datasetName 

    Dirname = "{}/src/nano/analysis/test/Batch/{}".format(os.environ['CMSSW_BASE'],jobName)
    if os.path.isdir(Dirname):
        print "ERROR: output directory already existing."
        sys.exit()
    else: os.makedirs(Dirname)

    files = np.array([])
    for f in open(fileList).readlines():
        f = f.strip()
        f = f.strip('\',"')
        if len(f) < 5: continue
        if '#' == f[0] or '.root' != f[-5:]: continue
        files = np.append(files,[f])
    nFiles = len(files) 
    nSection = int(ceil(1.0*nFiles/maxFiles))
    for section in range(nSection):
        begin = section*maxFiles
        end = min(begin + maxFiles, nFiles)
        FileNames = files[begin:end]
        FileNamesStr = " ".join(str(i) for i in FileNames)

        print "@@ Writing run script..."
        jds = "%s/submit.jds" %Dirname 
        fout = open(jds, "w")
        print>>fout, "# Job description file for condor job"
        print>>fout, """executable = {0}/src/nano/analysis/test/Batch/run.sh
universe   = vanilla

log = condor.log

getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = {1}/job_{2}.log
error = {1}/job_{2}.err
queue""" .format(os.environ['CMSSW_BASE'], Dirname, section)
        fout.close()

        subBatch = "condor_submit -batch-name {} -append 'arguments={} {} {} {}' {}".format(datasetName, analyser, SetDir, datasetName, FileNamesStr ,jds)    
        os.system(subBatch)

