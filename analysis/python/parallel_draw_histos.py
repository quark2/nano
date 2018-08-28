#!/usr/bin/env python


import ROOT, json, os, sys, copy
import time, datetime


def valToName(strVar): 
  return strVar.replace(".", "").replace("(", "").replace(")", "").replace(" ", "").replace("*", "")


def makeHisto(strName, strTitle, listBin):
  if listBin[ 0 ] >= 0: 
    if len(listBin) == 3: 
      return ROOT.TH1D(strName, strTitle, listBin[ 0 ], listBin[ 1 ], listBin[ 2 ])
    elif len(listBin) == 6: 
      return ROOT.TH2F(strName, strTitle, listBin[ 0 ], listBin[ 1 ], listBin[ 2 ], 
        listBin[ 3 ], listBin[ 4 ], listBin[ 5 ])
    elif len(listBin) == 9: 
      return ROOT.TH3F(strName, strTitle, listBin[ 0 ], listBin[ 1 ], listBin[ 2 ], 
        listBin[ 3 ], listBin[ 4 ], listBin[ 5 ], listBin[ 6 ], listBin[ 7 ], listBin[ 8 ])
  else: 
    if listBin[ 1 ] == "log1D": 
      import array
      
      nBin = listBin[ 2 ]
      fBinRatio = ( listBin[ 4 ] / listBin[ 3 ] ) ** ( 1.0 / nBin )
      
      listBinAct = []
      fBin = listBin[ 3 ]
      
      for i in range(nBin + 1): 
        listBinAct.append(fBin)
        fBin *= fBinRatio
      
      return ROOT.TH1D(strName, strTitle, nBin, array.array("d", listBinAct)) 


# Usage of the following class:
#   from [PATH].parallel_draw_histos import parallel_draw_histos
#   paraMain = parallel_draw_histos()
#   
#   paraMain.SetDirHistName("[NAME]")
#   paraMain.SetSrcPath("[PATH]")
#   paraMain.SetNameTree("[TREE NAME]")
#   paraMain.SetCut("[CUT]")
#   paraMain.SetWeight("[WEIGHT]")
#   paraMain.SetDatasets([LIST OF DATASETS])
#   paraMain.SetVars({DICT OF VARIABLES})
#   paraMain.SetPathDraw("[PATH]")
#   
#   paraMain.InitRun()
#   paraMain.RunOnCluster() or paraMain.RunOnWorknode()


class parallel_draw_histos:
  def __init__(self):
    self.strPathThis = os.path.join("%s/src/nano/analysis/python"%os.environ[ "CMSSW_BASE" ], 
      "parallel_draw_histos.py")
    
    self.strTypeArgOneRoot = "--oneroot"
    self.strTypeArgMerger = "--merger"
    
    self.nMaxProc = 20
  
  
  #### Get[...] area ####
  
  
  # strDirHist : The name of the directory in which the result will be
  def GetDirHistName(self): 
    return self.strDirHist
  
  
  # strSrcPath : The path (MUTE BE xrd protocol) in which the source ntuples is located
  def GetSrcPath(self): 
    return self.strSrcPath
  
  
  # strNameTree : The name of the tree in the ntuple
  def GetNameTree(self): 
    return self.strNameTree
  
  
  # strCut : Cut condition
  def GetCut(self): 
    return self.strCut
  
  
  # strWeight : Weight furmula
  def GetWeight(self): 
    return self.strWeight
  
  
  # listDatasets : List of datasets
  def GetDatasets(self): 
    return self.listDatasets
  
  
  # dicVars : Dictionary containing infos for variables and histograms which will be drawn
  def GetVars(self): 
    return self.dicVars
  
  
  # strPathDraw : The path for destination
  def GetPathDraw(self): 
    return self.strPathDraw
  
  
  # strFilenameDumpJSON : The path of dump JSON file
  def GetPathDumpJSON(self): 
    return self.strFilenameDumpJSON
  
  
  #### Set[...] area ####
  
  
  # strDirHist : The name of the directory in which the result will be
  def SetDirHistName(self, strName = ""): 
    self.strDirHist = strName if strName != "" else datetime.datetime.now().strftime("%y%m%d_%H%M%S")
  
  
  # strSrcPath : The path (MUTE BE xrd protocol) in which the source ntuples is located
  def SetSrcPath(self, strSrcPath): 
    self.strSrcPath = strSrcPath
  
  
  # strNameTree : The name of the tree in the ntuple
  def SetNameTree(self, strNameTree): 
    self.strNameTree = strNameTree
  
  
  # strCut : Cut condition
  def SetCut(self, strCut): 
    self.strCut = strCut
  
  
  # strWeight : Weight furmula
  def SetWeight(self, strWeight): 
    self.strWeight = strWeight
  
  
  # listDatasets : List of datasets
  def SetDatasets(self, listDatasets): 
    self.listDatasets = listDatasets
  
  
  # dicVars : Dictionary containing infos for variables and histograms which will be drawn
  def SetVars(self, dicVars): 
    self.dicVars = dicVars
  
  
  # strPathDraw : The path for destination
  def SetPathDraw(self, strPathDraw): 
    self.strPathDraw = strPathDraw
  
  
  # strFilenameDumpJSON : The path of dump JSON file
  def SetPathDumpJSON(self, strFilenameDumpJSON): 
    self.strFilenameDumpJSON = strFilenameDumpJSON
  
  
  def InitRun(self): 
    # Preparing to throwing jobs: making basic directories
    if not os.path.exists(os.path.join(self.GetPathDraw(), self.GetDirHistName())): 
      os.makedirs(os.path.join(self.GetPathDraw(), self.GetDirHistName()))
      os.makedirs(os.path.join(self.GetPathDraw(), self.GetDirHistName(), "logs"))
    
    # Preparing to throwing jobs: making a card (in JSON format)
    dicOut = {}
    
    dicOut[ "dataset" ] = []
    dicOut[ "rootfile" ] = ["root://cms-xrdr.sdfarm.kr:1094///xrd" + self.GetSrcPath()[ 7: ], 
      "dir_%(dataset)s", []]
    dicOut[ "treename" ] = self.GetNameTree()
    dicOut[ "cut" ] = self.GetCut()
    dicOut[ "weight" ] = self.GetWeight()
    dicOut[ "Vars" ] = self.GetVars()
    dicOut[ "res_dir" ] = os.path.join(self.GetPathDraw(), self.GetDirHistName())
    
    self.listDatasetOfFile = []
    self.dicNumFiles = {}
    self.nTotalNumFiles = 0
    
    # Preparing to throwing jobs: inserting lists of source ntuples into the card
    for strDataset in self.GetDatasets(): 
      listFiles = os.listdir(os.path.join(self.GetSrcPath(), "dir_" + strDataset))
      
      self.listDatasetOfFile.append(strDataset)
      self.dicNumFiles[ strDataset ] = len(listFiles)
      self.nTotalNumFiles += len(listFiles)
      
      strDirDataset = os.path.join(self.GetPathDraw(), self.GetDirHistName(), strDataset)
      if not os.path.exists(strDirDataset): 
        os.makedirs(strDirDataset)
      
      for strFilename in listFiles: 
        dicOut[ "dataset" ].append(strDataset)
        dicOut[ "rootfile" ][ 2 ].append(strFilename)
    
    self.SetPathDumpJSON(os.path.join(self.GetPathDraw(), self.GetDirHistName(), "config.json"))
    fWriteJSON = open(self.GetPathDumpJSON(), "w")
    json.dump(dicOut, fWriteJSON)
    fWriteJSON.close()
  
  
  def RunOnCluster(self): 
    # The card for condor submit
    strJDSTemplate  = "executable = " + self.strPathThis
    strJDSTemplate += """
universe   = vanilla
requirements = OpSysMajorVer == 6
arguments=%(flag)s %(json)s $(Process)

log = %(path)s/condor.log

getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = %(path)s/logs/job_$(Process).log
error = %(path)s/logs/job_$(Process).err
queue %(queue)i
    """
  
    # completing the submit card
    strJDS = strJDSTemplate%{"flag": self.strTypeArgOneRoot, 
      "json": self.GetPathDumpJSON(), "queue": self.nTotalNumFiles, "path": self.GetDirHistName()}
    
    # Submit!
    os.system("printf '%s' | condor_submit > /dev/null"%strJDS)
    
    print "All jobs have been thrown"
    
    dicMerge = self.dicNumFiles
    strDirDraw = os.path.join(self.GetPathDraw(), self.GetDirHistName())
    
    # Running codes for merging
    # Checking the progress (by counting the number of drawn root files) for each dataset, 
    # if it is done, then running the merger
    while len(dicMerge.keys()) > 0: 
      dicNext = {}
      
      for strDataset in dicMerge.keys(): 
        # Counting the result root files
        nNumDone = len([ s for s in os.listdir(strDirDraw) if strDataset in s and ".done" in s ])
        if nNumDone >= dicMerge[ strDataset ]: 
          # Merge!
          self.DoMerge(strDataset, self.dicNumFiles[ strDataset ])
          # If complete, it must be pulled out in the next dataset list
        else: 
          dicNext[ strDataset ] = dicMerge[ strDataset ]
      
      dicMerge = dicNext
      
    os.system("date")
    print self.GetDirHistName()
  
  
  def RunOnWorknode(self): 
    strDirHistFull = os.path.join(self.GetPathDraw(), self.GetDirHistName())
    nNumJob = 0
    
    for strDataset in self.listDatasetOfFile: 
      for iii in range(self.dicNumFiles[ strDataset ]): 
        while True:
          nDone = len([ s for s in os.listdir(strDirHistFull) if ".done" in s ])
          if nDone + self.nMaxProc > nNumJob: break
          time.sleep(1)
        
        os.system("python %s %s %s %i &"%(self.strPathThis, self.strTypeArgOneRoot, 
          self.GetPathDumpJSON(), nNumJob))
        nNumJob += 1
      
      # Merge!
      os.system("python %s %s %s %s %s %i &"%(self.strPathThis, self.strTypeArgMerger, 
        self.GetPathDraw(), self.strDirHist, strDataset, self.dicNumFiles[ strDataset ]))
    
    while True:
      nDone = len([ s for s in os.listdir(strDirHistFull) if ".root" in s ])
      #if nDone >= nNumJob: break
      if nDone >= len(self.dicNumFiles): break
      time.sleep(1)
    
    os.system("date")
    print self.GetDirHistName()
  
  
  def DoOneRoot(self, strPathJSON, nIdxJob): 
    fJSONQuery = open(strPathJSON, "r")
    dicMain = json.load(fJSONQuery)
    fJSONQuery.close()
    
    # Loading configurations
    strDataset = dicMain[ "dataset" ][ nIdxJob ].encode("ascii", "ignore")
    
    strWeight = dicMain[ "weight" ].encode("ascii", "ignore")
    strCut    = dicMain[ "cut" ].encode("ascii", "ignore")
    
    strPathMain = dicMain[ "rootfile" ][ 0 ].encode("ascii", "ignore")
    strDirName =  dicMain[ "rootfile" ][ 1 ].encode("ascii", "ignore")%{"dataset": strDataset}
    strNoFile =   dicMain[ "rootfile" ][ 2 ][ nIdxJob ].encode("ascii", "ignore")
    
    dicVar = dicMain[ "Vars" ]
    dicOutInfo = {}
    
    dicOutInfo[ "dataset" ] = strDataset
    dicOutInfo[ "ntuple" ] = os.path.join(strPathMain, strDirName, strNoFile)
    dicOutInfo[ "no_file" ] = strNoFile.split("_")[ 1 ].split(".")[ 0 ]
    
    dicOutInfo[ "weight" ] = strWeight
    dicOutInfo[ "cut" ] = strCut
    
    dicHist = {}
    
    # Loading the source ntuple root file
    strNameTree = dicMain[ "treename" ].encode("ascii", "ignore")
    fMain = ROOT.TFile.Open(dicOutInfo[ "ntuple" ])
    tree = fMain.Get(strNameTree)
    
    # Applying the cut
    # This way applies the cut only once for several histograms, so enhances the performances
    # Btw, if no entries in the existing tree, the following doesn't work, and also it is not needed
    if tree.GetEntries() > 0: 
      strNameAfterCut = strNameTree.replace("/", "") + "_aftercut"
      tree.Draw(">>" + strNameAfterCut, strCut, "entrylist")
      tree.SetEntryList(ROOT.gDirectory.Get(strNameAfterCut))
    
    # Drawing histograms
    for strVar in dicVar.keys(): 
      strNameVar = valToName(strVar)
      strNameHist  = "hist_%(dataset)s_%(var)s"%{"dataset": strDataset, "var": strNameVar}
      
      # Two cases: histogram drawn from the loaded tree, or histogram which exists already
      if type(dicVar[ strVar ]) is dict: 
        hTmp = makeHisto(strNameHist, strNameHist, dicVar[ strVar ][ "bin" ])
        tree.Project(strNameHist, strVar, strWeight)
        
        dicHist[ strVar ] = copy.deepcopy(hTmp)
        
        if dicHist[ strVar ].GetSumw2N() == 0:
          dicHist[ strVar ].Sumw2()
      else: 
        hTmp = copy.deepcopy(fMain.Get(strVar.encode("ascii", "ignore")))
        hTmp.SetName(strNameHist)
        
        dicHist[ strVar ] = hTmp
    
    # Now, ready to write
    strOutRoot = os.path.join(dicMain[ "res_dir" ], strDataset, 
      dicOutInfo[ "no_file" ] + ".root")
    fRootFile = ROOT.TFile(strOutRoot, "NEW")
    
    # Writing infos; the first is for string infos and the second is the histograms
    for strKey in dicOutInfo.keys(): ROOT.TNamed(strKey, dicOutInfo[ strKey ]).Write()
    for strVar in dicHist.keys(): dicHist[ strVar ].Write()
    
    # DO NOT FORGET!
    fRootFile.Write()
    fRootFile.Close()
    
    # Notifying that I'm done, in somewhat primitive way, but definite way
    os.system("touch %s"%(os.path.join(dicMain[ "res_dir" ], 
        strDataset + strNoFile.split("_")[ 1 ].split(".")[ 0 ] + ".done")))
  
  
  def DoMerge(self, strDataset, nNumJobs): 
    import shutil, time
    
    if nNumJobs <= 0: 
      print "Wrong arguments"
      sys.exit(1)
    
    strDirDraw = os.path.join(self.GetPathDraw(), self.GetDirHistName())
    
    while True:
      if len([ s for s in os.listdir(strDirDraw) if strDataset in s and ".done" in s ]) >= nNumJobs: break
      time.sleep(1)
    
    strDirDataset = os.path.join(self.GetPathDraw(), self.GetDirHistName(), strDataset)
    listRoot = [ s for s in os.listdir(strDirDataset) if ".root" in s ]
    
    fOrg = ROOT.TFile(os.path.join(strDirDataset, listRoot[ 0 ]))
    listRoot.pop(0)
    
    dicHisto = {}
    dicInfo  = {}
    i = 0
    
    while True:
      objkeyGet = fOrg.GetListOfKeys().At(i)
      if objkeyGet == None: break
      strNameObj = objkeyGet.GetName()
      
      if strNameObj.startswith("hist_"): 
        objGet = fOrg.Get(strNameObj)
        dicHisto[ strNameObj ] = copy.deepcopy(objGet)
        dicHisto[ strNameObj ].Sumw2(False)
      else: 
        dicInfo[ strNameObj ] = fOrg.Get(strNameObj).GetTitle()
      
      i += 1
    
    fOrg.Close()
    
    for strRootFile in listRoot: 
      fRead = ROOT.TFile(os.path.join(strDirDataset, strRootFile))
      
      for strNameHist in dicHisto.keys(): 
        histRead = fRead.Get(strNameHist)
        histRead.Sumw2(False)
        dicHisto[ strNameHist ].Add(histRead)
      
      fRead.Close()
    
    strOut = os.path.join(self.GetPathDraw(), self.GetDirHistName(), "hist_" + strDataset + ".root")
    fMerge = ROOT.TFile(strOut, "CREATE")
    
    for strKey in dicInfo.keys():  ROOT.TNamed(strKey, dicInfo[ strKey ]).Write()
    for strKey in dicHisto.keys(): dicHisto[ strKey ].Write()
    
    fMerge.Write()
    fMerge.Close()
    
    print "%s is done"%(strDataset)


if __name__ == "__main__": 
  paraMain = parallel_draw_histos()
  
  if len(sys.argv) < 2: 
    sys.stderr.write("Error: Too few arguments\n")
    sys.exit(1)
  
  if sys.argv[ 1 ] == paraMain.strTypeArgOneRoot: 
    if len(sys.argv) < 4: 
      sys.stderr.write("Error: Too few arguments\n")
      sys.exit(1)
    
    paraMain.DoOneRoot(sys.argv[ 2 ], int(sys.argv[ 3 ]))
  elif sys.argv[ 1 ] == paraMain.strTypeArgMerger: 
    if len(sys.argv) < 6: 
      sys.stderr.write("Error: Too few arguments\n")
      sys.exit(1)
    
    paraMain.SetPathDraw(sys.argv[ 2 ])
    paraMain.SetDirHistName(sys.argv[ 3 ])
    paraMain.DoMerge(sys.argv[ 4 ], int(sys.argv[ 5 ]))


