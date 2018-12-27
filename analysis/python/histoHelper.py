import math, array, ROOT, copy, CMS_lumi, tdrstyle
import PhysicsTools.PythonAnalysis.rootplot.core as rootplotcore
tdrstyle.setTDRStyle()
def defTH1(title, name, binning):
    if len(binning) == 3:
        hist = ROOT.TH1D(name, title, binning[0], binning[1], binning[2])
    else:
        hist = ROOT.TH1D(name, title, len(binning)-1, array.array('f', binning))
    return hist

def getTH1(title, binning, tree, plotvar, cut, scale = 0.):
    hist = defTH1(title, "name", binning)
    tree.Project("name", plotvar, cut)
    if hist.GetSumw2N() == 0:
        hist.Sumw2()
    if scale != 0:
        hist.Scale(scale)
    return copy.deepcopy(hist)

def makeTH1(filename, treename, title, binning, plotvar, cut, scale = 0.):
    tfile = ROOT.TFile(filename)
    tree  = tfile.Get(treename)
    
    hist = defTH1(title, "tmp", binning)
    for var in plotvar.split(','):
        hist.Add(getTH1(title, binning, tree, var, cut, scale))
        
    return copy.deepcopy(hist)

def getEntries(filename, treename):
    tfile = ROOT.TFile(filename)
    tree  = tfile.Get(treename)
    return tree.GetEntriesFast()

def getWeightedEntries(filename, treename, plotvar, weight):
    weighthist = makeTH1(filename, treename, '', [1, 0, 1], plotvar, weight)    
    return weighthist.Integral(-1,2)

def divide_canvas(canvas, ratio_fraction):
    margins = [ROOT.gStyle.GetPadTopMargin(), ROOT.gStyle.GetPadBottomMargin()]
    useable_height = 1 - (margins[0] + margins[1])
    canvas.Clear()
    pad = ROOT.TPad('mainPad', 'mainPad', 0., 0., 1., 1.)
    pad.SetFillStyle(4000)
    pad.Draw()
    pad.SetBottomMargin(margins[1] + ratio_fraction * useable_height)
    pad_ratio = ROOT.TPad('ratioPad', 'ratioPad', 0., 0., 1., 1.);
    pad_ratio.SetFillStyle(4000)
    pad_ratio.Draw()
    pad_ratio.SetTopMargin(margins[0] + (1 - ratio_fraction) * useable_height)
    return pad, pad_ratio

def makeCanvas(name, doRatio = False):
    H_ref = 600
    W_ref = 800
    if doRatio:
        H_ref = 800
    canv = ROOT.TCanvas(name,name,W_ref,H_ref)    
    return canv

def setMargins(canvas, doRatio = False):
    H_ref = 600
    W_ref = 800
    if doRatio:
        H_ref = 800
    W = W_ref
    H = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref 
    L = 0.12*W_ref
    R = 0.04*W_ref
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    #canvas.SetTopMargin( T/H )
    #canvas.SetBottomMargin( B/H )
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    return canvas

def setDefAxis(axis, title, offset, titlefont = 42, titlesize = 0.043, labelfont = 42, labelsize = 0.03):
    axis.SetTitle(title)
    axis.SetTitleOffset(offset)
    axis.SetTitleColor(1)
    axis.SetTitleFont(titlefont)
    axis.SetTitleSize(titlesize)
    axis.SetLabelColor(1)
    axis.SetLabelFont(labelfont)
    axis.SetLabelOffset(0.007)
    axis.SetLabelSize(labelsize)
    axis.SetAxisColor(1)
    axis.SetTickLength(0.03)
    axis.SetNdivisions(510)
    #axis.SetStripDecimals(True)
    #axis.SetPadTickX(1)
    #axis.SetPadTickY(1)

def setDefTH1Style(th1, x_name, y_name, titlefont = 42, titlesize = 0.043, labelfont = 42, labelsize = 0.03):
    setDefAxis(th1.GetYaxis(),y_name, 1.2, titlefont, titlesize, labelfont, labelsize)
    setDefAxis(th1.GetXaxis(),x_name, 1, titlefont, titlesize, labelfont, labelsize)
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.cd()
    return th1
    
def drawTH1(name, cmsLumi, mclist, data, x_name, y_name, doLog=False, doLogX=False, doRatio=True, ratioRange=0.45, legx=0.68, legfontsize=0.030, histSig=None):
    leg = ROOT.TLegend(legx,0.68,legx+0.2,0.91)
    leg.SetBorderSize(0)
    #leg.SetNColumns(2)
    leg.SetTextSize(legfontsize)
    leg.SetTextFont(42)
    leg.SetLineColor(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.AddEntry(data,"Data","lp")
    
    hs = ROOT.THStack("mcstack", "mcstack")    
    hratio = mclist[0].Clone("hratio")
    hratio.Reset()

    leghist = []
    for i, mc in enumerate(mclist):
        hnew = mc.Clone("hnew"+mc.GetName())
        hnew.Sumw2(False)
        hs.Add(hnew)
        hratio.Add(mc)
        inversed = mclist[len(mclist)-1-i]
        if not any(inversed.GetTitle() == s for s in leghist):
            leg.AddEntry(inversed, inversed.GetTitle(), "f")
            leghist.append(inversed.GetTitle())
                        
    #hratio = hs.GetStack().Last()
    hratio.Divide(data,hratio,1.,1.,"B")

    tdrstyle.setTDRStyle()

    setDefTH1Style(data, x_name, y_name)
    data.SetName('data')
    data.SetMaximum(data.GetMaximum()*1.8)
    if doLog:
        data.SetMaximum(data.GetMaximum()*100)
        #data.SetMinimum(10**-3)
    else:
        data.GetYaxis().SetTitleSize(0.04)
        data.GetYaxis().SetLabelSize(0.024)
        data.GetYaxis().SetTitleOffset(1.35)
        
    ratio_fraction = 0
    if doRatio:
        ratio_fraction = 0.3        
        data.GetXaxis().SetLabelSize(0)
        data.GetXaxis().SetTitleSize(0)
        setDefTH1Style(hratio, x_name, "Data/MC")
        hratio.GetYaxis().CenterTitle()
        hratio.GetYaxis().SetNdivisions(5)
            
    canv = makeCanvas(name, doRatio)
    pads=[canv]
    pads = rootplotcore.divide_canvas(canv, ratio_fraction)

    pads[0].cd()
    setMargins(pads[0], doRatio)
    if doLog:
        pads[0].SetLogy()
    
    if doLogX: 
        pads[ 0 ].SetLogx()
        pads[ 1 ].SetLogx()

    data.Draw()
    hs.Draw("same")
    data.Draw("esamex0")
    
    if histSig is not None: 
      histSig.Scale(data.Integral() / histSig.Integral())
      histSig.SetLineColor(2)
      
      histSig.Draw("same HIST")
      leg.AddEntry(histSig, "signal (scaled)", "l")
    
    leg.Draw("same")
    
    pads[0].Update()

    if doRatio:
        pads[1].cd()
        pads[1].SetGridy()
        setMargins(pads[1], doRatio)
        hratio.SetLineColor(1)
        hratio.Draw("e")
        hratio.SetMaximum(1.+ratioRange)
        hratio.SetMinimum(1.-ratioRange)
       # hratio.SetMaximum(2)
       # hratio.SetMinimum(0)

    for p in pads:
        p.RedrawAxis()
        p.Modified()
        p.Update()

    canv.cd()

    #iPos = 0 # in frame
    iPos = 11 # out frame
    if( iPos==0 ):
        cmsLumi.relPosX = 0.1
    cmsLumi.CMS_lumi(pads[0], 0, iPos)

    canv.Modified()
    canv.Update()
    return copy.deepcopy(canv)

def drellYanEstimation(mc_ee_in, mc_ee_out, mc_mm_in, mc_mm_out,
                       rd_ee_in, rd_mm_in, rd_em_in, kMM, kEE):
    #kMM = math.sqrt(rd_mm_in/rd_ee_in)/2.
    #kEE = math.sqrt(rd_ee_in/rd_mm_in)/2.

    rMC_mm = mc_mm_out/mc_mm_in
    rMC_ee = mc_ee_out/mc_ee_in
    print "rMC_mm  ", rMC_mm
    print "rMC_ee  ", rMC_ee
    
    nOutEst_mm = rMC_mm*(rd_mm_in - rd_em_in*kMM)
    nOutEst_ee = rMC_ee*(rd_ee_in - rd_em_in*kEE)
    return nOutEst_ee/mc_ee_out,nOutEst_mm/mc_mm_out

def findDataSet(name, datasets):
    for data in datasets:
        if data["name"] == name:
            return data
    return None

def adderrs(err1, err2, sign=1.):
    return math.sqrt(err1**2+sign*err2**2)

def table(mchistList, errList, signal_hist, signal_err):
    nums = {}
    errs = {}
    total = total_err = 0

    titles = list(set([mc.GetTitle() for mc in mchistList]))
    for t in titles:
        nums[t] = 0
        errs[t] = 0

    for i, mc in enumerate(mchistList):
        nbins = mc.GetSize()-2
        nums[mc.GetTitle()] += mc.Integral(0,nbins+1)
        errs[mc.GetTitle()] = adderrs(errs[mc.GetTitle()], errList[i])

        total += mc.Integral(0,nbins+1)
        total_err = adderrs(total_err, errList[i])
    
    nums['total'] = total
    errs['total'] = total_err

    bkg = total - signal_hist.Integral(0,signal_hist.GetSize()-1)
    bkg_err = adderrs(total_err, signal_err, -1)
    nums['bkg'] = bkg
    errs['bkg'] = bkg_err

    return nums, errs

def set_palette(name="", ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    gStyle.SetNumberContours(ncontours)

def overFlow(hist):
    nbins = hist.GetNbinsX()
    hist.SetBinContent(nbins, hist.GetBinContent(nbins)+hist.GetBinContent(nbins+1))
    hist.SetBinError(nbins, math.sqrt(hist.GetBinError(nbins)**2+hist.GetBinError(nbins+1)**2))
    hist.SetBinContent(1, hist.GetBinContent(1)+hist.GetBinContent(0))
    hist.SetBinError(1, math.sqrt(hist.GetBinError(1)**2+hist.GetBinError(0)**2))

def extraText(canv, position, content):
    canv.cd()
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.DrawLatex(position[0], position[1], content)
    canv.Update()

def sysUncertainty(filename, treename, binning, title, scale, cut, plotvar, h_nom, sysList):
    sysErrs_up = []
    sysErrs_dn = []
    for sys in sysList:
        if 'weight' not in sys:
            sysErrs_up.append(makeTH1(filename, "cattree/%s_u"%sys, title, binning, plotvar, cut, scale))
            sysErrs_dn.append(makeTH1(filename, "cattree/%s_d"%sys, title, binning, plotvar, cut, scale))
        else:
            sysErrs_up.append(makeTH1(filename, treename, title, binning, plotvar, cut.replace(sys,'%s_up'%sys), scale))
            sysErrs_dn.append(makeTH1(filename, treename, title, binning, plotvar, cut.replace(sys,'%s_dn'%sys), scale))

    for i in range(len(sysList)):
        sysErrs_up[i].Add(h_nom, -1)
        sysErrs_dn[i].Add(h_nom, -1)

    return sysErrs_up, sysErrs_dn


def drawRatioPlot(name, cmsLumi, mclist, data, x_name, y_name, doLog=False, doRatio=True, ratioRange=0.45, legx=0.68, legfontsize=0.030):
    leg = ROOT.TLegend(legx,0.68,legx+0.2,0.91)
    leg.SetBorderSize(0)
    #leg.SetNColumns(2)
    leg.SetTextSize(legfontsize)
    leg.SetTextFont(42)
    leg.SetLineColor(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.AddEntry(data,"Data","lp")
    
    hs = ROOT.THStack("mcstack", "mcstack")    
    hAllMC = mclist[0].Clone("hratio")
    hAllMC.Reset()

    leghist = []
    for i, mc in enumerate(mclist):
        hnew = mc.Clone("hnew"+mc.GetName())
        hnew.Sumw2(False)
        hs.Add(hnew)
        hAllMC.Add(mc)
        inversed = mclist[len(mclist)-1-i]
        if not any(inversed.GetTitle() == s for s in leghist):
            leg.AddEntry(inversed, inversed.GetTitle(), "f")
            leghist.append(inversed.GetTitle())
                        
    #hratio = hs.GetStack().Last()
    hratio = ROOT.TRatioPlot(data,hAllMC)
    hratio.SetSeparationMargin(0.0)

    tdrstyle.setTDRStyle()
    
    setDefTH1Style(hratio, x_name, y_name)

    canv = makeCanvas(name, False)
    hratio.Draw("")
    leg.Draw("same")

    cmsLumi.CMS_lumi(pads[0], 0, iPos)

    canv.Modified()
    canv.Update()
    return copy.deepcopy(canv)


class TH1drawer: 
  def __init__(self, cmsLumi = None, listMC = None, histRD = None, strXLabel = "", strYLabel = "", bDoLog = False, bDoRatio = True): 
    self.ratioMax = 1.8
    self.minPlot = None
    
    self.ratioRange = 0.45
    
    self.legx = 0.68
    self.legy = 0.91
    self.legwidth  = 0.20
    self.legheight = 0.23
    self.legfontsize = 0.030
    self.legfont = 42
    
    self.xtitlefont = 42
    self.xtitlesize = 0.043
    self.ytitlefont = 42
    self.ytitlesize = 0.043
    
    self.xlabelfont = 42
    self.xlabelsize = 0.03
    self.ylabelfont = 42
    self.ylabelsize = 0.03
    
    self.ytitleoffset = 1.35
    
    self.cmsLumi = cmsLumi
    self.mclist = listMC 
    
    self.data = histRD
    self.histSig = None
    
    self.x_name = strXLabel
    self.y_name = strYLabel
    
    self.doLog = bDoLog
    self.doLogX = False
    self.doRatio = bDoRatio
  
  # You must call at least once the Set methods in the followings before drawing
  
  def GetCMSLumi(self): return self.cmsLumi
  def SetCMSLumi(self, cmsLumi): self.cmsLumi = cmsLumi
  
  def GetMCPlot(self): return self.mclist
  def SetMCPlot(self, mclist): self.mclist = mclist
  
  def GetXaxis(self): return self.x_name
  def SetXaxis(self, x_name): self.x_name = x_name
  
  def GetYaxis(self): return self.y_name
  def SetYaxis(self, y_name): self.y_name = y_name
  
  # For unnecessary variables
  
  def GetDataPlot(self): return self.data
  def SetDataPlot(self, data): self.data = data
  
  def IsLog(self): return self.doLog
  def DoLog(self): self.doLog = True
  def UndoLog(self): self.doLog = False
  
  def IsLogX(self): return self.doLogX
  def DoLogX(self): self.doLogX = True
  def UndoLogX(self): self.doLogX = False
  
  def GetRatioMax(self): return self.ratioMax
  def SetRatioMax(self, ratioMax): self.ratioMax = ratioMax
  
  def GetMinimum(self): return self.minPlot
  def SetMinimum(self, minPlot): self.minPlot = minPlot
  
  def GetRatioRange(self): return self.ratioRange
  def SetRatioRange(self, ratioRange): self.ratioRange = ratioRange
  
  def GetXTitleFontSize(self): return self.xtitlesize
  def SetXTitleFontSize(self, xtitlesize): self.xtitlesize = xtitlesize
  
  def GetXTitleFont(self): return self.xtitlefont
  def SetXTitleFont(self, xtitlefont): self.xtitlefont = xtitlefont
  
  def GetYTitleFontSize(self): return self.ytitlesize
  def SetYTitleFontSize(self, ytitlesize): self.ytitlesize = ytitlesize
  
  def GetYTitleFont(self): return self.ytitlefont
  def SetYTitleFont(self, ytitlefont): self.ytitlefont = ytitlefont
  
  def GetXLabelFontSize(self): return self.xlabelsize
  def SetXLabelFontSize(self, xlabelsize): self.xlabelsize = xlabelsize
  
  def GetXLabelFont(self): return self.xlabelfont
  def SetXLabelFont(self, xlabelfont): self.xlabelfont = xlabelfont
  
  def GetYLabelFontSize(self): return self.ylabelsize
  def SetYLabelFontSize(self, ylabelsize): self.ylabelsize = ylabelsize
  
  def GetYLabelFont(self): return self.ylabelfont
  def SetYLabelFont(self, ylabelfont): self.ylabelfont = ylabelfont
  
  def GetYTitleOffset(self): return self.ytitleoffset
  def SetYTitleOffset(self, ytitleoffset): self.ytitleoffset = ytitleoffset
  
  def GetLegendX(self): return self.legx
  def SetLegendX(self, legx): self.legx = legx
  
  def GetLegendY(self): return self.legy
  def SetLegendY(self, legy): self.legy = legy
  
  def GetLegendW(self): return self.legwidth
  def SetLegendW(self, legwidth): self.legwidth = legwidth
  
  def GetLegendH(self): return self.legheight
  def SetLegendH(self, legheight): self.legheight = legheight
  
  def GetLegendFont(self): return self.legfont
  def SetLegendFont(self, legfont): self.legfont = legfont
  
  def GetLegendFontSize(self): return self.legfontsize
  def SetLegendFontSize(self, legfontsize): self.legfontsize = legfontsize
  
  
  # Functions for drawing; you can rewrite them in your daughter class
  
  
  def MakeLegend(self): 
    self.leg = ROOT.TLegend(self.legx, self.legy - self.legheight, self.legx + self.legwidth, self.legy)
    self.leg.SetBorderSize(0)
    #leg.SetNColumns(2)
    self.leg.SetTextSize(self.legfontsize)
    self.leg.SetTextFont(self.legfont)
    self.leg.SetLineColor(0)
    self.leg.SetFillColor(0)
    self.leg.SetFillStyle(0)
    
    if self.data is not None: self.leg.AddEntry(self.data, "Data", "lp")
  
  
  def ArrangeMCPlots(self): 
    self.hs = ROOT.THStack("mcstack", "mcstack")
    self.hratio = self.mclist[ 0 ].Clone("hratio")
    self.hratio.Reset()
    
    leghist = []
    for i, mc in enumerate(self.mclist):
      hnew = mc.Clone("hnew" + mc.GetName())
      hnew.Sumw2(False)
      
      self.hs.Add(hnew)
      self.hratio.Add(mc)
      
      inversed = self.mclist[ len(self.mclist) - 1 - i ]
      
      if not any(inversed.GetTitle() == s for s in leghist):
        self.leg.AddEntry(inversed, inversed.GetTitle(), "f")
        leghist.append(inversed.GetTitle())
    
    #self.hratio = self.hs.GetStack().Last()
    if self.doRatio: self.hratio.Divide(self.data, self.hratio, 1, 1, "B")
  
  
  def SetBaseHistogram(self): 
    if self.data is not None: 
      self.hBase = self.data
      self.strNameBase = "data"
      self.strDrawOption = ""
    else: 
      self.hBase = copy.deepcopy(self.hs.GetStack().Last())
      self.strNameBase = "data"
      self.strDrawOption = ""
  
  
  def DrawMCPlots(self, option = ""): 
    self.hs.Draw("same" + option)
  
  
  def DoDraw(self, name): 
    if self.cmsLumi is None: return None
    if self.mclist is None: return None
    
    if self.data is None: 
      self.doRatio = False
    
    self.MakeLegend()
    self.ArrangeMCPlots()
    
    tdrstyle.setTDRStyle()
    
    if self.data is not None: self.data.SetName("data")
    
    self.SetBaseHistogram()
    
    setDefTH1Style(self.hBase, self.x_name, self.y_name, 
      self.xtitlefont, self.xtitlesize, self.xlabelfont, self.xlabelsize)
    self.hBase.SetName(self.strNameBase)
    self.hBase.SetMaximum(self.hBase.GetMaximum() * self.ratioMax)
    
    if self.doLog:
      data.SetMaximum(self.hBase.GetMaximum() * 100)
      #data.SetMinimum(10 ** -3)
    else:
      self.hBase.GetYaxis().SetTitleSize(self.ytitlesize)
      self.hBase.GetYaxis().SetLabelSize(self.ylabelsize)
      self.hBase.GetYaxis().SetTitleOffset(self.ytitleoffset)
        
    ratio_fraction = 0
    if self.doRatio:
      ratio_fraction = 0.3
      self.hBase.GetXaxis().SetLabelSize(0)
      self.hBase.GetXaxis().SetTitleSize(0)
      setDefTH1Style(self.hratio, self.x_name, "Data/MC", 
        self.xtitlefont, self.xtitlesize, self.xlabelfont, self.xlabelsize)
      self.hratio.GetYaxis().CenterTitle()
      self.hratio.GetYaxis().SetNdivisions(5)
            
    canv = makeCanvas(name, self.doRatio)
    pads=[canv]
    pads = rootplotcore.divide_canvas(canv, ratio_fraction)
    
    pads[ 0 ].cd()
    setMargins(pads[ 0 ], self.doRatio)
    if self.doLog:
      pads[ 0 ].SetLogy()
    
    if self.doLogX: 
      pads[ 0 ].SetLogx()
      pads[ 1 ].SetLogx()
    
    if self.minPlot is not None: self.hBase.SetMinimum(self.minPlot)
    
    self.hBase.Draw(self.strDrawOption)
    self.DrawMCPlots()
    if self.data is not None: self.data.Draw("esamex0")
    
    if self.histSig is not None: 
      self.histSig.Scale(self.data.Integral() / self.histSig.Integral())
      self.histSig.SetLineColor(2)
      
      self.histSig.Draw("same HIST")
      if self.leg is not None: self.leg.AddEntry(self.histSig, "signal (scaled)", "l")
    
    if self.leg is not None: self.leg.Draw("same")
    
    pads[ 0 ].Update()
    
    if self.doRatio:
      pads[ 1 ].cd()
      pads[ 1 ].SetGridy()
      
      setMargins(pads[ 1 ], self.doRatio)
      self.hratio.SetLineColor(1)
      
      self.hratio.Draw("e")
      
      self.hratio.SetMaximum(1 + self.ratioRange)
      self.hratio.SetMinimum(1 - self.ratioRange)
      # hratio.SetMaximum(2)
      # hratio.SetMinimum(0)
    
    for p in pads:
      p.RedrawAxis()
      p.Modified()
      p.Update()
    
    canv.cd()
    
    #iPos = 0 # in frame
    iPos = 11 # out frame
    if iPos == 0:
      self.cmsLumi.relPosX = 0.1
    self.cmsLumi.CMS_lumi(pads[ 0 ], 0, iPos)
    
    canv.Modified()
    canv.Update()
    return copy.deepcopy(canv)


