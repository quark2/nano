import ROOT
import copy


f = ROOT.TFile.Open("res.root")
tree = f.Get("event")

fRes = ROOT.TFile.Open("hist.root", "recreate")

#histJEC   = ROOT.TH1F("hJEC",   "p_{T, nano} - p_{T, slimmedJetsNewJEC} for JEC; p_{T, nano} - p_{T, slimmedJetsNewJEC}; Events", 50, -0.000002, 0.000002)
#histJEC2D = ROOT.TH2F("hJEC2D", "p_{T, nano} - p_{T, slimmedJetsNewJEC} for JEC; p_{T, nano} - p_{T, slimmedJetsNewJEC}; p_{T, nano}", 50, 0, 1000, 50, -0.000002, 0.000002)

#histSmearPt   = ROOT.TH1F("hSmearPt",   "relative #Deltap_{T} for smearing; #Deltap_{T}/p_{T,PAT}; Events", 50, -0.000002, 0.000002)
#histSmearPt2D = ROOT.TH2F("hSmearPt2D", "relative #Deltap_{T} for smearing; #Deltap_{T}/p_{T,PAT}; p_{T, nano}", 50, 0, 1000, 50, -0.000002, 0.000002)

tree.Draw("jetPt - jetPtJEC >> hJEC")
tree.Draw("jetPt - jetPtJEC: jetPt >> hJEC2D")

histJEC   = ROOT.gDirectory.Get("hJEC")
histJEC2D = ROOT.gDirectory.Get("hJEC2D")

histJEC.SetTitle("p_{T, nano} - p_{T, slimmedJetsNewJEC} for JEC; p_{T, nano} - p_{T, slimmedJetsNewJEC}; Events")
histJEC2D.SetTitle("p_{T, nano} - p_{T, slimmedJetsNewJEC} for JEC; p_{T, nano} - p_{T, slimmedJetsNewJEC}; p_{T, nano}")

tree.Draw("( jetPt_Unc_an - jetPt_Unc_pp ) / jetPt_Unc_pp >> hSmearPt", "GenMatched == 1")
tree.Draw("( jetPt_Unc_an - jetPt_Unc_pp ) / jetPt_Unc_pp: jetPt_Unc_pp >> hSmearPt2D", "GenMatched == 1")

histSmearPt   = ROOT.gDirectory.Get("hSmearPt")
histSmearPt2D = ROOT.gDirectory.Get("hSmearPt2D")

histSmearPt.SetTitle("relative #Deltap_{T} for smearing; #Deltap_{T}/p_{T,PAT}; Events")
histSmearPt2D.SetTitle("relative #Deltap_{T} for smearing; #Deltap_{T}/p_{T,PAT}; p_{T, nano}")

histJEC.SetBinContent(1,  histJEC.GetBinContent(0)  + histJEC.GetBinContent(-1))
histJEC.SetBinContent(50, histJEC.GetBinContent(50) + histJEC.GetBinContent(51))

histSmearPt.SetBinContent(1,  histSmearPt.GetBinContent(0)  + histSmearPt.GetBinContent(-1))
histSmearPt.SetBinContent(50, histSmearPt.GetBinContent(50) + histSmearPt.GetBinContent(51))

histJEC.Write()
histJEC2D.Write()

histSmearPt.Write()
histSmearPt2D.Write()

fRes.Write()


