import ROOT
import glob

chain = ROOT.TChain("h2000")

files = glob.glob("../data/Kpi/mc/*0.root")
print(files)

for file in files:
    chain.Add(file)

c = ROOT.TCanvas("myCanvasName","The Canvas Title",800,600)

chain.Draw("thetab", "evmcflag!=1")
c.Draw()

