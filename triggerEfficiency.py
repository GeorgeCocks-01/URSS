import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plotter(key, xMin, xMax, xVar, units = ""):
  withoutTrigger, bin_edges = np.histogram(simQCD.loc[:, key], range = [xMin, xMax], bins = 50)
  withTrigger, bin_edges = np.histogram(simQCDtrigger.loc[:, key], bins = bin_edges)

  efficiency = 100*withTrigger/withoutTrigger
  w = abs(bin_edges[1:] - bin_edges[:-1])
  plt.bar(bin_edges[:-1], efficiency, width = w, align = "edge")
  plt.title("Trigger efficiency of QCD simulation data")
  plt.xlabel(xVar + " " + units)
  plt.ylabel("Trigger Efficiency (%)")
  plt.savefig("img/trigger" + xVar)
  plt.close()

with uproot.open("/storage/epp2/phshgg/Public/DVTuples__v24g/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root:WpNoMuID/DecayTree") as file:
  simQCD = file.arrays(["mu_PTSUMCONE040", "mu_MIPCHI2DV", "mu_TRCHI2DOF", "mu_TRUEID", "mu_ISMUON", "mu_ETA", "mu_PT", "mu_HcalE", "mu_EcalE", "mu_PHI", "mu_L0Global_TOS"], library = "pd")

simQCD[["mu_PTSUMCONE040", "mu_PT", "mu_HcalE", "mu_EcalE"]] = simQCD[["mu_PTSUMCONE040", "mu_PT", "mu_HcalE", "mu_EcalE"]]/1000
simQCD.loc[:, "HcalET"] = simQCD["mu_HcalE"]/np.cosh(simQCD["mu_ETA"])
simQCD.loc[:, "EcalET"] = simQCD["mu_EcalE"]/np.cosh(simQCD["mu_ETA"])


simQCDtrigger = simQCD[simQCD["mu_L0Global_TOS"] == True]

#Plotting pT
plotter("mu_PT", 20, 60, "pT", "(GeV)")

#Plotting ETA
plotter("mu_ETA", 1.5, 5, "ETA")

#Plotting isolation graph
plotter("mu_PTSUMCONE040", 0, 125, "Isolation", "(GeV)")

#Plotting ipchi2 graph
plotter("mu_MIPCHI2DV", 0, 15, "IPCHI2")

#Plotting trchi2 graph
plotter("mu_TRCHI2DOF", 0, 3, "TRCHI2")

#Plotting Hcal ET
plotter("HcalET", 0, 60, "HcalET", "(GeV)")

#Plotting Ecal ET
plotter("EcalET", 0, 10, "EcalETasdfasdfadfa", "(GeV)")

#Plotting PHI
plotter("mu_PHI", -4, 4, "PHI")
