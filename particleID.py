import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plotter(data, xMin, xMax, xVar, units = ""):
  fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
  fig.suptitle(xVar + " for different particles")

  for i, j in particleDict.items():
    particleData = data.loc[simQCD["mu_TRUEID"] == i]
    ax1.hist(particleData, bins = 100, histtype = "step", label = j, range = [xMin, xMax])
    ax2.hist(particleData, bins = 100, histtype = "step", label = j, range = [xMin, xMax], density = True)
  ax1.set(xlabel = xVar + " " + units, ylabel = "Counts")
  ax2.set(xlabel = xVar + " " + units, ylabel = "Normalised Counts")
  ax1.legend(loc = "upper right", fontsize = "x-small")
  plt.savefig("img/particles" + xVar)

  plt.close()


with uproot.open("/storage/epp2/phshgg/Public/DVTuples__v24g/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root:WpNoMuID/DecayTree") as file:
  simQCD = file.arrays(["mu_PTSUMCONE040", "mu_MIPCHI2DV", "mu_TRCHI2DOF", "mu_TRUEID", "mu_ISMUON", "mu_ETA", "mu_PT", "mu_HcalE", "mu_EcalE", "mu_PHI"], library = "pd")

simQCD = simQCD.loc[simQCD["mu_ISMUON"] == False]
simQCD[["mu_PTSUMCONE040", "mu_PT", "mu_HcalE", "mu_EcalE"]] = simQCD[["mu_PTSUMCONE040", "mu_PT", "mu_HcalE", "mu_EcalE"]]/1000
simQCD.loc[:, "HcalET"] = simQCD["mu_HcalE"]/np.cosh(simQCD["mu_ETA"])
simQCD.loc[:, "EcalET"] = simQCD["mu_EcalE"]/np.cosh(simQCD["mu_ETA"])


particleDict = {
  2212: "p",
  321: "K+",
  211: "pi+",
  0: "Type 0",
  -11: "e+",
  -13: "mu+",
  # -15: "tau+",
  # -211: "pi-",
  # -321: "K-",
  # -2212: "anti-p",
}


#Plotting isolation graph
plotter(simQCD["mu_PTSUMCONE040"], 0, 125, "Isolation", "(GeV)")
#Plotting log10(Isolation)
plotter(np.log10(simQCD["mu_PTSUMCONE040"]), -0.5, 2.5, "log(Isolation)")

#Plotting ipchi2 graph
plotter(simQCD["mu_MIPCHI2DV"], 0, 15, "IPCHI2")
#Plotting log10(IPCHI2)
plotter(np.log10(simQCD["mu_MIPCHI2DV"]), -3, 4, "log(IPCHI2)")

#Plotting trchi2 graph
plotter(simQCD["mu_TRCHI2DOF"], 0, 3, "TRCHI2")

#Plotting ETA
plotter(simQCD["mu_ETA"], 1.5, 5, "ETA")

#Plotting PT
plotter(simQCD["mu_PT"], 20, 60, "PT", "(GeV)")

#Plotting Hcal ET
plotter(simQCD["HcalET"], 0, 60, "HcalET", "(GeV)")

#Plotting Ecal ET
plotter(simQCD["EcalET"], 0, 15, "EcalET", "(GeV)")

#Plotting PHI
plotter(simQCD["mu_PHI"], -4, 4, "PHI")
