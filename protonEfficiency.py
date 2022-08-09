import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def f(x, a, b, c, d):
  return a*x**3 + b*x**2 + c*x + d

with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root:WpNoMuID/DecayTree") as file:
  qcd = file.arrays(["mu_PT", "mu_ISMUON", "mu_TRUEID", "mu_ETA"], library = "pd")


qcd["mu_PT"] = qcd["mu_PT"]/1000
protons = qcd["mu_TRUEID"] == 2212


trueCounts, trueBins, _ = plt.hist(qcd.loc[(protons) & (qcd["mu_ISMUON"] == True), "mu_PT"], bins = 15, range = [20, 70], label = "Protons with ISMUON = True", histtype = "step")
allCounts, allBins, _ = plt.hist(qcd.loc[protons, "mu_PT"], bins = 15, range = [20, 70], label = "All protons", histtype = "step")

plt.legend()
plt.xlabel("pT (GeV)")
plt.ylabel("Normalised Counts")
plt.title("QCD simulation comparison of all protons to those with ISMUON = True")
plt.savefig("img/protonEfficiency/protonComparison.png")
plt.close()

efficiency = trueCounts/allCounts


w = abs(trueBins[1:] - trueBins[:-1] - 0.2)
c = (trueBins[1:] + trueBins[:-1])/2
params, _ = curve_fit(f, c, efficiency*100)

plt.bar(c, efficiency*100, align = "center", width = w, label = "Efficiency data")
plt.plot(c, f(c, *params), color = "r", label = "Efficiency fit")
plt.legend()
plt.xlabel("pT (GeV)")
plt.ylabel("Efficiency (%)")
plt.title("Proton Efficiency")
plt.savefig("img/protonEfficiency/protonEfficiency.png")
plt.close()

