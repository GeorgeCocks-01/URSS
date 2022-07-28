import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import stats

with uproot.open("/tmp/13TeV_2018_34_Up_EW.root:WpNoMuID/DecayTree") as file:
  realQCD = file.arrays(["mu_PT", "mu_TRCHI2DOF", "mu_ISMUON", "mu_ETA"], library = "pd")

realQCD["mu_PT"] = realQCD["mu_PT"]/1000
muonTrue = realQCD["mu_ISMUON"] == True
realQCD.loc[:, "P"] = realQCD["mu_PT"]*np.cosh(realQCD["mu_ETA"])

#2D histogram comparing TRCHI2 to pT for any ISMUON
plt.hist2d(realQCD["P"], realQCD["mu_TRCHI2DOF"], bins = 200, norm = LogNorm(), range = [[40, 1500], [0, 4]])

average = stats.linregress(realQCD["P"], realQCD["mu_TRCHI2DOF"])
plt.plot(realQCD["P"], average.intercept + average.slope*realQCD["P"], label = "Average TRCHI2", color = "r")

plt.colorbar()
plt.ylabel("TRCHI2")
plt.xlabel("P (GeV)")
plt.title("Comparing TRCHI2 to pT for real QCD background")
plt.savefig("img/TRCHI2vPTreal.png")
plt.close()


#2D histogram comparing TRCHI2 to pT for ISMUON == True
plt.hist2d(realQCD.loc[muonTrue, "P"], realQCD.loc[muonTrue, "mu_TRCHI2DOF"], bins = 200, norm = LogNorm(), range = [[40, 800], [0, 2.5]])
plt.colorbar()

averageTrue = stats.linregress(realQCD.loc[muonTrue, "P"], realQCD.loc[muonTrue, "mu_TRCHI2DOF"])
plt.plot(realQCD.loc[muonTrue, "P"], averageTrue.intercept + averageTrue.slope*realQCD.loc[muonTrue, "P"], label = "Average TRCHI2", color = "r")

plt.ylabel("TRCHI2")
plt.xlabel("P (GeV)")
plt.title("Comparing TRCHI2 to pT for ISMUON = True real QCD background")
plt.savefig("img/TRCHI2vPT_TRUEreal.png")
plt.close()
