import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import stats

with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root:WpNoMuID/DecayTree") as file:
  simQCD = file.arrays(["mu_PT", "mu_TRCHI2DOF", "mu_ISMUON", "mu_ETA"], library = "pd")

simQCD["mu_PT"] = simQCD["mu_PT"]/1000
muonTrue = simQCD["mu_ISMUON"] == True
simQCD.loc[:, "P"] = simQCD["mu_PT"]*np.cosh(simQCD["mu_ETA"])

#2D histogram comparing TRCHI2 to pT for any ISMUON
plt.hist2d(simQCD["P"], simQCD["mu_TRCHI2DOF"], bins = 200, norm = LogNorm(), range = [[40, 2500], [0, 4]])

average = stats.linregress(simQCD["P"], simQCD["mu_TRCHI2DOF"])
plt.plot(simQCD["P"], average.intercept + average.slope*simQCD["P"], label = "Average TRCHI2", color = "r")

plt.colorbar()
plt.ylabel("TRCHI2")
plt.xlabel("P (GeV)")
plt.title("Comparing TRCHI2 to pT for simulated QCD background")
plt.savefig("img/TRCHI2vPTsim.png")
plt.close()


#2D histogram comparing TRCHI2 to pT for ISMUON == True
plt.hist2d(simQCD.loc[muonTrue, "P"], simQCD.loc[muonTrue, "mu_TRCHI2DOF"], bins = 200, norm = LogNorm(), range = [[40, 1000], [0, 2.5]])
plt.colorbar()

averageTrue = stats.linregress(simQCD.loc[muonTrue, "P"], simQCD.loc[muonTrue, "mu_TRCHI2DOF"])
plt.plot(simQCD.loc[muonTrue, "P"], averageTrue.intercept + averageTrue.slope*simQCD.loc[muonTrue, "P"], label = "Average TRCHI2", color = "r")

plt.ylabel("TRCHI2")
plt.xlabel("P (GeV)")
plt.title("Comparing TRCHI2 to pT for ISMUON = True simulated QCD background")
plt.savefig("img/TRCHI2vPT_TRUEsim.png")
plt.close()
