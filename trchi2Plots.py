import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import stats
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000


def plotter(data, xType, yType, muonTF, yLim):
  muon = data["mu_ISMUON"] == muonTF
  plt.hist2d(data.loc[muon, str(xType)], data.loc[muon, str(yType)], bins = 150, norm = LogNorm(), range = [[40, 1400], [0, yLim]])

  average = stats.linregress(data.loc[muon, str(xType)], data.loc[muon, str(yType)])
  plt.plot(data.loc[muon, str(xType)], average.intercept + average.slope*data.loc[muon, str(xType)], color = "r")
  plt.colorbar()
  if str(yType) == "mu_PTSUMCONE040":
    plt.ylabel("Isolation (GeV)")
  else:
    plt.ylabel("TRCHI2")
  plt.xlabel(str(xType) + " (GeV)")
  return average.intercept, average.slope


with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root:WpNoMuID/DecayTree") as file:
  simQCD = file.arrays(["mu_PT", "mu_PTSUMCONE040", "mu_TRCHI2DOF", "mu_ISMUON", "mu_ETA"], library = "pd")
with uproot.open("/tmp/13TeV_2018_34_Up_EW.root:WpNoMuID/DecayTree") as file:
  realQCD = file.arrays(["mu_PT","mu_PTSUMCONE040", "mu_TRCHI2DOF", "mu_ISMUON", "mu_ETA"], library = "pd")

simQCD["mu_PT"] = simQCD["mu_PT"]/1000
realQCD["mu_PT"] = realQCD["mu_PT"]/1000
simQCD["mu_PTSUMCONE040"] = simQCD["mu_PTSUMCONE040"]/1000
realQCD["mu_PTSUMCONE040"] = realQCD["mu_PTSUMCONE040"]/1000
simQCD.loc[:, "P"] = simQCD["mu_PT"]*np.cosh(simQCD["mu_ETA"])
realQCD.loc[:, "P"] = realQCD["mu_PT"]*np.cosh(realQCD["mu_ETA"])
x = np.linspace(40, 1000, 10000)


#####REAL DATA#####
#ISMUON FALSE, PT
falseInterceptPT, falseSlopePT = plotter(realQCD, "P", "mu_TRCHI2DOF", False, 4)
plt.title("Real QCD background with ISMUON = False")
plt.savefig("img/TRCHI2vPreal.png")
plt.close()

#ISMUON TRUE, PT
trueInterceptPT, trueSlopePT = plotter(realQCD, "P", "mu_TRCHI2DOF", True, 4)
plt.title("Real QCD background with ISMUON = True")
plt.savefig("img/TRCHI2vP_TRUEreal.png")
plt.close()

#ISMUON FALSE, ISO
falseInterceptISO, falseSlopeISO = plotter(realQCD, "P", "mu_PTSUMCONE040", False, 200)
plt.title("Real QCD background with ISMUON = False")
plt.savefig("img/ISOvPreal.png")
plt.close()

#ISMUON TRUE, ISO
trueInterceptISO, trueSlopeISO = plotter(realQCD, "P", "mu_PTSUMCONE040", True, 200)
plt.title("Real QCD background with ISMUON = True")
plt.savefig("img/ISOvP_TRUEreal.png")
plt.close()
############
#Comparing the average TRCHI2 fit of ISMUON = True AND ISMUON = False, PT
ratio = (falseInterceptPT + falseSlopePT*x)/(trueInterceptPT + trueSlopePT*x)
plt.plot(x, ratio)
plt.xlabel("pT (GeV)")
plt.title("Ratio of ISMUON False to ISMUON True average for real data")
plt.savefig("img/trchi2PRatioreal.png")
plt.close()

#Comparing the average TRCHI2 fit of ISMUON = True AND ISMUON = False, ISO
ratio = (falseInterceptISO + falseSlopePT*x)/(trueInterceptISO + trueSlopeISO*x)
plt.plot(x, ratio)
plt.xlabel("Isolation (GeV)")
plt.title("Ratio of ISMUON False to ISMUON True average for real data")
plt.savefig("img/ISOPRatioreal.png")
plt.close()




#####SIM DATA#####
#ISMUON FALSE, PT
falseInterceptPT, falseSlopePT = plotter(simQCD, "P", "mu_TRCHI2DOF", False, 4)
plt.title("Simulated QCD background with ISMUON = False")
plt.savefig("img/TRCHI2vPsim.png")
plt.close()

#ISMUON TRUE, PT
trueInterceptPT, trueSlopePT = plotter(simQCD, "P", "mu_TRCHI2DOF", True, 4)
plt.title("Simulated QCD background with ISMUON = True")
plt.savefig("img/TRCHI2vP_TRUEsim.png")
plt.close()

#ISMUON FALSE, ISO
falseInterceptISO, falseSlopeISO = plotter(simQCD, "P", "mu_PTSUMCONE040", False, 200)
plt.title("Simulated QCD background with ISMUON = False")
plt.savefig("img/ISOvPsim.png")
plt.close()

#ISMUON TRUE, ISO
trueInterceptISO, trueSlopeISO = plotter(simQCD, "P", "mu_PTSUMCONE040", True, 200)
plt.title("Simulated QCD background with ISMUON = True")
plt.savefig("img/ISOvP_TRUEsim.png")
plt.close()
###############
#Comparing the average TRCHI2 fit of ISMUON = True AND ISMUON = False, PT
ratio = (falseInterceptPT + falseSlopePT*x)/(trueInterceptPT + trueSlopePT*x)
plt.plot(x, ratio)
plt.xlabel("pT (GeV)")
plt.title("Ratio of ISMUON False to ISMUON True average for sim data")
plt.savefig("img/trchi2PRatiosim.png")
plt.close()

#Comparing the average TRCHI2 fit of ISMUON = True AND ISMUON = False, ISO
ratio = (falseInterceptISO + falseSlopeISO*x)/(trueInterceptISO + trueSlopeISO*x)
plt.plot(x, ratio)
plt.xlabel("Isolation (GeV)")
plt.title("Ratio of ISMUON False to ISMUON True average for sim data")
plt.savefig("img/ISOPRatiosim.png")
plt.close()
