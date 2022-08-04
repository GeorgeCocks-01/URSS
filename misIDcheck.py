import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit

c = 3*10**8
pionMass = 0.13957
pionTau = 2.6033*10**-8
pionbfrac = 0.999877
kaonMass = 0.493677
kaonTau = 1.238*10**-8
kaonbfrac = 0.6356

def getFracErr(data, data_misID, rangeMax):
  hist, bins = np.histogram(data, bins = 50, range = [55, rangeMax])
  frac = np.histogram(data_misID, bins = 50, range = [55, rangeMax])[0]/hist
  err = np.sqrt((frac*(1 - frac))/hist)
  return np.array([frac, err]).T, bins

def integral(p, length, mass, tau, bfrac):
  return 100*bfrac*(1 - np.exp(-(length*mass)/(c*p*tau)))

def compositeIntegral(p, prefac1, prefac2, length = None, mass = None, tau = None, bfrac = None):
  if mass == None or tau == None or bfrac == None or length == None:
    return prefac1*p + prefac2
  else:
    return 100*bfrac*(1 - np.exp(-(length*mass)/(c*p*tau))) + prefac1*p + prefac2

with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root:WpNoMuID/DecayTree") as file:
  qcd = file.arrays(["mu_PT", "mu_ISMUON", "mu_TRUEID", "mu_ETA"], library = "pd")

#Fit the QCD simulation data to the theoretical pion and kaon plots (above) to give the length
qcd["mu_PT"] = qcd["mu_PT"]/1000
qcd = qcd.loc[qcd["mu_PT"] > 20]
qcd["P"] = qcd["mu_PT"]*np.cosh(qcd["mu_ETA"])
kaons = qcd["mu_TRUEID"] == 321
pions = qcd["mu_TRUEID"] == 211
kaons_misID = (kaons & (qcd["mu_ISMUON"] == True))
pions_misID = (pions & (qcd["mu_ISMUON"] == True))

#Create dataframe for data to plot
plottingData = pd.DataFrame(np.zeros((50, 4)), columns = ["kaons_frac", "pions_frac", "kaons_err", "pions_err"], dtype = float)

plottingData[["kaons_frac", "kaons_err"]], kaons_bins = getFracErr(qcd.loc[kaons, "P"], qcd.loc[kaons_misID, "P"], 200)
plottingData[["pions_frac", "pions_err"]], pions_bins = getFracErr(qcd.loc[pions, "P"], qcd.loc[pions_misID, "P"], 200)

plt.errorbar(x = kaons_bins[:-1], y = plottingData["kaons_frac"]*100, label = "Kaon data", yerr = plottingData["kaons_err"]*100)
plt.errorbar(x = pions_bins[:-1], y = plottingData["pions_frac"]*100, label = "Pion data", yerr = plottingData["pions_err"]*100)

#Finding the fit parameters from the data
kaonLength, _ = curve_fit(lambda p, length: integral(p, length, kaonMass, kaonTau, kaonbfrac), kaons_bins[:-1], plottingData["kaons_frac"]*100)
pionLength, _ = curve_fit(lambda p, length: integral(p, length, pionMass, pionTau, pionbfrac), pions_bins[:-1], plottingData["pions_frac"]*100)

print("Kaon fit length:", round(kaonLength[0], 2))
print("Pion fit length:", round(pionLength[0], 2))

pRange = np.linspace(55, 200, 1000)

plt.plot(pRange, integral(pRange, kaonLength, kaonMass, kaonTau, kaonbfrac), label = "Kaon Fit")
plt.plot(pRange, integral(pRange, pionLength, pionMass, pionTau, pionbfrac), label = "Pion Fit")
plt.legend()
plt.xlabel("P (GeV)")
plt.ylabel("Percent Mislabeled")
plt.title("Fitting QCD simulation misID data to get length")
plt.savefig("img/misIDfit.png")
plt.close()


######################################
#Getting data for protons
protons = qcd["mu_TRUEID"] == 2212
protons_misID = (protons & (qcd["mu_ISMUON"] == True))

#Fitting the data on the log scale using the composite function of decay + log terms
plottingData[["kaons_frac", "kaons_err"]], kaons_bins = getFracErr(qcd.loc[kaons, "P"], qcd.loc[kaons_misID, "P"], 1000)
plottingData[["pions_frac", "pions_err"]], pions_bins = getFracErr(qcd.loc[pions, "P"], qcd.loc[pions_misID, "P"], 1000)
plottingData[["protons_frac", "protons_err"]], protons_bins = getFracErr(qcd.loc[protons, "P"], qcd.loc[protons_misID, "P"], 1000)

kaonParams, kaonCov = curve_fit(lambda p, kaonPrefac1, kaonPrefac2, kaonLength: compositeIntegral(p, kaonPrefac1, kaonPrefac2, kaonLength, kaonMass, kaonTau, kaonbfrac), kaons_bins[:-1], plottingData["kaons_frac"]*100, bounds = ((-np.inf, -np.inf, -np.inf), (np.inf, np.inf, 100)))

pionParams, pionCov = curve_fit(lambda p, pionPrefac1, pionPrefac2, pionLength: compositeIntegral(p, pionPrefac1, pionPrefac2, pionLength, pionMass, pionTau, pionbfrac), pions_bins[:-1], plottingData["pions_frac"]*100, bounds = ((-np.inf, -np.inf, -np.inf), (np.inf, np.inf, 100)))

protonParams, protonCov = curve_fit(lambda p, protonPrefac1, protonPrefac2: compositeIntegral(p, protonPrefac1, protonPrefac2), protons_bins[:-1], plottingData["protons_frac"]*100)

# print("kaon parameters:", kaonParams)
# print("pion parameters:", pionParams)
# print("proton parameters:", protonParams)

pRange = np.linspace(55, 1000, 10000)
plt.errorbar(x = kaons_bins[:-1], y = plottingData["kaons_frac"]*100, label = "Kaon data", yerr = plottingData["kaons_err"]*100)
plt.errorbar(x = pions_bins[:-1], y = plottingData["pions_frac"]*100, label = "Pion data", yerr = plottingData["pions_err"]*100)
plt.errorbar(x = protons_bins[:-1], y = plottingData["protons_frac"]*100, label = "Proton data", yerr = plottingData["protons_err"]*100)
plt.xscale("log")
plt.plot(pRange, compositeIntegral(pRange, *kaonParams, kaonMass, kaonTau, kaonbfrac), label = "Kaon Fit")
plt.plot(pRange, compositeIntegral(pRange, *pionParams, pionMass, pionTau, pionbfrac), label = "Pion Fit")
plt.plot(pRange, compositeIntegral(pRange, *protonParams), label = "Proton Fit", color = "k")

plt.legend(loc = "upper center")
plt.ylabel("Percent mislabeled")
plt.xlabel("P (GeV)")
plt.title("Fitting QCD simulation misID data with punch through to get length")
plt.savefig("img/misIDcompositefit.png")
plt.close()


################################################
#Comparing 1/P weights to weights with integrals
qcd.loc[kaons, "weights"] = compositeIntegral(qcd.loc[kaons, "mu_PT"], *kaonParams, kaonMass, kaonTau, kaonbfrac)
qcd.loc[pions, "weights"] = compositeIntegral(qcd.loc[pions, "mu_PT"], *pionParams, pionMass, pionTau, pionbfrac)
qcd.loc[protons, "weights"] = compositeIntegral(qcd.loc[protons, "mu_PT"], *protonParams)

plt.hist(qcd.loc[kaons, "mu_PT"], bins = 50, label = "Kaons with 1/P weights", weights = 1/qcd.loc[kaons, "P"], range = [20, 70], histtype = "step", density = True, color = "b")
plt.hist(qcd.loc[pions, "mu_PT"], bins = 50, label = "Pions with 1/P weights", weights = 1/qcd.loc[pions, "P"], range = [20, 70], histtype = "step", density = True, color = "r")
plt.hist(qcd.loc[protons, "mu_PT"], bins = 50, label = "Protons with no weight", range = [20, 70], histtype = "step", density = True, color = "g")
plt.hist(qcd.loc[kaons, "mu_PT"], bins = 50, label = "Kaons with Composite weights", weights = qcd.loc[kaons, "weights"], range = [20, 70], histtype = "step", density = True, linestyle = "--", color = "b")
plt.hist(qcd.loc[pions, "mu_PT"], bins = 50, label = "Pions with Composite weights", weights = qcd.loc[pions, "weights"], range = [20, 70], histtype = "step", density = True, linestyle = "--", color = "r")
plt.hist(qcd.loc[protons, "mu_PT"], bins = 50, label = "Protons with Composite weights", weights = qcd.loc[protons, "weights"], range = [20, 70], histtype = "step", density = True, linestyle = "--", color = "g")

plt.legend()
plt.ylabel("Normalised Counts")
plt.xlabel("pT (GeV)")
plt.title("Comparing 1/P weights to Composite function weights")
plt.savefig("img/misIDWeightsComparison.png")
plt.close()
