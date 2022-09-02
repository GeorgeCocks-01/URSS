import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

c = 3*10**8
pionMass = 0.13957
pionTau = 2.6033*10**-8
pionbfrac = 0.999877
kaonMass = 0.493677
kaonTau = 1.238*10**-8
kaonbfrac = 0.6356

def getFracErr(data, data_misID, rangeMin, rangeMax):
  hist, bins = np.histogram(data, bins = 50, range = [rangeMin, rangeMax])
  frac = np.histogram(data_misID, bins = 50, range = [rangeMin, rangeMax])[0]/hist
  err = np.sqrt((frac*(1 - frac))/hist)
  frac = np.nan_to_num(frac, nan=0, posinf=0, neginf=0)
  return np.array([frac, err]).T, bins


def compositeIntegral(p, prefac1, prefac2, length = None, mass = None, tau = None, bfrac = None):
  if mass == None or tau == None or bfrac == None or length == None:
    return prefac1*p + prefac2
  else:
    return 100*bfrac*(1 - np.exp(-(length*mass)/(c*p*tau))) + prefac1*p + prefac2


with uproot.open("../gcocks-background-simulation/mW/data/tuples/Wp_QcdBgdPt18GeV_13TeV_SmearingAll_AlignCorrOn_InitialMomScaleCorrOn.root:DecayTree") as file:
  qcd = file.arrays(["mu_PT", "p", "mu_TRUEID", "mu_ISMUON"], aliases={"mu_PT": "mu_pt", "p": "mu_P", "mu_TRUEID": "mu_trueID", "mu_ISMUON": "mu_realIsMuon"}, library = "pd")
with uproot.open("../gcocks-background-simulation/mW/data/tuples/Wm_QcdBgdPt18GeV_13TeV_SmearingAll_AlignCorrOn_InitialMomScaleCorrOn.root:DecayTree") as file:
  qcd2 = file.arrays(["p", "mu_TRUEID", "mu_ISMUON"], aliases={"p": "mu_P", "mu_TRUEID": "mu_trueID", "mu_ISMUON": "mu_realIsMuon"}, library = "pd")

qcd = qcd.append(qcd2, ignore_index = True)
rangeMin = qcd["p"].min()
rangeMax = qcd["p"].max()
print("range: ", rangeMin, "to ", rangeMax)

kaons = qcd["mu_TRUEID"] == 321
pions = qcd["mu_TRUEID"] == 211
kaons_misID = (kaons & (qcd["mu_ISMUON"] == True))
pions_misID = (pions & (qcd["mu_ISMUON"] == True))

#Create dataframe for data to plot
plottingData = pd.DataFrame(np.zeros((50, 4)), columns = ["kaons_frac", "pions_frac", "kaons_err", "pions_err"], dtype = float)

#Getting data for protons
protons = qcd["mu_TRUEID"] == 2212
protons_misID = (protons & (qcd["mu_ISMUON"] == True))

#Fitting the data on the log scale using the composite function of decay + log terms
plottingData[["kaons_frac", "kaons_err"]], kaons_bins = getFracErr(qcd.loc[kaons, "p"], qcd.loc[kaons_misID, "p"], rangeMin, rangeMax)
plottingData[["pions_frac", "pions_err"]], pions_bins = getFracErr(qcd.loc[pions, "p"], qcd.loc[pions_misID, "p"], rangeMin, rangeMax)
plottingData[["protons_frac", "protons_err"]], protons_bins = getFracErr(qcd.loc[protons, "p"], qcd.loc[protons_misID, "p"], rangeMin, rangeMax)

kaonParams, kaonCov = curve_fit(lambda p, kaonPrefac1, kaonPrefac2, kaonLength: compositeIntegral(p, kaonPrefac1, kaonPrefac2, kaonLength, kaonMass, kaonTau, kaonbfrac), kaons_bins[:-1], plottingData["kaons_frac"]*100, bounds = ((-np.inf, -np.inf, -np.inf), (np.inf, np.inf, 100)))

pionParams, pionCov = curve_fit(lambda p, pionPrefac1, pionPrefac2, pionLength: compositeIntegral(p, pionPrefac1, pionPrefac2, pionLength, pionMass, pionTau, pionbfrac), pions_bins[:-1], plottingData["pions_frac"]*100, bounds = ((-np.inf, -np.inf, -np.inf), (np.inf, np.inf, 100)))

protonParams, protonCov = curve_fit(lambda p, protonPrefac1, protonPrefac2: compositeIntegral(p, protonPrefac1, protonPrefac2), protons_bins[:-1], plottingData["protons_frac"]*100)

print("kaon parameters:", kaonParams)
print("pion parameters:", pionParams)
print("proton parameters:", protonParams)

pRange = np.linspace(rangeMin, rangeMax, 10000)
plt.errorbar(x = kaons_bins[:-1], y = plottingData["kaons_frac"]*100, label = "Kaon data", yerr = plottingData["kaons_err"]*100)
plt.errorbar(x = pions_bins[:-1], y = plottingData["pions_frac"]*100, label = "Pion data", yerr = plottingData["pions_err"]*100)
plt.errorbar(x = protons_bins[:-1], y = plottingData["protons_frac"]*100, label = "Proton data", yerr = plottingData["protons_err"]*100)
plt.xscale("log")

plt.plot(pRange, compositeIntegral(pRange, *kaonParams, kaonMass, kaonTau, kaonbfrac), label = "Kaon Fit")
plt.plot(pRange, compositeIntegral(pRange, *pionParams, pionMass, pionTau, pionbfrac), label = "Pion Fit")
plt.plot(pRange, compositeIntegral(pRange, *protonParams), label = "Proton Fit", color = "k")

plt.legend(loc = "upper center")
plt.ylabel("Percent mislabeled")
plt.xlabel("p (GeV)")
plt.title("Fitting QCD simulation misID data with punch through to get length")
plt.savefig("img/test/misIDcheck.png")
plt.close()
