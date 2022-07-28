import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

c = 3*10**8
pionMass = 0.13957
pionTau = 2.6033*10**-8
pionbfrac = 0.999877
kaonMass = 0.493677
kaonTau = 1.238*10**-8
kaonbfrac = 0.6356
kaonpibfrac = 0.03352

with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root:WpNoMuID/DecayTree") as file:
  qcd = file.arrays(["mu_PT", "mu_ISMUON", "mu_TRUEID", "mu_ETA"], library = "pd")
with uproot.open("/tmp/13TeV_2018_34_Up_EW.root:WpNoMuID/DecayTree") as file:
  realQCD = file.arrays(["mu_PT", "mu_ISMUON", "mu_ETA"], library = "pd")
with uproot.open("/tmp/13TeV_2018_34_Up_EW.root:WpIso/DecayTree") as file:
  realW = file.arrays(["mu_PT", "mu_ISMUON"], library = "pd")

qcd.loc[:, "mu_PT"] = qcd["mu_PT"]/1000
qcd.loc[:, "Pweight"] = np.zeros(len(qcd))
qcd.loc[:, "weight"] = np.zeros(len(qcd))
qcd.loc[:, "P"] = qcd["mu_PT"]*np.cosh(qcd["mu_ETA"])

pions = qcd["mu_TRUEID"] == 211
kaons = qcd["mu_TRUEID"] == 321
protons = qcd["mu_TRUEID"] == 2212

#1/P weights for pions and kaons
qcd.loc[pions, "Pweight"] = 1/qcd.loc[pions, "P"]
qcd.loc[kaons, "Pweight"] = 1/qcd.loc[kaons, "P"]

#Equation weights for pions and kaons
qcd.loc[pions, "weight"] = pionbfrac*pionMass/(qcd.loc[pions, "P"]*pionTau)
qcd.loc[kaons, "weight"] = (kaonbfrac + kaonpibfrac)*kaonMass/(qcd.loc[kaons, "P"]*kaonTau)



#No weight plots
plt.hist(qcd.loc[pions, "mu_PT"], bins = 50, histtype = "step", label = "Pions no weight", density = True, color = "r", range = [20, 60])
plt.hist(qcd.loc[kaons, "mu_PT"], bins = 50, histtype = "step", label = "Kaons no weight", density = True, color = "b", range = [20, 60])
plt.hist(qcd.loc[protons, "mu_PT"], bins = 50, histtype = "step", label = "Protons no weight", density = True, color = "g", range = [20, 60])
#1/P weight plots
plt.hist(qcd.loc[pions, "mu_PT"], bins = 50, histtype = "step", label = "Pions 1/P weights", weights = qcd.loc[pions, "Pweight"], density = True, color = "r", linestyle = "dashed", range = [20, 60])
plt.hist(qcd.loc[kaons, "mu_PT"], bins = 50, histtype = "step", label = "Kaons 1/P weights", weights = qcd.loc[kaons, "Pweight"], density = True, color = "b", linestyle = "dashed", range = [20, 60])
#Equation plots
plt.hist(qcd.loc[pions, "mu_PT"], bins = 50, histtype = "step", label = "Pions m/(tau*p) weights", weights = qcd.loc[pions, "weight"], density = True, color = "r", linestyle = ":", range = [20, 60])
plt.hist(qcd.loc[kaons, "mu_PT"], bins = 50, histtype = "step", label = "Kaons m/(tau*p) weights", weights = qcd.loc[kaons, "weight"], density = True, color = "b", linestyle = ":", range = [20, 60])
#ISMUON = True plot
plt.hist(qcd.loc[(protons) & (qcd["mu_ISMUON"] == True), "mu_PT"], bins = 50, histtype = "step", label = "Protons with ISMUON True", density = True, color = "g", linestyle = ":", range = [20, 60])

plt.legend()
plt.xlabel("PT (GeV)")
plt.ylabel("Normalised Counts")
plt.title("Comparing weights for particle types")
plt.savefig("img/weightscomparison.png")
plt.close()

##########################
#Plotting with log scale

#No weight log plots
plt.hist(qcd.loc[pions, "mu_PT"], bins = 50, histtype = "step", label = "Pions no weight", density = True, color = "r", range = [20, 200], log = True)
plt.hist(qcd.loc[kaons, "mu_PT"], bins = 50, histtype = "step", label = "Kaons no weight", density = True, color = "b", range = [20, 200], log = True)
plt.hist(qcd.loc[protons, "mu_PT"], bins = 50, histtype = "step", label = "Protons no weight", density = True, color = "g", range = [20, 200], log = True)
#1/P weight log plots
plt.hist(qcd.loc[pions, "mu_PT"], bins = 50, histtype = "step", label = "Pions 1/P weights", weights = qcd.loc[pions, "Pweight"], density = True, color = "r", linestyle = "dashed", range = [20, 200], log = True)
plt.hist(qcd.loc[kaons, "mu_PT"], bins = 50, histtype = "step", label = "Kaons 1/P weights", weights = qcd.loc[kaons, "Pweight"], density = True, color = "b", linestyle = "dashed", range = [20, 200], log = True)
#Equation log plots
plt.hist(qcd.loc[pions, "mu_PT"], bins = 50, histtype = "step", label = "Pions m/(tau*p) weights", weights = qcd.loc[pions, "weight"], density = True, color = "r", linestyle = ":", range = [20, 200], log = True)
plt.hist(qcd.loc[kaons, "mu_PT"], bins = 50, histtype = "step", label = "Kaons m/(tau*p) weights", weights = qcd.loc[kaons, "weight"], density = True, color = "b", linestyle = ":", range = [20, 200], log = True)
#ISMUON = True log plot
plt.hist(qcd.loc[(protons) & (qcd["mu_ISMUON"] == True), "mu_PT"], bins = 50, histtype = "step", label = "Protons with ISMUON True", density = True, color = "g", linestyle = ":", range = [20, 200], log = True)

plt.legend()
plt.xlabel("PT (GeV)")
plt.ylabel("Log(Normalised Counts)")
plt.title("Comparing weights for particle types")
plt.savefig("img/weightscomparisonLog.png")
plt.close()


###############################
#Requiring ismuon == true and comparing to real QCD data

realQCD.loc[:, "mu_PT"] = realQCD["mu_PT"]/1000
realQCD.loc[:, "Pweight"] = np.zeros(len(realQCD))
realQCD.loc[:, "P"] = realQCD["mu_PT"]*np.cosh(realQCD["mu_ETA"])
realQCD["Pweight"] = 1/realQCD["P"]


plt.hist(qcd.loc[protons, "mu_PT"], bins = 50, histtype = "step", label = "Protons no weight", range = [20, 65], density = True)
plt.hist(qcd.loc[(protons) & (qcd["mu_ISMUON"] == True), "mu_PT"], bins = 50, histtype = "step", label = "Protons with ISMUON True", range = [20, 65], density = True)
plt.hist(realQCD["mu_PT"], bins = 50, histtype = "step", label = "Real QCD with 1/P weights", weights = realQCD["Pweight"], range = [20, 65], density = True)
plt.hist(realQCD.loc[realQCD["mu_ISMUON"] == True, "mu_PT"], bins = 50, histtype = "step", label = "Real QCD with ISMUON True", range = [20, 65], density = True)


plt.legend()
plt.xlabel("PT (GeV)")
plt.ylabel("Normalised Counts")
plt.title("Comparing weights protons against real QCD data")
plt.savefig("img/weightscomparisonMuonTrue.png")
plt.close()


###############################
#Comparing the QCD simulation to the real W data
realW.loc[:, "mu_PT"] = realW["mu_PT"]/1000
realW["invPT"] = 1/realW["mu_PT"]

qcd["invPT"] = 1/qcd["mu_PT"]


qcd.loc[:, "weight2m"] = np.zeros(len(qcd))
qcd.loc[:, "weight15m"] = np.zeros(len(qcd))
#2m weights for pions and kaons
qcd.loc[pions, "weight2m"] = pionbfrac*pionMass*2/(qcd.loc[pions, "P"]*pionTau*c)
qcd.loc[kaons, "weight2m"] = (kaonbfrac + kaonpibfrac)*kaonMass*2/(qcd.loc[kaons, "P"]*kaonTau*c)

#15m weights for pions and kaons
qcd.loc[pions, "weight15m"] = pionbfrac*pionMass*15/(qcd.loc[pions, "P"]*pionTau*c)
qcd.loc[kaons, "weight15m"] = (kaonbfrac + kaonpibfrac)*kaonMass*15/(qcd.loc[kaons, "P"]*kaonTau*c)

#Weights of 1 for protons
qcd.loc[protons, "weight"] = np.ones(len(qcd.loc[protons]))


histData, bins, _ = plt.hist(
  [(qcd.loc[(kaons) | (pions), "invPT"]),
  qcd.loc[(protons) & (qcd["mu_ISMUON"] == True), "invPT"]],
  bins = 50,
  stacked = True,
  label = ["Kaons and Pions d = 15m weights", "Protons with ISMUON True"],
  weights = [qcd.loc[(kaons) | (pions), "weight15m"],
    qcd.loc[(protons) & (qcd["mu_ISMUON"] == True), "weight"]],
  density = True,
  histtype = "step",
  range = [1/22, 1/20],
  color = ["r", "b"]
)

plt.hist(realW["invPT"], bins = 50, histtype = "step", label = "Real W with ISMUON True", density = True, color = "g", range = [1/22, 1/20])

histData = histData[1, :]
result = stats.linregress(bins[:-1], histData)

plt.plot(bins[:-1], result.intercept + result.slope*bins[:-1], label = "Kaons, Pions and Protons fit", color = "k", linestyle = "--")

plt.legend(loc = "lower center")
plt.xlabel("1/pT (GeV^-1)")
plt.ylabel("Normalised Counts")
plt.title("Comparing real W data to simulation")
plt.savefig("img/weightscomparisonW.png")
plt.close()


#################################
#Kaon and Pion weights with length included to compare to protons with ISMUON = True


#2m weight plots
plt.hist(
  [qcd.loc[pions, "mu_PT"],
  qcd.loc[kaons, "mu_PT"],
  qcd.loc[(protons) & (qcd["mu_ISMUON"] == True), "mu_PT"]],
  bins = 50,
  stacked = True,
  label = ["Pions d = 2m weights", "Kaons d = 2m weights", "Protons with ISMUON True"],
  weights = [
    qcd.loc[pions, "weight2m"],
    qcd.loc[kaons, "weight2m"],
    qcd.loc[(protons) & (qcd["mu_ISMUON"] == True), "weight"]],
  density = True,
  range = [20, 55],
  histtype = "step",
  color = ["limegreen", "r", "b"]
)
plt.hist(
  [qcd.loc[pions, "mu_PT"],
  qcd.loc[kaons, "mu_PT"],
  qcd.loc[(protons) & (qcd["mu_ISMUON"] == True), "mu_PT"]],
  bins = 50,
  stacked = True,
  label = ["Pions d = 15m weights", "Kaons d = 15m weights", "Protons with ISMUON True"],
  weights = [
    qcd.loc[pions, "weight15m"],
    qcd.loc[kaons, "weight15m"],
    qcd.loc[(protons) & (qcd["mu_ISMUON"] == True), "weight"]],
  density = True,
  range = [20, 55],
  color = ["k", "c", "m"],
  alpha = 0.5
)


plt.legend()
plt.xlabel("PT (GeV)")
plt.ylabel("Normalised Counts")
plt.title("Comparing protons to pions and kaons at different decay lengths")
plt.savefig("img/weightscomparisonLengths.png")
plt.close()
