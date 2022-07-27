import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
  realW = file.arrays(["mu_PT", "mu_ISMUON", "mu_ETA"], library = "pd")


qcd.loc[:, "mu_PT"] = qcd["mu_PT"]/1000
qcd.loc[:, "Pweight"] = np.zeros(len(qcd))
qcd.loc[:, "weight"] = np.zeros(len(qcd))
qcd.loc[:, "P"] = qcd["mu_PT"]*np.cosh(qcd["mu_ETA"])

#1/P weights for pions and kaons
qcd.loc[qcd["mu_TRUEID"] == 211, "Pweight"] = 1/qcd.loc[qcd["mu_TRUEID"] == 211, "P"]
qcd.loc[qcd["mu_TRUEID"] == 321, "Pweight"] = 1/qcd.loc[qcd["mu_TRUEID"] == 321, "P"]

#Equation weights for pions and kaons
qcd.loc[qcd["mu_TRUEID"] == 211, "weight"] = pionbfrac*pionMass/(qcd.loc[qcd["mu_TRUEID"] == 211, "P"]*pionTau)
qcd.loc[qcd["mu_TRUEID"] == 321, "weight"] = (kaonbfrac + kaonpibfrac)*kaonMass/(qcd.loc[qcd["mu_TRUEID"] == 321, "P"]*kaonTau)



#No weight plots
plt.hist(qcd.loc[qcd["mu_TRUEID"] == 211, "mu_PT"], bins = 50, histtype = "step", label = "Pions no weight", range = [20, 60], density = True)
plt.hist(qcd.loc[qcd["mu_TRUEID"] == 321, "mu_PT"], bins = 50, histtype = "step", label = "Kaons no weight", range = [20, 60], density = True)
plt.hist(qcd.loc[qcd["mu_TRUEID"] == 2212, "mu_PT"], bins = 50, histtype = "step", label = "Protons no weight", range = [20, 60], density = True)
# #1/P weight plots
plt.hist(qcd.loc[qcd["mu_TRUEID"] == 211, "mu_PT"], bins = 50, histtype = "step", label = "Pions 1/P weights", weights = qcd.loc[qcd["mu_TRUEID"] == 211, "Pweight"], range = [20, 60], density = True)
plt.hist(qcd.loc[qcd["mu_TRUEID"] == 321, "mu_PT"], bins = 50, histtype = "step", label = "Kaons 1/P weights", weights = qcd.loc[qcd["mu_TRUEID"] == 321, "Pweight"], range = [20, 60], density = True)
#Equation plots
plt.hist(qcd.loc[qcd["mu_TRUEID"] == 211, "mu_PT"], bins = 50, histtype = "step", label = "Pions m/tau*p weights", weights = qcd.loc[qcd["mu_TRUEID"] == 211, "weight"], range = [20, 60], density = True)
plt.hist(qcd.loc[qcd["mu_TRUEID"] == 321, "mu_PT"], bins = 50, histtype = "step", label = "Kaons m/tau*p weights", weights = qcd.loc[qcd["mu_TRUEID"] == 321, "weight"], range = [20, 60], density = True)


plt.legend()
plt.xlabel("PT (GeV)")
plt.ylabel("Normalised Counts")
plt.title("Comparing weights for particle types")
plt.savefig("img/weightscomparison.png")
plt.close()


###############################
#Requiring ismuon == true and comparing to real QCD data

qcdT = qcd.loc[qcd["mu_ISMUON"] == True]
realQCD.loc[:, "mu_PT"] = realQCD["mu_PT"]/1000
realQCD.loc[:, "Pweight"] = np.zeros(len(realQCD))
realQCD.loc[:, "P"] = realQCD["mu_PT"]*np.cosh(realQCD["mu_ETA"])
realQCD["Pweight"] = 1/realQCD["P"]
realQCDT = realQCD.loc[realQCD["mu_ISMUON"] == True]


plt.hist(qcdT.loc[qcdT["mu_TRUEID"] == 2212, "mu_PT"], bins = 50, histtype = "step", label = "Protons no weight", range = [20, 75], density = True)
plt.hist(realQCD["mu_PT"], bins = 50, histtype = "step", label = "Real QCD with 1/P weights", weights = realQCD["Pweight"], range = [20, 75], density = True)
plt.hist(realQCDT["mu_PT"], bins = 50, histtype = "step", label = "Real QCD with ISMUON True", range = [20, 75], density = True)


plt.legend()
plt.xlabel("PT (GeV)")
plt.ylabel("Normalised Counts")
plt.title("Comparing weights protons against real QCD data")
plt.savefig("img/weightscomparisonMuonTrue.png")
plt.close()



# realW.loc[:, "mu_PT"] = realW["mu_PT"]/1000
# realW.loc[:, "P"] = realW["mu_PT"]*np.cosh(realW["mu_ETA"])
# realW["Pweight"] = 1/realW["P"]
# realWT = realW.loc[realW["mu_ISMUON"] == True]
