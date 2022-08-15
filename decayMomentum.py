import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

c = 3*10**8
pionMass = 0.13957
pionTau = 2.6033*10**-8
pionbfrac = 0.999877
kaonMass = 0.493677
kaonTau = 1.238*10**-8
kaonbfrac = 0.6356
muonMass = 0.105658

with uproot.open("/storage/epp2/phshgg/Public/DVTuples__v24g/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root:WpNoMuID/DecayTree") as file:
  simQCD = file.arrays(["mu_PT", "mu_ETA", "mu_TRUEID"], library = "pd")

simQCD = simQCD.loc[(simQCD["mu_TRUEID"] == 321) | (simQCD["mu_TRUEID"] == 211)]
simQCD["mu_PT"] = simQCD["mu_PT"]/1000
simQCD = simQCD.loc[simQCD["mu_PT"] > 20]
simQCD["mu_P"] = simQCD["mu_PT"]*np.cosh(simQCD["mu_ETA"])

kaons = simQCD["mu_TRUEID"] == 321
pions = simQCD["mu_TRUEID"] == 211


kaontest = 13*kaonMass/(20*c) # maybe change the minimum p here from 20 to a maximum p?? (1000 ish)
piontest = 13*pionMass/(20*c)
kaonBound = np.exp(-(13*kaonMass)/(20*c*kaonTau))
pionBound = np.exp(-(13*pionMass)/(20*c*pionTau))

print(kaonBound, pionBound)
print(kaontest, piontest)

#rng generator for exponential
gen = np.random.default_rng()
#random lifetime for kaons and pions
simQCD.loc[:, "lifetime"] = np.zeros(len(simQCD))
# simQCD.loc[kaons, "lifetime"] = gen.exponential(kaonTau, len(simQCD.loc[kaons, "lifetime"]))
# simQCD.loc[pions, "lifetime"] = gen.exponential(pionTau, len(simQCD.loc[pions, "lifetime"]))
simQCD.loc[kaons, "lifetime"] = -kaonTau*np.log(np.random.uniform(low = 0, high = kaonBound, size = len(simQCD.loc[kaons, "lifetime"])))
simQCD.loc[pions, "lifetime"] = -pionTau*np.log(np.random.uniform(low = 0, high = pionBound, size = len(simQCD.loc[pions, "lifetime"]))) # or low = kaon/pionBound and high = 1?

# #decay length from random lifetimes
simQCD.loc[:, "decayLength"] = np.zeros(len(simQCD))
kMom = simQCD.loc[kaons, "mu_P"]
pMom = simQCD.loc[pions, "mu_P"]
simQCD.loc[kaons, "decayLength"] = c*simQCD.loc[kaons, "lifetime"]*kMom/kaonMass #np.sqrt(kaonMass**2 + kMom**2/c**2)
simQCD.loc[pions, "decayLength"] = c*simQCD.loc[pions, "lifetime"]*pMom/pionMass #np.sqrt(pionMass**2 + pMom**2/c**2)
# print(simQCD.loc[kaons, "decayLength"])

print(-c*kaonTau*np.log(kaonBound)*200/kaonMass)

#Only include hadrons which decay before the calorimeter (13m)
simQCD = simQCD.loc[simQCD["decayLength"] < 13]

print(len(simQCD.loc[kaons, "decayLength"]))

#change only P values which have decay length < 5m (i.e. before magnet)
simQCD.loc[:, "adjusted_P"] = simQCD["mu_P"]
simQCD.loc[kaons & (simQCD["decayLength"] < 5), "adjusted_P"] = simQCD.loc[kaons &(simQCD["decayLength"] < 5), "mu_P"]*kaonMass/muonMass
simQCD.loc[pions & (simQCD["decayLength"] < 5), "adjusted_P"] = simQCD.loc[pions &(simQCD["decayLength"] < 5), "mu_P"]*pionMass/muonMass

#################
#Plotting p graph
plt.hist(simQCD["mu_P"], bins = 50, label = "Pions and Kaons", histtype = "step")
plt.hist(simQCD["adjusted_P"], bins = 50, label = "Pions and Kaons with Decay Model", histtype = "step")
plt.legend()
plt.title("Comparing the p difference from the addition of a decay model")
plt.xlabel("p (GeV)")
plt.ylabel("Counts")
plt.savefig("img/decayMomentum/decayModelP.png")
plt.close()

#adjust PT values, likewise to P earlier
simQCD.loc[:, "adjusted_PT"] = simQCD["adjusted_P"]/np.cosh(simQCD["mu_ETA"])

#Plotting PT graph
plt.hist(simQCD["mu_PT"], bins = 50, label = "Pions and Kaons", histtype = "step")
plt.hist(simQCD["adjusted_PT"], bins = 50, label = "Pions and Kaons with Decay Model", histtype = "step")
plt.legend()
plt.title("Comparing the pT difference from the addition of a decay model")
plt.xlabel("pT (GeV)")
plt.ylabel("Counts")
plt.savefig("img/decayMomentum/decayModelPT.png")
plt.close()
