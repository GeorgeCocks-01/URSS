import uproot
import numpy as np
import matplotlib.pyplot as plt

with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root") as file:
  PTSUMCONE040 = file["WpNoMuID/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000
  IPCHI2 = file["WpNoMuID/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  TRCHI2 = file["WpNoMuID/DecayTree/mu_TRCHI2DOF"].array(library = "np")
  trueID = file["WpNoMuID/DecayTree/mu_TRUEID"].array(library = "np")
  isMuon = file["WpNoMuID/DecayTree/mu_ISMUON"].array(library = "np")

# particles = np.unique(trueID)

PTSUMCONE040 = PTSUMCONE040[isMuon == False]
IPCHI2 = IPCHI2[isMuon == False]
TRCHI2 = TRCHI2[isMuon == False]
trueID = trueID[isMuon == False]

particleDict = {
  2212: "p",
  321: "K+",
  211: "pi+",
  0: "Unknown",
  -11: "e+",
  -13: "mu+",
  # -15: "tau+",
  # -211: "pi-",
  # -321: "K-",
  # -2212: "anti-p",
}

# particles = np.intersect1d(particles, np.array([*particleDict.keys()]))


#Plotting the isolation graph
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
fig.suptitle("Isolation for different particles")

for i, j in particleDict.items():
  particle_PTSUMCONE040 = PTSUMCONE040[trueID == i]
  ax1.hist(particle_PTSUMCONE040, bins = 100, histtype = "step", label = j, range = [0, 150])
  ax2.hist(particle_PTSUMCONE040, bins = 100, histtype = "step", label = j, range = [0, 150], density = True)

ax1.set(xlabel = "Isolation (GeV)", ylabel = "Counts")
ax2.set(xlabel = "Isolation (GeV)", ylabel = "Normalised Counts")
ax1.legend(loc = "upper right", fontsize = "x-small")
plt.savefig("img/particlesISO")

plt.close()


#Plotting the ipchi2 graph
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
fig.suptitle("IPCHI2 for different particles")

for i, j in particleDict.items():
  particle_IPCHI2 = IPCHI2[trueID == i]
  ax1.hist(particle_IPCHI2, bins = 100, histtype = "step", label = j, range = [0, 12])
  ax2.hist(particle_IPCHI2, bins = 100, histtype = "step", label = j, range = [0, 20], density = True)

ax1.set(xlabel = "IPCHI2", ylabel = "Counts")
ax2.set(xlabel = "IPCHI2", ylabel = "Normalised Counts")
ax1.legend(loc = "upper right", fontsize = "x-small")
plt.savefig("img/particlesIPCHI2")

plt.close()


#Plotting the trchi2 graph
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
fig.suptitle("TRCHI2 for different particles")

for i, j in particleDict.items():
  particle_TRCHI2 = TRCHI2[trueID == i]
  ax1.hist(particle_TRCHI2, bins = 100, histtype = "step", label = j)
  ax2.hist(particle_TRCHI2, bins = 100, histtype = "step", label = j, density = True)

ax1.set(xlabel = "TRCHI2", ylabel = "Counts")
ax2.set(xlabel = "TRCHI2", ylabel = "Normalised Counts")
ax1.legend(loc = "upper right", fontsize = "x-small")
plt.savefig("img/particlesTRCHI2")
