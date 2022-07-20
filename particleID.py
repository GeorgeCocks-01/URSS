import uproot
import numpy as np
import matplotlib.pyplot as plt

with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root") as file:
  PTSUMCONE040 = file["WpNoMuID/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000
  IPCHI2 = file["WpNoMuID/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  TRCHI2 = file["WpNoMuID/DecayTree/mu_TRCHI2DOF"].array(library = "np")
  trueID = file["WpNoMuID/DecayTree/mu_TRUEID"].array(library = "np")

particles = np.unique(trueID)

particleDict = {
  4232: "Xi_c+",
  4122: "Lambda_c+",
  3222: "Sigma+",
  2212: "p",
  521: "B+",
  431: "D_s+",
  411: "D+",
  321: "K+",
  211: "Pi+",
  0: "0",
  -11: "e+",
  -13: "mu+",
  -15: "tau+",
  -211: "Pi-",
  -321: "K-",
  -411: "D-",
  -431: "D_s-",
  -521: "B-",
  -2212: "p-",
  -3112: "Anti-Sigma-",
  -3312: "Xi+",
  -3334: "Omega+",
  -4122: "Lambda_c-",
  -4232: "Xi_c-"
}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
fig.suptitle("Isolation for different particles")

#Plotting the isolation graph
for i in range(0, len(particles)):
  print(particles[i])
  particle_PTSUMCONE040 = PTSUMCONE040[trueID == particles[i]]
  ax1.hist(particle_PTSUMCONE040, bins = 100, histtype = "step", label = particleDict[particles[i]], range = [0, 150])
  ax2.hist(particle_PTSUMCONE040, bins = 100, histtype = "step", label = particleDict[particles[i]], range = [0, 150], density = True)

ax1.set(xlabel = "Isolation (GeV)", ylabel = "Counts")
ax2.set(xlabel = "Isolation (GeV)", ylabel = "Normalised Counts")
ax1.legend(loc = "right", fontsize = "x-small")
# plt.title("Isolation for different particles")
# plt.ylabel("Counts")
# plt.xlabel("Isolation (GeV)")
plt.savefig("img/particlesISO")

plt.figure()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
fig.suptitle("IPCHI2 for different particles")

#Plotting the ipchi2 graph
for i in range(0, len(particles)):
  particle_IPCHI2 = IPCHI2[trueID == particles[i]]
  ax1.hist(particle_IPCHI2, bins = 100, histtype = "step", label = particleDict[particles[i]], range = [0, 12])
  ax2.hist(particle_IPCHI2, bins = 100, histtype = "step", label = particleDict[particles[i]], range = [0, 4000], density = True)

ax1.set(xlabel = "IPCHI2", ylabel = "Counts")
ax2.set(xlabel = "IPCHI2", ylabel = "Normalised Counts")
ax1.legend(loc = "right", fontsize = "x-small")
# plt.title("IPCHI2 for different particles")
# plt.ylabel("Counts")
# plt.xlabel("IPCHI2")
plt.savefig("img/particlesIPCHI2")

plt.figure()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
fig.suptitle("TRCHI2 for different particles")

#Plotting the trchi2 graph
for i in range(0, len(particles)):
  particle_TRCHI2 = TRCHI2[trueID == particles[i]]
  ax1.hist(particle_TRCHI2, bins = 100, histtype = "step", label = particleDict[particles[i]])
  ax2.hist(particle_TRCHI2, bins = 100, histtype = "step", label = particleDict[particles[i]], density = True)

ax1.set(xlabel = "TRCHI2", ylabel = "Counts")
ax2.set(xlabel = "TRCHI2", ylabel = "Normalised Counts")
ax1.legend(loc = "right", fontsize = "x-small")
# plt.title("TRCHI2 for different particles")
# plt.ylabel("Counts")
# plt.xlabel("TRCHI2")
plt.savefig("img/particlesTRCHI2")
