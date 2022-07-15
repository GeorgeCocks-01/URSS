import uproot
import numpy as np
import matplotlib.pyplot as plt

with uproot.open("/tmp/13TeV_2018_34_Up_W_Sim09k.root") as file:
  mu_PT = file["WpIso/DecayTree/mu_PT"].array(library = "np")/1000
  mu_PTSUMCONE040 = file["WpIso/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000

#plots the transverse momentum of the muons from simulated W decays
plt.hist(mu_PT, bins = 200, range = [10, 125], histtype = 'step')
plt.title("pT distribution of $\mu$ for simulated event")
plt.ylabel("Counts")
plt.xlabel("pT (GeV)")
plt.savefig("img/simWmu_PT")

plt.figure()


#plots the transverse momentum of all particles in a cone around muons
plt.hist(mu_PTSUMCONE040, bins = 200, range = [0, 12], histtype = 'step')
# plt.ylim((None, 10000))
plt.title("pT distribution of particles in a cone around $\mu$ for simulated event")
plt.ylabel("Counts")
plt.xlabel("pT (GeV)")
plt.savefig("img/simWmu_PTSUMCONE040")

plt.figure()