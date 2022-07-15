import uproot
import numpy as np
import matplotlib.pyplot as plt

with uproot.open("/tmp/13TeV_2018_34_Up_EW.root") as file:
  mu_PT = file["WpNoMuID/DecayTree/mu_PT"].array(library = "np")/1000
  mu_PTSUMCONE040 = file["WpNoMuID/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000

#plots the transverse momentum of the muons from real data with enriched background
plt.hist(mu_PT, bins = 200, range = [10, 60], histtype = 'step')
plt.title("pT distribution of $\mu$ with enriched background")
plt.ylabel("Counts")
plt.xlabel("pT (GeV)")
plt.savefig("img/enrichedmu_PT")

plt.figure()


#plots the transverse momentum of all particles in a cone around muons
plt.hist(mu_PTSUMCONE040, bins = 200, range = [-2, 140], histtype = 'step')
plt.title("pT distribution of particles in a cone around $\mu$ with enriched background")
plt.ylabel("Counts")
plt.xlabel("pT (GeV)")
plt.savefig("img/enrichedmu_PTSUMCONE040")

