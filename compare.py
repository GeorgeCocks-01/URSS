#compares the enriched background real data from enrichedW.py to the simulated QCD background data from simulatedQCDbackground.py

import uproot
import numpy as np
import matplotlib.pyplot as plt

with uproot.open("/tmp/13TeV_2018_34_Up_EW.root") as file:
  real_PT = file["WpNoMuID/DecayTree/mu_PT"].array(library = "np")/1000
  real_PTSUMCONE040 = file["WpNoMuID/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000
with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root") as file:
  QCD_PT = file["WpNoMuID/DecayTree/mu_PT"].array(library = "np")/1000
  QCD_PTSUMCONE040 = file["WpNoMuID/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000

#plots the transverse momentum of the muons. Range starts at 20GeV to mitigate the artefact in the simulated QCD data
plt.hist(real_PT, bins = 200, range = [20, 60], histtype = 'step', density = True, label = "Real pT")
plt.hist(QCD_PT, bins = 200, range = [20, 60], histtype = 'step', density = True, label = "Simulated QCD pT")
plt.legend()
plt.title("muon pT of background rich data vs. simulated QCD data")
plt.ylabel("Normalised Counts")
plt.xlabel("pT (GeV)")
plt.savefig("img/compare_PT")

plt.figure()


#plots the transverse momentum of all particles in a cone around muons
plt.hist(real_PTSUMCONE040, bins = 200, range = [-2, 140], histtype = 'step', density = True, label = "Real pT")
plt.hist(QCD_PTSUMCONE040, bins = 200, range = [-2, 140], histtype = 'step', density = True, label = "Simulated QCD pT")
plt.legend()
plt.title("isolation of background rich data vs. simulated QCD data")
plt.ylabel("Normalised Counts")
plt.xlabel("pT (GeV)")
plt.savefig("img/compare_PTSUMCONE040")
