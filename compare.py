#compares the enriched background real data from enrichedW.py to the simulated QCD background data from simulatedQCDbackground.py

import uproot
import numpy as np
import matplotlib.pyplot as plt

with uproot.open("/tmp/13TeV_2018_34_Up_EW.root") as file:
  real_PT = file["WpNoMuID/DecayTree/mu_PT"].array(library = "np")/1000
  real_PTSUMCONE040 = file["WpNoMuID/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000
  isMuon = file["WpNoMuID/DecayTree/mu_ISMUON"].array(library = "np")
with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root") as file:
  QCD_PT = file["WpNoMuID/DecayTree/mu_PT"].array(library = "np")/1000
  QCD_PTSUMCONE040 = file["WpNoMuID/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000

print(sum(isMuon)*100/len(isMuon), "% of all interactions are muons")

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
plt.xlabel("isolation (GeV)")
plt.savefig("img/compare_PTSUMCONE040")

plt.figure()

#subplots for the isolation each at different small ranges of muon pT
start = float(input("Enter the start value for small pT range: "))

fig, ax = plt.subplots(2, 2)
fig.suptitle("Isolation comparison at different values of muon pT")

for i in range(0,2):
  for j in range(0,2):
    cropped_real_PTSUMCONE040 = real_PTSUMCONE040[(real_PT > start) & (real_PT < start + 5)]
    cropped_QCD_PTSUMCONE040 = QCD_PTSUMCONE040[(QCD_PT > start) & (QCD_PT < start + 5)]

    ax[i, j].hist(cropped_real_PTSUMCONE040, bins = 100, range = [-2, 140], histtype = 'step', density = True, label = "Real pT")
    ax[i, j].hist(cropped_QCD_PTSUMCONE040, bins = 100, range = [-2, 140], histtype = 'step', density = True, label = "Simulated QCD pT")
    title = str(start) + " < pT < " + str(start + 5)
    ax[i, j].set_title(title, fontsize = "small")
    plt.setp(ax[-1, :], xlabel='isolation')
    plt.setp(ax[:, 0], ylabel='Normalised Counts')

    start += 5

plt.subplots_adjust(hspace = 0.3, wspace = 0.3)
plt.legend(fontsize = "x-small")
plt.savefig("img/compare")

start -= 20

#find % of interactions which are muons (between 0-40 isolation)
cropped_isMuon = isMuon[start < real_PT]
adaptedreal_PT = real_PT[start < real_PT]
cropped_isMuon = cropped_isMuon[adaptedreal_PT < (start + 5)]

realRange = real_PTSUMCONE040[start < real_PT]
adaptedreal_PT = real_PT[start < real_PT]
realRange = realRange[adaptedreal_PT < (start + 5)]

cropped_isMuon = cropped_isMuon[realRange < 40]

print(sum(cropped_isMuon)*100/len(cropped_isMuon), "% of interactions", start, "-", start+5,"pT are muons")
