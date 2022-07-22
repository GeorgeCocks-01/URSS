#Calculates fractions of Kaons, Pions, Protons and type 0 that have isMuon == True
import uproot
import numpy as np
import matplotlib.pyplot as plt

with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root") as file:
  PT = file["WpNoMuID/DecayTree/mu_PT"].array(library = "np")/1000
  isMuon = file["WpNoMuID/DecayTree/mu_ISMUON"].array(library = "np")
  trueID = file["WpNoMuID/DecayTree/mu_TRUEID"].array(library = "np")
  ETA = file["WpNoMuID/DecayTree/mu_ETA"].array(library = "np")

P = PT[PT > 20]*np.cosh(ETA[PT > 20])
trueID = trueID[PT > 20]
isMuon = isMuon[PT > 20]

kaons_P = P[(trueID == 321)]
pions_P = P[(trueID == 211)]
protons_P = P[(trueID == 2212)]
type0_P = P[(trueID == 0)]

kaons_misID = P[(isMuon == True) & (trueID == 321)]
pions_misID = P[(isMuon == True) & (trueID == 211)]
protons_misID = P[(isMuon == True) & (trueID == 2212)]
type0_misID = P[(isMuon == True) & (trueID == 0)]

kaons_hist, kaons_bins = np.histogram(kaons_P, bins = 50, range = [55, 1000])
pions_hist, pions_bins = np.histogram(pions_P, bins = 50, range = [55, 1000])
protons_hist, protons_bins = np.histogram(protons_P, bins = 50, range = [55, 1000])
type0_hist, type0_bins = np.histogram(type0_P, bins = 50, range = [55, 1000])


kaons_frac = np.histogram(kaons_misID, bins = 50, range = [55, 1000])[0]/kaons_hist
pions_frac = np.histogram(pions_misID, bins = 50, range = [55, 1000])[0]/pions_hist
protons_frac = np.histogram(protons_misID, bins = 50, range = [55, 1000])[0]/protons_hist
type0_frac = np.histogram(type0_misID, bins = 50, range = [55, 1000])[0]/type0_hist


kaons_err = np.sqrt((kaons_frac*(1 - kaons_frac))/kaons_hist)
pions_err = np.sqrt((pions_frac*(1 - pions_frac))/pions_hist)
protons_err = np.sqrt((protons_frac*(1 - protons_frac))/protons_hist)
type0_err = np.sqrt((type0_frac*(1 - type0_frac))/type0_hist)

plt.xscale("log")
plt.errorbar(x = kaons_bins[:-1], y = kaons_frac*100, label = "Kaons", yerr = kaons_err*100)
plt.errorbar(x = pions_bins[:-1], y = pions_frac*100, label = "Pions", yerr = pions_err*100)
plt.errorbar(x = protons_bins[:-1], y = protons_frac*100, label = "Protons", yerr = protons_err*100)
plt.errorbar(x = type0_bins[:-1], y = type0_frac*100, label = "type 0", yerr = type0_err*100)



plt.legend(loc = "upper center")
plt.xlabel("P (GeV)")
plt.ylabel("Percent Mislabeled")
plt.title("Misidentified Particles in QCD Simulation data")
plt.savefig("img/misID.png")
