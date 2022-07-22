import uproot
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit

c = 3*10**8
pionMass = 0.13957
pionTau = 2.6033*10**-8
pionbfrac = 0.999877
kaonMass = 0.493677
kaonTau = 1.238*10**-8
kaonbfrac = 0.6356

def f(L, p, mass, tau):
  return np.exp(-(L*mass)/(p*c*tau))

def fintergral(p, mass, tau, bfrac, length):
  return 100*bfrac*integrate.quad(f, 0, length, args = (p, mass, tau))[0]/integrate.quad(f, 0, np.inf, args = (p, mass, tau))[0]

def makeArray(p, mass, tau, bfrac, length):
  arr = np.array([])
  for i in p:
    arr = np.append(arr, fintergral(i, mass, tau, bfrac, length))
  return arr

#Plotting graph for Pion, using length = 5
pRange = np.linspace(20, 100, 1000)
plt.plot(pRange, makeArray(pRange, pionMass, pionTau, pionbfrac, 5), label = "Pion")
plt.xlabel("p (GeV)")
plt.ylabel("f(p)")
plt.title("Theoretical misID for pions and kaons")

#Plotting graph for Kaon, using length = 5
plt.plot(pRange, makeArray(pRange, kaonMass, kaonTau, kaonbfrac, 5), label = "Kaon")
plt.legend()

plt.savefig("img/misIDcheck.png")
plt.close()

#Fit the QCD simulation data to the theoretical pion and kaon plots (above) to give the length
with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root") as file:
  PT = file["WpNoMuID/DecayTree/mu_PT"].array(library = "np")/1000
  isMuon = file["WpNoMuID/DecayTree/mu_ISMUON"].array(library = "np")
  trueID = file["WpNoMuID/DecayTree/mu_TRUEID"].array(library = "np")
  ETA = file["WpNoMuID/DecayTree/mu_ETA"].array(library = "np")

#code taken from misID.py
Ptot = PT[PT > 20]*np.cosh(ETA[PT > 20])
trueID = trueID[PT > 20]
isMuon = isMuon[PT > 20]
kaons_P = Ptot[(trueID == 321)]
pions_P = Ptot[(trueID == 211)]
kaons_misID = Ptot[(isMuon == True) & (trueID == 321)]
pions_misID = Ptot[(isMuon == True) & (trueID == 211)]
kaons_hist, kaons_bins = np.histogram(kaons_P, bins = 50, range = [55, 200])
pions_hist, pions_bins = np.histogram(pions_P, bins = 50, range = [55, 200])
kaons_frac = np.histogram(kaons_misID, bins = 50, range = [55, 200])[0]/kaons_hist
pions_frac = np.histogram(pions_misID, bins = 50, range = [55, 200])[0]/pions_hist
kaons_err = np.sqrt((kaons_frac*(1 - kaons_frac))/kaons_hist)
pions_err = np.sqrt((pions_frac*(1 - pions_frac))/pions_hist)

plt.errorbar(x = kaons_bins[:-1], y = kaons_frac*100, label = "Kaon Data", yerr = kaons_err*100)
plt.errorbar(x = pions_bins[:-1], y = pions_frac*100, label = "Pion Data", yerr = pions_err*100)
plt.xlabel("P (GeV)")
plt.ylabel("Percent Mislabeled")
plt.title("Fitting QCD simulation misID data to get length")
##########################

kaonParams, kaonCov = curve_fit(lambda p, length: makeArray(p, kaonMass, kaonTau, kaonbfrac, length), kaons_bins[:-1], kaons_frac*100)
pionParams, pionCov = curve_fit(lambda p, length: makeArray(p, pionMass, pionTau, pionbfrac, length), pions_bins[:-1], pions_frac*100)

print(kaonParams[0], np.sqrt(np.diag(kaonCov)))
print(pionParams[0], np.sqrt(np.diag(pionCov)))

pRange = np.linspace(55, 200, 1000)
plt.plot(pRange, makeArray(pRange, kaonMass, kaonTau, kaonbfrac, kaonParams[0]), label = "Kaon Fit")
plt.plot(pRange, makeArray(pRange, pionMass, pionTau, pionbfrac, pionParams[0]), label = "Pion Fit")
plt.legend()
plt.savefig("img/misIDfit.png")
