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
  return (mass*100*bfrac*integrate.quad(f, 0, length, args = (p, mass, tau))[0])/(p*c*tau)

def makeArray(p, mass, tau, bfrac, length):
  arr = np.array([])
  for i in p:
    arr = np.append(arr, fintergral(i, mass, tau, bfrac, length))
  return arr


def fgrowintergral(p, mass, tau, bfrac, length, prefac1, prefac2):
  return 100*(prefac1*mass*bfrac*integrate.quad(f, 0, length, args = (p, mass, tau))[0])/(p*c*tau) + prefac2*np.log10(p)

def makegrowArray(p, mass, tau, bfrac, length, prefac1, prefac2):
  arr = np.array([])
  for i in p:
    arr = np.append(arr, fgrowintergral(i, mass, tau, bfrac, length, prefac1, prefac2))
  return arr

#Plotting graph for Pion, using length = 5
pRange = np.linspace(55, 200, 1000)
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

#Finding the fit parameters from the data
kaonParams, kaonCov = curve_fit(lambda p, length: makeArray(p, kaonMass, kaonTau, kaonbfrac, length), kaons_bins[:-1], kaons_frac*100)
pionParams, pionCov = curve_fit(lambda p, length: makeArray(p, pionMass, pionTau, pionbfrac, length), pions_bins[:-1], pions_frac*100)

kaonLength = kaonParams[0]
pionLength = pionParams[0]

print("Kaon fit length:", round(kaonLength, 2), "and standard deviation:", round(np.sqrt(np.diag(kaonCov))[0], 2))
print("Pion fit length:", round(pionLength, 2), "and standard deviation:", round(np.sqrt(np.diag(pionCov))[0], 2))

pRange = np.linspace(55, 200, 1000)
plt.plot(pRange, makeArray(pRange, kaonMass, kaonTau, kaonbfrac, kaonLength), label = "Kaon Fit")
plt.plot(pRange, makeArray(pRange, pionMass, pionTau, pionbfrac, pionLength), label = "Pion Fit")
plt.legend()
plt.savefig("img/misIDfit.png")
plt.close()




#Now using a composite function of decay term and log growth term
#Plotting graph for Pion, using length = 5
pRange = np.linspace(55, 1000, 10000)
plt.xscale("log")
plt.plot(pRange, makegrowArray(pRange, pionMass, pionTau, pionbfrac, 5, 0.7, 0.3), label = "Pion")
plt.xlabel("p (GeV)")
plt.ylabel("f(p)")
plt.title("Theoretical misID for pions and kaons")

#Plotting graph for Kaon, using length = 5
plt.plot(pRange, makegrowArray(pRange, kaonMass, kaonTau, kaonbfrac, 5, 0.85, 0.25), label = "Kaon")
plt.legend()

plt.savefig("img/misIDcheckcomposite.png")
plt.close()




#Fitting the data on the log scale using the composite function of decay + log terms
kaons_hist, kaons_bins = np.histogram(kaons_P, bins = 50, range = [55, 1000])
pions_hist, pions_bins = np.histogram(pions_P, bins = 50, range = [55, 1000])
kaons_frac = np.histogram(kaons_misID, bins = 50, range = [55, 1000])[0]/kaons_hist
pions_frac = np.histogram(pions_misID, bins = 50, range = [55, 1000])[0]/pions_hist
kaons_err = np.sqrt((kaons_frac*(1 - kaons_frac))/kaons_hist)
pions_err = np.sqrt((pions_frac*(1 - pions_frac))/pions_hist)


kaonParams, kaonCov = curve_fit(lambda p, prefac1, prefac2: makegrowArray(p, kaonMass, kaonTau, kaonbfrac, kaonLength, prefac1, prefac2), kaons_bins[:-1], kaons_frac*100)
pionParams, pionCov = curve_fit(lambda p, prefac1, prefac2: makegrowArray(p, pionMass, pionTau, pionbfrac, pionLength, prefac1, prefac2), pions_bins[:-1], pions_frac*100)

print(kaonParams)
print(pionParams)

pRange = np.linspace(55, 1000, 10000)
plt.errorbar(x = kaons_bins[:-1], y = kaons_frac*100, label = "Kaons", yerr = kaons_err*100)
plt.errorbar(x = pions_bins[:-1], y = pions_frac*100, label = "Pions", yerr = pions_err*100)
plt.xscale("log")
plt.plot(pRange, makegrowArray(pRange, kaonMass, kaonTau, kaonbfrac, kaonLength, kaonParams[0], kaonParams[1]), label = "Kaon Fit")
plt.plot(pRange, makegrowArray(pRange, pionMass, pionTau, pionbfrac, pionLength, pionParams[0], pionParams[1]), label = "Pion Fit")

plt.legend()
plt.xlabel("Percent mislabeled")
plt.ylabel("P (GeV)")
plt.title("Fitting QCD simulation misID data with punch through to get length")
plt.savefig("img/misIDcompositefit.png")
plt.close()
