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

def getFracErr(data, data_misID, rangeMax):
  hist, bins = np.histogram(data, bins = 50, range = [55, rangeMax])
  frac = np.histogram(data_misID, bins = 50, range = [55, rangeMax])[0]/hist
  err = np.sqrt((frac*(1 - frac))/hist)
  return frac, err, bins

def integral(p, mass, tau, bfrac, length):
  return 100*bfrac*(1 - np.exp(-(length*mass)/(c*p*tau)))

def compositeIntegral(p, prefac1, prefac2, length = None, mass = None, tau = None, bfrac = None):
  if mass == None or tau == None or bfrac == None or length == None:
    return prefac1*p + prefac2
  else:
    return 100*bfrac*(1 - np.exp(-(length*mass)/(c*p*tau))) + prefac1*p + prefac2

# #Plotting graph for Pion, using length = 5
# pRange = np.linspace(55, 200, 1000)
# plt.plot(pRange, integral(pRange, pionMass, pionTau, pionbfrac, 5), label = "Pion")
# plt.xlabel("p (GeV)")
# plt.ylabel("f(p)")
# plt.title("Theoretical misID for pions and kaons")

# #Plotting graph for Kaon, using length = 5
# plt.plot(pRange, integral(pRange, kaonMass, kaonTau, kaonbfrac, 5), label = "Kaon")
# plt.legend()

# plt.savefig("img/misID/misIDcheck.png")
# plt.close()

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

print(np.shape(kaons_P))
print(np.shape(pions_P))
kaons_frac, kaons_err, kaons_bins = getFracErr(kaons_P, kaons_misID, 200)
pions_frac, pions_err, pions_bins = getFracErr(pions_P, pions_misID, 200)

plt.errorbar(x = kaons_bins[:-1], y = kaons_frac*100, label = "Kaon Data", yerr = kaons_err*100)
plt.errorbar(x = pions_bins[:-1], y = pions_frac*100, label = "Pion Data", yerr = pions_err*100)
plt.xlabel("P (GeV)")
plt.ylabel("Percent Mislabeled")
plt.title("Fitting QCD simulation misID data to get length")
##########################

#Finding the fit parameters from the data
kaonParams, kaonCov = curve_fit(lambda p, length: integral(p, kaonMass, kaonTau, kaonbfrac, length), kaons_bins[:-1], kaons_frac*100)
pionParams, pionCov = curve_fit(lambda p, length: integral(p, pionMass, pionTau, pionbfrac, length), pions_bins[:-1], pions_frac*100)

kaonLength = kaonParams[0]
pionLength = pionParams[0]

print("Kaon fit length:", round(kaonLength, 2), "and standard deviation:", round(np.sqrt(np.diag(kaonCov))[0], 2))
print("Pion fit length:", round(pionLength, 2), "and standard deviation:", round(np.sqrt(np.diag(pionCov))[0], 2))

pRange = np.linspace(55, 200, 1000)
plt.plot(pRange, integral(pRange, kaonMass, kaonTau, kaonbfrac, kaonLength), label = "Kaon Fit")
plt.plot(pRange, integral(pRange, pionMass, pionTau, pionbfrac, pionLength), label = "Pion Fit")
plt.legend()
plt.savefig("img/misID/misIDfit.png")
plt.close()




# #Now using a composite function of decay term and log growth term
# #Plotting graph for Pion, using length = 5
# pRange = np.linspace(55, 1000, 10000)
# plt.xscale("log")
# plt.plot(pRange, compositeIntegral(pRange, pionMass, pionTau, pionbfrac, pionLength, 0.7, 0.3), label = "Pion")
# plt.xlabel("p (GeV)")
# plt.ylabel("f(p)")
# plt.title("Theoretical misID for pions and kaons")

# #Plotting graph for Kaon, using length = 5
# plt.plot(pRange, compositeIntegral(pRange, kaonMass, kaonTau, kaonbfrac, kaonLength, 0.85, 0.25), label = "Kaon")
# plt.legend()

# plt.savefig("img/misID/misIDcheckcomposite.png")
# plt.close()




#Fitting the data on the log scale using the composite function of decay + log terms
kaons_frac, kaons_err, kaons_bins = getFracErr(kaons_P, kaons_misID, 1000)
pions_frac, pions_err, pions_bins = getFracErr(pions_P, pions_misID, 1000)


#Getting data for protons
protons_P = Ptot[(trueID == 2212)]
protons_misID = Ptot[(isMuon == True) & (trueID == 2212)]
protons_frac, protons_err, protons_bins = getFracErr(protons_P, protons_misID, 1000)

kaonParams, kaonCov = curve_fit(lambda p, kaonPrefac1, kaonPrefac2, kaonLength: compositeIntegral(p, kaonPrefac1, kaonPrefac2, kaonLength, kaonMass, kaonTau, kaonbfrac), kaons_bins[:-1], kaons_frac*100, bounds = ((-np.inf, -np.inf, -np.inf), (np.inf, np.inf, 100)))
pionParams, pionCov = curve_fit(lambda p, pionPrefac1, pionPrefac2, pionLength: compositeIntegral(p, pionPrefac1, pionPrefac2, pionLength, pionMass, pionTau, pionbfrac), pions_bins[:-1], pions_frac*100, bounds = ((-np.inf, -np.inf, -np.inf), (np.inf, np.inf, 100)))
protonParams, protonCov = curve_fit(lambda p, protonPrefac1, protonPrefac2: compositeIntegral(p, protonPrefac1, protonPrefac2), protons_bins[:-1], protons_frac*100)

print("kaon parameters:", kaonParams)
print("pion parameters:", pionParams)
print("proton parameters:", protonParams)

pRange = np.linspace(55, 1000, 10000)
plt.errorbar(x = kaons_bins[:-1], y = kaons_frac*100, label = "Kaon data", yerr = kaons_err*100)
plt.errorbar(x = pions_bins[:-1], y = pions_frac*100, label = "Pion data", yerr = pions_err*100)
plt.errorbar(x = protons_bins[:-1], y = protons_frac*100, label = "Proton data", yerr = protons_err*100)
plt.xscale("log")
plt.plot(pRange, compositeIntegral(pRange, *kaonParams, kaonMass, kaonTau, kaonbfrac), label = "Kaon Fit")
plt.plot(pRange, compositeIntegral(pRange, *pionParams, pionMass, pionTau, pionbfrac), label = "Pion Fit")
plt.plot(pRange, compositeIntegral(pRange, *protonParams), label = "Proton Fit", color = "k")

plt.legend(loc = "upper center")
plt.ylabel("Percent mislabeled")
plt.xlabel("P (GeV)")
plt.title("Fitting QCD simulation misID data with punch through to get length")
plt.savefig("img/misID/misIDcompositefit.png")
plt.close()



################################################
#Comparing 1/P weights to weights with integrals
PT = PT[PT > 20]

kaonIntegralWeight = compositeIntegral(PT[trueID == 321], *kaonParams, kaonMass, kaonTau, kaonbfrac)
pionIntegralWeight = compositeIntegral(PT[(trueID == 211)], *pionParams, pionMass, pionTau, pionbfrac)
protonIntegralWeight = compositeIntegral(PT[(trueID == 2212)], *protonParams)

plt.hist(PT[(trueID == 321)], bins = 50, label = "Kaons with 1/P weights", weights = 1/(Ptot[(trueID == 321)]), range = [20, 70], histtype = "step", density = True, color = "b")
plt.hist(PT[(trueID == 211)], bins = 50, label = "Pions with 1/P weights", weights = 1/(Ptot[(trueID == 211)]), range = [20, 70], histtype = "step", density = True, color = "r")
plt.hist(PT[(trueID == 2212)], bins = 50, label = "Protons with no weight", range = [20, 70], histtype = "step", density = True, color = "g")
plt.hist(PT[(trueID == 321)], bins = 50, label = "Kaons with Composite weights", weights = kaonIntegralWeight, range = [20, 70], histtype = "step", density = True, linestyle = "--", color = "b")
plt.hist(PT[(trueID == 211)], bins = 50, label = "Pions with Composite weights", weights = pionIntegralWeight, range = [20, 70], histtype = "step", density = True, linestyle = "--", color = "r")
plt.hist(PT[(trueID == 2212)], bins = 50, label = "Protons with Composite weights", weights = protonIntegralWeight, range = [20, 70], histtype = "step", density = True, linestyle = "--", color = "g")

plt.legend()
plt.ylabel("Normalised Counts")
plt.xlabel("pT (GeV)")
plt.title("Comparing 1/P weights to Composite function weights")
plt.savefig("img/misID/misIDWeightsComparison")
plt.close()
