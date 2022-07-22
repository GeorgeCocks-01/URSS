import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

c = 3*10**8
pionMass = 0.13957
pionTau = 2.6033*10**-8
pionbfrac = 0.999877
kaonMass = 0.493677
kaonTau = 1.238*10**-8
kaonbfrac = 0.6356

def f(L, p, mass, tau):
  return np.exp(-(L*mass)/(p*c*tau))

def fintergral(p, mass, tau, bfrac):
  return bfrac*integrate.quad(f, 0, 5, args = (p, mass, tau))[0]/integrate.quad(f, 0, np.inf, args = (p, mass, tau))[0]

def makeArray(p, mass, tau, bfrac):
  arr = np.array([])
  for i in p:
    arr = np.append(arr, fintergral(i, mass, tau, bfrac))
  return arr

#Plotting graph for Pion
p = np.linspace(1, 100, 1000)
plt.plot(p, makeArray(p, pionMass, pionTau, pionbfrac), label = "Pion")
plt.xlabel("p (GeV)")
plt.ylabel("f(p)")
plt.title("Theoretical misID for pions and kaons")

#plotting for Kaon
plt.plot(p, makeArray(p, kaonMass, kaonTau, kaonbfrac), label = "Kaon")
plt.legend()

plt.savefig("img/misIDcheck.png")
