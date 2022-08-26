import numpy as np
import matplotlib.pyplot as plt

c = 3*10**8
pionMass = 0.13957
pionTau = 2.6033*10**-8
kaonMass = 0.493677
kaonTau = 1.238*10**-8


p = np.linspace(55, 1000, 10000)

pionRatio = 0.5*(15*kaonMass/(c*kaonTau*p))
kaonRatio = 0.5*(15*pionMass/(c*pionTau*p))

plt.plot(p, pionRatio, label = "pions")
plt.plot(p, kaonRatio, label = "kaons")
plt.legend()
plt.savefig("img/taylorTest.png")
plt.close()
