import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

with uproot.open("/tmp/13TeV_2018_34_Up_EW.root") as file:
  mu_PT = file["WpIso/DecayTree/mu_PT"].array(library = "np")/1000
  mu_PTSUMCONE040 = file["WpIso/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000

#plots the transverse momentum of the muons from possible W decays
plt.hist(mu_PT, bins = 200, range = [0, 170], histtype = 'step')
plt.title("pT distribution of $\mu$")
plt.ylabel("Counts")
plt.xlabel("pT (GeV)")
plt.savefig("img/Wmu_PT")

plt.figure()


#plots the transverse momentum of all particles in a cone around muons
plt.hist(mu_PTSUMCONE040, bins = 200, range = [0, 170], histtype = 'step')
plt.title("pT distribution of particles in a cone around $\mu$")
plt.ylabel("Counts")
plt.xlabel("pT (GeV)")
plt.savefig("img/Wmu_PTSUMCONE040")

plt.figure()

#plots the momentum of mu vs the log10 of the cone momentum (heat plot)
log_cone = np.log10(np.clip(mu_PTSUMCONE040, 0.1, None))

plt.hist2d(mu_PT, log_cone, bins=200, norm=LogNorm(), range = [[20,120], [-1.3, 4]])
plt.colorbar()
plt.ylabel("Log10(pT particles in cone (GeV))")
plt.xlabel("muon pT (GeV)")
plt.savefig("img/muonVScone")

plt.figure()

#plots the transverse momentum of just mu, but only for values of cone momentum < 5GeV
sub5mu_PT = mu_PT[mu_PTSUMCONE040 < 5.0]

plt.hist(sub5mu_PT, bins = 200, range = [20, 140], histtype = 'step')
plt.title("pT distribution of $\mu$ for isolation < 5GeV")
plt.ylabel("Counts")
plt.xlabel("pT (GeV)")
plt.savefig("img/1D")
