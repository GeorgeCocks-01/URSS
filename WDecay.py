import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

with uproot.open("/storage/epp2/phshgg/Public/DVTuples__v24g/13TeV_2018_34_Up_EW.root") as file:
  mu_PT = file["WpIso/DecayTree/mu_PT"].array(library = "np")/1000
  mu_PTSUMCONE040 = file["WpIso/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000

mu_PT = mu_PT[mu_PT >= 20]

#plots the transverse momentum of the muons from possible W decays
plt.hist(mu_PT, bins = 200, range = [20, 60], histtype = 'step')
plt.title("pT distribution of $\mu$")
plt.ylabel("Counts")
plt.xlabel("pT (GeV)")
plt.savefig("img/WDecay/Wmu_PT")

plt.close()

plt.hist(1/mu_PT, bins = 200, histtype = 'step')
plt.title("1/pT distribution of $\mu$")
plt.ylabel("Counts")
plt.xlabel("1/pT (1/GeV)")
plt.savefig("img/WDecay/test.png")

plt.close()


#plots the transverse momentum of all particles in a cone around muons
plt.hist(mu_PTSUMCONE040, bins = 200, range = [0, 125], histtype = 'step')
plt.title("pT distribution of particles in a cone around $\mu$")
plt.ylabel("Counts")
plt.xlabel("isolation (GeV)")
plt.savefig("img/WDecay/Wmu_PTSUMCONE040")

plt.close()

#plots the momentum of mu vs the log10 of the cone momentum (heat plot)
log_cone = np.log10(np.clip(mu_PTSUMCONE040, 0.1, None))

plt.hist2d(mu_PT, log_cone, bins=200, norm=LogNorm(), range = [[20,120], [-1.3, 4]])
plt.colorbar()
plt.ylabel("Log10(Isolation (GeV))")
plt.xlabel("muon pT (GeV)")
plt.savefig("img/WDecay/WmuonVScone")

plt.close()

#plots the transverse momentum of just mu, but only for values of cone momentum < 5GeV
sub5mu_PT = mu_PT[mu_PTSUMCONE040 < 5.0]

plt.hist(sub5mu_PT, bins = 200, range = [20, 80], histtype = 'step')
plt.title("pT distribution of $\mu$ for isolation < 5GeV")
plt.ylabel("Counts")
plt.xlabel("pT (GeV)")
plt.savefig("img/WDecay/W1D")
