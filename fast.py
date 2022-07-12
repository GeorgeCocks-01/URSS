import ROOT
import uproot
import numpy as np
import matplotlib.pyplot as plt

restMass = 0.105658

with uproot.open("/tmp/13TeV__2018__magnet_down_data__Z_candidates.root") as file:
  mup_PT = file["DecayTree"]["mup_PT"].array(library="np")/1000
  mup_ETA = file["DecayTree"]["mup_ETA"].array(library="np")
  mup_PHI = file["DecayTree"]["mup_PHI"].array(library="np")
  mum_PT = file["DecayTree"]["mum_PT"].array(library="np")/1000
  mum_ETA = file["DecayTree"]["mum_ETA"].array(library="np")
  mum_PHI = file["DecayTree"]["mum_PHI"].array(library="np")

mup_theta = 2*np.arctan(np.exp(-mup_ETA))
mum_theta = 2*np.arctan(np.exp(-mum_ETA))

mup_r = mup_PT/np.sin(mup_theta)
mum_r = mum_PT/np.sin(mum_theta)

mup_x = mup_r*np.sin(mup_theta)*np.cos(mup_PHI)
mup_y = mup_r*np.sin(mup_theta)*np.sin(mup_PHI)
mup_z = mup_r*np.cos(mup_theta)
mum_x = mum_r*np.sin(mum_theta)*np.cos(mum_PHI)
mum_y = mum_r*np.sin(mum_theta)*np.sin(mum_PHI)
mum_z = mum_r*np.cos(mum_theta)

mup_E = np.sqrt(restMass**2 + mup_x**2 + mup_y**2 + mup_z**2)
mum_E = np.sqrt(restMass**2 + mum_x**2 + mum_y**2 + mum_z**2)

invMass = np.sqrt((mup_E + mum_E)**2 - (mup_x + mum_x)**2 - (mup_y + mum_y)**2 - (mup_z + mum_z)**2)

plt.hist(invMass, bins = 200, range = [0, 200], histtype='step')
plt.title("Invariant mass distribution from $\mu^+ \mu^-$")
plt.ylabel("Counts")
plt.xlabel("Invariant mass (GeV)")
plt.savefig("zpeak_f.png")
