import ROOT
import matplotlib.pyplot as plt
import numpy as np
tree = ROOT.TChain("DecayTree")

tree.Add("/tmp/13TeV__2018__magnet_down_data__Z_candidates.root")

restMass = 0.105658

fig = plt.figure()

invMassArray = np.array([])

for entry in tree:
  mup_PT = entry.mup_PT/1000
  mup_ETA = entry.mup_ETA
  mup_PHI = entry.mup_PHI

  mum_PT = entry.mum_PT/1000
  mum_ETA = entry.mum_ETA
  mum_PHI = entry.mum_PHI

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

  invMassArray = np.append(invMassArray, np.sqrt((mup_E + mum_E)**2 - (mup_x + mum_x)**2 - (mup_y + mum_y)**2 - (mup_z + mum_z)**2))


plt.hist(invMassArray, bins = 200, range = [0, 200], histtype='step')
plt.title("Invariant mass distribution from $\mu^+ \mu^-$")
plt.ylabel("Counts")
plt.xlabel("Invariant mass (GeV)")
plt.savefig("img/zpeak_s.png")
