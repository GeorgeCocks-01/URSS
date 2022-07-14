#Creates plots of pT, eta and phi for the mu+ and mu- separately
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


fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(10,4))
fig.suptitle('$\mu^+$')

ax1.hist(mup_PT, bins = 200, range = [0, 200], histtype='step')
ax1.set(xlabel='pT (GeV)', ylabel='Counts')
ax2.hist(mup_ETA, bins = 200, range = [0, 6], histtype='step')
ax2.set(xlabel='$\eta$')
ax3.hist(mup_PHI, bins = 200, range = [0, 4], histtype='step')
ax3.set(xlabel='$\phi$')

plt.subplots_adjust(wspace = 0.4)
plt.savefig("mup.png")



fig, (ax4, ax5, ax6) = plt.subplots(1,3, figsize=(10,4))
fig.suptitle('$\mu^-$')

ax4.hist(mum_PT, bins = 200, range = [0, 200], histtype='step')
ax4.set(xlabel='pT (GeV)', ylabel='Counts')
ax5.hist(mum_ETA, bins = 200, range = [0, 6], histtype='step')
ax5.set(xlabel='$\eta$')
ax6.hist(mum_PHI, bins = 200, range = [0, 4], histtype='step')
ax6.set(xlabel='$\phi$')

plt.subplots_adjust(wspace = 0.4)
plt.savefig("mum.png")

# Plots for the mass, pT and rapidity for the {mu+ mu-} system

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

dimuon_E = np.sqrt(restMass**2 + mup_x**2 + mup_y**2 + mup_z**2) + np.sqrt(restMass**2 + mum_x**2 + mum_y**2 + mum_z**2)
dimuon_x = mup_x + mum_x
dimuon_y = mup_y + mum_y
dimuon_z = mup_z + mum_z

invMass = np.sqrt((dimuon_E)**2 - (dimuon_x)**2 - (dimuon_y)**2 - (dimuon_z)**2)

total_PT = np.sqrt((dimuon_x)**2 + (dimuon_y)**2)

y = 0.5*np.log((dimuon_E + dimuon_z)/(dimuon_E - dimuon_z))

fig, (ax7, ax8, ax9) = plt.subplots(1,3, figsize=(10,4))
fig.suptitle('$Di-\mu$ system')

ax7.hist(invMass, bins = 200, range = [30, 155], histtype='step')
ax7.set(xlabel="Invariant Mass (GeV)", ylabel="Counts")
ax8.hist(total_PT, bins = 200, range = [0, 200], histtype='step')
ax8.set(xlabel="Total pT (GeV)")
ax9.hist(y, bins = 200, range = [1, 5], histtype='step')
ax9.set(xlabel="Rapidity")

plt.subplots_adjust(wspace = 0.4)
plt.savefig("dimuon.png")
