import uproot
import numpy as np
import matplotlib.pyplot as plt

with uproot.open("/tmp/13TeV_2018_34_Up_EW.root") as file:
  real_PT = file["WpNoMuID/DecayTree/mu_PT"].array(library = "np")/1000
  real_PTSUMCONE040 = file["WpNoMuID/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000
  isRealMuon = file["WpNoMuID/DecayTree/mu_ISMUON"].array(library = "np")
  realIPCHI2 = file["WpNoMuID/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  realTRCHI2 = file["WpNoMuID/DecayTree/mu_TRCHI2DOF"].array(library = "np")
  realETA = file["WpNoMuID/DecayTree/mu_ETA"].array(library = "np")
with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root") as file:
  QCD_PT = file["WpNoMuID/DecayTree/mu_PT"].array(library = "np")/1000
  QCD_PTSUMCONE040 = file["WpNoMuID/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000
  isSimMuon = file["WpNoMuID/DecayTree/mu_ISMUON"].array(library = "np")
  simIPCHI2 = file["WpNoMuID/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  simTRCHI2 = file["WpNoMuID/DecayTree/mu_TRCHI2DOF"].array(library = "np")
  QCDETA = file["WpNoMuID/DecayTree/mu_ETA"].array(library = "np")
  QCDtrueID = file["WpNoMuID/DecayTree/mu_TRUEID"].array(library = "np")


#subplots for log10(isolation) each at different small ranges of muon pT (code from compare.py)
start = float(input("Enter the start value for small muon pT range: "))

ipLim = 4
trLim = 2

logQCD_PTSUMCONE040 = np.log10(np.clip(QCD_PTSUMCONE040, 0.1, None))
logreal_PTSUMCONE040 = np.log10(np.clip(real_PTSUMCONE040, 0.1, None))

fig, ax = plt.subplots(2, 2)
fig.suptitle("Log10(Isolation) comparison at for IPCHI2 < 4 and TRCHI2 < 2")

for i in range(0,2):
  for j in range(0,2):
    cropped_real_PTSUMCONE040 = logreal_PTSUMCONE040[(real_PT > start) & (real_PT < start + 5) & (realIPCHI2 < ipLim) & (realTRCHI2 < trLim) & (realETA < 4.5) & (realETA > 2)]
    cropped_QCD_PTSUMCONE040 = logQCD_PTSUMCONE040[(QCD_PT > start) & (QCD_PT < start + 5) & (simIPCHI2 < ipLim) & (simTRCHI2 < trLim) & (QCDETA < 4.5) & (QCDETA > 2) & (QCDtrueID != 0 )]


    ax[i, j].hist(cropped_real_PTSUMCONE040, bins = 50, histtype = 'step', label = "Real", density = True)
    ax[i, j].hist(cropped_QCD_PTSUMCONE040, bins = 50, histtype = 'step', label = "Simulated QCD", density = True)
    title = str(start) + " < pT < " + str(start + 5)
    ax[i, j].set_title(title, fontsize = "small")
    plt.setp(ax[-1, :], xlabel='log10(isolation (GeV))')
    plt.setp(ax[:, 0], ylabel='Normalised Counts')

    start += 5

plt.subplots_adjust(hspace = 0.3, wspace = 0.3)
plt.legend(fontsize = "x-small", loc = "upper left")
plt.savefig("img/muoninvestigationLog.png")
