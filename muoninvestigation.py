import uproot
import numpy as np
import matplotlib.pyplot as plt


with uproot.open("/tmp/13TeV_2018_34_Up_EW.root") as file:
  real_PT = file["WpNoMuID/DecayTree/mu_PT"].array(library = "np")/1000
  real_PTSUMCONE040 = file["WpNoMuID/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000
  isRealMuon = file["WpNoMuID/DecayTree/mu_ISMUON"].array(library = "np")
  realIPCHI2 = file["WpNoMuID/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  realTRCHI2 = file["WpNoMuID/DecayTree/mu_TRCHI2DOF"].array(library = "np")
with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root") as file:
  QCD_PT = file["WpNoMuID/DecayTree/mu_PT"].array(library = "np")/1000
  QCD_PTSUMCONE040 = file["WpNoMuID/DecayTree/mu_PTSUMCONE040"].array(library = "np")/1000
  isSimMuon = file["WpNoMuID/DecayTree/mu_ISMUON"].array(library = "np")
  simIPCHI2 = file["WpNoMuID/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  simTRCHI2 = file["WpNoMuID/DecayTree/mu_TRCHI2DOF"].array(library = "np")


#subplots for the isolation each at different small ranges of muon pT (code from compare.py)
start = float(input("Enter the start value for small pT range: "))

ipLimit = 4
trLim = 2

fig, ax = plt.subplots(2, 2)
fig.suptitle("Isolation comparison at different values of muon pT for ")

for i in range(0,2):
  for j in range(0,2):
    realRange = real_PTSUMCONE040[start < real_PT]
    adaptedreal_PT = real_PT[start < real_PT]
    cropped_isRealMuon = isRealMuon[start < real_PT]
    cropped_realIPCHI2 = realIPCHI2[start < real_PT]
    cropped_realTRCHI2 = realTRCHI2[start < real_PT]
    realRange = realRange[adaptedreal_PT < (start + 5)]
    cropped_isRealMuon = cropped_isRealMuon[adaptedreal_PT < (start + 5)]
    cropped_realIPCHI2 = cropped_realIPCHI2[adaptedreal_PT < (start + 5)]
    cropped_realTRCHI2 = cropped_realTRCHI2[adaptedreal_PT < (start + 5)]

    QCDRange = QCD_PTSUMCONE040[start < QCD_PT]
    adaptedQCD_PT = QCD_PT[start < QCD_PT]
    cropped_isSimMuon = isSimMuon[start < QCD_PT]
    cropped_simIPCHI2 = simIPCHI2[start < QCD_PT]
    cropped_simTRCHI2 = simTRCHI2[start < QCD_PT]
    QCDRange = QCDRange[adaptedQCD_PT < (start + 5)]
    cropped_isSimMuon = cropped_isSimMuon[adaptedQCD_PT < (start + 5)]
    cropped_simIPCHI2 = cropped_simIPCHI2[adaptedQCD_PT < (start + 5)]
    cropped_simTRCHI2 = cropped_simTRCHI2[adaptedQCD_PT < (start + 5)]

    #require isRealMuon = True
    # realRange = realRange[cropped_isRealMuon == True]
    # QCDRange = QCDRange[cropped_isSimMuon == True]

    #require IPCHI2 of to be below ipLimit
    realRange = realRange[cropped_realIPCHI2 < ipLimit]
    cropped_realTRCHI2 = cropped_realTRCHI2[cropped_realIPCHI2 < ipLimit]
    QCDRange = QCDRange[cropped_simIPCHI2 < ipLimit]
    cropped_simTRCHI2 = cropped_simTRCHI2[cropped_simIPCHI2 < ipLimit]

    #require TRCHI2 to be below trLim
    realRange = realRange[cropped_realTRCHI2 < trLim]
    QCDRange = QCDRange[cropped_simTRCHI2 < trLim]

    ax[i, j].hist(realRange, bins = 200, range = [-2, 140], histtype = 'step', density = True, label = "Real pT")
    ax[i, j].hist(QCDRange, bins = 200, range = [-2, 140], histtype = 'step', density = True, label = "Simulated QCD pT")
    title = str(start) + " < pT < " + str(start + 5)
    ax[i, j].set_title(title, fontsize = "small")
    plt.setp(ax[-1, :], xlabel='isolation')
    plt.setp(ax[:, 0], ylabel='Normalised Counts')

    start += 5

plt.subplots_adjust(hspace = 0.3, wspace = 0.3)
plt.legend(fontsize = "x-small")
plt.savefig("img/isMuonTrue")
