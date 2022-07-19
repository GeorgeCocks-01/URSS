import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


with uproot.open("/tmp/13TeV_2018_34_Up_EW.root") as file:
  realIPCHI2 = file["WpNoMuID/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  realTRCHI2 = file["WpNoMuID/DecayTree/mu_TRCHI2DOF"].array(library = "np")
  isRealMuon = file["WpNoMuID/DecayTree/mu_ISMUON"].array(library = "np")
with uproot.open("/tmp/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root") as file:
  QCDsimIPCHI2 = file["WpNoMuID/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  QCDsimTRCHI2 = file["WpNoMuID/DecayTree/mu_TRCHI2DOF"].array(library = "np")
  isQCDMuon = file["WpNoMuID/DecayTree/mu_ISMUON"].array(library = "np")
with uproot.open("/tmp/13TeV_2018_34_Up_W_Sim09k.root") as file:
  WsimIPCHI2 = file["WpIso/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  WsimTRCHI2 = file["WpIso/DecayTree/mu_TRCHI2DOF"].array(library = "np")
  isWMuon = file["WpIso/DecayTree/mu_ISMUON"].array(library = "np")

realIPCHI2 = realIPCHI2[isRealMuon == True]
realTRCHI2 = realTRCHI2[isRealMuon == True]
QCDsimIPCHI2 = QCDsimIPCHI2[isQCDMuon == True]
QCDsimTRCHI2 = QCDsimTRCHI2[isQCDMuon == True]
WsimIPCHI2 = WsimIPCHI2[isWMuon == True]
WsimTRCHI2 = WsimTRCHI2[isWMuon == True]

#2D plot the real IPCHI2 and TRCHI2 data
plt.hist2d(realIPCHI2, realTRCHI2, bins=200, norm=LogNorm(), range = [[0, 17], [0, 2.5]])
plt.colorbar()
plt.ylabel("TRCHI2")
plt.xlabel("IPCHI2")
plt.title("Real W data")
plt.savefig("img/real_trVSipCHI")

plt.figure()

#2D plot the IPCHI2 and TRCHI2 data from the QCD background simulations
plt.hist2d(QCDsimIPCHI2, QCDsimTRCHI2, bins=200, norm=LogNorm(), range = [[0, 17], [0, 2.5]])
plt.colorbar()
plt.ylabel("TRCHI2")
plt.xlabel("IPCHI2")
plt.title("Simulation QCD data")
plt.savefig("img/simQCD_trVSipCHI")

plt.figure()

#2D plot the IPCHI2 and TRCHI2 data from the W -> mu nu simulation
plt.hist2d(WsimIPCHI2, WsimTRCHI2, bins=200, norm=LogNorm(), range = [[0, 14], [0, 2]])
plt.colorbar()
plt.ylabel("TRCHI2")
plt.xlabel("IPCHI2")
plt.title("Simulation W data")
plt.savefig("img/simW_trVSipCHI")

plt.figure()




#Plot of real, QCD sim and W sim IPCHI2
plt.hist(realIPCHI2, bins = 200, range = [0, 15], histtype = "step", label = "Real", density = True)
plt.hist(QCDsimIPCHI2, bins = 200, range = [0, 15], histtype = "step", label = "QCD sim", density = True)
plt.hist(WsimIPCHI2, bins = 200, range = [0, 15], histtype = "step", label = "W sim", density = True)
plt.legend()
plt.xlabel("IPCHI2")
plt.ylabel("Normalised Counts")
plt.title("IPCHI2 comparison")
plt.savefig("img/ipchi2")

plt.figure()

#Plot of real, QCD sim and W sim TRCHI2
plt.hist(realTRCHI2, bins = 200, range = [0, 3.5], histtype = "step", label = "Real", density = True)
plt.hist(QCDsimTRCHI2, bins = 200, range = [0, 3.5], histtype = "step", label = "QCD sim", density = True)
plt.hist(WsimTRCHI2, bins = 200, range = [0, 3.5], histtype = "step", label = "W sim", density = True)
plt.legend()
plt.xlabel("TRCHI2")
plt.ylabel("Normalised Counts")
plt.title("TRCHI2 comparison")
plt.savefig("img/trchi2")

plt.figure()

#Log plots of the real, QCD sim and W sim IPCHI2 data
plt.hist(realIPCHI2, bins = 200, range = [0, 40000], histtype = "step", label = "Real", log = True, density = True)
plt.hist(QCDsimIPCHI2, bins = 200, range = [0, 40000], histtype = "step", label = "QCD sim", log = True, density = True)
plt.hist(WsimIPCHI2, bins = 200, range = [0, 40000], histtype = "step", label = "W sim", log = True, density = True)
plt.legend()
plt.xlabel("IPCHI2")
plt.ylabel("Log(Normalised Counts)")
plt.title("IPCHI2 comparison")
plt.savefig("img/logipchi2")
