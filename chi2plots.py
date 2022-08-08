import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

with uproot.open("/storage/epp2/phshgg/Public/DVTuples__v24g/13TeV_2018_34_Up_EW.root") as file:
  realQCDIPCHI2 = file["WpNoMuID/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  realQCDTRCHI2 = file["WpNoMuID/DecayTree/mu_TRCHI2DOF"].array(library = "np")
  isRealQCDMuon = file["WpNoMuID/DecayTree/mu_ISMUON"].array(library = "np")

  realWIPCHI2 = file["WpIso/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  realWTRCHI2 = file["WpIso/DecayTree/mu_TRCHI2DOF"].array(library = "np")
  isRealWMuon = file["WpIso/DecayTree/mu_ISMUON"].array(library = "np")
with uproot.open("/storage/epp2/phshgg/Public/DVTuples__v24g/13TeV_2017_29r2_Up_QcdBgdPt18GeV_Sim09k.root") as file:
  QCDsimIPCHI2 = file["WpNoMuID/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  QCDsimTRCHI2 = file["WpNoMuID/DecayTree/mu_TRCHI2DOF"].array(library = "np")
  isQCDMuon = file["WpNoMuID/DecayTree/mu_ISMUON"].array(library = "np")
  QCDtrueID = file["WpNoMuID/DecayTree/mu_TRUEID"].array(library = "np")
with uproot.open("/storage/epp2/phshgg/Public/DVTuples__v24g/13TeV_2018_34_Up_W_Sim09k.root") as file:
  WsimIPCHI2 = file["WpIso/DecayTree/mu_MIPCHI2DV"].array(library = "np")
  WsimTRCHI2 = file["WpIso/DecayTree/mu_TRCHI2DOF"].array(library = "np")
  isWMuon = file["WpIso/DecayTree/mu_ISMUON"].array(library = "np")

realQCDIPCHI2 = realQCDIPCHI2[isRealQCDMuon == False]
realQCDTRCHI2 = realQCDTRCHI2[isRealQCDMuon == False]
realWIPCHI2 = realWIPCHI2[isRealWMuon == True]
realWTRCHI2 = realWTRCHI2[isRealWMuon == True]
QCDsimIPCHI2 = QCDsimIPCHI2[isQCDMuon == False]
QCDsimTRCHI2 = QCDsimTRCHI2[isQCDMuon == False]
QCDtrueID = QCDtrueID[isQCDMuon == False]
WsimIPCHI2 = WsimIPCHI2[isWMuon == True]
WsimTRCHI2 = WsimTRCHI2[isWMuon == True]

# #2D plot the real IPCHI2 and TRCHI2 data
# plt.hist2d(realQCDIPCHI2, realQCDTRCHI2, bins=200, norm=LogNorm(), range = [[0, 17], [0, 2.5]])
# plt.colorbar()
# plt.ylabel("TRCHI2")
# plt.xlabel("IPCHI2")
# plt.title("Real W data")
# plt.savefig("img/real_trVSipCHI")

# plt.close()

# #2D plot the IPCHI2 and TRCHI2 data from the QCD background simulations
# plt.hist2d(QCDsimIPCHI2, QCDsimTRCHI2, bins=200, norm=LogNorm(), range = [[0, 17], [0, 2.5]])
# plt.colorbar()
# plt.ylabel("TRCHI2")
# plt.xlabel("IPCHI2")
# plt.title("Simulation QCD data")
# plt.savefig("img/simQCD_trVSipCHI")

# plt.close()

# #2D plot the IPCHI2 and TRCHI2 data from the W -> mu nu simulation
# plt.hist2d(WsimIPCHI2, WsimTRCHI2, bins=200, norm=LogNorm(), range = [[0, 14], [0, 2]])
# plt.colorbar()
# plt.ylabel("TRCHI2")
# plt.xlabel("IPCHI2")
# plt.title("Simulation W data")
# plt.savefig("img/simW_trVSipCHI")

# plt.close()




#Plot of real, QCD sim and W sim IPCHI2
plt.hist(realQCDIPCHI2, bins = 200, range = [0, 15], histtype = "step", label = "Real QCD", density = True)
plt.hist(realWIPCHI2, bins = 200, range = [0, 15], histtype = "step", label = "Real W", density = True)
plt.hist(QCDsimIPCHI2, bins = 200, range = [0, 15], histtype = "step", label = "QCD sim", density = True)
plt.hist(WsimIPCHI2, bins = 200, range = [0, 15], histtype = "step", label = "W sim", density = True)
plt.legend()
plt.xlabel("IPCHI2")
plt.ylabel("Normalised Counts")
plt.title("IPCHI2 comparison")
plt.savefig("img/ipchi2")

plt.close()

#Plot of real, QCD sim and W sim TRCHI2, first with and then without type 0 particles
plt.hist(realQCDTRCHI2, bins = 200, range = [0, 3.5], histtype = "step", label = "Real QCD", density = True)
plt.hist(realWTRCHI2, bins = 200, range = [0, 3.5], histtype = "step", label = "Real W", density = True)
plt.hist(QCDsimTRCHI2, bins = 200, range = [0, 3.5], histtype = "step", label = "QCD sim", density = True)
plt.hist(WsimTRCHI2, bins = 200, range = [0, 3.5], histtype = "step", label = "W sim", density = True)
plt.legend()
plt.vlines(1, 0, 1.75)
plt.xlabel("TRCHI2")
plt.ylabel("Normalised Counts")
plt.title("TRCHI2 comparison")
plt.savefig("img/trchi2")

plt.close()

QCDsimTRCHI2 = QCDsimTRCHI2[QCDtrueID != 0] #removes "unknown" (type = 0) particles from qcd sim TRCHI2 data
plt.hist(realQCDTRCHI2, bins = 200, range = [0, 3.5], histtype = "step", label = "Real QCD", density = True)
plt.hist(realWTRCHI2, bins = 200, range = [0, 3.5], histtype = "step", label = "Real W", density = True)
plt.hist(QCDsimTRCHI2, bins = 200, range = [0, 3.5], histtype = "step", label = "QCD sim", density = True)
plt.hist(WsimTRCHI2, bins = 200, range = [0, 3.5], histtype = "step", label = "W sim", density = True)
plt.legend()
plt.vlines(1, 0, 1.75)
plt.xlabel("TRCHI2")
plt.ylabel("Normalised Counts")
plt.title("TRCHI2 comparison without type 0 particles")
plt.savefig("img/trchi2_new")

plt.close()

#Log plots of the real, QCD sim and W sim IPCHI2 data
plt.hist(realQCDIPCHI2, bins = 200, range = [0, 40000], histtype = "step", label = "Real QCD", log = True, density = True)
plt.hist(realWIPCHI2, bins = 200, range = [0, 40000], histtype = "step", label = "Real W", log = True, density = True)
plt.hist(QCDsimIPCHI2, bins = 200, range = [0, 40000], histtype = "step", label = "QCD sim", log = True, density = True)
plt.hist(WsimIPCHI2, bins = 200, range = [0, 40000], histtype = "step", label = "W sim", log = True, density = True)
plt.legend()
plt.xlabel("IPCHI2")
plt.ylabel("Normalised Counts")
plt.title("IPCHI2 comparison")
plt.savefig("img/ipchi2logscale")
plt.close()

#Plot of the real, QCD sim and W sim log(IPCHI2) data
plt.hist(np.log10(realQCDIPCHI2), bins = 200, range = [-3.5, 4], histtype = "step", label = "Real QCD", density = True)
plt.hist(np.log10(realWIPCHI2), bins = 200, range = [-3.5, 4], histtype = "step", label = "Real W", density = True)
plt.hist(np.log10(QCDsimIPCHI2), bins = 200, range = [-3.5, 4], histtype = "step", label = "QCD sim", density = True)
plt.hist(np.log10(WsimIPCHI2), bins = 200, range = [-3.5, 4], histtype = "step", label = "W sim", density = True)
plt.legend()
plt.xlabel("log10(IPCHI2)")
plt.ylabel("Normalised Counts")
plt.title("IPCHI2 comparison")
plt.savefig("img/ipchi2log10")
plt.close()