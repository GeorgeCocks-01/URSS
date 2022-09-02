import uproot
import numpy as np

with uproot.open("../gcocks-background-simulation/mW/data/tuples/Wm_QcdBgdPt18GeV_13TeV.root") as file:
  Wm = file["DecayTree/mu_P"].array(library = "np")
with uproot.open("../gcocks-background-simulation/mW/data/tuples/Wp_QcdBgdPt18GeV_13TeV.root") as file:
  Wp = file["DecayTree/mu_P"].array(library = "np")

print("Wm: ", len(Wm))

print("Wp: ", len(Wp))

print(Wp)
print(Wm)
