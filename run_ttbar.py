# ----------------------------------------------------------------------
# Import Statements
# ----------------------------------------------------------------------
from numpythia import Pythia, hepmc_write, hepmc_read
from numpythia import STATUS, HAS_END_VERTEX, ABS_PDG_ID
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyjet import cluster
# ----------------------------------------------------------------------
# Function Definitions and Global Variables
# ----------------------------------------------------------------------
# njets, the number of jets (associated with each cluster?)
njets = 0

def calc_ET(E, C_1, C_2):
   # Eis the energy, C_1 and C_2 are component values
   # of a triange to determine the angle theta
   # E_T = Esin(\theta)
   return E*np.sin(np.arctan(C_2/C_1))

def calc_jetmass(jet, ep=True):
   if ep == True:
      # Jet constituent arrays are in the form of E, px, py, pz
      # sum the four vectors
      a = 0
      b = 0
      c = 0
      d = 0
      if (len(jet)) == 1:
         return np.sqrt(jet[0][0]**2
                      - jet[0][1]**2
                      - jet[0][2]**2 
                      - jet[0][3]**2)
      for tup in jet:
         a += tup[0]
         b += tup[1]
         c += tup[2]
         d += tup[3]
      summed = np.array([a, b, c, d])
      # calc the invariant mass
      return np.sqrt(summed[0]**2
                     - summed[1]**2
                     - summed[2]**2
                     - summed[3]**2
                    )
   else:
      # arrays are in the form of pT, eta, phi, mass
      print("ep=False is unimplemented\n")

def calc_effr(jet_mass, px, py):
   return 2*jet_mass / np.sqrt(px**2 + py**2)

def calc_eta(C_1, C_2):
   # Again, C stands for components
   theta = np.arctan(C_2/C_1)
   return -1*np.log(C_2/(2*C_1))

def calc_phi(C_1, C_2):
   # C stands for componenets
   return np.arctan(C_2/C_1)

def delR(true_w, w_jet):
   print("Function delR unimplemneted")
   pass
   
def nConsts(jet):
   #return the number of constituents in a jet
   return len(jet)
# ---------------------------------------------------------------------
# Generate Events
# ---------------------------------------------------------------------
pythia = Pythia('ttbar.cmnd', random_state=1)
selection = ((STATUS == 1) & ~HAS_END_VERTEX)
unclustered_particles = list()
# ---------------------------------------------------------------------
# Generate Jets
# ---------------------------------------------------------------------
for event in pythia(events=100):
   vectors = event.all(selection)
   sequence = cluster(np.array(vectors), R=0.4, p=-1, ep=True)
   jets = sequence.inclusive_jets()
   unclustered_particles.append(sequence.unclustered_particles())
   for jet in jets:
      const_arry = jet.constituents_array(ep=True)
      print(const_arry)
      print("Jet Mass = ", calc_jetmass(const_arry))
      print("nConsts = ", nConsts(const_arry))


