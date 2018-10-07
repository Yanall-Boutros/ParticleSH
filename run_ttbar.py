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
# -----------------------------------------------------------------------
# Generate Events
# -----------------------------------------------------------------------
pythia = Pythia('ttbar.cmnd', random_state=1)
selection = ((STATUS == 1) & ~HAS_END_VERTEX)
unclustered_particles = list()
# -----------------------------------------------------------------------
# Generate Jets and Histogram Data
# -----------------------------------------------------------------------
event_data = list()
for event in pythia(events=1):
   vectors = event.all(selection)
   sequence = cluster(vectors, R=0.4, p=-1, ep=True)
   jets = sequence.inclusive_jets()
   unclustered_particles.append(sequence.unclustered_particles())
   # position i of jets_data links to the list containing various bins
   # data for a jet
   jets_data = list()
   for jet in jets:
      jet_data = list()
      jet_data.append(jet.mass)
      jet_data.append(jet.eta)
      jet_data.append(jet.phi)
      jet_data.append(jet.e)
      jet_data.append(jet.et)
      jet_data.append(jet.pt)
      jet_data.append(int(len(jet.constituents_array())))
      jet_data.append(jet.mass*2/(np.sqrt(jet.px**2 + jet.py**2)))
      jets_data.append(np.array(jet_data))
   jets_data.append(len(jets))
   event_data.append(np.array(jets_data))
# Data Structure Artchitecture Logic: event_data is the set of all jets
# data. The cardinality of event_data = the number of events. Each member
# in event_data is also a set. Let each member {jets} in the set of event
# data also be a set, with the cardinality equal to the number of jets.
# Each jet is a member in the set of jets and is an ordered set which
# contains various physical properties of the jet 
# -----------------------------------------------------------------------
# Plot Data
# -----------------------------------------------------------------------
# Create a histogram of counts per event with respect to the event number
# for each physical property

# Plot of number of counts of mass of jet in the jets of that event
jets_njets = list()
c_v_index_to_name = {
   0 : "Mass"
   1 : "Eta"
   2 : "Phi"
   3 : "Energy"
   4 : "Transverse Energy"
   5 : "Transverse Momentum"
   6 : "Number of constituents"
   7 : "Eff r"
}
c_v_index_to_step = {
   0 :  
   1 e
   2 : 
   3 : 
   4 : 
   5 : 
   6 : 
   7 : 
}
c_v_index_to_range = dict()
for jets_elem in event_data:
   jets_properties = list()
   i = 0
   for jet_elem in jets_elem[:-1]:
      for physical_property in jet_elem:
         if (len(jets_properties) == 0):
            jets_properties.append([physical_property])
            i += 1
            continue
         if i < 8:
            jets_properties.append([physical_property])
         else: jets_properties[i % 8].append(physical_property)
         i += 1
   for column_variable_index in range(len(jets_properties)):
      
      plt.hist(jets_properties[column_variable_index],
                              100,
                              range=[0, 768])
      plt.title(str(column_variable_index))
      plt.xlabel("Jet Variable Property")
      plt.ylabel("Coutnts per event")
      plt.savefig(str(column_variable_index)+".pdf")
      plt.figure()
