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
def get_min_max(array):
   return (int(min(array)), int(max(array))+1)
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
for event in pythia(events=10):
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
   0 : "Mass",
   1 : "Eta",
   2 : "Phi",
   3 : "Energy",
   4 : "Transverse Energy",
   5 : "Transverse Momentum",
   6 : "Number of constituents",
   7 : "Eff r",
}
c_v_index_to_units = {
   0 : "Mass units",
   1 : "Eta units",
   2 : "Phi units",
   3 : "Energy units",
   4 : "Transverse Energy units",
   5 : "Transverse Momentum units",
   6 : "Number of constituents units",
   7 : "Eff r units",
}
c_v_index_to_step = {
   0 : 100,
   1 : 100,
   2 : 100,
   3 : 100,
   4 : 100,
   5 : 100,
   6 : 100,
   7 : 100,
}
# To Do, remove index_to_range dict and replace with a dynamic method
# of determining the range
c_v_index_to_range = {
   0 : [0, 15],
   1 : [-10, 10],
   2 : [-5, 5],
   3 : [0, 1500],
   4 : [0, 85],
   5 : [0, 85],
   6 : [0, 20],
   7 : [0, 20],
}
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
                        c_v_index_to_step[column_variable_index],
                        range=c_v_index_to_range[column_variable_index])
      plt.title(c_v_index_to_name[column_variable_index])
      plt.xlabel(c_v_index_to_name[column_variable_index] + " in " +
                              c_v_index_to_units[column_variable_index])
      plt.ylabel("Counts per event")
      plt.savefig(c_v_index_to_name[column_variable_index]+".png")
      plt.figure()
