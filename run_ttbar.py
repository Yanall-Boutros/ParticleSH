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
   event_data.append(np.array(jets_data).T)
# event_data is an n x m x o matrix, where n is the event, m is an
# observable of the jet (See line 53), and o is the data entry.
# len(event_data) = number of numpythia events
# len(event_data[0]) = number of observables
# len(event_data[0][0]) = number of jets
# -----------------------------------------------------------------------
# Plot Data
# -----------------------------------------------------------------------
# Create a histogram of counts per event with respect to the event number
# for each physical property

# Plot of number of counts of mass of jet in the jets of that event
njets_in_event = list()
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
# Function to determine the range of an array for histogram
def get_min_max(array):
   return [int(min(array)), int(max(array))+6]

for n in range(len(event_data)):
   njets_in_event.append(len(event_data[n][0]))
   for m in range(len(event_data[0])):
      r = get_min_max(event_data[n][m])
      plt.hist(event_data[n][m],
                        c_v_index_to_step[m],
                        range=r)
      plt.title(c_v_index_to_name[m])
      plt.xlabel(c_v_index_to_name[m] + " in " +
                              c_v_index_to_units[m])
      plt.ylabel("Counts per event")
      plt.savefig(c_v_index_to_name[m]+".png")
      plt.figure()
njets_in_event = np.array(njets_in_event)
plt.hist(np.array(range(10)), njets_in_event)
plt.title("Number of jets in each event")
plt.xlabel("Event Number")
plt.ylabel("Counts per Event")
plt.savefig("Number of Jets by Event.png")
