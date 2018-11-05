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
event_data = []
jets_data = []
# Change lists to using 
for event in pythia(events=100):
   vectors = event.all(selection)
   sequence = cluster(vectors, R=0.4, p=-1, ep=True)
   jets = sequence.inclusive_jets()
   unclustered_particles.append(sequence.unclustered_particles())
   for i, jet in enumerate(jets):
      data = (
              jet.mass, jet.eta, jet.phi, jet.pt,
              len(jet.constituents_array()), 2*jet.mass/jet.pt
             )
      # Append Leading Jet information
      if i == 0:
         event_data.append(data)
      jets_data.append(data)
columns = [
           ("Mass", "GeV"), ("Eta", ""), ("Phi", ""),
           ("Transverse Momentum", "GeV"),
           ("Number of constituents", ""), ("Eff Radius", "")
          ]
nbins = [100]*len(columns)
event_data = np.array(event_data)
jets_data = np.array(jets_data)

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
# Aggregate data list
agg_data = [[], [], [], [], [], [], [], []]
# Function to determine the range of an array for histogram
def get_min_max(array):
   return [int(min(array)), int(max(array))+6]
for i,(name,units) in enumerate(columns):
   fig, ax = plt.subplots()
   ax.hist(event_data[:,i])
   ax.set_xlabel('{0:s} [{1:s}]'.format(name, units))
   ax.set_ylabel('Events')
   fig.savefig('{0:s}_event.png'.format(name))

   fig, ax = plt.subplots()
   ax.hist(jets_data[:,i])
   ax.set_xlabel('{0:s} [{1:s}]'.format(name, units))
   ax.set_ylabel('Jets')
   fig.savefig('{0:s}_jets.png'.format(name))
