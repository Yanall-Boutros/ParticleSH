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
def is_massless_and_isolated(jet):
   # Returns true if a jet has nconsts = 1 and has a pdgid equal to that
   # of a photon or a gluon, (todo: When to also discard muons)
   if len(jet.constituents_array()) == 1: 
      if np.abs(jet.info['pdgid']) == 21 or np.abs(jet.info['pdgid']) == 22:
         return True
      return np.abs(jet.info['pdgid']) == 13 and 2*jet.mass/jet.pt > 0.4
event_data = []
jets_data = []
discarded_data = [] 
for event in pythia(events=10000):
   vectors = event.all(selection)
   sequence = cluster(vectors, R=0.4, p=-1, ep=True)
   jets = sequence.inclusive_jets()
   unclustered_particles.append(sequence.unclustered_particles())
   for i, jet in enumerate(jets):
      data = (
              jet.mass, jet.eta, jet.phi, jet.pt,
              len(jet.constituents_array()), 2*jet.mass/jet.pt
             )
      if is_massless_and_isolated(jet):
         discarded_data.append(data)
      else:
         # Append Leading Jet information
         if i == 0:
            event_data.append(data)
         jets_data.append(data)
columns = [
           ("Mass", "GeV"), ("Eta", ""), ("Phi", ""),
           ("Transverse Momentum", "GeV"),
           ("Number of constituents", ""), ("Eff Radius", "")
          ]
leading_ranges = [
                  (-1,100), (-5,5), (-4,4), (-5,800), (-1,80),
                  (-0.25, 0.75)
                 ]
agg_ranges = [
              (-1,100), (-15,15), (-4,4), (-5,100), (-1,40), (-0.25,1.25)
             ]
nbins = [100]*len(columns)
event_data = np.array(event_data)
jets_data = np.array(jets_data)
# -----------------------------------------------------------------------
# Plot Data
# -----------------------------------------------------------------------
# Create a histogram of counts per event with respect to the event number
# for each physical property
# Plot of number of counts of mass of jet in the jets of that event
for i,(name,units) in enumerate(columns):
   # Plot the same data but in the 0-10GeV range for Transverse Momentum
   # and Mass. As well as plot the entire data but with 5GeV increments
   if i == 0 or i == 3:
      # length of range divided by 5 = num of bins
      b_l = int((leading_ranges[i][1] - leading_ranges[i][0])/5)+1
      b_a = int((agg_ranges[i][1] - agg_ranges[i][0])/5)+1
      # plot the 10GeV range for for leading jet data
      fig, ax = plt.subplots()
      r = (-1, 10) # Set the range to 10 GeV
      ag_tit = "Aggregate Data 0-10Gev"
      lead_tit = "Leading Jet Data 0-10GeV"
      if i == 3: # Unless we're plotting the Momentum
         r=(-1, 100) # then set the range to -1 to 100
         ag_tit = "Aggregate Data 0-100Gev" # and adjust the titles
         lead_tit = "Leading Jet Data 0-100GeV" # accordingly
      ax.hist(event_data[:,i], bins=nbins[i], range=r)
      plt.title(lead_tit)
      ax.set_xlabel('{0:s} [{1:s}]'.format(name, units))
      ax.set_ylabel('Jets')
      plt.savefig('{0:s}_event10GeV.png'.format(name))

      # plot the 10GeV range for aggregate data
      fig, ax = plt.subplots()
      ax.hist(jets_data[:,i], bins=nbins[i], range=r)
      plt.title(ag_tit)
      ax.set_xlabel('{0:s} [{1:s}]'.format(name, units))
      ax.set_ylabel('Jets')
      plt.savefig('{0:s}_jets10GeV.png'.format(name))
      
      # plot the standard ranges.
      # Plot the data from each leading jet in each event
      fig, ax = plt.subplots()
      ax.hist(event_data[:,i], bins=nbins[i], range=leading_ranges[i])
      ax.set_xlabel('{0:s} [{1:s}]'.format(name, units))
      ax.set_ylabel('Events')
      plt.title("Leading Jet Data")
      plt.savefig('{0:s}_event.png'.format(name))
      
      # Plot the data from all the jets in all events
      fig, ax = plt.subplots()
      ax.hist(jets_data[:,i], bins=nbins[i], range=agg_ranges[i])
      plt.title("Aggregate Data")
      ax.set_xlabel('{0:s} [{1:s}]'.format(name, units))
      ax.set_ylabel('Jets')
      plt.savefig('{0:s}_jets.png'.format(name))
      
      # Plots with increments being 5GeV
      # Plot the data from each leading jet in each event
      fig, ax = plt.subplots()
      ax.hist(event_data[:,i], bins=b_l, range=leading_ranges[i])
      ax.set_xlabel('{0:s} [{1:s}]'.format(name, units))
      ax.set_ylabel('Events')
      plt.title("Leading Jet Data (bin width = 5 GeV)")
      plt.savefig('{0:s}_event_by_5.png'.format(name))
      
      # Plot the data from all the jets in all events
      fig, ax = plt.subplots()
      ax.hist(jets_data[:,i], bins=b_a, range=agg_ranges[i])
      plt.title("Aggregate Data (bin width = 5 GeV)")
      ax.set_xlabel('{0:s} [{1:s}]'.format(name, units))
      ax.set_ylabel('Jets')
      plt.savefig('{0:s}_jets_by_5.png'.format(name))
   else: 
       # Plot the data from each leading jet in each event
       fig, ax = plt.subplots()
       ax.hist(event_data[:,i], bins=nbins[i], range=leading_ranges[i])
       ax.set_xlabel('{0:s} [{1:s}]'.format(name, units))
       ax.set_ylabel('Events')
       plt.title("Leading Jet Data")
       plt.savefig('{0:s}_event.png'.format(name))
       # Plot the data from all the jets in all events
       fig, ax = plt.subplots()
       ax.hist(jets_data[:,i], bins=nbins[i], range=agg_ranges[i])
       plt.title("Aggregate Data")
       ax.set_xlabel('{0:s} [{1:s}]'.format(name, units))
       ax.set_ylabel('Jets')
       plt.savefig('{0:s}_jets.png'.format(name))
