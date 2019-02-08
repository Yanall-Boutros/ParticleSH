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
from mpl_toolkits.mplot3d import Axes3D
def print_debug_data(dd):
    for tup in dd:
        j = tup[1]
        print("Mass = ", j[0], "\tpT = ", j[3], "\tEffR = ", j[5])
def allequs(jet):
    mag = (jet.px**2 + jet.py**2 + jet.pz**2)**0.5
    eng = jet.constituents_array(ep=True)[0][0]
    if mag == eng and eng == jet.e:
        return True
    return False
# -----------------------------------------------------------------------
# Generate Events
# -----------------------------------------------------------------------
pythia = Pythia('ttbar.cmnd', random_state=1)
selection = ((STATUS == 1) & ~HAS_END_VERTEX)
unclustered_particles = list()
# -----------------------------------------------------------------------
# Generate Jets and Histogram Data
# -----------------------------------------------------------------------
def is_massless_or_isolated(jet):
   # Returns true if a jet has nconsts = 1 and has a pdgid equal to that
   # of a photon or a gluon
   if len(jet.constituents_array()) == 1: 
      if np.abs(jet.info['pdgid']) == 21 or np.abs(jet.info['pdgid']) == 22:
         return True
      # if a muon is outside of the radius of the jet, discard it
      if np.abs(jet.info['pdgid']) == 13:
         if 2*jet.mass/jet.pt > 0.4: return True
   # Remove Jets with too high an eta
   if np.abs(jet.eta) > 5.0:
      return True
   # Remove any jets less than an arbitrary near zero mass
   if jet.mass < 0.4:
      return True
   return False
debug_data = []
leading_data = []
jets_data = []
discarded_data = [] 
# -----------------------------------------------------------------------
# Main Loop for Storing Jet Data
# -----------------------------------------------------------------------
for event in pythia(events=1000):
   vectors = event.all(selection)
   sequence = cluster(vectors, R=0.4, p=-1, ep=True)
   jets = sequence.inclusive_jets()
   unclustered_particles.append(sequence.unclustered_particles())
   for i, jet in enumerate(jets):
      data = (
              jet.eta, jet.phi, jet.e,
             )
      if is_massless_or_isolated(jet):
         discarded_data.append(jet)
      else:
         # Append Leading Jet information
         if i == 0:
            leading_data.append(data)
         jets_data.append(data)
# -----------------------------------------------------------------------
# Plot 3D-histograms
# -----------------------------------------------------------------------
leading_xdata = np.array(leading_data[0])
leading_ydata = np.array(leading_data[1])
leading_zdata = np.array(leading_data[2])
jet_xdata = np.array(jets_data[0])
jet_ydata = np.array(jets_data[1])
jet_zdata = np.array(jets_data[2])

fig = plt.figure()
plt.hist2d(leading_xdata, leading_ydata, range=[(-5,5), (0, 2*np.pi)])
plt.xlabel("$\eta$")
plt.ylabel("$\phi$")
plt.title("$T\overline{T}$")
plt.savefig("Test.pdf")
