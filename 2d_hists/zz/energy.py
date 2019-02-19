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
from matplotlib.colors import ListedColormap
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
# Load Colormap info
# -----------------------------------------------------------------------
chans = open("colormap.txt", 'r') # Credit to Ian Heywood for the 
for chan in chans: exec(chan) # Python implementation of CubeHelix
                              # Green, D. A., 2011, `A colour scheme for
                              # the display of astronomical intensity
                              # images', Bulletin of the Astronomical
                              #  Society of India, 39, 289.
                              # (2011BASI...39..289G at ADS.) 
ctab = []
for i in range(len(r_chan)):
    ctab.append([r_chan[i], g_chan[i], b_chan[i]])
cubehelixcm = ListedColormap(ctab,name='cubehelix', N=None)
# -----------------------------------------------------------------------
# Generate Events
# -----------------------------------------------------------------------
pythia = Pythia('zz.cmnd', random_state=1)
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
for event in pythia(events=10000):
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
plt.figure()
plt.hist2d(leading_xdata, leading_ydata,
    range=[(-5,5),(-1*np.pi, np.pi)],
    bins=(50, 50), cmap='cubehelix')
#cax = ax.imshow((leading_xdata, leading_ydata), cmap=
plt.xlabel("$\eta$")
plt.ylabel("$\phi$")
plt.title("$ZZ$")
plt.colorbar()
plt.savefig("LeadingZZ.pdf")
plt.figure()
plt.hist2d(jet_xdata, jet_ydata,
    range=[(-5,5),(-1*np.pi, np.pi)],
    bins=(50, 50), cmap='cubehelix')
plt.xlabel("$\eta$")
plt.ylabel("$\phi$")
plt.title("$ZZ$")
plt.colorbar()
plt.savefig("AggregateZZ.pdf")
