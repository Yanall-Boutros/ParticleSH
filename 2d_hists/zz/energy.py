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
id = 'pdgid'
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
        if np.abs(jet.info[id]) == 21 or np.abs(jet.info[id]) == 22:
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
def return_particle_data(jet):
    # return the array containing all the eta, phi, and energies of the
    # particles in a jets constituent array
    vals = []
    has_eta_phi = jet.constituents_array()
    eta = []
    phi = []
    m = []
    pt = []
    for i in range(len(has_eta_phi)):
        pt.append(has_eta_phi[i][0])
        eta.append(has_eta_phi[i][1])
        phi.append(has_eta_phi[i][2])
        m.append(has_eta_phi[i][3])
    m = np.array(m)
    pt = np.array(pt)
    e = (pt**2 + m**2)**0.5
    return [eta, phi, e]
debug_data = []
leading_data = []
leading_particle_eta = []
leading_particle_phi = []
leading_particle_energy = []
jets_data = []
jets_particle_eta = []
jets_particle_phi = []
jets_particle_energy = []
discarded_data = [] 
# -----------------------------------------------------------------------
# Main Loop for Storing Jet Data
# -----------------------------------------------------------------------
a = 0
for event in pythia(events=10):
    vectors = event.all(selection)
    sequence = cluster(vectors, R=0.4, p=-1, ep=True)
    jets = sequence.inclusive_jets()
    unclustered_particles.append(sequence.unclustered_particles())
    part_data = []
    for i, jet in enumerate(jets):
        data = (
                jet.eta, jet.phi, jet.e
               )
        part_data = return_particle_data(jet)
        if is_massless_or_isolated(jet):
            discarded_data.append(jet)
        else:
            jets_particle_eta.extend(part_data[0])
            jets_particle_phi.extend(part_data[1])
            jets_particle_energy.extend(part_data[2])
    plt.figure()
    plt.hist2d(jets_particle_eta, jets_particle_phi,
           weights=jets_particle_energy,
           range=[(-5,5),(-1*np.pi,np.pi)],
           bins=(5,5), cmap='cubehelix')
    plt.xlabel("$\eta$")
    plt.ylabel("$\phi$")
    plt.title("$ZZ$")
    plt.colorbar()
    plt.savefig("Jets_Particles_ZZ"+str(a)+".png")
    a += 1
