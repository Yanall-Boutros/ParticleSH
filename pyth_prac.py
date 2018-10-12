# ----------------------------------------------------------------------
# Import Statements
# ----------------------------------------------------------------------
import numpy as np
import matplotlib
import pandas as pd
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpythia import Pythia, hepmc_write, hepmc_read
from numpythia import STATUS, HAS_END_VERTEX, ABS_PDG_ID
from numpythia.testcmnd import get_cmnd
from numpy.testing import assert_array_equal
# ----------------------------------------------------------------------
# Constant Definitions
# ----------------------------------------------------------------------
pythia = Pythia(get_cmnd('w'), random_state=1)
pdg_to_part = {
   1  : 'd',
   2  : 'u',
   3  : 's',
   4  : 'c',
   5  : 'b',
   6  : 't',
   7  : 'b\'',
   8  : 't\'',
   11 : 'e^-',  
   12 : 'v_e', 
   13 : 'mu^-',  
   14 : 'v_mu', 
   15 : 't^-', 
   16 : 'v_t', 
   17 : 't\'^-',  
   18 : 'v_{t\'}', 
   21 : 'g', 
   22 : 'gamma', 
   23 : 'Z^0', 
   24 : 'W^+', 
   25 : 'h^0/H_1^0', 
   32 : 'Z\'/Z_2^0', 
   33 : 'Z\'\'/Z_3^0', 
   34 : 'W\'/W_2+', 
   35 : 'H0/H_2^0', 
   36 : 'A0/H_3^0', 
   37 : 'H+', 
}
part_to_pdg = {
  'd'         : 1,
  'u'         : 2, 
  's'         : 3, 
  'c'         : 4, 
  'b'         : 5, 
  't'         : 6, 
  'b\''       : 7, 
  't\''       : 8, 
  'e^-'       : 11,
  'v_e'       : 12,
  'mu^-'      : 13,
  'v_mu'      : 14,
  't^-'       : 15,
  'v_t'       : 16,
  't\'^-'     : 17,
  'v_{t\'}'   : 18,
  'g'         : 21,
  'gamma'     : 22,
  'Z^0'       : 23,
  'W^+'       : 24,
  'h^0/H_1^0' : 25,
  'Z\'/Z_2^0' : 32,
  'Z\'/Z_3^0' : 33,
  'W\'/W_2+'  : 34,
  'H0/H_2^0'  : 35,
  'A0/H_3^0'  : 36,
  'H+'        : 37,
}
# ----------------------------------------------------------------------
# Function Definitions
# ----------------------------------------------------------------------
def gen_arrays(selection):
   print("Selection = ", selection)
   for event in hepmc_write('events.hepmc', pythia(events=1)):
      array1 = event.all(selection)
   for event in hepmc_read('events.hepmc'):
      array2 = event.all(selection)
   return np.array(array1), np.array(array2)
def set_particle_selection(params):
   # sets the filter for the selection variable based on an array of
   # string or int input vars where a string represents the particulate
   # name and the int represents the pdg_id 
   # convert to int if necessary
   pdgs = list()
   for elem in params:
      if type(elem) is str: pdgs.append(part_to_pdg[elem])
      else: pdgs.append(elem)
   selection = ((STATUS == 1) & ~HAS_END_VERTEX)
   for pdg in pdgs: selection = selection & (ABS_PDG_ID != int(pdg))
   return selection
def calc_mom(x, y, z):
   return x**2 + y**2 + z**2
def calc_mass(fvec):
   i_mass = fvec[0]**2
   p_sq = calc_mom(fvec[1], fvec[2], fvec[3])
   return ((i_mass + p_sq)**0.5)
# ----------------------------------------------------------------------
# Main Function
# ----------------------------------------------------------------------
inp = [12, 14, 16]
inv_mass = list()
selection = set_particle_selection(inp)
a1, a2 = gen_arrays(selection)
for elem in a1:
   print(elem)
   if elem[9] == 211 or elem[9] == '211':
      inv_mass.append(calc_mass(elem))
inv_mass = np.array(inv_mass)
plt.hist(inv_mass, 200, range=[0, 400])
plt.title("Histogram of Invariant Mass")
plt.xlabel("Invariant Mass [GeV]")
plt.ylabel("Counts per Event")
plt.savefig("h.png")
