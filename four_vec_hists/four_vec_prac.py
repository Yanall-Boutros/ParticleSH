# ----------------------------------------------------------------------
# Import Statements
# ----------------------------------------------------------------------
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
# ----------------------------------------------------------------------
# Function Definitions
# ----------------------------------------------------------------------
def get_event_number(line):
   if line[0] != 'E': return -1
   for i in range(len(line)):
      if line[i].isdigit():
         return int(line[i:])
def return_values(line):
   if line[0] != 'E':
      trimmed = line.split(" ")
      trimmed = (float(trimmed[0]), float(trimmed[2]), float(trimmed[4]),
                float(trimmed[6]))
      return trimmed
def gaussian(data):
   std_dev = np.std(vecs)
   avg = np.average(vecs)
   x = np.linspace(-3000, 3000, 6001)
   coeff = (2.0*np.pi*std_dev**2)**-0.5
   return coeff*np.exp((-1*((x-avg)**2))/(2*std_dev**2))
def calc_mom(vec):
   mom = 0
   for elem in vec[:3]:
      mom += (elem)**2
   return np.sqrt(mom)
def calc_mass(vec):
   i_mass = vec[3] * vec[3]
   p = calc_mom(vec)
   return ((i_mass + (p**2)) ** 0.5)
# ----------------------------------------------------------------------
# Data Initialization and Calculations
# ----------------------------------------------------------------------
with open("events2.out") as f:
   content = f.readlines()
content = [x.strip() for x in content]
# Go through the file line by line. If there's a next event, then add
# the two events together and append to the data
vecs = list()
for i in range(len(content)):
   line = content[i]
   if i < len(content)-1:
      next_line = content[i+1]
   if len(line) > 2 and len(next_line) > 2:
      if line[0].isdigit() or line[0] == '-':
         if next_line[0].isdigit() or next_line[0] == '-':
            vec_a = return_values(line)
            vec_b = return_values(next_line)
            vec_sum = np.add(np.array(vec_a), np.array(vec_b))
            vecs.append(calc_mass(vec_sum))
vecs = np.array(vecs)
# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------
x = np.linspace(-3000, 3000, 6001)
y = gaussian(vecs)
plt.hist(vecs, 600, range=[0, 1000])
plt.title("Hist. of Inv. Mass"), plt.xlabel("Mass [GeV]"), plt.ylabel("Events")
plt.savefig("hist.pdf")
plt.figure()
plt.plot(x, y)
plt.title("Gaussian Distribution")
plt.xlim(0, 1000)
plt.xlabel("Event Number")
plt.ylabel("Mass")
plt.savefig("gass.pdf")
