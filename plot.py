#! ../plot_robert/plot.py
#--------------------------------------------------------------------------------
# plot_robert.pro: This program is used to plot the potential temperature anomaly
# of 2d-robert dry bubble test.
#
# This program is written by Huazhi Ge.
#--------------------------------------------------------------------------------

# Modules
import os
import glob
from os import listdir
from os.path import isfile, join
import numpy as np
import matplotlib.pyplot as plt
import moviepy.editor as mpy

from netCDF4 import MFDataset
from netCDF4 import Dataset

# Constants
Rd = 287
Cp = 1004
Ps = 1E5
dt = 10
CFL = 0.4

# Path for output files.
figurepath  = "./figures/"
animatepath = "./animation/"

# Read the name of the block files from the 2d-bryan directory.
dirpath = "../Explicit/"
files = [f for f in listdir(dirpath) if isfile(join(dirpath, f))]

os.chdir("./figures/")
os.system("rm -f *.png")
os.chdir("../animation")
os.system("rm -f *.gif *.mp4")
os.chdir("../")

for file in files:
   if ("straka-" in file) and (".nc" in file):
      filename = dirpath + file
      data = MFDataset(filename)
      info = Dataset(filename, "r", format="NETCDF3_CLASSIC")
      print(info)
      time = data.variables["time"][:]
      x1  = data.variables["x1"][:]
      x2  = data.variables["x2"][:]
      theta = data.variables["theta"][:]-299.85

      # Plot figures
      level = np.linspace(-17, 1, num = 19)
      # Name the plots
      nfile = 0
      for itime in time:
         plotname = figurepath + "Figure" + str(nfile).zfill(4)
         # Plot the figure with the plots of density, pressure and veolocity.
         plt.figure(nfile)
         plt.figure(figsize = (20, 8))
         axes = plt.gca()
         axes.set_xlim([np.min(x2), 20000])
         axes.set_ylim([np.min(x1), 4000])
         title = "Straka Bubble Test: Potential Temperature at " + str(itime) + \
                 "s, CFL = " + str(CFL) #+ "(" + \
                 #str(CFL*(25600./x2.size)/(6400./x1.size)) + ")"
         plt.title(title, fontsize=20)
         plt.xticks([0, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000])
         plt.yticks([0, 1000, 2000, 3000, 4000])
         fig = plt.contour(x2, x1, theta[nfile,0,:,:].transpose(), \
               levels = level, colors="black", linestyles="solid")
         plt.xlabel("X (Horizontal) m", fontsize=20)
         plt.ylabel("Y (Vertical) m", fontsize=20)
         plt.text( 15000, 3000, r'max $\theta$ = ' + "{:6.3f}".format(np.max(theta)) \
                   + "\n" + r'min $\theta$ = ' + "{:6.3f}".format(np.min(theta)), size = 20, \
                   ha="center", va="center" )
         snapshot = plt.gcf()
         snapshot.savefig(plotname, dpi=100)
         plt.close()
         nfile += 1        
 
      data.close()

# Create animation from the figures.
gifname = "Straka"
figurelist = glob.glob(figurepath + "*.png")
figurelist.sort()
gif = mpy.ImageSequenceClip(figurelist, fps = 15)
gif.write_gif(animatepath + '{}.gif'.format(gifname), fps = 15)
   
   
