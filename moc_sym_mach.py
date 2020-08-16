import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

from numpy.random import uniform, seed
import pandas as pd


M_unit = 0.4
M_exit = 2.4

boundaries=np.linspace(M_unit, M_exit, 11, endpoint=True)




mach = pd.read_csv("point_folder/machcontour25.csv", delim_whitespace=True)

x = mach['x']
y = mach['y']
z = mach['Ma']
# define grid.
xi = np.linspace(x.min(), x.max(), 1000)
yi = np.linspace(y.min(), y.max(), 1000)

# grid the data.
zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')

npts = 200
clev = np.linspace(M_unit,M_exit,npts)

plt.figure(figsize=(10,9))

CS = plt.contourf(xi,yi,zi,clev, cmap=plt.cm.jet, extend='both') #
CS1 = plt.contourf(xi,-yi,zi,clev, cmap=plt.cm.jet, extend='both') #

for c in CS.collections:
    c.set_edgecolor("face")
for c in CS1.collections:
    c.set_edgecolor("face")
    
plt.colorbar(ticks=boundaries, label="Mach Number", orientation="horizontal") 
plt.xlabel("Nozzle length (m)",fontsize=16) #r"Nozzle length  $\frac{x}{x0}$",fontsize=16
plt.ylabel("Nozzle radius (m)",fontsize=16) #r"Nozzle radius  $\frac{y}{y0}$",fontsize=16
# plt.title('' % npts)
#plt.show()

file_format = "png"
file_name = "moc_sym."+file_format
plt.savefig(("plots/"+file_name), format=file_format, dpi=600) #, 
print("Image file ",file_name," has been generated!")