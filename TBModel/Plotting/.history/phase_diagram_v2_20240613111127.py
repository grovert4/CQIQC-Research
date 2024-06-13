import h5py as h5
import yaml as yml
from pathlib import Path
import os,sys
import numpy as np
import matplotlib.pyplot as plt
loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Bilayer_Data/"
t1 = -1.0
os.chdir("/home/andrewhardy/Documents/Graduate/Codes/Skyrmion/TBModel/Plotting")
os.getcwd()

filename = "05.20.2024_Bilayer" # J = 0, N = 2
#filename = "05.15.2024_Bilayer"  # J = 1, N = 2
filename = "05.14.2024_Bilayer" # J = 2, N = 2
#filename = "05.12.2024_Bilayer" # J = 4, N = 2

#filename = "05.10.2024_Bilayer"  # J = 1, N = 3
#filename = "05.16.2024_Bilayer"  # J = 4, N = 3 DNE?
filename = "05.23.2024_Bilayer"  
filename = "05.01-0.5.2024_Bilayer"  

filename = "05.02-0.25.2024_Bilayer"  
#filename = "05.02-0.5.2024_Bilayer"  
filename = "05.03-0.66.2024_Bilayer"  



plt.style.use("lake.mplstyle")
plt.rcParams.update({"text.usetex": True})

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"] + plt.rcParams["font.serif"]
plt.rcParams.update({"text.usetex": True})
params = yml.safe_load(Path(f"../Input/{filename}.yml").read_text())
U_array = np.linspace(params["U_min"], params["U_max"], params["U_length"])
filling_arr = np.linspace(params["filling_min"], params["filling_max"], params["filling_length"]) / (params["filling_max"] *2 )
J_array = np.linspace(params["J_min"], params["J_max"], params["J_length"])

Uniform_Status = False
polarization = np.zeros((params["U_length"], params["J_length"]))
energy = np.zeros((params["U_length"], params["J_length"]))
conduct = np.zeros((params["U_length"], params["J_length"]))
for (ind_J,J) in enumerate(J_array):
    for (ind_u, U_var) in enumerate(U_array):
        print(J)
        fileName = loc + f"Last_Itr_{filename}_J={round(J, 3)}_U={round(U_var, 2)}.jld2"
        TBResults = h5.File(fileName, 'r')
        conduct[ind_u, ind_J] =  np.mean(TBResults["Chern Fill"])
        energy[ind_u, ind_J] = TBResults["MFT_Energy"][-1]
        if Uniform_Status == True:
            polarization[ind_u, ind_J] = np.abs(TBResults["Expectations"][3:])
        else:
            temp_up = np.array(TBResults["ssf_up"])
            temp_dn = np.array(TBResults["ssf_dn"])
            ssf_up = temp_up['re'] + 1j*temp_up['im']
            ssf_dn = temp_dn['re'] + 1j*temp_dn['im']

            polarization[ind_u, ind_J] = np.abs(np.max(ssf_up)-np.max(ssf_dn))
        if polarization[ind_u, ind_J] < 0.05
            polarization[ind_u, ind_J] = 0.0
            conduct[ind_u, ind_J] = 0.0
        # need to make this python compatible 
U_array_flat = U_array.repeat(len(J_array))
J_array_flat = np.tile(J_array, len(U_array))
conduct_flat = conduct.flatten()
polarization_flat = polarization.flatten()
energy_flat = energy.flatten()
# Create the scatter plot
plt.scatter(J_array_flat, U_array_flat,c=conduct_flat, cmap='viridis', vmax = 4)
plt.colorbar(label=r'$\sigma_{xy}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$J$')

plt.show()
plt.scatter(J_array_flat, U_array_flat,c=polarization_flat, cmap='viridis')
plt.colorbar(label=r'$P$')
plt.ylabel(r'$U$')
plt.xlabel(r'$J$')

plt.show()
plt.scatter(J_array_flat, U_array_flat,c=energy_flat, cmap='viridis')
plt.colorbar(label=r'$E$')
plt.ylabel(r'$U$')
plt.xlabel(r'$J$')

plt.show()
fig = plt.figure(figsize=(8, 8))
plt.imshow(conduct, aspect='auto', cmap='PRGn',vmin = -2, vmax=2, origin='lower',
           extent=[J_array.min(), J_array.max(), U_array.min(), U_array.max()])

plt.colorbar(label=r'$\sigma_{xy}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$J$')
plt.ylim(U_array.min(), min(U_array.max(), 7))

plt.savefig("Plots/Bilayer_Conductivity_Extended.pdf")
plt.show()
fig = plt.figure(figsize=(8, 8))
plt.imshow(np.abs(conduct), aspect='auto', cmap='Blues',vmin = 0, vmax=2, origin='lower',
           extent=[J_array.min(), J_array.max(), U_array.min(), U_array.max()])

plt.colorbar(label=r'$\sigma_{xy}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$J$')
plt.ylim(U_array.min(), min(U_array.max(), 7))

plt.savefig("Plots/Bilayer_Conductivity_Extended.pdf")
plt.show()

fig = plt.figure(figsize=(8, 8))
plt.imshow(polarization, aspect='auto', cmap='viridis', origin='lower',
           extent=[J_array.min(), J_array.max(), U_array.min(), U_array.max()])
plt.colorbar(label=r'$N(k)_{max}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$J$')
plt.ylim(U_array.min(), min(U_array.max(), 7))

plt.savefig("Plots/Bilayer_Polarization_Extended.pdf")

plt.show()
fig, ax = plt.subplots(figsize=(8, 8))

# Main plot
im = ax.imshow(polarization, aspect='auto', cmap='PuBuGn', origin='lower',
               extent=[J_array.min(), J_array.max(), U_array.min(), U_array.max()])
ax.set_ylabel(r'$U$')
ax.set_xlabel(r'$J$')
ax.set_ylim(U_array.min(), min(U_array.max(), 7))
plt.colorbar(im, ax=ax, label=r'$N(k)_{max}$')

# Inset plot
axins = inset_axes(ax, width="45%", height="45%", loc='lower self',bbox_to_anchor=(0., 0.05, 1.0, 1.05), bbox_transform=ax.transAxes)
im_ins = axins.imshow(conduct, aspect='auto', cmap='PRGn',vmin = -2, vmax=2, origin='lower',
                      extent=[J_array.min(), J_array.max(), U_array.min(), U_array.max()])
cax = inset_axes(axins,
                 width="5%",  # width = 5% of parent_bbox width
                 height="100%",  # height : 100%
                 loc='right',
                 bbox_to_anchor=(0.5, 0., 1, 1),
                 bbox_transform=ax.transAxes,
                 borderpad=0,
                 )

plt.colorbar(im_ins, cax=cax, label=r'$\sigma_{xy}$')

plt.savefig("Plots/Bilayer_Polarization_Extended.pdf")

plt.show()