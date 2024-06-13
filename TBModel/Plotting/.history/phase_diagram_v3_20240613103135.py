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


#filename = "05.02-0.5.2024_Bilayer"  
filename = "05.03-0.4.2024_Bilayer"  
filename = "05.04-0.4.2024_Bilayer"  
filename = "05.04-0.66.2024_Bilayer"  

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


plt.style.use("lake.mplstyle")
plt.rcParams.update({"text.usetex": True})

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"] + plt.rcParams["font.serif"]
plt.rcParams.update({"text.usetex": True})
params = yml.safe_load(Path(f"../Input/{filename}.yml").read_text())
U_array = np.linspace(params["U_min"], params["U_max"], params["U_length"])
filling_arr = np.linspace(params["filling_min"], params["filling_max"], params["filling_length"]) / (params["filling_max"] *2 )
V_array = np.linspace(params["V_min"], params["V_max"], params["V_length"])

Uniform_Status = False
polarization = np.zeros((params["U_length"], params["V_length"]))
energy = np.zeros((params["U_length"], params["V_length"]))
conduct = np.zeros((params["U_length"], params["V_length"]))
for (ind_V,V) in enumerate(V_array):
    for (ind_u, U_var) in enumerate(U_array):
        print(V)
        fileName = loc + f"Last_Itr_{filename}_V={round(V, 3)}_U={round(U_var, 2)}.jld2"
        print(fileName)
        try:
            TBResults = h5.File(fileName, 'r')
            conduct[ind_u, ind_V] =  np.mean(TBResults["Chern Fill"])
            energy[ind_u, ind_V] = TBResults["MFT_Energy"][-1]
            if Uniform_Status == True:
                polarization[ind_u, ind_V] = np.abs(TBResults["Expectations"][3:])
            else:
                temp_up = np.array(TBResults["ssf_up"])
                temp_dn = np.array(TBResults["ssf_dn"])
                ssf_up = temp_up['re'] + 1j*temp_up['im']
                ssf_dn = temp_dn['re'] + 1j*temp_dn['im']

                polarization[ind_u, ind_V] = np.abs(np.max(ssf_up)-np.max(ssf_dn))
        except:
            continue

        # need to make this python compatible 
U_array_flat = U_array.repeat(len(V_array))
V_array_flat = np.tile(V_array, len(U_array))
conduct_flat = conduct.flatten()
polarization_flat = polarization.flatten()
energy_flat = energy.flatten()
# Create the scatter plot
plt.scatter(V_array_flat, U_array_flat,c=conduct_flat, cmap='viridis', vmax = 4)
plt.colorbar(label=r'$\sigma_{xy}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')

plt.show()
plt.scatter(V_array_flat, U_array_flat,c=polarization_flat, cmap='viridis')
plt.colorbar(label=r'$P$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')

plt.show()
plt.scatter(V_array_flat, U_array_flat,c=energy_flat, cmap='viridis')
plt.colorbar(label=r'$E$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')

plt.show()
fig = plt.figure(figsize=(8, 8))
plt.imshow(conduct, aspect='auto', cmap='PRGn',vmin = -1.5, vmax=1.5, origin='lower',
           extent=[V_array.min(), V_array.max(), U_array.min(), U_array.max()])

plt.colorbar(label=r'$\sigma_{xy}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')
plt.ylim(U_array.min(), min(U_array.max(), 7))

plt.savefig("Plots/Bilayer_Conductivity_Extended.pdf")
plt.show()
fig = plt.figure(figsize=(8, 8))
plt.imshow(polarization, aspect='auto', cmap='viridis', origin='lower',
           extent=[V_array.min(), V_array.max(), U_array.min(), U_array.max()])
plt.colorbar(label=r'$N(k)_{max}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$V$')
plt.ylim(U_array.min(), min(U_array.max(), 7))

plt.savefig("Plots/Bilayer_Polarization_Extended.pdf")

plt.show()

fig, ax = plt.subplots(figsize=(8, 8))

# Main plot
im = ax.imshow(polarization, aspect='auto', cmap='viridis', origin='lower',
               extent=[V_array.min(), V_array.max(), U_array.min(), U_array.max()])
ax.set_ylabel(r'$U$')
ax.set_xlabel(r'$V$')
ax.set_ylim(U_array.min(), min(U_array.max(), 7))
plt.colorbar(im, ax=ax, label=r'$N(k)_{max}$')

# Inset plot
axins = inset_axes(ax, width="50%", height="50%", loc='lower right')
im_ins = axins.imshow(conduct, aspect='auto', cmap='PRGn',vmin = -2, vmax=2, origin='lower',
                      extent=[V_array.min(), V_array.max(), U_array.min(), U_array.max()])
cax = inset_axes(axins,
                 width="5%",  # width = 5% of parent_bbox width
                 height="100%",  # height : 100%
                 loc='right',
                 bbox_to_anchor=(1.05, 0., 1, 1),
                 bbox_transform=axins.transAxes,
                 borderpad=0,
                 )

plt.colorbar(im_ins, cax=cax, label=r'$\sigma_{xy}$')

plt.savefig("Plots/Bilayer_Polarization_Extended.pdf")

plt.show()



# Uniform_Status = False
# polarization = np.zeros((len(U_array), len(filling_arr)))
# energy_2 = np.zeros((len(U_array), len(filling_arr)))
# conduct = np.zeros((len(U_array), len(filling_arr)))
# for (ind_V, filling) in enumerate(filling_arr):
#     for (ind_u, U_var) in enumerate(U_array):
#         if Uniform_Status == True:
#             fileName = loc + f"Last_Itr/Last_Itr_{filename}_UNIFORM_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
#         else:
#             fileName = loc + f"Last_Itr/Last_Itr_{filename}_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
#         TBResults = h5.File(fileName, 'r')
#         conduct[ind_u, ind_V] =  np.abs(TBResults["Chern Fill"])
#         energy_2[ind_u, ind_V] = TBResults["MFT_Energy"][-1]
#         if Uniform_Status == True:
#             polarization[ind_u, ind_V] = np.abs(TBResults["Expectations"][0])
#         else:
#             polarization[ind_u, ind_V] = np.mean(np.abs(TBResults["Expectations"]))


# conduct_flat = conduct.flatten()
# polarization_flat = polarization.flatten()
# energy_2_flat = energy_2.flatten()
# # Create the scatter plot
# plt.scatter(filling_arr_flat, U_array_flat,c=conduct_flat, cmap='viridis', vmax =1 )
# plt.colorbar(label=r'$\sigma_{xy}$')
# plt.ylabel(r'$U$')
# plt.xlabel(r'$n$')

# plt.show()
# plt.scatter(filling_arr_flat, U_array_flat,c=polarization_flat, cmap='viridis')
# plt.colorbar(label=r'$P$')
# plt.ylabel(r'$U$')
# plt.xlabel(r'$n$')

# plt.show()
# plt.scatter(filling_arr_flat, U_array_flat,c=energy_2_flat, cmap='viridis')
# plt.colorbar(label=r'$E$')
# plt.ylabel(r'$U$')
# plt.xlabel(r'$n$')

# plt.show()

# plt.scatter(filling_arr_flat, U_array_flat,c=(energy_flat-energy_2_flat), cmap='viridis')
# plt.colorbar(label=r'$\Delta E$')
# plt.ylabel(r'$U$')
# plt.xlabel(r'$n$')

# plt.show()

