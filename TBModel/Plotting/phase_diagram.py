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



plt.style.use("lake.mplstyle")
plt.rcParams.update({"text.usetex": True})

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"] + plt.rcParams["font.serif"]
plt.rcParams.update({"text.usetex": True})
params = yml.safe_load(Path(f"../Input/{filename}.yml").read_text())
U_array = np.linspace(params["U_min"], params["U_max"], params["U_length"])
filling_arr = np.linspace(params["filling_min"], params["filling_max"], params["filling_length"]) / (params["filling_max"] *2 )
filling = filling_arr[1]
Uniform_Status = False
polarization = np.zeros((len(U_array), len(filling_arr)))
energy = np.zeros((len(U_array), len(filling_arr)))
conduct = np.zeros((len(U_array), len(filling_arr)))
for (ind_f, filling) in enumerate(filling_arr[1:]):
    for (ind_u, U_var) in enumerate(U_array):
        if Uniform_Status == True:
            fileName = loc + f"Last_Itr_{filename}_UNIFORM_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
        else:
            fileName = loc + f"Last_Itr_{filename}_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
        TBResults = h5.File(fileName, 'r')
        conduct[ind_u, ind_f] =  np.mean(TBResults["Chern Fill"])
        energy[ind_u, ind_f] = TBResults["MFT_Energy"][-1]
        if Uniform_Status == True:
            polarization[ind_u, ind_f] = np.abs(TBResults["Expectations"][3:])
        else:
            polarization[ind_u, ind_f] = np.abs(TBResults["Expectations"])[5]#np.mean(np.abs(TBResults["Expectations"]))


        # need to make this python compatible 
U_array_flat = U_array.repeat(len(filling_arr))
filling_arr_flat = np.tile(filling_arr, len(U_array))
conduct_flat = conduct.flatten()
polarization_flat = polarization.flatten()
energy_flat = energy.flatten()
# Create the scatter plot
plt.scatter(filling_arr_flat, U_array_flat,c=conduct_flat, cmap='viridis', vmax = 4)
plt.colorbar(label=r'$\sigma_{xy}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')

plt.show()
plt.scatter(filling_arr_flat, U_array_flat,c=polarization_flat, cmap='viridis')
plt.colorbar(label=r'$P$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')

plt.show()
plt.scatter(filling_arr_flat, U_array_flat,c=energy_flat, cmap='viridis')
plt.colorbar(label=r'$E$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')

plt.show()
fig = plt.figure(figsize=(8, 8))
plt.imshow(conduct, aspect='auto', cmap='PRGn',vmin = -1.5, vmax=1.5, origin='lower',
           extent=[filling_arr.min(), filling_arr.max(), U_array.min(), U_array.max()])

plt.colorbar(label=r'$\sigma_{xy}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')
plt.ylim(U_array.min(), min(U_array.max(), 7))

plt.savefig("Plots/Bilayer_Conductivity_Extended.pdf")
plt.show()
fig = plt.figure(figsize=(8, 8))
plt.imshow(polarization, aspect='auto', cmap='viridis', origin='lower',
           extent=[filling_arr.min(), filling_arr.max(), U_array.min(), U_array.max()])
plt.colorbar(label=r'$P$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')
plt.ylim(U_array.min(), min(U_array.max(), 7))

plt.savefig("Plots/Bilayer_Polarization_Extended.pdf")

plt.show()

# Uniform_Status = False
# polarization = np.zeros((len(U_array), len(filling_arr)))
# energy_2 = np.zeros((len(U_array), len(filling_arr)))
# conduct = np.zeros((len(U_array), len(filling_arr)))
# for (ind_f, filling) in enumerate(filling_arr):
#     for (ind_u, U_var) in enumerate(U_array):
#         if Uniform_Status == True:
#             fileName = loc + f"Last_Itr/Last_Itr_{filename}_UNIFORM_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
#         else:
#             fileName = loc + f"Last_Itr/Last_Itr_{filename}_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
#         TBResults = h5.File(fileName, 'r')
#         conduct[ind_u, ind_f] =  np.abs(TBResults["Chern Fill"])
#         energy_2[ind_u, ind_f] = TBResults["MFT_Energy"][-1]
#         if Uniform_Status == True:
#             polarization[ind_u, ind_f] = np.abs(TBResults["Expectations"][0])
#         else:
#             polarization[ind_u, ind_f] = np.mean(np.abs(TBResults["Expectations"]))


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

