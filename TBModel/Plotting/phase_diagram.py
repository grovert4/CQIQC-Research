import h5py as h5
import yaml as yml
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Bilayer_Data/"
t1 = -1.0
filename = "02.15.2024_Bilayer"
plt.style.use("lake.mplstyle")
plt.rcParams.update({"text.usetex": True})

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"] + plt.rcParams["font.serif"]
plt.rcParams.update({"text.usetex": True})
params = yml.safe_load(Path(f"../Input/{filename}.yml").read_text())
U_array = np.linspace(params["U_min"], params["U_max"], params["U_length"])
filling_arr = np.linspace(params["filling_min"], params["filling_max"], params["filling_length"]) / 48
filling = filling_arr[6]
Uniform_Status = True
polarization = np.zeros((len(U_array), len(filling_arr)))
energy = np.zeros((len(U_array), len(filling_arr)))
conduct = np.zeros((len(U_array), len(filling_arr)))
for (ind_f, filling) in enumerate(filling_arr):
    for (ind_u, U_var) in enumerate(U_array):
        if Uniform_Status == True:
            fileName = loc + f"Last_Itr/Last_Itr_{filename}_UNIFORM_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
        else:
            fileName = loc + f"Last_Itr/Last_Itr_{filename}_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
        TBResults = h5.File(fileName, 'r')
        conduct[ind_u, ind_f] =  np.abs(TBResults["Chern Fill"])
        energy[ind_u, ind_f] = TBResults["MFT_Energy"][-1]
        if Uniform_Status == True:
            polarization[ind_u, ind_f] = np.abs(TBResults["Outputs"][0])
        else:
            polarization[ind_u, ind_f] = np.mean(np.abs(TBResults["Outputs"]))


        # need to make this python compatible 
U_array_flat = U_array.repeat(len(filling_arr))
filling_arr_flat = np.tile(filling_arr, len(U_array))
conduct_flat = conduct.flatten()
polarization_flat = polarization.flatten()
energy_flat = energy.flatten()
# Create the scatter plot
plt.scatter(filling_arr_flat, U_array_flat,c=conduct_flat, cmap='viridis', vmax = 3)
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

Uniform_Status = False
polarization = np.zeros((len(U_array), len(filling_arr)))
energy_2 = np.zeros((len(U_array), len(filling_arr)))
conduct = np.zeros((len(U_array), len(filling_arr)))
for (ind_f, filling) in enumerate(filling_arr):
    for (ind_u, U_var) in enumerate(U_array):
        if Uniform_Status == True:
            fileName = loc + f"Last_Itr/Last_Itr_{filename}_UNIFORM_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
        else:
            fileName = loc + f"Last_Itr/Last_Itr_{filename}_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
        TBResults = h5.File(fileName, 'r')
        conduct[ind_u, ind_f] =  np.abs(TBResults["Chern Fill"])
        energy_2[ind_u, ind_f] = TBResults["MFT_Energy"][-1]
        if Uniform_Status == True:
            polarization[ind_u, ind_f] = np.abs(TBResults["Outputs"][0])
        else:
            polarization[ind_u, ind_f] = np.mean(np.abs(TBResults["Outputs"]))


conduct_flat = conduct.flatten()
polarization_flat = polarization.flatten()
energy_2_flat = energy_2.flatten()
# Create the scatter plot
plt.scatter(filling_arr_flat, U_array_flat,c=conduct_flat, cmap='viridis', vmax =1 )
plt.colorbar(label=r'$\sigma_{xy}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')

plt.show()
plt.scatter(filling_arr_flat, U_array_flat,c=polarization_flat, cmap='viridis')
plt.colorbar(label=r'$P$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')

plt.show()
plt.scatter(filling_arr_flat, U_array_flat,c=energy_2_flat, cmap='viridis')
plt.colorbar(label=r'$E$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')

plt.show()

plt.scatter(filling_arr_flat, U_array_flat,c=(energy_flat-energy_2_flat), cmap='viridis')
plt.colorbar(label=r'$\Delta E$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')

plt.show()

plt.imshow(conduct, aspect='auto', cmap='viridis', vmax=1, origin='lower',
           extent=[filling_arr.min(), filling_arr.max(), U_array.min(), U_array.max()])

plt.colorbar(label=r'$\sigma_{xy}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')
plt.show()
plt.imshow(polarization, aspect='auto', cmap='viridis', origin='lower',
           extent=[filling_arr.min(), filling_arr.max(), U_array.min(), U_array.max()])
plt.colorbar(label=r'$P$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')
plt.show()