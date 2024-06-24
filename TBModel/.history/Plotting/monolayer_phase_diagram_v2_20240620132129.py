import h5py as h5
import yaml as yml
from pathlib import Path
import os,sys
import numpy as np
import matplotlib.pyplot as plt
loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Monolayer_Data/"
t1 = -1.0
filename = "02.05.2024_Monolayer" # J = 4, N = 2
#filename = "04.05.2024_Monolayer_E" # J = 4, N = 2
#filename = "05.05.2024_Monolayer_4" # J = 4, N = 2

filename = "05.01-0.5.2024_Monolayer"
filename = "06.19.2024_Monolayer"


os.chdir("/home/andrewhardy/Documents/Graduate/Codes/Skyrmion/TBModel/Plotting")
os.getcwd()
sys.path.append(os.getcwd())
plt.style.use("lake.mplstyle")
plt.rcParams.update({"text.usetex": True})

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"] + plt.rcParams["font.serif"]
plt.rcParams.update({"text.usetex": True})
params = yml.safe_load(Path(f"../Input/{filename}.yml").read_text())
U_array = np.linspace(params["U_min"], params["U_max"], params["U_length"])
filling_arr = np.linspace(params["filling_min"], params["filling_max"], params["filling_length"]) / (params["filling_max"]*2)
#J_array = np.linspace(params["J_min"], params["J_max"], params["J_length"])
filling_array = ((12+np.linspace(params["filling_min"], params["filling_max"], params["filling_length"])) / 24 )

polarization = np.zeros((len(U_array), len(J_array)))
energy = np.zeros((len(U_array), len(J_array)))
conduct = np.zeros((len(U_array), len(J_array)))
#for (ind_J, J) in enumerate(J_array):
for (ind_n,filling) in enumerate(filling_array):
    for (ind_u, U_var) in enumerate(U_array):
        #fileName = loc + f"Last_Itr_{filename}_J={round(J, 3)}_U={round(U_var, 2)}.jld2"
        fileName = loc + f"Last_Itr_{filename}_n={round(filling, 3)}_U={round(U_var, 2)}.jld2"

        TBResults = h5.File(fileName, 'r')
        conduct[ind_u, ind_J] =  np.mean(TBResults["Chern Fill"])
        energy[ind_u, ind_J] = TBResults["MFT_Energy"][-1]
        temp = np.array(TBResults["ssf"])
        ssf = temp['re'] + 1j*temp['im']
        polarization[ind_u, ind_J] = np.abs(np.max(ssf))

U_array_flat = U_array.repeat(len(J_array))
J_array_flat = np.tile(J_array, len(U_array))
conduct_flat = conduct.flatten()
polarization_flat = polarization.flatten()
energy_flat = energy.flatten()
# Create the scatter plot
# plt.scatter(J_array_flat, U_array_flat,c=conduct_flat, cmap='viridis', vmax = 3)
# plt.colorbar(label=r'$\sigma_{xy}$')
# plt.ylabel(r'$U$')
# plt.xlabel(r'$n$')

# plt.show()
# plt.scatter(J_array_flat, U_array_flat,c=polarization_flat, cmap='viridis')
# plt.colorbar(label=r'$P$')
# plt.ylabel(r'$U$')
# plt.xlabel(r'$n$')

# plt.show()
# plt.scatter(J_array_flat, U_array_flat,c=energy_flat, cmap='viridis')
# plt.colorbar(label=r'$E$')
# plt.ylabel(r'$U$')
# plt.xlabel(r'$n$')

# plt.show()

# Uniform_Status = False
# polarization = np.zeros((len(U_array), len(J_array)))
# energy_2 = np.zeros((len(U_array), len(J_array)))
# conduct = np.zeros((len(U_array), len(J_array)))
# for (ind_J, filling) in enumerate(J_array):
#     for (ind_u, U_var) in enumerate(U_array):
#         if Uniform_Status == True:
#             fileName = loc + f"Last_Itr/Last_Itr_{filename}_UNIFORM_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
#         else:
#             fileName = loc + f"Last_Itr/Last_Itr_{filename}_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
#         TBResults = h5.File(fileName, 'r')
#         conduct[ind_u, ind_J] =  np.abs(TBResults["Chern Fill"])
#         energy_2[ind_u, ind_J] = TBResults["MFT_Energy"][-1]
#         if Uniform_Status == True:
#             polarization[ind_u, ind_J] = np.abs(TBResults["Expectations"][0])
#         else:
#             polarization[ind_u, ind_J] = np.mean(np.abs(TBResults["Expectations"]))


conduct_flat = conduct.flatten()
polarization_flat = polarization.flatten()
# Create the scatter plot
plt.scatter(J_array_flat, U_array_flat,c=conduct_flat, cmap='viridis', vmax =1 )
plt.colorbar(label=r'$\sigma_{xy}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')

plt.show()
plt.scatter(J_array_flat, U_array_flat,c=polarization_flat, cmap='viridis')
plt.colorbar(label=r'$P$')
plt.ylabel(r'$U$')
plt.xlabel(r'$n$')

plt.show()


fig = plt.figure(figsize=(8, 8))

plt.imshow(np.round(conduct,3), aspect='auto', cmap='PRGn',vmin = -6, vmax=6, origin='lower',
           extent=[J_array.min(), J_array.max(), U_array.min(), U_array.max()])

plt.colorbar(label=r'$\sigma_{xy}$')
plt.ylabel(r'$U_1$')
plt.xlabel(r'$n$')
plt.ylim(U_array.min(), min(U_array.max(), 5))

plt.savefig("Plots/Monolayer_Conductivity.pdf")

plt.show()
fig = plt.figure(figsize=(8, 8))

plt.imshow(polarization, aspect='auto', cmap='viridis', origin='lower',
           extent=[J_array.min(), J_array.max(), U_array.min(), U_array.max()])
plt.colorbar(label=r'$P_1$')
plt.ylabel(r'$U_1$')
plt.xlabel(r'$n$')
#plt.ylim(U_array.min(), min(U_array.max(), 5))

plt.savefig("Plots/Monolayer_Polarization.pdf")

plt.show()