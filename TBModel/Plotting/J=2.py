import h5py as h5
import yaml as yml
from pathlib import Path
import os,sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap

loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Bilayer_Data/"
t1 = -1.0
os.chdir("/home/andrewhardy/Documents/Graduate/Codes/Skyrmion/TBModel/Plotting")
os.getcwd()

 
filename = "06.27-27.2024_Bilayer"  
filename = "07.19.2024_Bilayer"  
#filename = "07.13-29.2024_Bilayer"
#filename = "07.15-27.2024_Bilayer"
#filename = "07.20-25.2024_Bilayer"
filename = "07.21-25.2024_Bilayer"
#filename = "07.19-29.2024_Bilayer"
filename = "07.24-25.2024_Bilayer"
filename = "07.31-25.2024_Bilayer"
#filename = "07.31_4-25.2024_Bilayer"
filename = "08.02-25.2024_Bilayer"
filename = "09.03-25.2024_Bilayer"

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


plt.style.use("lake.mplstyle")
plt.rcParams.update({"text.usetex": True})

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"] + plt.rcParams["font.serif"]
plt.rcParams.update({"text.usetex": True})
params = yml.safe_load(Path(f"../Input/{filename}.yml").read_text())
U_array = np.linspace(params["U_min"], params["U_max"], params["U_length"])
#filling_arr = np.linspace(params["filling_min"], params["filling_max"], params["filling_length"]) / (params["filling_max"] *2 )
V_array = np.linspace(params["V_min"], params["V_max"], params["V_length"])
V_array = np.linspace(params["V_min"], params["V_max"], params["V_length"])

Uniform_Status = False
polarization = np.zeros((params["U_length"], params["V_length"]))
energy = np.zeros((params["U_length"], params["V_length"]))
conduct = np.zeros((params["U_length"], params["V_length"]))
gap = np.zeros((params["U_length"], params["V_length"]))
ssf = np.zeros((params["U_length"], params["V_length"]))

for (ind_V,V) in enumerate(V_array):
    for (ind_u, U_var) in enumerate(U_array):
        #print(V)
        fileName = loc + f"Last_Itr_{filename}_V={round(V, 3)}_U={round(U_var, 2)}.jld2"
        #print(fileName)
        try:
            TBResults = h5.File(fileName, 'r')
            conduct[ind_u, ind_V] =  np.mean(TBResults["Chern Fill"])
            energy[ind_u, ind_V] = TBResults["MFT_Energy"][-1]
            gap[ind_u, ind_V]    = float(TBResults[TBResults["Bands"][1]][30]-TBResults[TBResults["Bands"][1]][28])#np.mean(TBResults["Gap"])np.abs(TBResults["Gap"])#

            temp_up = np.array(TBResults["ssf_up"])
            temp_dn = np.array(TBResults["ssf_dn"])
            ssf_up = temp_up['re'] + 1j*temp_up['im']
            ssf_dn = temp_dn['re'] + 1j*temp_dn['im']
            #polarization[ind_u, ind_V] = np.max(np.abs(ssf_up-ssf_dn))
            polarization[ind_u, ind_V] = np.abs(np.sum(TBResults["Expectations"][-12:]) - np.sum(TBResults["Expectations"][-24:-12]))
            ssf[ind_u,ind_V] = np.max(np.abs(ssf_up-ssf_dn))

            # if V < 0.3:
            #     conduct[ind_u, ind_V] = 0.000
            # if polarization[ind_u, ind_V] < 0.05:
            #     polarization[ind_u, ind_V] = 0.000
            #print("everything OK?")# - ssf_dn

        except:
            print(fileName)
            if ssf[ind_u,ind_V] < 0.001:
                print("IN HERE")
                ssf[ind_u,ind_V] =  ssf[ind_u,ind_V] = (ssf[ind_u,ind_V-1])
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
plt.xlabel(r'$V$')

plt.show()
plt.scatter(V_array_flat, U_array_flat,c=polarization_flat, cmap='viridis')
plt.colorbar(label=r'$P$')
plt.ylabel(r'$U$')
plt.xlabel(r'$V$')

plt.show()
plt.scatter(V_array_flat, U_array_flat,c=energy_flat, cmap='viridis')
plt.colorbar(label=r'$E$')
plt.ylabel(r'$U$')
plt.xlabel(r'$V$')

plt.show()
fig = plt.figure(figsize=(8, 8))
# plt.imshow(conduct, aspect='auto', cmap='PRGn',vmin = -1.5, vmax=1.5, origin='lower',
#            extent=[V_array.min(), V_array.max(), U_array.min(), U_array.max()])
plt.imshow(np.abs(conduct), aspect='auto', cmap='PuBuGn',vmin = 0, vmax=2, origin='lower',
           extent=[V_array.min(), V_array.max(), U_array.min(), U_array.max()])

plt.colorbar(label=r'$\sigma_{xy}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$V$')
plt.ylim(U_array.min(), min(U_array.max(), 7))

plt.savefig("Plots/Bilayer_Conductivity_Extended.pdf")
plt.show()
fig = plt.figure(figsize=(8, 8))
plt.imshow(polarization, aspect='auto', cmap='viridis', origin='lower', vmax = 0.01,
           extent=[V_array.min(), V_array.max(), U_array.min(), U_array.max()])
plt.colorbar(label=r'$N(k)_{max}$')
plt.ylabel(r'$U$')
plt.xlabel(r'$V$')
plt.ylim(U_array.min(), min(U_array.max(), 7))

plt.savefig("Plots/Bilayer_Polarization_Extended.pdf")

plt.show()

fig = plt.figure(figsize=(8, 8))
plt.imshow(gap, aspect='auto', cmap='viridis', origin='lower',
           extent=[V_array.min(), V_array.max(), U_array.min(), U_array.max()])
plt.colorbar(label=r'$\Delta$')
plt.ylabel(r'$U$')
plt.xlabel(r'$V$')

plt.savefig("Plots/Bilayer_Polarization_Extended.pdf")

plt.show()
fig = plt.figure(figsize=(18, 8))
gs = fig.add_gridspec(1, 5, width_ratios=[1,1, 0.05, 0.2, 0.05], wspace = 0.05)

ax1 = fig.add_subplot(gs[0, 0])
cax1 = fig.add_subplot(gs[0, 2])
ax2 = fig.add_subplot(gs[0, 1])
cax2 = fig.add_subplot(gs[0, 4])

# Main plot (left subplot)
Vidx = np.searchsorted(V_array, 1.5)
Uidx = np.searchsorted(U_array, 1.5)

X, Y = np.meshgrid(V_array[15:Vidx], U_array[:Uidx])
blues = plt.get_cmap("Blues")

# Create a new colormap with white at the bottom
colors = blues(np.linspace(0, 1, 256))
colors[0] = [1, 1, 1, 1]  # Set the first color to white
new_cmap = LinearSegmentedColormap.from_list("new_Blues", colors)

im1 = ax1.pcolormesh(X, Y, ssf[:Uidx,15:Vidx], cmap=new_cmap, vmin=0, vmax=0.521, shading='auto')
ax1.set_ylabel(r'$U$', fontsize=30)
ax1.set_xlabel(r'$V$', fontsize=30)
ax1.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
ax1.tick_params(axis='both', which='major', labelsize=25)
ax1.tick_params(axis='both', which='minor', labelsize=20)
ax1.yaxis.set_major_locator(MaxNLocator(nbins=5))  # Set maximum number of y-axis ticks

# Second plot (right subplot)
im2 = ax2.pcolormesh(X, Y, np.abs(conduct[:Uidx,15:Vidx]), cmap='Purples', vmin=0, vmax=1, shading='auto')
ax2.set_xlabel(r'$V$', fontsize=30)
ax2.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
ax2.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
ax2.tick_params(axis='both', which='major', labelsize=25)
ax2.tick_params(axis='both', which='minor', labelsize=20)
ax2.set_yticklabels([])  # Remove y-tick labels for the second plot

# Colorbars
cbar1 = plt.colorbar(im1, cax=cax1, label=r'$\Delta n$')
cbar1.ax.yaxis.label.set_size(30)
cbar2 = plt.colorbar(im2, cax=cax2, label=r'$\sigma_{xy}$')
cbar2.ax.yaxis.label.set_size(30)
cbar1.ax.tick_params(labelsize=25)  # Increase colorbar tick label size

cbar2.ax.tick_params(labelsize=25)  # Increase colorbar tick label size

# Rasterize
ax1.set_rasterized(True)
ax2.set_rasterized(True)

plt.savefig("Plots/"+filename+"_U_V.pdf", format='pdf', bbox_inches="tight")
plt.show()

fig, ax = plt.subplots(figsize=(8, 8))
ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
# Main plot
Vidx = np.searchsorted(V_array, 2.0)
Uidx = np.searchsorted(U_array, 1.5)

X, Y = np.meshgrid(V_array[:Vidx], U_array[:Uidx])

# im = ax.imshow(gap[:,:Vidx], aspect='auto', cmap='PuBuGn', origin='lower',
#                extent=[V_array.min(), V_array[Vidx], U_array.min(), U_array.max()])
blues = plt.get_cmap("Blues")

# Create a new colormap with white at the bottom
colors = blues(np.linspace(0, 1, 256))
colors[0] = [1, 1, 1, 1]  # Set the first color to white
new_cmap = LinearSegmentedColormap.from_list("new_Blues", colors)

im = ax.pcolormesh(X, Y, ssf[:Uidx,:Vidx], cmap=new_cmap, vmin=0, vmax=0.521, shading='auto')

ax.set_ylabel(r'$U$')
ax.set_xlabel(r'$V$')
# Inset plot
axins = inset_axes(ax, width="49%", height="49%", loc='lower left', bbox_to_anchor=(0.08, 0.05, 1.08, 1.05 ),bbox_transform=ax.transAxes)
#im_ins = axins.imshow(np.abs(conduct[:,:Vidx]), aspect='auto', cmap='PuBuGn',vmin = 0, vmax=2, origin='lower',
#                      extent=[V_array.min(), V_array[Vidx], U_array.min(), U_array.max()])
im_ins = axins.pcolormesh(X, Y, np.abs(conduct[:Uidx,:Vidx]), cmap='Purples', vmin=0, vmax=1, shading='auto')
#shading='gouraud'
axins.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
axins.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
cax_1 = inset_axes(ax,
                 width="5%",  # width = 5% of parent_bbox width
                 height="100%",  # height : 100%
                 loc='right',
                 bbox_to_anchor=(0.15, 0., 1, 1),
                 bbox_transform=ax.transAxes,
                 borderpad=0,
                 )
cax_2 = inset_axes(ax,
width="5%",  # width = 5% of parent_bbox width
height="100%",  # height : 100%
loc='right',
bbox_to_anchor=(0.4, 0., 1, 1),
bbox_transform=ax.transAxes,
borderpad=0,
)
ax.set_rasterized(True)

axins.set_rasterized(True)

cax_1.set_rasterized(True)
cax_2.set_rasterized(True)

plt.colorbar(im_ins, cax=cax_2, label=r'$\sigma_{xy}$')
plt.colorbar(im, cax=cax_1, label=r'$\Delta n(Q)$')

#plt.savefig("Plots/"+filename+"_U_V.pdf", format = 'pdf')

plt.show()
u_idx = 9
v_idx = 25
from matplotlib import cycler
colors = plt.cm.Set1(np.linspace(0, 1, 9))  # "Set1" has 9 distinct colors
plt.rc('axes', prop_cycle=(cycler('color', colors)))
fig, ax = plt.subplots(figsize=(8, 8))
ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
p = polarization[u_idx,v_idx:]#/np.max(polarization[0,:])
s = ssf[u_idx,v_idx:]
print(U_array[u_idx])
#for i in range(len(p)):

#     if (1-p[i]) < 0.05:
#         try: 
#             p[i] = 1/2*(p[i-1]+p[i+1])
#         except:
#             p[i] = 1
#     if np.abs(conduct[u_idx,i]) > 2:
#         conduct[0,i] = conduct[0,i+1]*1.1
#     if np.abs(conduct[u_idx,i]) > 0.1:
#         gap[0,i] = 0.0
#plt.plot(V_array, gap[u_idx,:]/np.max(gap[u_idx,:]), linewidth = 2, markersize = 10,marker = "^", label = r"$\Delta$")
plt.plot(V_array[v_idx:], np.abs(conduct[u_idx,v_idx:]), linewidth = 2.5, markersize = 8,marker = "o", label = r"$|\sigma_{xy}|$", color = colors[3], markeredgecolor='black', markeredgewidth=0.65)
plt.plot(V_array[v_idx:], 100*p/12, linewidth = 2.5, markersize = 8,marker = "X", label = r"$100 \times \Delta n$",color = colors[1], markeredgecolor='black', markeredgewidth=0.65)
plt.plot(V_array[v_idx:], s, linewidth = 2.5, markersize = 8,marker = "^", label = r"$\Delta n(Q)$", color = colors[2], markeredgecolor='black', markeredgewidth=0.65)

plt.legend(fontsize = 20)
plt.xlabel(r"$V$", fontsize = 20)
plt.savefig("Plots/"+filename+"_U.pdf", format = 'pdf')
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

