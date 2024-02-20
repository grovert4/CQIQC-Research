import h5py as h5
import yaml as yml
from pathlib import Path
import numpy as np
loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Bilayer_Data/"
t1 = -1.0
filename = "02.15.2024_Bilayer"
Uniform_Status = True

params = yml.safe_load(Path(f"../Input/{filename}.yml").read_text())
U_array = np.linspace(params["U_min"], params["U_max"], params["U_length"])
filling_arr = np.linspace(params["filling_min"], params["filling_max"], params["filling_length"]) / 48
filling = filling_arr[6]
for (ind, U_var) in enumerate(U_array):
    if Uniform_Status == True:
        fileName = loc + f"Last_Itr/Last_Itr_{filename}_UNIFORM_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
    else:
        fileName = loc + f"Last_Itr/Last_Itr_{filename}_UNIFORM_p={round(filling, 3)}_U={round(U_var, 2)}_t1={round(t1, 2)}.jld2"
    TBResults = h5.File(fileName, 'r')
        # need to make this python compatible 