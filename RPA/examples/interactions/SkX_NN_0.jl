using JLD2, TightBindingToolkit, LinearAlgebra

#####* Define the honeycomb lattice unit vectors
params = Dict()
SkXSize = get!(params, "SkXSize", 2)
SkX = get!(params, "SkX", "Bloch")
SkX = "Bloch"
a1 = SkXSize / 2 * [-3.0, sqrt(3)]
a2 = SkXSize / 2 * [3.0, sqrt(3)]
l1 = [1.0, 0]
l2 = [-0.5, sqrt(3) / 2]
UC = UnitCell([a1, a2], 2, 2)
##Parameters
n = get!(params, "n", 20)
kSize = 6 * n + 3
t = get!(params, "t", 1.0)
jh = get!(params, "jh", 0.0)
U = get!(params, "U", 0.0)
##### Thermodynamic parameters
filling = get!(params, "filling", 12.5/24)
T = get!(params, "T", 0.0)
t1 = -t
t1Param = Param(t1, 2)
jhParam = Param(jh, 2)
HoppingParams = [t1Param, jhParam]
su2spin = SpinMats(1 // 2)

##Adding inner-hexagon structure
for j = 0:(SkXSize-1)
    for i = 0:(SkXSize*3-1)
        AddBasisSite!(UC, i .* l1 + j .* l2)
    end
end
AddIsotropicBonds!(t1Param, UC, 1.0, su2spin[4], "t1", checkOffsetRange=1)
#####* Adding a parameter which tracks the nearest neighbor Heisenberg interaction in the spin-spin basis
const V  =   +1.0
VParam   =   Param(V, 2)
AddIsotropicBonds!(VParam, UC , 1.0,
    [1.0;;] ,
    "NN density-density repulsion")


values = Dict()

params = [VParam]
values["NN_repulsive_density-density"] = [1.0]
#values["NN_attractive_density-density"] = [0.0, 0.0, -1.0]
#####* Saving the unit cell in a JLD2 file
file_name = "/scratch/a/aparamek/andykh/Data/Monolayer_Data/RPA/SkX_NN_0.jld2"
#save(file_name, Dict("parameters" => params))
save(file_name, Dict("parameters" => params, "values" => values))
