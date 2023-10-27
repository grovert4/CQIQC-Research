using Plots, LinearAlgebra, ColorSchemes
using MeanFieldToolkit, TightBindingToolkit
loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Bilayer_Data"

##Triangular Lattice 

const a1 = [-3.0, sqrt(3)]
const a2 = [3.0, sqrt(3)]

const l1 = [1.0, 0]
const l2 = [-0.5, sqrt(3) / 2]
UC = UnitCell([a1, a2], 4)

##Parameters
const n = 5
const kSize = 6 * n + 3
const t = 1.0
const t_inter = -0.3
const jh = -1.0
const U = 1.0
const t_density = -0.8
U_array = collect(LinRange(0.0, 10.0, 21))
SpinVec = SpinMats(1 // 2)
##### Thermodynamic parameters
const T = 0.001
const stat = -1
const mixingAlpha = 0.5
const ep = 1.0

tinter_param = Param(t_inter, 2)

const t_density = -0.8
t1 = -t
t1Param = Param(t1, 2)
jhParam = Param(jh, 2)
tdParam = Param(t_density, 2)
tiParam = Param(t_inter, 2)
HoppingParams = [t1Param, tdParam, tiParam]

su2spin = SpinMats(1 // 2)
su4spin = SpinMats(3 // 2)

##Adding inner-hexagon structure  
for j = 1:2
    for i = -1:4
        AddBasisSite!(UC, i .* l1 + j .* l2)
    end
end

##Istrotropic bonds
AddIsotropicBonds!(t1Param, UC, 1.0, su4spin[4], "t1", checkOffsetRange=1)
AddIsotropicBonds!(tiParam, UC, 0.0, 2 * kron(su2spin[1], su2spin[4]), "interlayer")
AddIsotropicBonds!(tdParam, UC, 0.0, 2 * kron(su2spin[3], su2spin[4]), "imbalance")

##Functions that will be useful for adding anisotropic bonds
weiss0(v) = [sin(pi * (1 - norm(v) / 2)) * v[1] / norm(v), sin(pi * (1 - norm(v) / 2)) * v[2] / norm(v), cos(pi * (1 - norm(v) / 2))]
weiss1(v) = [sin(pi * (1 - norm(v) / 2)) * v[1] / norm(v), sin(pi * (1 - norm(v) / 2)) * v[2] / norm(v), -cos(pi * (1 - norm(v) / 2))]
sigmav(i, j) = 2 .* [su2spin[1][i, j], su2spin[2][i, j], su2spin[3][i, j]]

s11 = sigmav(1, 1)
s12 = sigmav(1, 2)
s21 = sigmav(2, 1)
s22 = sigmav(2, 2)
intermat(s1, s2) = [dot(s1, s11) dot(s1, s12) 0 0; dot(s1, s21) dot(s1, s22) 0 0; 0 0 dot(s2, s11) dot(s2, s12); 0 0 dot(s2, s21) dot(s2, s22)]



CreateUnitCell!(UC, HoppingParams)


##Creating BZ and Hamiltonian Model
bz = BZ(kSize)
FillBZ!(bz, UC)
kSize = 6 * 12 + 3
bz = BZ(kSize, 2)
FillBZ!(bz, UC)
path = CombinedBZPath(bz, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]]; nearest=true)
# Adding MFT Parameters
HoppingParams = [t1Param]
n_up = real.(kron([1.0 0.0; 0.0 0.0], su2spin[4]))
n_down = real.(kron([0.0 0.0; 0.0 1.0], su2spin[4]))
Hubbard = DensityToPartonCoupling(n_up, n_down)
UParam = Param(1.0, 4)
AddIsotropicBonds!(UParam, UC, 0.0, Hubbard, "Hubbard Interaction") # Do I need to add this to all sites?
##Adding anisotropic bonds and normalizing if needed
for (ind, bas) in enumerate(UC.basis)
    if 1 < norm(bas) < 2
        mat = intermat(normalize(weiss0(bas) + weiss0(-bas)), normalize(weiss1(bas) + weiss1(-bas)))
    else
        closest = [bas, bas - a1, bas - a2]
        clv = closest[findmin(x -> norm(x), closest)[2]]
        mat = intermat(replace!(weiss0(clv), NaN => 0.0), replace!(weiss1(clv), NaN => 0.0))
    end
    AddAnisotropicBond!(jhParam, UC, ind, ind, [0, 0], mat, 0.0, "Hunds")
end

##Creating BZ and Hamiltonian Model
for U_var in U_array
    Density = []
    UParam.value = [U_var]

    for (ind, bas) in enumerate(UC.basis)
        push!(Density, Param(1.0, 2))
        AddAnisotropicBond!(Density[ind], UC, ind, ind, [0, 0], kron(su2spin[3], su2spin[4]), 0.0, "Dens-" * string(ind))
    end

    ChiParams = vcat(Density)
    ChiParams = Vector{Param{2,Float64}}(ChiParams)


    path = CombinedBZPath(bz, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]]; nearest=true)
    H = Hamiltonian(UC, bz)
    DiagonalizeHamiltonian!(H)
    filling = 0.5
    Mdl = Model(UC, bz, H; filling=filling, T=T) # Does T matter, don't I want 0 T, or is that technically impossible? 
    SolveModel!(Mdl; get_gap=true)
    mft = TightBindingMFT(Mdl, ChiParams, [UParam], IntraQuarticToHopping)
    fileName = loc * "/Monolayer=$(round(filling, digits=3))_U=$(round(U_var, digits=2))_t1=$(round(t1, digits=2)).jld2"
    SolveMFT!(mft, fileName; max_iter=100, tol=1e-4)
end












# ##Plotting the unit cell
# plot_UC = Plot_UnitCell!(UC);

# # JH_range = collect(-1.0:-1.0)
# #tinter_range = collect(-0.0:-0.2:-1.0)
# #t_density_range = collect(0.0:0.2:1.0)
# for tdensity_val in t_density_range
#     push!(t_density_param.value, tdensity_val)
#     ModifyUnitCell!(UC, [t_density_param])

#     for tinter_val in tinter_range

#         push!(tinter_param.value, tinter_val)
#         ModifyUnitCell!(UC, [tinter_param])

#         global H = Hamiltonian(UC, bz)
#         DiagonalizeHamiltonian!(H)
#         global Mdl = Model(UC, bz, H; filling=1 / 24)
#         SolveModel!(Mdl; get_gap=true)
#         global bands = Plot_Band_Structure!(Mdl, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]], collect(1:12), labels=["G", "K1", "M2"], plot_legend=false)
#         display(bands)
#         # savefig(bands,"Jh  v t-inter bands/Jh = $jh_val, t-inter = $tinter_val.png")

#         push!(combined_Cnums, ChernNumber(H, collect(1:2)))
#         push!(gaps, Mdl.gap)
#         push!(mus, Mdl.mu)
#         println(ChernNumber(H, collect(1:2)))
#     end
# end