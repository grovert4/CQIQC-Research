using Plots, TightBindingToolkit, LinearAlgebra, ColorSchemes
##Triangular Lattice
SkXSize = get!(params, "SkXSize", 2)
SkX = get!(params, "SkX", "Bloch")
SkX = "Neel"
a1 = SkXSize / 2 * [-3.0, sqrt(3)]
a2 = SkXSize / 2 * [3.0, sqrt(3)]
l1 = [1.0, 0]
l2 = [-0.5, sqrt(3) / 2]
UC = UnitCell([a1, a2], 2, 2)
##Parameters
n = get!(params, "n", 20)
kSize = 6 * n + 3
t = get!(params, "t", 1.0)
jh = get!(params, "jh", -1.0)
U = get!(params, "U", 0.0)
##### Thermodynamic parameters
filling = get!(params, "filling", 0.5)
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
##Functions that will be useful for adding anisotropic bonds
weiss_neel(v) = [sin(pi * (norm(v) / (SkXSize))) * v[1] / norm(v), sin(pi * (norm(v) / (SkXSize))) * v[2] / norm(v), cos(pi * (norm(v) / (SkXSize)))]
weiss_bloch(v) = [sin(pi * (norm(v) / (SkXSize))) * v[2] / norm(v), sin(pi * (norm(v) / (SkXSize))) * -v[1] / norm(v), cos(pi * (norm(v) / (SkXSize)))]
weiss = Dict("Neel" => weiss_neel, "Bloch" => weiss_bloch)
sigmav(i, j) = 2 .* [su2spin[1][i, j], su2spin[2][i, j], su2spin[3][i, j]]
s11 = sigmav(1, 1)
s12 = sigmav(1, 2)
s21 = sigmav(2, 1)
s22 = sigmav(2, 2)

intermat(s) = [dot(s, s11) dot(s, s12); dot(s, s21) dot(s, s22)]

##Adding anisotropic bonds and normalizing if needed

bz = BZ(kSize)
FillBZ!(bz, UC)
##Adding anisotropic bonds and normalizing if needed
# Adding MFT Parameters
n_up = [1.0 0.0; 0.0 0.0]
n_down = [0.0 0.0; 0.0 1.0]
Hubbard = DensityToPartonCoupling(n_up, n_down)
UParam = Param(0.0, 4)
Sz = []
tParam = Param(1.0, 2)
AddIsotropicBonds!(UParam, UC, 1.0, Hubbard, "Hubbard Interaction", checkOffsetRange=1) # Do I need to add this to all sites?
for (ind, bas) in enumerate(UC.basis)
    closest = [bas, bas - a1, bas - a2, bas - a1 - a2, bas + a1, bas + a2, bas + a1 + a2, bas + a1 - a2, bas - a1 + a2]
    minimal = findmin(x -> norm(x), closest)[2]
    if (SkXSize - 1) < norm(closest[minimal]) < SkXSize
        mat = intermat(normalize(weiss[SkX](closest[minimal]) + weiss[SkX](-closest[minimal])))
    else
        spn = weiss[SkX](closest[minimal])
        replace!(spn, NaN => 0.0)
        mat = intermat(normalize(spn))
    end
    AddAnisotropicBond!(jhParam, UC, ind, ind, [0, 0], mat, 0.0, "Hunds")
end
CreateUnitCell!(UC, HoppingParams)
AddIsotropicBonds!(tParam, UC, 1.0, su2spin[4], "s_H") # Am I not double counting the hopping ??
for (ind, bas) in enumerate(UC.basis)
    push!(Sz, Param(1.0, 2))
    AddAnisotropicBond!(Sz[ind], UC, ind, ind, [0, 0], su2spin[3], 0.0, "Sz-" * string(ind))
    #AddAnisotropicBond!(Nd[ind], UC, ind, ind, [0, 0], n_down, 0.0, "Ndown-" * string(ind))
end

ChiParams = vcat(tParam, Sz)
ChiParams = Vector{Param{2,Float64}}(ChiParams)
##Creating BZ and Hamiltonian Model
bz = BZ(kSize)
FillBZ!(bz, UC)
H = Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)

p = Plot_Fields!(UC; use_lookup=true, site_size=1.0,
    field_thickness=1.0, range=2, field_opacity=0.9, scale=0.5,
    cmp=:viridis)

display(p)
# p.legend = false##Plotting the unit cell
# plot_UC = Plot_UnitCell!(UC);
# display(plot_UC)

#Creating BZ and Hamiltonian Model
kSize = 6 * 5 + 3
bz = BZ(kSize)
FillBZ!(bz, UC)
path = CombinedBZPath(bz, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]]; nearest=true)
H = Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)
Mdl = Model(UC, bz, H; filling=0.5)
SolveModel!(Mdl; get_gap=true)

#Plotting the band structure
bands = Plot_Band_Structure!(Mdl, [bz.HighSymPoints["G"], bz.HighSymPoints["M2"], bz.HighSymPoints["M3"]], labels=["G", "M2", "M3"], plot_legend=false);
plot!(bands, legend=false);
display(bands)

#Calculating Chern Numbers for bands
# for i in 1:2*length(UC.basis)
#     c = ChernNumber(H, [i])
#     println(round(c))
# end
# println("Chern Number for first 12 bands: ", ChernNumber(H, collect(1:12)))
