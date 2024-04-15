using Plots, TightBindingToolkit, LinearAlgebra, ColorSchemes
SkXSize = get!(params, "SkXSize", 3)
SkXSize = 2
SkX = get!(params, "SkX", "Bloch")
a1 = SkXSize / 2 * [-3.0, sqrt(3)]
a2 = SkXSize / 2 * [3.0, sqrt(3)]
UC = UnitCell([a1, a2], 2)

##Triangular Lattice
l1 = [1.0, 0]
l2 = [-0.5, sqrt(3) / 2]

##Parameters
t = 1.0
pspin = SpinMats(1 // 2)
JH = -1.0


##Adding inner-hexagon structure
for j = 0:(SkXSize-1)
    for i = 0:(SkXSize*3-1)
        AddBasisSite!(UC, i .* l1 + j .* l2)
    end
end

##Istrotropic bonds
AddIsotropicBonds!(UC, 1.0, -t * pspin[4], "iso", checkOffsetRange=1)

##Functions that will be useful for adding anisotropic bonds

weiss_neel(v) = [sin(pi * (norm(v) / (SkXSize))) * v[1] / norm(v), sin(pi * (norm(v) / (SkXSize))) * v[2] / norm(v), cos(pi * (norm(v) / (SkXSize)))]
weiss_bloch(v) = [sin(pi * (norm(v) / (SkXSize))) * v[1] / norm(v), sin(pi * (norm(v) / (SkXSize))) * v[2] / norm(v), cos(pi * (norm(v) / (SkXSize)))]
weiss = Dict("Neel" => weiss_neel, "Bloch" => weiss_bloch)
sigmav(i, j) = 2 .* [pspin[1][i, j], pspin[2][i, j], pspin[3][i, j]]
s11 = sigmav(1, 1)
s12 = sigmav(1, 2)
s21 = sigmav(2, 1)
s22 = sigmav(2, 2)

intermat(s) = [dot(s, s11) dot(s, s12); dot(s, s21) dot(s, s22)]

##Adding anisotropic bonds and normalizing if needed
for (ind, bas) in enumerate(UC.basis)
    closest = [bas, bas - a1, bas - a2, bas - a1 - a2, bas + a1, bas + a2, bas + a1 + a2, bas + a1 - a2, bas - a1 + a2]
    minimal = findmin(x -> norm(x), closest)[2]
    if (SkXSize - 1) < norm(closest[minimal]) < SkXSize
        mat = intermat(normalize(weiss[SkX](closest[minimal]) + weiss[SkX](-closest[minimal])))
    else

        # println(findmin(x -> norm(x), closest))
        spn = weiss[SkX](closest[minimal])
        replace!(spn, NaN => 0.0)
        mat = intermat(normalize(spn))
    end
    AddAnisotropicBond!(UC, ind, ind, [0, 0], -JH * mat, 0.0, "interaction")
end

p = Plot_Fields!(UC; use_lookup=true, site_size=1.0,
    field_thickness=1.0, range=2, field_opacity=0.9, scale=0.5,
    cmp=:viridis)


# p.legend = false##Plotting the unit cell
# plot_UC = Plot_UnitCell!(UC);
# display(plot_UC)

#Creating BZ and Hamiltonian Model
kSize = 6 * 10 + 3
bz = BZ(kSize)
FillBZ!(bz, UC)
path = CombinedBZPath(bz, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]]; nearest=true)
H = Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)
Mdl = Model(UC, bz, H; filling=0.5)
SolveModel!(Mdl; get_gap=true)

#Plotting the band structure
bands = Plot_Band_Structure!(Mdl, [bz.HighSymPoints["G"], bz.HighSymPoints["M1"], bz.HighSymPoints["M2"]], labels=["G", "M2", "M3"], plot_legend=false);
plot!(bands, legend=false);
display(bands)

#Calculating Chern Numbers for bands
# for i in 1:2*length(UC.basis)
#     c = ChernNumber(H, [i])
#     println(round(c))
# end
# println("Chern Number for first 12 bands: ", ChernNumber(H, collect(1:12)))
