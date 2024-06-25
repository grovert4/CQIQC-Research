using Plots, TightBindingToolkit, LinearAlgebra, ColorSchemes, MeanFieldToolkit
using LaTeXStrings
# using PyCall

# pyplot() # or pgfplotsx()

# rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
# mpl = pyimport("matplotlib")

# mpl.style.use("./TBModel/Plotting/lake.mplstyle")
# rcParams["text.usetex"] = true

# rcParams["font.family"] = "serif"


function KuboChern(Ham::Hamiltonian, bz::BZ, mu::Float64)

    Vx = conj.(permutedims.(Ham.states)) .* Ham.velocity[1] .* Ham.states
    Vy = conj.(permutedims.(Ham.states)) .* Ham.velocity[2] .* Ham.states

    chern = 0.0 + im * 0.0
    for k in eachindex(Ham.bands)
        Es = Ham.bands[k]
        vx = Vx[k]
        vy = Vy[k]

        ind = searchsortedfirst(Es, mu)
        if ind == 1 || ind == length(Es)
            continue
        else
            for i in 1:ind-1
                for j in ind:length(Es)
                    chern += (vx[i, j] * vy[j, i] - vx[j, i] * vy[i, j]) / ((Es[j] - Es[i])^2)
                end
            end
        end

    end

    b1 = [bz.basis[1]; 0.0]
    b2 = [bz.basis[2]; 0.0]
    bzUnitArea = cross(b1, b2)[3] / (4 * pi^2)

    return imag(chern) * bzUnitArea * 2 * pi / length(Ham.bands)

end
##Triangular Lattice
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

Mdl = Model(UC, bz, H; filling=filling)
SolveModel!(Mdl; get_gap=true)

#Plotting the band structure
bands = Plot_Band_Structure!(Mdl, [bz.HighSymPoints["G"], bz.HighSymPoints["M1"], bz.HighSymPoints["K1"]], labels=["G", "M", "K"], plot_legend=false);
plot!(bands, legend=false);
display(bands)
savefig(bands, "TBModel/Plotting/Plots/Monolayer_Band_Structure_$(SkXSize)_$(jh).pdf")
filling_arr = LinRange(0.0001, 0.5, 60)
#mu_arr = LinRange(-7.5, -1, 100)
mu_arr = LinRange(-7.5, 0, 150)

c_fill = zeros(length(mu_arr))
ssf_plot = plot(framestyle=:box, aspect_ratio=:equal, xlabel=L"k_x", ylabel=L"k_y", grid=false, legend=:outerright)
heatmap!(berry, color=:viridis)  # Adjust the size as needed
display(ssf_plot)
# #for (i, filling) in enumerate(filling_arr)
# for (i, mu) in enumerate(mu_arr)
# DiagonalizeHamiltonian!(H)
# GetVelocity!(H, bz)
# #Mdl = Model(UC, bz, H; filling=filling)
# #SolveModel!(Mdl; get_gap=true)
# c_fill[i] = KuboChern(H, bz, mu)
# end
# p = plot(mu_arr, -1*c_fill,  ylims=(-5, 15), xlabel=L"\mu", ylabel=L"\sigma_{xy}", label=L"\sigma_{xy}", legend=false)
# display(p)
# savefig(p, "TBModel/Plotting/Plots/Monolayer_Hall_$(SkXSize)_$(jh).pdf")

#Calculating Chern Numbers for bands
# for i in 1:2*length(UC.basis)
#     c = ChernNumber(H, [i])
#     println(round(c))
# end
# println("Chern Number for first 12 bands: ", ChernNumber(H, collect(1:12)))
