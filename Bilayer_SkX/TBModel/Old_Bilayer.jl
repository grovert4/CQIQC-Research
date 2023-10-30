using Plots, TightBindingToolkit, LinearAlgebra, ColorSchemes


#Skyrmion Lattice
# a1 = [3.0 , sqrt(3)]
# a2 = [3.0 , -sqrt(3)


a1 = [-3.0, sqrt(3)]
a2 = [3.0, sqrt(3)]
UC = UnitCell([a1, a2], 4)

##Triangular Lattice 
l1 = [1.0, 0]
l2 = [-0.5, sqrt(3) / 2]

##Parameters
t = -1.0
t_param = Param(t, 2)

tinter = -0.3
tinter_param = Param(tinter, 2)

t_density_param = Param(t_density, 2)

jh = -1.0
jh_param = Param(jh, 2)

su2spin = SpinMats(1 // 2)
su4spin = SpinMats(3 // 2)

##Adding inner-hexagon structure  
for j = 1:2
    for i = -1:4
        AddBasisSite!(UC, i .* l1 + j .* l2)
    end
end

##Istrotropic bonds
AddIsotropicBonds!(t_param, UC, 1.0, su4spin[4], "iso", checkOffsetRange=1)
AddIsotropicBonds!(tinter_param, UC, 0.0, 2 * kron(su2spin[1], su2spin[4]), "isoHop")
AddIsotropicBonds!(t_density_param, UC, 0.0, 2 * kron(su2spin[3], su2spin[4]), "imbalance")

##Functions that will be useful for adding anisotropic bonds
tau0(v) = [sin(pi * (1 - norm(v) / 2)) * v[1] / norm(v), sin(pi * (1 - norm(v) / 2)) * v[2] / norm(v), cos(pi * (1 - norm(v) / 2))]
tau1(v) = [sin(pi * (1 - norm(v) / 2)) * v[1] / norm(v), sin(pi * (1 - norm(v) / 2)) * v[2] / norm(v), -cos(pi * (1 - norm(v) / 2))]
sigmav(i, j) = 2 .* [su2spin[1][i, j], su2spin[2][i, j], su2spin[3][i, j]]

s11 = sigmav(1, 1)
s12 = sigmav(1, 2)
s21 = sigmav(2, 1)
s22 = sigmav(2, 2)
intermat(s1, s2) = [dot(s1, s11) dot(s1, s12) 0 0; dot(s1, s21) dot(s1, s22) 0 0; 0 0 dot(s2, s11) dot(s2, s12); 0 0 dot(s2, s21) dot(s2, s22)]

##Adding anisotropic bonds and normalizing if needed
for (ind, bas) in enumerate(UC.basis)
    if 1 < norm(bas) < 2
        mat = intermat(normalize(tau0(bas) + tau0(-bas)), normalize(tau1(bas) + tau1(-bas)))
    else
        closest = [bas, bas - a1, bas - a2]
        clv = closest[findmin(x -> norm(x), closest)[2]]
        mat = intermat(replace!(tau0(clv), NaN => 0.0), replace!(tau1(clv), NaN => 0.0))
    end
    AddAnisotropicBond!(jh_param, UC, ind, ind, [0, 0], mat, 0.0, "interaction")
end

CreateUnitCell!(UC, [t_param, tinter_param, jh_param, t_density_param])

##Plotting the unit cell
plot_UC = Plot_UnitCell!(UC);

# JH_range = collect(-1.0:-1.0)
tinter_range = collect(-0.0:-0.2:-1.0)
t_density_range = collect(0.0:0.2:1.0)
gaps = []
mus = []
combined_Cnums = []
individual_Cnums = []

##Creating BZ and Hamiltonian Model
kSize = 6 * 12 + 3
bz = BZ(kSize, 2)
FillBZ!(bz, UC)
path = CombinedBZPath(bz, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]]; nearest=true)


for tdensity_val in t_density_range
    push!(t_density_param.value, tdensity_val)
    ModifyUnitCell!(UC, [t_density_param])

    for tinter_val in tinter_range

        push!(tinter_param.value, tinter_val)
        ModifyUnitCell!(UC, [tinter_param])

        global H = Hamiltonian(UC, bz)
        DiagonalizeHamiltonian!(H)
        global Mdl = Model(UC, bz, H; filling=1 / 24)
        SolveModel!(Mdl; get_gap=true)
        global bands = Plot_Band_Structure!(Mdl, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]], collect(1:12), labels=["G", "K1", "M2"], plot_legend=false)
        display(bands)
        # savefig(bands,"Jh  v t-inter bands/Jh = $jh_val, t-inter = $tinter_val.png")

        push!(combined_Cnums, ChernNumber(H, collect(1:2)))
        push!(gaps, Mdl.gap)
        push!(mus, Mdl.mu)
        println(ChernNumber(H, collect(1:2)))
    end
end