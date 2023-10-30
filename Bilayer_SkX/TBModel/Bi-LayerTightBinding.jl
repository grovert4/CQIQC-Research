using Plots, LinearAlgebra, ColorSchemes
using MeanFieldToolkit, TightBindingToolkit, FixedPointToolkit
loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Bilayer_Data"

##Triangular Lattice 

const a1 = [-3.0, sqrt(3)]
const a2 = [3.0, sqrt(3)]

const l1 = [1.0, 0]
const l2 = [-0.5, sqrt(3) / 2]
UC = UnitCell([a1, a2], 4)
##Parameters
const n = 10
const kSize = 6 * n + 3
const t = 1.0
const t_inter = 0.0
const jh = -1.0
const U = 1.0
const t_density = 0
U_array = collect(LinRange(0.0, 7.0, 12))
SpinVec = SpinMats(1 // 2)
##### Thermodynamic parameters
filling = 0.5
const T = 0.001

tinter_param = Param(t_inter, 2)
t1 = -t
t1Param = Param(t1, 2)
jhParam = Param(jh, 2)
tdParam = Param(t_density, 2)
tiParam = Param(t_inter, 2)
HoppingParams = [t1Param, tdParam, tiParam, jhParam]

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
path = CombinedBZPath(bz, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]]; nearest=true)

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
path = CombinedBZPath(bz, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]]; nearest=true)
H = Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)
Mdl = Model(UC, bz, H; filling=filling, T=T) # Does T matter, don't I want 0 T, or is that technically impossible? 
SolveModel!(Mdl; get_gap=true)
for i in 1:2*length(UC.basis)
    c = ChernNumber(H, [i])
    println(round(c))
end
