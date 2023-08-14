using Plots, TightBindingToolkit, LinearAlgebra

##Skyrmion Lattice
a1 = [3.0 , sqrt(3)]
a2 = [3.0 , -sqrt(3)]
UC = UnitCell([a1 , a2] , 2) 

##Triangular Lattice 
l1 = [1.0, 0]
l2 = [-0.5, sqrt(3)/2]

##Parameters
t = 0.5
jh = 2.5
pspin = SpinMats(1//2)


##Adding inner-hexagon structure  
for i = -2:2
    for j = -2:2
        if norm(i .* l1 + j .* l2) <= 1.0
            AddBasisSite!(UC, i .* l1 + j .* l2)
        end
    end
end

##Adding boundary points of outer-hexagon 
AddBasisSite!(UC, 2 * l2)
AddBasisSite!(UC, 2 * l2 + l1)
AddBasisSite!(UC, 2 * l2 + l1)
AddBasisSite!(UC, 2 * l2 + 2 * l1)
AddBasisSite!(UC, -l2 + 2 * l1)

##Istrotropic bonds
AddIsotropicBonds!(UC, 1.0, -t * pspin[4], "iso", checkOffsetRange = 1)

##Functions that will be useful for adding anisotropic bonds
tau(v) = [sin(pi * (1 - norm(v)/2)) * cos(atan(v[2]/v[1])), cos(pi * (1 - norm(v)/2 )), sin(pi * (1 - norm(v)/2 )) * sin(atan(v[2]/v[1]))]
sigmav(i,j) = [pspin[1][i,j], pspin[2][i,j], pspin[3][i,j]]

s11 = sigmav(1,1)
s12 = sigmav(1,2)
s21 = sigmav(2,1)
s22 = sigmav(2,2)
intermat(s) = [dot(s, s11) dot(s, s12);dot(s, s21) dot(s, s22)]

##Adding anisotropic bonds and normalizing if needed
for (ind, bas) in enumerate(UC.basis)
    if 1 < norm(bas) < 2
        mat = intermat(normalize(tau(bas) + tau(-bas))) 
    else 
        mat = intermat(tau(bas))
    end
    replace!(mat, NaN + NaN * im=> 0.0 + 0.0 * im)
    AddAnisotropicBond!(UC, ind, ind, [0,0], -jh * mat, 0.0, "interaction")
end

##Plotting the unit cell
plot_UC = Plot_UnitCell!(UC);
# display(plot_UC)

##Creating BZ and Hamiltonian Model
const kSize = 6 * 15 + 3  
bz = BZ(kSize)
FillBZ!(bz, UC)
path = CombinedBZPath(bz, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M1"]] ; nearest=true)
H = Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)
Mdl = Model(UC, bz, H)
SolveModel!(Mdl)

##Plotting the band structure
bands = Plot_Band_Structure!(Mdl, path, labels = ["G", "K", "M"]);
plot!(bands, legend = false);
display(bands)

##Calculating Chern Numbers for bands
totC = 0
for i in 1:length(UC.basis)
    c = ChernNumber(H, [i])
    global totC += c
    println("Chern Number is: ", c)
end
println("Sum = ", totC)
