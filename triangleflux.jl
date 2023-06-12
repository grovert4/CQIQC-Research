using TightBindingToolkit, Plots, LinearAlgebra

a1 = [ 1/2 , sqrt(3)/2]
a2 = [2.0, 0 ]
UC = UnitCell( [a1 , a2] , 1)
b1 = [ 0.0 , 0.0 ]
b2 = [ 1.0 , 0.0 ]


addBasisSite!(UC , b1)
addBasisSite!(UC , b2)


const t1  =  -1.0


addAnisotropicBond!(UC, 1, 1, [1,0], t1, 1.0, "a1")
addAnisotropicBond!(UC, 1, 2, [0,0], t1, 1.0, "a1")
addAnisotropicBond!(UC, 2, 1, [0,1], t1, 1.0, "a1")
addAnisotropicBond!(UC, 2, 2, [1,0], t1*exp(im * pi), 1.0, "b1")
addAnisotropicBond!(UC, 2, 1, [1,0], t1*exp(im * pi/2), 1.0, "c1")
addAnisotropicBond!(UC, 2, 1, [-1,1], t1*exp(im * pi/2), 1.0, "d1")



kSize   =   6 * 20 + 3   ##### a Monkhorst grid of size N =    6n+3 covers all the High-symemtry points.
bz =  BZ(kSize)
fillBZ!(bz, UC)

H = Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)

"""
Filling the model with fermions at half-filling, with the default temperature T=0.0.
"""
triangularFlux = Model(UC, bz, H ; filling=0.5)
SolveModel!(triangularFlux)
chern   =   ChernNumber(H, [1])
println(chern)
# push!(cherns, chern)


# addIsotropicBonds!(UC , NNdistance, t1 , "t1")



fluxs = range(-1*pi,1*pi,50)
cherns = []

for f in fluxs
    RemoveBonds!(UC, "b1")
    RemoveBonds!(UC, "c1")
    RemoveBonds!(UC, "d1")
    addAnisotropicBond!(UC, 2, 2, [1,0], t1*exp(im * 2 * f), 1.0, "b1")
    addAnisotropicBond!(UC, 2, 1, [1,0], t1*exp(im * f), 1.0, "c1")
    addAnisotropicBond!(UC, 2, 1, [-1,1], t1*exp(im * f), 1.0, "d1")
    kSize   =   6 * 15 + 3   ##### a Monkhorst grid of size N =    6n+3 covers all the High-symemtry points.
    bz =  BZ(kSize)
    fillBZ!(bz, UC)

    H = Hamiltonian(UC, bz)
    DiagonalizeHamiltonian!(H)

    """
    Filling the model with fermions at half-filling, with the default temperature T=0.0.
    """
    triangularFlux = Model(UC, bz, H ; filling=0.5)
    SolveModel!(triangularFlux)
    chern   =   ChernNumber(H, [1])
    push!(cherns, chern)
end

plot(fluxs, cherns)

