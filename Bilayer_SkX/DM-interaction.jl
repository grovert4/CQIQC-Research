using SpinMC, LinearAlgebra, Plots
include("functions.jl")

#Unit Cell Construction
a1 = (1.0 , 0.0, 0.0)  #-
a2 = (1/2 , sqrt(3)/2, 0.0)  #/
a3 = (0.0, 0.0, 2.0)  #|
UC = UnitCell(a1 , a2, a3) 

b1 = addBasisSite!(UC, (0.0, 0.0, 0.0)) ##layer A (z = 0)
b2 = addBasisSite!(UC, (0.0, 0.0, 1.0)) ##layer A (z = 1)


#Parameters
D = 0.5
J_ll = 0.1

D1v = D .* [0, -1, 0]
D2v = D .* [sqrt(3)/2, 1/2, 0] 
D3v = D .* [-sqrt(3)/2, 1/2, 0] 
D1ex = [0 -D1v[3] D1v[2]; D1v[3] 0 -D1v[1]; -D1v[2] D1v[1] 0]
D2ex = [0 -D2v[3] D2v[2]; D2v[3] 0 -D2v[1]; -D2v[2] D2v[1] 0]
D3ex = [0 -D3v[3] D3v[2]; D3v[3] 0 -D3v[1]; -D3v[2] D3v[1] 0]


##Same layer interactions
for i in 1:length(UC.basis)
    ##DM interaction
    addInteraction!(UC, i, i, (-1)^(i-1) * D1ex, (1,0,0))
    addInteraction!(UC, i, i, (-1)^(i-1) * D2ex, (-1,1,0))
    addInteraction!(UC, i, i, (-1)^(i-1) * D3ex, (0,-1,0))
end


L = (24, 24, 1)
global lattice = Lattice(UC, L)

thermalizationSweeps = 25000
measurementSweeps = 75000
m = MonteCarlo(lattice, 1000.0, thermalizationSweeps, measurementSweeps, reportInterval = 50000);
run!(m)
