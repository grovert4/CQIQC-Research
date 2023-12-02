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
H = 0.5
A_ion = 0.2
J_ll = 0.1
Sz = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]

##Same layer interactions
for i in 1:length(UC.basis)
    #Magnetic field
    setField!(UC, i, [-H,0,0])

    # ##onsite anisotropy
    setInteractionOnsite!(UC, i, A_ion * S_z)
end


L = (24, 24, 1)
global lattice = Lattice(UC, L)

thermalizationSweeps = 25000
measurementSweeps = 75000
m = MonteCarlo(lattice, 1000.0, thermalizationSweeps, measurementSweeps, reportInterval = 50000);
run!(m)
