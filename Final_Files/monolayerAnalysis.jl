using Plots, LinearAlgebra, SpinMC_more_more, HDF5, JSON, LaTeXStrings
include("functions.jl")
include("plottingFunctions.jl")

inputFile = JSON.parsefile("inputParametersMonoLayer.json")

##Importing parameters from input files
Hs = collect(range(inputFile["H_min"], inputFile["H_max"] ,length = inputFile["H_length"])) 
J2s = collect(range(inputFile["J2_min"], inputFile["J2_max"] ,length = inputFile["J2_length"]))

##Setting up temporary UnitCell and Lattice
a1 = (1.0 , 0.0)
a2 = (-1/2 , sqrt(3)/2)
UCtemp = UnitCell(a1 , a2)
b1 = addBasisSite!(UCtemp, (0.0, 0.0)) 
tempLattice = Lattice(UCtemp, (inputFile["System_Size"], inputFile["System_Size"]))
vertex = getVertex(tempLattice)

phases = zeros(length(Hs), length(J2s)) ##Keeps track of whether or not there is a non-zero skyrmion number
exactskxnums = zeros(length(Hs), length(J2s)) ##Keep track of exact number of skyrmions
energies = zeros(length(Hs), length(J2s)) ## Keeps track of energy per site of system
sztrack = zeros(length(Hs), length(J2s)) ##Keeps track of average Sz of system

J1 = 1
D = 0.25
A_ion = 0.2

# 0 == Fully Polarized
# 1 == Spin Spiral
# 2 == Skyrmion
# 3 == Error/Recompute for params

for (j2idx, j2) in enumerate(J2s)
   global J2 = round(j2,sigdigits=5)
   
   for (hidx,h) in enumerate(Hs)
   
      global H = round(h,sigdigits=5)
      updateSpins!("C:/Users/tanma/Monolayer_Runs_Take2/13.10.2024-Monolayer_H=$H,J2=$J2.h5", tempLattice) ##Updates temprorary lattice
      skxNum = round(getSkyrmionNumber(0, tempLattice, vertex),digits=2)
      if skxNum == 0.0
         phases[hidx,j2idx] = 0
      elseif skxNum != 0.0 
         phases[hidx,j2idx] = 1 
      end
      exactskxnums[hidx,j2idx] = skxNum
      sztrack[hidx, j2idx] = avgSz(tempLattice)
      energies[hidx, j2idx] = readEnergy("C:/Users/tanma/Monolayer_Runs_Take2/13.10.2024-Monolayer_H=$H,J2=$J2.h5")[1]
   end
end

#Generate random lattices to see what spin configuration looks like for them
#randomLattices()

plotMiniPhase(J2s, Hs, phases; flipx=true, flipy=false, xlab=L"J_2", ylab=L"H", title="Phase Diagram for Monolayer")
