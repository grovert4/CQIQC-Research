using Plots, LinearAlgebra, SpinMC_more_more, HDF5, JSON, LaTeXStrings
include("functions.jl")
include("plottingFunctions.jl")

inputFile = JSON.parsefile("inputParametersBiLayer.json")

Jperps = collect(range(inputFile["Jperp_min"], inputFile["Jperp_max"] ,length = inputFile["Jperp_length"]))
J2s = collect(range(inputFile["J2_min"], inputFile["J2_max"] ,length = inputFile["J2_length"]))

a1 = (1.0 , 0.0, 0.0)
a2 = (-1/2 , sqrt(3)/2, 0.0)
a3 = (0.0, 0.0, 2.0)
UCtemp = UnitCell(a1 , a2, a3)
b1 = addBasisSite!(UCtemp, (0.0, 0.0, 0.0)) 
b2 = addBasisSite!(UCtemp, (0.0, 0.0, 1.0)) 

J1 = 1.0
D = 0.25
A_ion = 0.2
J2 = -0.25

tempLattice = Lattice(UCtemp, (inputFile["System_Size"],inputFile["System_Size"], 1))
vertex = getVertex(tempLattice)

phases = zeros(length(Jperps), length(J2s)) ##Keeps track of whether or not there is antiferro skyrmion or not
skxnums0 = zeros(length(Jperps), length(J2s)) ##Keep track of exact number of skyrmions in layer 0
skxnums1 = zeros(length(Jperps), length(J2s)) ##Keep track of exact number of skyrmions in layer 1
energies = zeros(length(Jperps), length(J2s)) ##Keeps track of energy per site of system
sztrack = zeros(length(Jperps), length(J2s)) ##Keeps track of average Sz of system

Id = [1 0 0; 0 1 0; 0 0 1]
Sz = [0 0 0; 0 0 0; 0 0 1]

#0 == 0 Skymrions both layer
#1 == Different SkxNumms in layers 
#2 == Same Skyrmion Nums in layers
#3 == Opposite Skyrmions in layers 
#4 == Something else

for (j2idx, j2) in enumerate(J2s)
   global J2 = round(j2,sigdigits=5)
   for (jpidx,jperp) in enumerate(Jperps)
      global J_ll = round(jperp,sigdigits=5)
            
      updateSpins!("C:/Users/tanma/bigboybilayer/firstphasediagram/01.05.2024-Bilayer-decreasingfield_Jperp=$J_ll,J2=$J2.h5", tempLattice)

      skxNum0 = round(getSkyrmionNumber(0, tempLattice, vertex),digits=2)
      skxNum1 = round(getSkyrmionNumber(1, tempLattice, vertex),digits=2)
      if skxNum0 == 0.0 && skxNum1 == 0.0
         phases[jpidx, j2idx] = 0
      elseif skxNum0 > 0 && skxNum1 < 0
         phases[jpidx, j2idx] = 3
      elseif skxNum0 == skxNum1
         phases[jpidx, j2idx] = 2
      elseif skxNum0 != skxNum1 
         phases[jpidx, j2idx] = 1 
      else
         phases[jpidx, j2idx] = 4
      end

      skxnums0[jpidx,j2idx] = skxNum0
      skxnums1[jpidx,j2idx] = skxNum1
      sztrack[jpidx, j2idx] = avgSz(tempLattice)
      energies[jpidx, j2idx] = readEnergy("C:/Users/tanma/bigboybilayer/firstphasediagram/01.05.2024-Bilayer-decreasingfield_Jperp=$J_ll,J2=$J2.h5")[1]
   end
end

#Generate random lattices to see what spin configuration looks like for them
#randomLattices()

plotMiniPhase(J2s, Jperps, phases; flipx=true, flipy=true, xlab=L"J_2", ylab=L"J_{||}", title="Phase Diagram for Bilayer")


