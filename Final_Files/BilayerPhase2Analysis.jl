using Plots, LinearAlgebra, SpinMC_more_more, HDF5, JSON, LaTeXStrings, PyCall
include("functions.jl")
include("plottingFunctions.jl")

inputFile = JSON.parsefile("inputParametersBiLayerphase2.json")

Jperps = collect(range(inputFile["Jperp_min"], inputFile["Jperp_max"] ,length = inputFile["Jperp_length"]))
Hs = collect(range(inputFile["H_min"], inputFile["H_max"] ,length = inputFile["H_length"]))

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

phases = zeros(length(Hs), length(Jperps))
skxnums0 = zeros(length(Hs), length(Jperps))
skxnums1 = zeros(length(Hs), length(Jperps))
energies = zeros(length(Hs), length(Jperps))
sztrack = zeros(length(Hs), length(Jperps))

#0 == 0 Skymrions both layer
#1 == Different SkxNumms in layers 
#2 == Same Skyrmion Nums in layers
#3 == Opposite Skyrmions in layers 
#4 == Something else

for (hidx, h) in enumerate(Hs)
   global H = round(h,sigdigits=5)
   for (jpidx,jperp) in enumerate(Jperps)
      global J_ll = round(jperp,sigdigits=5)
         updateSpins!("C:/Users/tanma/bigboybilayer/secondphasediagram/01.05.2024-Bilayer-decreasingfield-secondphase_Jperp=$J_ll,H=$H.h5", tempLattice)
   
         skxNum0 = round(getSkyrmionNumber(0, tempLattice, vertex),digits=4)
         skxNum1 = round(getSkyrmionNumber(1, tempLattice, vertex),digits=4)
         if skxNum0 == 0.0 && skxNum1 == 0.0
            phases[hidx, jpidx] = 0
         elseif skxNum0 == skxNum1
            phases[hidx, jpidx] = 2
         elseif skxNum0 > 0 && skxNum1 < 0 
            phases[hidx, jpidx] = 3
         elseif (skxNum0 != skxNum1) && (skxNum0 >= 0) && (skxNum1 >= 0) 
            phases[hidx, jpidx] = 1 
         else
            phases[hidx, jpidx] = 4
         end
         skxnums0[hidx, jpidx] = skxNum0
         skxnums1[hidx, jpidx] = skxNum1
         sztrack[jpidx, jpidx] = avgSz(tempLattice)
         energies[jpidx, jpidx] = readEnergy("C:/Users/tanma/bigboybilayer/secondphasediagram/01.05.2024-Bilayer-decreasingfield-secondphase_Jperp=$J_ll,H=$H.h5")[1]
   end
end

#Generate random lattices to see what spin configuration looks like for them
#randomLattices()

plotMiniPhase(Jperps, Hs, phases; flipx=true, flipy=false, xlab=L"J_{||}", ylab=L"H", title="Phase Diagram for Bilayer")
