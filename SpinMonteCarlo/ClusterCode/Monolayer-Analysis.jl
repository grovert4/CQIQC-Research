using Plots, LinearAlgebra, SpinMC_more_more, HDF5, JSON, LaTeXStrings
include("functions.jl")
include("plottingFunctions.jl")

inputFile = JSON.parsefile("./inputParametersMonoLayer.json")

Hs = collect(range(inputFile["H_min"], inputFile["H_max"] ,length = inputFile["H_length"]))
J2s = collect(range(inputFile["J2_min"], inputFile["J2_max"] ,length = inputFile["J2_length"]))

a1 = (1.0 , 0.0)
a2 = (-1/2 , sqrt(3)/2)
UCtemp = UnitCell(a1 , a2)
b1 = addBasisSite!(UCtemp, (0.0, 0.0)) 
tempLattice = Lattice(UCtemp, (inputFile["System_Size"], inputFile["System_Size"]))
vertex = getVertex(tempLattice)

phases = zeros(length(Hs), length(J2s))

##0 == Fully Polarized
##1 == Spin Spiral
##2 == Skyrmion
##3 == Error/Recompute for params
global c = 0  
for (j2idx, j2) in enumerate(J2s)
   global J2 = round(j2,sigdigits=5)
   
   for (hidx,h) in enumerate(Hs)
      global c += 1
   
      global H = round(h,sigdigits=5)
      if true
         println("H = $H, J2 = $J2")
         # plotReadSpins("Official-Cluster-Run-1-MonoLayer/H=$H,J2=$J2", tempLattice, vertex)
         updateSpins!("D:/CQIQC_Data/MonoLayer_Runs/08.17.2023_Monolayer_H=$H,J2=$J2.h5", tempLattice)
         skxNum = round(getSkyrmionNumber(0, tempLattice, vertex),digits=2)
         sZavg = avgSz(tempLattice)
         if round(sZavg,digits=2) >= 0.8
            phases[hidx,j2idx] = 0
         elseif skxNum != 0.0 
            phases[hidx,j2idx] = 2 
         else
            phases[hidx,j2idx] = 1
         end
      else
         phases[hidx,j2idx] = 3
         println(c)
      end
   end
end

pDfigure = PyPlot.figure()
phaseDiagram = pcolor(Hs, J2s, phases)
display(phaseDiagram)
