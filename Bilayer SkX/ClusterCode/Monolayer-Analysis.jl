function updateSpins(file, lat)
   file = h5open(file)["mc"]
   energy = read(file["observables"]["energyDensity"]["mean"])
   std = read(file["observables"]["energyDensity"]["error"])
   sites = parse.(Int64,collect(keys(read(file["lattice"]["spins"]))))
   spins = collect(values(read(file["lattice"]["spins"])))

   sorted = sortperm(sites)
   lat.spins = reshape(vcat(spins[sorted]...),(3,lat.length))  
end

function plotReadSpins(file, lat, vertex)
   updateSpins(file,lat)
   PyPlot.isjulia_display[] = false  
   xpos,ypos = [],[]

   for i in 1:length(lat.sitePositions)
       push!(xpos, lat.sitePositions[i][1])
       push!(ypos, lat.sitePositions[i][2])
   end
   xspin1,yspin1,zspin1 = [],[],[]
   for i in 1:length(lat.spins[1,:])
       push!(xspin1,lat.spins[1,i])
       push!(yspin1,lat.spins[2,i])
       push!(zspin1,lat.spins[3,i])
   end
   len1 = sqrt.(xspin1.^2  + yspin1.^2)
   max1 = maximum(len1)
   min1 = minimum(len1) 
   s1 = ((len1 .- min1) .* ((1 - 0.4)/(max1 - min1)) .+ 0.4) ./ len1
   xspin1 = @. xspin1 .* s1
   yspin1 = @. yspin1 .* s1

   PyPlot.close("all")

   PyPlot.rc("xtick", labelsize=8)
   PyPlot.rc("ytick", labelsize=8)

   sizer = (lat.size[1],lat.size[2]) 

   PyPlot.isjulia_display[] = false
   f1 = PyPlot.figure(dpi=250)
   title(L"J_1 = %$J1, J_2 = %$J2, H = %$H, D = %$D, A = %$A_ion")
   pcolor(reshape(xpos,sizer),reshape(ypos,sizer), reshape(zspin1,sizer),vmin = -1, vmax = 1,cmap = "bwr", alpha=0.6) 
   colorbar(shrink = 0.585, pad=0.01)
   p1 = PyPlot.quiver(xpos,ypos, xspin1,yspin1, width = 0.00325, headlength = 4,headaxislength = 3.5,angles="xy", scale_units="xy", scale=1)
   ax1 = PyPlot.gca()            
   ax1.set_aspect("equal","box")
   f1.tight_layout(pad=-2.5)

   if length(vertex) != 0
      skxnum0 = round(getSkyrmionNumber(0, lat, vertex), sigdigits=3)
      chirality0 = round(getScalarChirality(0,lat,vertex),sigdigits=3)
      energy = round(energy,sigdigits=4)
      std = round(std,sigdigits=3)
      # sf1 = StructureFactor(mc, 0)
      # maxvalsidx1 = partialsortperm(vec(sf1), 1:3,rev=true)
      # maxvals1 = round.(partialsort(vec(sf1), 1:3,rev=true),sigdigits=3)
      # z11 = maxvals1[1]
      # z12 = maxvals1[2]
      # z13 = maxvals1[3]
      # kx = collect(range(-5pi/3,5pi/3,length=256))
      # ky = collect(range(-5pi/3,5pi/3, length=256))
      # qxs1 = round.(kx[Int64.(ceil.(maxvalsidx1 ./ 256))],sigdigits=3)
      # x11 = qxs1[1]
      # x12 = qxs1[2]
      # x13 = qxs1[3] 
      # qys1 = round.(ky[(maxvalsidx1 .% 256)],sigdigits=3) 
      # y11 = qys1[1]
      # y12 = qys1[2]
      # y13 = qys1[3]

      annotate(L"SkX Number: %$skxnum0, Chirality: %$chirality0, Energy: %$energy \pm %$std,", xy=(0.5, -0.1), xycoords="axes fraction", ha = "center", fontsize=9)
      # annotate(L"SF Peak: (%$x11, %$y11): %$z11", xy=(0.5, -0.15), xycoords="axes fraction", ha = "center", fontsize=9)
  end

   display(f1)

end

function avgSz(lat)
   avgSz = 0
   for i in 1:length(lat.sitePositions)
      avgSz += lat.spins[3,i]
   end
   return avgSz/length(lat.sitePositions)
end

using Plots, LinearAlgebra, SpinMC_more_more, HDF5, JSON
include("functions.jl")
include("plottingFunctions.jl")

inputFile = JSON.parsefile("./Input_Files/"*ARGS[1]*".json")

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

for (j2idx, j2) in enumerate(J2s)
   global J2 = round(j2,sigdigits=5)
   
   for (hidx,h) in enumerate(Hs)
   
      global H = round(h,sigdigits=5)
   
      # plotReadSpins("Official-Cluster-Run-1-MonoLayer/H=$H,J2=$J2", tempLattice, vertex)
      tempLattice = updateSpins("Official-Cluster-Run-1-MonoLayer/H=$H,J2=$J2", tempLattice)
      skxNum = round(getSkyrmionNumber(0, tempLattice, vertex),digits=2)
      avgSz = avgSz(tempLattice)
      if round(avgSz,digits=2) >= 0.8
         phases[hidx,j2idx] = 0
      elseif skxNum != 0.0 
         phases[hidx,j2idx] = 2 
      else
         phases[hidx,j2idx] = 1
      end
   end
end

pDfigure = PyPlot.figure()
phaseDiagram = pcolor(Hs, J2s, phases)

