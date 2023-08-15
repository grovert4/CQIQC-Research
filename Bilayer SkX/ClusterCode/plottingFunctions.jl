using SpinMC_more_more, LinearAlgebra, Plots, PyPlot
include("functions.jl")


function plotSpins(lat::Lattice, bg::Bool = true)
   PyPlot.isjulia_display[] = false
   if length(lat.unitcell.basis) == 2
       xpos,ypos = [],[]
       for i in 1:2:length(lat.sitePositions)
           push!(xpos, lat.sitePositions[i][1])
           push!(ypos, lat.sitePositions[i][2])
       end
       xspin1,yspin1,zspin1 = [],[],[]
       xspin2,yspin2, zspin2 = [],[],[]
       for i in 1:2:length(lat.spins[1,:])
           # push!(xspin1,lat.spins[1,i]/ (sqrt(lat.spins[1,i]^2 + lat.spins[2,i]^2)))
           # push!(yspin1,lat.spins[2,i]/ (sqrt(lat.spins[1,i]^2 + lat.spins[2,i]^2)))

           push!(xspin1,lat.spins[1,i])
           push!(yspin1,lat.spins[2,i])
           push!(zspin1,lat.spins[3,i])
       end
       for i in 2:2:length(lat.spins[1,:])
           # push!(xspin2,lat.spins[1,i]/ (sqrt(lat.spins[1,i]^2 + lat.spins[2,i]^2)))
           # push!(yspin2,lat.spins[2,i]/ (sqrt(lat.spins[1,i]^2 + lat.spins[2,i]^2)))

           push!(xspin2,lat.spins[1,i])
           push!(yspin2,lat.spins[2,i])
           push!(zspin2,lat.spins[3,i])
       end
       len1 = sqrt.(xspin1.^2  + yspin1.^2)
       max1 = maximum(len1)
       min1 = minimum(len1) 
       s1 = ((len1 .- min1) .* ((1 - 0.4)/(max1 - min1)) .+ 0.4) ./ len1
       
       len2 = sqrt.(xspin2.^2  + yspin2.^2)
       max2 = maximum(len2)
       min2 = minimum(len2)
       s2 = ((len2 .- min2) .* ((1 - 0.4)/(max2 - min2)) .+ 0.4) ./ len2

       # s2 = len2 .* 0.75 .* m2 .+ 0.25 .* m2
       # # s2 = 1 ./ (1 .+ log.(m2 ./ len2))
       xspin1 = @. xspin1 .* s1
       yspin1 = @. yspin1 .* s1
       xspin2 = @. xspin2 .* s2
       yspin2 = @. yspin2 .* s2

       PyPlot.close("all")

       PyPlot.rc("xtick", labelsize=8)
       PyPlot.rc("ytick", labelsize=8)

       PyPlot.isjulia_display[] = false
       f1 = PyPlot.figure(dpi=250)
       if bg == true
           pcolor(reshape(xpos,(lat.size[1],lat.size[1])),reshape(ypos,(lat.size[1],lat.size[1])), reshape(zspin1,(lat.size[1],lat.size[1])),vmin = -1, vmax = 1,cmap = "bwr", alpha=0.6)
           colorbar(shrink = 0.585, pad=0.01)
       end
       title(L"J_1 = %$J1, J_2 = %$J2, J_{||} = %$J_ll, H = %$H, D = %$D, A = %$A_ion")
       p1 = PyPlot.quiver(xpos,ypos, xspin1,yspin1, width = 0.00325, headlength = 4,headaxislength = 3.5,angles="xy", scale_units="xy", scale=1)
       ax1 = PyPlot.gca()            
       ax1.set_aspect("equal","box")
       f1.tight_layout(pad=-2.5)
       display(f1)
       f2 = PyPlot.figure(dpi=250)
       if bg == true
           pcolor(reshape(xpos,(lat.size[1],lat.size[1])),reshape(ypos,(lat.size[1],lat.size[1])),reshape(zspin2,(lat.size[1],lat.size[1])), cmap = "bwr", alpha=0.6, vmin = -1, vmax = 1)
           colorbar(shrink = 0.585, pad=0.01)
       end
       title(L"J_1 = %$J1, J_2 = %$J2, J_{||} = %$J_ll, H = %$H, D = %$D, A = %$A_ion")
       p2 = PyPlot.quiver(xpos,ypos, xspin2,yspin2, width = 0.00325, headlength = 4,headaxislength = 3.5,angles="xy", scale_units="xy", scale=1)
       ax2 = PyPlot.gca()
       ax2.set_aspect("equal","box")
       f2.tight_layout(pad=-2.5)
       display(f2)


   elseif length(lat.unitcell.basis) == 1
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
       title(L"J_1 = %$J1, J_2 = %$J2, J_3 = %$J3, J_{||} = %$J_ll, H = %$H, D = %$D, A = %$A_ion")
       pcolor(reshape(xpos,sizer),reshape(ypos,sizer), reshape(zspin1,sizer),vmin = -1, vmax = 1,cmap = "bwr", alpha=0.6) 
       colorbar(shrink = 0.585, pad=0.01)
       p1 = PyPlot.quiver(xpos,ypos, xspin1,yspin1, width = 0.00325, headlength = 4,headaxislength = 3.5,angles="xy", scale_units="xy", scale=1)
       ax1 = PyPlot.gca()            
       ax1.set_aspect("equal","box")
       f1.tight_layout(pad=-2.5)
       return f1
   end
end

function plotSquareSpins(lat::Lattice)
   xpos,ypos = [],[]
   scalMat = inv([1 1/2;0 sqrt(3)/2])
   scalMat2 = inv([1 0.5 0; 0 sqrt(3)/2 0; 0 0 1])
   # scalMat = [1 0; 0 1]
   for i in 1:2:length(lat.sitePositions)
       push!(xpos, sum(scalMat * collect(lat.sitePositions[i][1:2])))
       push!(ypos, (scalMat * collect(lat.sitePositions[i][1:2]))[2])
   end
   xspin1,yspin1,zspin1 = [],[],[]
   xspin2,yspin2, zspin2 = [],[],[]
   for i in 1:2:length(lat.spins[1,:])
       push!(xspin1,normalize(inv(scalMat2) * lat.spins[1:3,i])[1])
       push!(yspin1,normalize(inv(scalMat2) * lat.spins[1:3,i])[2])
       # push!(zspin1,lat.spins[3,i])
   end
   for i in 2:2:length(lat.spins[1,:])
       push!(xspin2,lat.spins[1,i])
       push!(yspin2,lat.spins[2,i])
       push!(zspin2,lat.spins[3,i])
   end 
   # println(xpos)
   # println(ypos)
   p1 = Plots.quiver(xpos,ypos,quiver = (xspin1,yspin1), legend = false, aspect_ratio = :equal); 
   p2 = Plots.quiver(xpos,ypos,quiver = (xspin2,yspin2), legend = false, aspect_ratio = :equal);
   limx = (min(xlims(p1)[1],xlims(p2)[1]), max(xlims(p1)[2],xlims(p2)[2]))
   limy =  (min(ylims(p1)[1],ylims(p2)[1]), max(ylims(p1)[2],ylims(p2)[2]))
   p1 = Plots.quiver(xpos,ypos,quiver = (xspin1,yspin1), legend = false, aspect_ratio = :equal, xlims=limx,ylims=limy)
   display(p1)
   scatter!(xpos, ypos, markersize = 1,legend = false, aspect_ratio = :equal, xlims = limx, ylims = limy)
   p2 = Plots.quiver(xpos,ypos,quiver = (xspin2,yspin2), legend = false, aspect_ratio = :equal, xlims=limx,ylims=limy)
   display(p2)
end

function plotChirality(centers,chirals)
   Plots.scatter(eachrow(hcat(centers...))..., markersize = 5.0, legend = true, zcolor = chirals, colormap=:bwr,clims=(-1,1))
end

function plotStructureFactor(mc::MonteCarlo, layer::Int64)
  N = 256
  correlation = mean(mc.observables.correlation) # The correlation is measured with respect to the spin, i.e. the i-th entry is the correlation dot(S_1,S_i). 
  kx = collect(range(-5pi/3,5pi/3,length=N))
  ky = collect(range(-5pi/3,5pi/3,length=N))
  structurefactor = StructureFactor(mc, layer)
  
  # Plot result
  hmap = heatmap(kx,ky,structurefactor,aspect_ratio=1,xrange=(-5*pi/3,5*pi/3),yrange=(-5*pi/3,5*pi/3),xlabel="kx", ylabel="ky",background_color=:black,clims=(0,1))
  xs = [4.0pi/3.0, 2.0pi/3.0, -2.0pi/3.0, -4.0pi/3.0, -2.0pi/3.0, 2.0pi/3.0, 4.0pi/3.0]
  ys = [0.0, 2.0pi/sqrt(3.0), 2.0pi/sqrt(3.0), 0.0, -2.0pi/sqrt(3.0), -2.0pi/sqrt(3.0), 0.0]
  plot!(xs, ys, legend=false)  
end


function subPlotsSpins(lat, bg::Bool=true)
   PyPlot.isjulia_display[] = false
   xpos,ypos = [],[]
   for i in 1:2:length(lat.sitePositions)
      push!(xpos, lat.sitePositions[i][1])
      push!(ypos, lat.sitePositions[i][2])
   end
   xspin1,yspin1,zspin1 = [],[],[]
   xspin2,yspin2, zspin2 = [],[],[]
   for i in 1:2:length(lat.spins[1,:])

      push!(xspin1,lat.spins[1,i])
      push!(yspin1,lat.spins[2,i])
      push!(zspin1,lat.spins[3,i])
   end
   for i in 2:2:length(lat.spins[1,:])
      # push!(xspin2,lat.spins[1,i]/ (sqrt(lat.spins[1,i]^2 + lat.spins[2,i]^2)))
      # push!(yspin2,lat.spins[2,i]/ (sqrt(lat.spins[1,i]^2 + lat.spins[2,i]^2)))

      push!(xspin2,lat.spins[1,i])
      push!(yspin2,lat.spins[2,i])
      push!(zspin2,lat.spins[3,i])
   end
   len1 = sqrt.(xspin1.^2  + yspin1.^2)
   max1 = maximum(len1)
   min1 = minimum(len1) 
   s1 = ((len1 .- min1) .* ((1 - 0.4)/(max1 - min1)) .+ 0.4) ./ len1
   
   len2 = sqrt.(xspin2.^2  + yspin2.^2)
   max2 = maximum(len2)
   min2 = minimum(len2)
   s2 = ((len2 .- min2) .* ((1 - 0.4)/(max2 - min2)) .+ 0.4) ./ len2

   # s2 = len2 .* 0.75 .* m2 .+ 0.25 .* m2
   # # s2 = 1 ./ (1 .+ log.(m2 ./ len2))
   xspin1 = @. xspin1 .* s1
   yspin1 = @. yspin1 .* s1
   xspin2 = @. xspin2 .* s2
   yspin2 = @. yspin2 .* s2

   PyPlot.close("all")

   PyPlot.rc("xtick", labelsize=8)
   PyPlot.rc("ytick", labelsize=8)

   PyPlot.isjulia_display[] = false
   
   f1 = PyPlot.figure(L"J_1 = %$J1, J_2 = %$J2, J_3 = %$J3, J_{||} = %$J_ll, H = %$H, D = %$D, A = %$A_ion", dpi=450)
   f1.suptitle(L"J_1 = %$J1, J_2 = %$J2, J_3 = %$J3, J_{||} = %$J_ll, H = %$H, D = %$D, A = %$A_ion", fontsize=11.5)
   
   PyPlot.subplot(121)
   PyPlot.title("Layer 0", fontsize=10)
   ax1 = PyPlot.gca()            

   if bg == true
      pcolor(reshape(xpos,(lat.size[1],lat.size[1])),reshape(ypos,(lat.size[1],lat.size[1])), reshape(zspin1,(lat.size[1],lat.size[1])),vmin = -1, vmax = 1,cmap = "bwr", alpha=0.6)
   end
   p1 = PyPlot.quiver(xpos,ypos, xspin1,yspin1, width = 0.00325, headlength = 4,headaxislength = 3.5,angles="xy", scale_units="xy", scale=1)
   ax1.set_aspect("equal","box")
   
   PyPlot.subplot(122)
   PyPlot.title("Layer 1", fontsize=10)
   if bg == true
      pcolor(reshape(xpos,(lat.size[1],lat.size[1])),reshape(ypos,(lat.size[1],lat.size[1])),reshape(zspin2,(lat.size[1],lat.size[1])), cmap = "bwr", alpha=0.6, vmin = -1, vmax = 1)
   end
   p2 = PyPlot.quiver(xpos,ypos, xspin2,yspin2, width = 0.00325, headlength = 4,headaxislength = 3.5,angles="xy", scale_units="xy", scale=1)
   ax2 = PyPlot.gca()
   ax2.set_yticklabels([])
   ax2.set_yticks([])
   ax2.set_aspect("equal","box")
   # colorbar(shrink = 0.585, pad=0.01    
   f1.tight_layout(pad=0.0)
   f1.subplots_adjust(top=1.42)
   display(f1) 

end

function DetailedMonoPlot(mc, lat, vertex = none, bg = true)
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
      energy = round(means(mc.observables.energy)[1],sigdigits=4)
      std = round(std_errors(mc.observables.energy)[1],sigdigits=3)
      sf1 = StructureFactor(mc, 0)
      maxvalsidx1 = partialsortperm(vec(sf1), 1:3,rev=true)
      maxvals1 = round.(partialsort(vec(sf1), 1:3,rev=true),sigdigits=3)
      z11 = maxvals1[1]
      z12 = maxvals1[2]
      z13 = maxvals1[3]
      kx = collect(range(-5pi/3,5pi/3,length=256))
      ky = collect(range(-5pi/3,5pi/3, length=256))
      qxs1 = round.(kx[Int64.(ceil.(maxvalsidx1 ./ 256))],sigdigits=3)
      x11 = qxs1[1]
      x12 = qxs1[2]
      x13 = qxs1[3] 
      qys1 = round.(ky[(maxvalsidx1 .% 256)],sigdigits=3) 
      y11 = qys1[1]
      y12 = qys1[2]
      y13 = qys1[3]

      annotate(L"SkX Number: %$skxnum0, Chirality: %$chirality0, Energy: %$energy \pm %$std,", xy=(0.5, -0.1), xycoords="axes fraction", ha = "center", fontsize=9)
      annotate(L"SF Peak: (%$x11, %$y11): %$z11", xy=(0.5, -0.15), xycoords="axes fraction", ha = "center", fontsize=9)
  end

   display(f1)
end

function DetailedsubPlotsSpins(mc, lat, vertex = none, bg=true)
   PyPlot.isjulia_display[] = false
   xpos,ypos = [],[]
   for i in 1:2:length(lat.sitePositions)
      push!(xpos, lat.sitePositions[i][1])
      push!(ypos, lat.sitePositions[i][2])
   end
   xspin1,yspin1,zspin1 = [],[],[]
   xspin2,yspin2, zspin2 = [],[],[]
   for i in 1:2:length(lat.spins[1,:])

      push!(xspin1,lat.spins[1,i])
      push!(yspin1,lat.spins[2,i])
      push!(zspin1,lat.spins[3,i])
   end
   for i in 2:2:length(lat.spins[1,:])
      # push!(xspin2,lat.spins[1,i]/ (sqrt(lat.spins[1,i]^2 + lat.spins[2,i]^2)))
      # push!(yspin2,lat.spins[2,i]/ (sqrt(lat.spins[1,i]^2 + lat.spins[2,i]^2)))

      push!(xspin2,lat.spins[1,i])
      push!(yspin2,lat.spins[2,i])
      push!(zspin2,lat.spins[3,i])
   end
   len1 = sqrt.(xspin1.^2  + yspin1.^2)
   max1 = maximum(len1)
   min1 = minimum(len1) 
   s1 = ((len1 .- min1) .* ((1 - 0.4)/(max1 - min1)) .+ 0.4) ./ len1
   
   len2 = sqrt.(xspin2.^2  + yspin2.^2)
   max2 = maximum(len2)
   min2 = minimum(len2)
   s2 = ((len2 .- min2) .* ((1 - 0.4)/(max2 - min2)) .+ 0.4) ./ len2

   # s2 = len2 .* 0.75 .* m2 .+ 0.25 .* m2
   # # s2 = 1 ./ (1 .+ log.(m2 ./ len2))
   xspin1 = @. xspin1 .* s1
   yspin1 = @. yspin1 .* s1
   xspin2 = @. xspin2 .* s2
   yspin2 = @. yspin2 .* s2

   PyPlot.close("all")

   PyPlot.rc("xtick", labelsize=8)
   PyPlot.rc("ytick", labelsize=8)

   PyPlot.isjulia_display[] = false
   
   f1 = PyPlot.figure(L"J_1 = %$J1, J_2 = %$J2, J_{||} = %$J_ll, H = %$H, D = %$D, A = %$A_ion", dpi=450)
   f1.suptitle(L"J_1 = %$J1, J_2 = %$J2, J_{||} = %$J_ll, H = %$H, D = %$D, A = %$A_ion", fontsize=11.5)
   
   PyPlot.subplot(121)
   PyPlot.title("Layer 0", fontsize=10)
   ax1 = PyPlot.gca()            

   if bg == true
      pcolor(reshape(xpos,(lat.size[1],lat.size[1])),reshape(ypos,(lat.size[1],lat.size[1])), reshape(zspin1,(lat.size[1],lat.size[1])),vmin = -1, vmax = 1,cmap = "bwr", alpha=0.6)
   end
   p1 = PyPlot.quiver(xpos,ypos, xspin1,yspin1, width = 0.00325, headlength = 4,headaxislength = 3.5,angles="xy", scale_units="xy", scale=1)
   ax1.set_aspect("equal","box")
   
   PyPlot.subplot(122)
   PyPlot.title("Layer 1", fontsize=10)
   if bg == true
      pcolor(reshape(xpos,(lat.size[1],lat.size[1])),reshape(ypos,(lat.size[1],lat.size[1])),reshape(zspin2,(lat.size[1],lat.size[1])), cmap = "bwr", alpha=0.6, vmin = -1, vmax = 1)
   end
   p2 = PyPlot.quiver(xpos,ypos, xspin2,yspin2, width = 0.00325, headlength = 4,headaxislength = 3.5,angles="xy", scale_units="xy", scale=1)
   ax2 = PyPlot.gca()
   ax2.set_yticklabels([])
   ax2.set_yticks([])
   ax2.set_aspect("equal","box")
   # colorbar(shrink = 0.585, pad=0.01)    
   f1.tight_layout(pad=0.0)
   f1.subplots_adjust(top=1.42)
 

   if length(vertex) != 0
      skxnum0 = round(getSkyrmionNumber(0, lat, vertex), sigdigits=3)
      skxnum1 = round(getSkyrmionNumber(1, lat, vertex), sigdigits=3)
      chirality0 = round(getScalarChirality(0,lat,vertex),sigdigits=3)
      chirality1 = round(getScalarChirality(1,lat,vertex), sigdigits=3)
      energy = round(means(mc.observables.energy)[1],sigdigits=4)
      std = round(std_errors(mc.observables.energy)[1],sigdigits=3)
      sf1 = StructureFactor(mc, 0)
      sf2 = StructureFactor(mc, 1)
      maxvalsidx1 = partialsortperm(vec(sf1), 1:3, rev = true)
      maxvals1 = round.(partialsort(vec(sf1), 1:3, rev = true),sigdigits=3)
      z11 = maxvals1[1]
      z12 = maxvals1[2]
      z13 = maxvals1[3]
      kx = collect(range(-5pi/3,5pi/3,length=256))
      ky = collect(range(-5pi/3,5pi/3, length=256))
      qxs1 = round.(kx[Int64.(ceil.(maxvalsidx1 ./ 256))],sigdigits=3)
      x11 = qxs1[1]
      x12 = qxs1[2]
      x13 = qxs1[3] 
      qys1 = round.(ky[Int64.(ceil.(maxvalsidx1 .% 256))],sigdigits=3)
      y11 = qys1[1]
      y12 = qys1[2]
      y13 = qys1[3]
      maxvalsidx2 = partialsortperm(vec(sf2), 1:3, rev = true)
      maxvals2 = round.(partialsort(vec(sf2), 1:3, rev = true),sigdigits=3)
      z21 = maxvals2[1]
      z22 = maxvals2[2]
      z23 = maxvals2[3]
      qxs2 = round.(kx[Int64.(ceil.(maxvalsidx2 ./ 256))],sigdigits=3)
      x21 = qxs2[1]
      x22 = qxs2[2]
      x23 = qxs2[3] 

      qys2 = round.(ky[Int64.(ceil.(maxvalsidx1 .% 256))],sigdigits=3)
      y21 = qys2[1]
      y22 = qys2[2]
      y23 = qys2[3]

      annotate(L"SkX Number: %$skxnum0 / %$skxnum1, Chirality: %$chirality0 / %$ chirality1, Energy: %$energy \pm %$std,", xy=(-0.05, -0.2), xycoords="axes fraction", ha = "center", fontsize=7)
      annotate(L"SF Peak: (%$x11, %$y11) %$z11 / (%$x21 , %$y21) %$z21", xy=(-0.05, -0.27), xycoords="axes fraction", ha = "center", fontsize=7)
   end
   
   display(f1) 

end


function plotSx(lat::Lattice)
    if length(lat.unitcell.basis) == 2
        PyPlot.isjulia_display[] = false
        xpos,ypos = [],[]
        for i in 1:2:length(lat.sitePositions)
            push!(xpos, lat.sitePositions[i][1])
            push!(ypos, lat.sitePositions[i][2])
        end
        xspin1,yspin1,zspin1 = [],[],[]
        xspin2,yspin2, zspin2 = [],[],[]
        for i in 1:2:length(lat.spins[1,:])
            push!(xspin1,lat.spins[1,i])
            push!(yspin1,lat.spins[2,i])
            # push!(zspin1,lat.spins[3,i])
        end
        for i in 2:2:length(lat.spins[1,:])
            push!(xspin2,lat.spins[1,i])
            push!(yspin2,lat.spins[2,i])
            # push!(zspin2,lat.spins[3,i])
        end
        # len1 = sqrt.(xspin1.^2  + yspin1.^2)
        # max1 = maximum(len1)
        # min1 = minimum(len1) 
        # s1 = ((len1 .- min1) .* ((1 - 0.4)/(max1 - min1)) .+ 0.4) ./ len1
        
        # len2 = sqrt.(xspin2.^2  + yspin2.^2)
        # max2 = maximum(len2)
        # min2 = minimum(len2)
        # s2 = ((len2 .- min2) .* ((1 - 0.4)/(max2 - min2)) .+ 0.4) ./ len2

        # xspin1 = @. xspin1 .* s1
        # yspin1 = @. yspin1 .* s1
        # xspin2 = @. xspin2 .* s2
        # yspin2 = @. yspin2 .* s2

        PyPlot.close("all")

        PyPlot.rc("xtick", labelsize=8)
        PyPlot.rc("ytick", labelsize=8)

        PyPlot.isjulia_display[] = false
        
        f1 = PyPlot.figure(L"J_1 = %$J1, J_2 = %$J2, J_3 = %$J3, J_{||} = %$J_ll, H = %$H, D = %$D, A = %$A_ion", dpi=450)
        f1.suptitle(L"S_x,  J_1 = %$J1, J_2 = %$J2, J_3 = %$J3, J_{||} = %$J_ll, H = %$H, D = %$D, A = %$A_ion", fontsize=11.5)
        
        PyPlot.subplot(121)
        PyPlot.title("Layer 0", fontsize=10)
        ax1 = PyPlot.gca()            

        pcolor(reshape(xpos,(lat.size[1],lat.size[1])),reshape(ypos,(lat.size[1],lat.size[1])), reshape(xspin1,(lat.size[1],lat.size[1])),vmin = -1, vmax = 1,cmap = "bwr", alpha=0.6)
        # p1 = PyPlot.quiver(xpos,ypos, xspin1,yspin1, width = 0.00325, headlength = 4,headaxislength = 3.5,angles="xy", scale_units="xy", scale=1)
        ax1.set_aspect("equal","box")
        
        PyPlot.subplot(122)
        PyPlot.title("Layer 1", fontsize=10)
        pcolor(reshape(xpos,(lat.size[1],lat.size[1])),reshape(ypos,(lat.size[1],lat.size[1])),reshape(xspin2,(lat.size[1],lat.size[1])), cmap = "bwr", alpha=0.6, vmin = -1, vmax = 1)
        # p2 = PyPlot.quiver(xpos,ypos, xspin2,yspin2, width = 0.00325, headlength = 4,headaxislength = 3.5,angles="xy", scale_units="xy", scale=1)
        ax2 = PyPlot.gca()
        ax2.set_yticklabels([])
        ax2.set_yticks([])
        ax2.set_aspect("equal","box")
        # colorbar(shrink = 0.585, pad=0.01)    
        f1.tight_layout(pad=0.0)
        f1.subplots_adjust(top=1.42)
        display(f1) 
    elseif length(lat.unitcell.basis) == 1
        xpos,ypos = [],[]

        for i in 1:length(lat.sitePositions)
            push!(xpos, lat.sitePositions[i][1])
            push!(ypos, lat.sitePositions[i][2])
        end
        xspin1,yspin1,zspin1 = [],[],[]
        for i in 1:length(lat.spins[1,:])
            push!(xspin1,lat.spins[1,i])
            push!(yspin1,lat.spins[2,i])
            # push!(zspin1,lat.spins[3,i])
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
        title(L"J_1 = %$J1, J_2 = %$J2, J_3 = %$J3, J_{||} = %$J_ll, H = %$H, D = %$D, A = %$A_ion")
        pcolor(reshape(xpos,sizer),reshape(ypos,sizer), reshape(xspin1,sizer),vmin = -1, vmax = 1,cmap = "bwr", alpha=0.6) 
        colorbar(shrink = 0.585, pad=0.01)
        # p1 = PyPlot.quiver(xpos,ypos, xspin1,yspin1, width = 0.00325, headlength = 4,headaxislength = 3.5,angles="xy", scale_units="xy", scale=1)
        ax1 = PyPlot.gca()            
        ax1.set_aspect("equal","box")
        f1.tight_layout(pad=-2.5)
        return f1
    end
end
