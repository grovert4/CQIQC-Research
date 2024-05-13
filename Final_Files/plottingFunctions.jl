using SpinMC_more_more, LinearAlgebra, Plots, PyPlot
include("functions.jl")

##Plot local scalar chiralities (need centers using getCenters() and chirals using localChiralities())
function plotChirality(centers,chirals, scale = false)
    if scale == true
        m = 1.5 * max(abs.(findmax(chirals)[1]), abs.(findmin(chirals)[1]))
    else
        m = 1
    end
    Plots.scatter(eachrow(hcat(centers...))..., markersize = 5.0, legend = true, zcolor = chirals, colormap=:bwr,clims=(-m,m))
end

##Plot Structure Factor using <S(0) . S(r)> for lattice 
function plot_file_StructureFactor(lat::Lattice{D,N}, layer::Int64)
   N = 256
   correlation = getCorrelation(lat) # The correlation is measured with respect to the spin, i.e. the i-th entry is the correlation dot(S_1,S_i). 
   kx = collect(range(-5pi/3,5pi/3,length=N))
   ky = collect(range(-5pi/3,5pi/3,length=N))
   structurefactor = readStructureFactor(lat, layer)
   
   # Plot result
   hmap = heatmap(kx,ky,structurefactor,aspect_ratio=1,xrange=(-5*pi/3,5*pi/3),yrange=(-5*pi/3,5*pi/3),xlabel="kx", ylabel="ky",background_color=:black, clims=(0,1))
   xs = [4.0pi/3.0, 2.0pi/3.0, -2.0pi/3.0, -4.0pi/3.0, -2.0pi/3.0, 2.0pi/3.0, 4.0pi/3.0]
   ys = [0.0, 2.0pi/sqrt(3.0), 2.0pi/sqrt(3.0), 0.0, -2.0pi/sqrt(3.0), -2.0pi/sqrt(3.0), 0.0]
   plot!(xs, ys, legend=false)  
end

##Plot Monolayer Spins for given input file
function plotMonoReadSpins(file, lat::Lattice{D,N}, vertex, tit=true)
    updateSpins!(file,lat)
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
    if tit
        title(L"J_1 = %$J1, J_2 = %$J2, H = %$H, D = %$D, A = %$A_ion")
    end
    pcolor(reshape(xpos,sizer),reshape(ypos,sizer), reshape(zspin1,sizer),vmin = -1, vmax = 1,cmap = "bwr", alpha=0.6) 
    colorbar(shrink = 0.585, pad=0.01)
    p1 = PyPlot.quiver(xpos,ypos, xspin1,yspin1, width = 0.00325, headlength = 4,headaxislength = 3.5,angles="xy", scale_units="xy", scale=1)
    ax1 = PyPlot.gca()            
    ax1.set_aspect("equal","box")
    f1.tight_layout(pad=-2.5)
 
    if length(vertex) != 0
       skxnum0 = round(getSkyrmionNumber(0, lat, vertex), sigdigits=5)
       chirality0 = round(getScalarChirality(0,lat,vertex),sigdigits=5)
       energy, std = readEnergy(file)
       energy = round(energy, sigdigits=4)
       std = round(std, sigdigits=3)
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

##Plot Bilayer Spins for given input file
function plotBiReadSpins(file, lat::Lattice{D,N}, vertex = none, bg=true)
   updateSpins!(file,lat)
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
      energy, std = round.(readEnergy(file), sigdigits=3)
      skxnum0 = round(getSkyrmionNumber(0, lat, vertex), sigdigits=3)
      skxnum1 = round(getSkyrmionNumber(1, lat, vertex), sigdigits=3)
      chirality0 = round(getScalarChirality(0,lat,vertex),sigdigits=3)
      chirality1 = round(getScalarChirality(1,lat,vertex), sigdigits=3)
      annotate(L"SkX Number: %$skxnum0 / %$skxnum1, Chirality: %$chirality0 / %$ chirality1, Energy: %$energy \pm %$std,", xy=(-0.05, -0.2), xycoords="axes fraction", ha = "center", fontsize=7)
   end

   display(f1)  
end

##Plot <S(r)> Correlations have to fix this
function plot_file_Correlations(lat::Lattice{D,N}, layer::Int64,index::Int64; sf)
    N = 256
    # correlation = getCorrelation(lat) # The correlation is measured with respect to the spin, i.e. the i-th entry is the correlation dot(S_1,S_i). 
    kx = collect(range(-5pi/3,5pi/3,length=N))
    ky = collect(range(-5pi/3,5pi/3,length=N))
    structurefactor = nothing
    if length(sf) == 0
        structurefactor = abs.(readCorrelations(lat, layer)[index])
    else
        structurefactor = abs.(sf)
    end
    # Plot result
    hmap = heatmap(kx,ky,structurefactor,aspect_ratio=1,xrange=(-5*pi/3,5*pi/3),yrange=(-5*pi/3,5*pi/3),xlabel="kx", ylabel="ky",background_color=:black)
    xs = [4.0pi/3.0, 2.0pi/3.0, -2.0pi/3.0, -4.0pi/3.0, -2.0pi/3.0, 2.0pi/3.0, 4.0pi/3.0]
    ys = [0.0, 2.0pi/sqrt(3.0), 2.0pi/sqrt(3.0), 0.0, -2.0pi/sqrt(3.0), -2.0pi/sqrt(3.0), 0.0]
    plot!(xs, ys, legend=false)  
end
