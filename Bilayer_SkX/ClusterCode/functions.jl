#using Pkg
#Pkg.add(url="https://github.com/grovert4/SpinMC_more_more.jl")
using SpinMC_more_more, LinearAlgebra #, Plots, ColorSchemes, PyPlot

function getSkyrmionNumber(layer,lat,vertex)
    Q = 0
    for i in vertex
        Si = collect(getSpin(lat, i[1] + layer))
        Sj = collect(getSpin(lat, i[2] + layer))
        Sk = collect(getSpin(lat, i[3] + layer))            
        # rho = sqrt(2 * (1 + dot(Si,Sj)) * (1 + dot(Sj, Sk)) * (1 + dot(Sk, Si)))
        Q += atan( dot(Si, cross(Sj,Sk)) / (1 + dot(Si, Sj) + dot(Sj, Sk) + dot(Sk, Si) ) )
    end 
    return Q/(2 * pi)  

end

function getScalarChirality(layer, lat, vertex)
    chi = 0
    for i in vertex
        Si = collect(getSpin(lat, i[1] + layer))
        Sj = collect(getSpin(lat, i[2] + layer))
        Sk = collect(getSpin(lat, i[3] + layer))            
        chi += dot(Si, cross(Sj,Sk))/((size(lat)[1])^2)
    end    
    return chi
end

function localChiralities(layer,lat, vertex)
    chiralsA = []
    chiralsB = []
    for i in vertex
        Si = collect(getSpin(lat, i[1] + layer))
        Sj = collect(getSpin(lat, i[2] + layer))
        Sk = collect(getSpin(lat, i[3] + layer))            
        if layer == 0
            push!(chiralsA, dot(Si, cross(Sj,Sk)))
        else
            push!(chiralsB, dot(Si, cross(Sj,Sk)))
        end
    end
    if layer == 0
        return chiralsA
    else
        return chiralsB
    end
end

function getCenters(lat)
    centerpos = []
    if length(lat.unitcell.basis) == 2
        for (ind, i) in enumerate(lat.sitePositions)
            if i[3] == 0
                if (i .+ a1 in lat.sitePositions) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + 50)
                    push!(centerpos, collect(i .+ a1 ./ 2 .+ (0,sqrt(3)/6,0))[1:2])
                end
                if (findmin(x -> norm(i .+ a2 .- x),lat.sitePositions)[2] == ind + 2) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + 50)
                    push!(centerpos, collect(i .+ (0,sqrt(3)/3,0))[1:2])
            
                end
            end
        end
    elseif length(lat.unitcell.basis) == 1
        for (ind, i) in enumerate(lat.sitePositions)
            if (i .+ a1 in lat.sitePositions) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + 25)
                push!(centerpos, collect(i .+ a1 ./ 2 .+ (0,sqrt(3)/6))[1:2])
            end
            if (findmin(x -> norm(i .+ a2 .- x),lat.sitePositions)[2] == ind + 1) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + 25)
                push!(centerpos, collect(i .+ (0,sqrt(3)/3))[1:2])
            end
        end
    end

    return centerpos
end

function getVertex(lat)
    if size(lat)[1] != 24
        return "Lattice size is not 24x24"
    else
        vertex = []
        if length(lat.unitcell.basis) == 2
            for (ind, i) in enumerate(lat.sitePositions)
                if i[3] == 0
                    if (i .+ a1 in lat.sitePositions) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + 50)
                        push!(vertex, [ind, ind + 48, ind + 50])
                    end
                    if (findmin(x -> norm(i .+ a2 .- x),lat.sitePositions)[2] == ind + 2) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + 50)
                        push!(vertex, [ind, ind + 50, ind+2])
                    end
                end
            end
            return vertex
        elseif length(lat.unitcell.basis) == 1
            for (ind, i) in enumerate(lat.sitePositions)
                if (i .+ a1 in lat.sitePositions) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + 25)
                    push!(vertex, [ind, ind + 24, ind + 25])
                end
                if (findmin(x -> norm(i .+ a2 .- x),lat.sitePositions)[2] == ind + 1) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + 25)
                    push!(vertex, [ind, ind + 25, ind+1])
                end
            end
            return vertex
        end
    end
end



function runAnneal(t0,tf,lat,thermSweeps,MeasureSweeps, coolRate, H, J2, outfile)
   betas = [1/t0:coolRate:1/tf]
   monte = nothing
   for (ind,beta) in enumerate(betas) 
      thermalizationSweeps = thermSweeps
      measurementSweeps = 0
      if ind == 1
        m = MonteCarlo(lat, beta, thermalizationSweeps, measurementSweeps, reportInterval = 50000, rewrite = true);
        run!(m)
      else
        if ind == length(ts)
           thermalizationSweeps = 0
           measurementSweeps = MeasureSweeps
        end
        m = MonteCarlo(monte.lattice, beta, thermalizationSweeps, measurementSweeps, reportInterval = 50000, rewrite = false);

        if ind == length(ts)
           h = round(H,sigdigits=3)
           j2 = round(J2,sigdigits=3)
           run!(m, outfile = outfile)
        else
           run!(m)
        end  
      end
      monte = deepcopy(m)
   end
   return monte
end
