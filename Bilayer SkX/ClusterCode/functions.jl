using SpinMC_more_more, LinearAlgebra, Plots, ColorSchemes, PyPlot
ioff()

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
            if (i .+ a1 in lat.sitePositions) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + numBasis * (size + 1))
               push!(centerpos, collect(i .+ a1 ./ 2 .+ (0,sqrt(3)/6))[1:2])
            end
            if (findmin(x -> norm(i .+ a2 .- x),lat.sitePositions)[2] == ind + numBasis) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + numBasis * (size + 1))
               push!(centerpos, collect(i .+ (0,sqrt(3)/3))[1:2])
            end
        end
    end

    return centerpos
end

function getVertex(lat)
   vertex = []
   size = lat.size[1]
   numBasis = length(lat.unitcell.basis)

   if numBasis == 2   
      for (ind, i) in enumerate(lat.sitePositions)
         if i[3] == 0
            if (i .+ a1 in lat.sitePositions) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + numBasis * (size + 1))
               push!(vertex, [ind, ind + size * numBasis, ind + numBasis * (size + 1)])
            end
            if (findmin(x -> norm(i .+ a2 .- x),lat.sitePositions)[2] == ind + numBasis) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + numBasis * (size + 1))
               push!(vertex, [ind, ind + numBasis * (size + 1), ind + numBasis])
            end
         end
      end
   elseif numBasis == 1
      for (ind, i) in enumerate(lat.sitePositions)
         if (i .+ a1 in lat.sitePositions) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + numBasis * (size + 1))
            push!(vertex, [ind, ind + size * numBasis, ind + numBasis * (size + 1)])
         end
         if (findmin(x -> norm(i .+ a2 .- x),lat.sitePositions)[2] == ind + numBasis) && (findmin(x -> norm(i .+ a2 .+ a1 .- x),lat.sitePositions)[2] == ind + numBasis * (size + 1))
            push!(vertex, [ind, ind + numBasis * (size + 1), ind + numBasis])
         end
      end
   end
   return vertex
end

function StructureFactor(mc, layer)
   N = 256
   correlation = mean(mc.observables.correlation) # The correlation is measured with respect to the first spin, i.e. the i-th entry is the correlation dot(S_1,S_i). 
   kx = collect(range(-5pi/3,5pi/3,length=N))
   ky = collect(range(-5pi/3,5pi/3,length=N))
   structurefactor = zeros(N,N)
   for i in 1:N
      for j in 1:N
         z = 0.0
         numBasis = length(mc.lattice.unitcell.basis)
         # Compute Fourier transformation at momentum (kx, ky). The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
         for k in (layer + 1):numBasis:length(mc.lattice)
            z += cos(dot((kx[i],ky[j]),getSitePosition(mc.lattice,k)[1:2])) * correlation[k,layer + 1]
         end
         structurefactor[j,i] =  numBasis * z / length(mc.lattice)
      end
   end
   
   return structurefactor 
end


function runAnneal(t0,tf,lat,thermSweeps,MeasureSweeps, coolRate, outfile=nothing)
   ts = [t0 * coolRate^t for t in -500:5000 if t0 >= t0 * coolRate^t >= tf]
   monte = nothing
   for (ind,temp) in enumerate(ts) 
      thermalizationSweeps = thermSweeps
      measurementSweeps = 0
   
      if ind == 1
            m = MonteCarlo(lat, 1/temp, thermalizationSweeps, measurementSweeps, reportInterval = MeasureSweeps, rewrite = true);
            run!(m, disableOutput = true)
      else
            if ind == length(ts)
               thermalizationSweeps = 0
               measurementSweeps = MeasureSweeps
            end
            m = MonteCarlo(monte.lattice, 1/temp, thermalizationSweeps, measurementSweeps, reportInterval = MeasureSweeps, rewrite = false);

            if ind != length(ts)
               run!(m, disableOutput = true)
            else
               run!(m, outfile = outfile)
            end  
      end
      monte = deepcopy(m);
   end
   return monte
end
