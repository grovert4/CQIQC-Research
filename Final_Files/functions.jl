using SpinMC_more_more, LinearAlgebra, HDF5

##Convention layer 0 is bottom layer (starting with index 1) and layer 1 is top layer starting with index 2 

##Getting the Skyrmions number for a given layer, lattice, and input vertices(use getVertex() below for this)
function getSkyrmionNumber(layer::Int64, lat::Lattice, vertex)
    Q = 0
    for i in vertex
        Si = collect(getSpin(lat, i[1] + layer))
        Sj = collect(getSpin(lat, i[2] + layer))
        Sk = collect(getSpin(lat, i[3] + layer))            
        Q += atan( dot(Si, cross(Sj,Sk)) / (1 + dot(Si, Sj) + dot(Sj, Sk) + dot(Sk, Si) ) )
    end 
    return Q/(2 * pi)  

end

##Get Skyrmion number locally for different vertices 
function localSkxNum(layer::Int64,lat::Lattice,vertex)
   Q = 0
   skxs = []
   for i in vertex
      Si = collect(getSpin(lat, i[1] + layer))
      Sj = collect(getSpin(lat, i[2] + layer))
      Sk = collect(getSpin(lat, i[3] + layer))            
      push!(skxs, atan( dot(Si, cross(Sj,Sk)) / (1 + dot(Si, Sj) + dot(Sj, Sk) + dot(Sk, Si) ) ) / (2 * pi))
   end 
   return skxs  

end

##Get the scalar chirality for a given layer, lattice, and input vertices(use getVertex() below for this)
function getScalarChirality(layer::Int64, lat::Lattice, vertex)
    chi = 0
    for i in vertex
        Si = collect(getSpin(lat, i[1] + layer))
        Sj = collect(getSpin(lat, i[2] + layer))
        Sk = collect(getSpin(lat, i[3] + layer))            
        chi += dot(Si, cross(Sj,Sk))/((size(lat)[1])^2)
    end    
    return chi
end

##Get the scalar chirality locally for different vertices
function localChiralities(layer::Int64, lat::Lattice, vertex)
    chirals = []
    for i in vertex
        Si = collect(getSpin(lat, i[1] + layer))
        Sj = collect(getSpin(lat, i[2] + layer))
        Sk = collect(getSpin(lat, i[3] + layer))            
        push!(chirals, dot(Si, cross(Sj,Sk)))
    end
    return chirals
end

##Get the positions of centers of the lattice of a layer
function getCenters(lat::Lattice)
   centerpos = []
   numBasis = length(lat.unitcell.basis)
   size = lat.size[1]
   if numBasis == 2
      sitea12(ind) = checkbounds(Bool, lat.sitePositions, ind + numBasis * size) ? ind + numBasis*size  : ind - size * (size - 1) * numBasis 
      sitea22(ind) = (((ind + 1) % (numBasis * size)) != 0) ? ind + numBasis : (ind - numBasis * (size - 1))
      sitea1pa22(ind) = sitea22(sitea12(ind))

      for (ind, i) in enumerate(lat.sitePositions)
         if i[3] == 0
               if (sitea12(ind) == ind + numBasis * size)  && (sitea1pa22(ind) == ind + 50)
                  push!(centerpos, collect(i .+ a1 ./ 2 .+ (0,sqrt(3)/6,0))[1:2])
               end
               if sitea22(ind) == ind + numBasis * size && (sitea1pa22(ind) == ind + 50)
                  push!(centerpos, collect(i .+ (0,sqrt(3)/3,0))[1:2])
               end
         end
      end
   elseif numBasis == 1
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

##Get the indeices lattice sites which are vertices the lattice of a layer
function getVertex(lat::Lattice)
   vertex = []
   size = lat.size[1]
   numBasis = length(lat.unitcell.basis)

   if numBasis == 2
      sitea12(ind) = checkbounds(Bool, lat.sitePositions, ind + numBasis * size) ? ind + numBasis*size  : ind - size * (size - 1) * numBasis 
      sitea22(ind) = (((ind + 1) % (numBasis * size)) != 0) ? ind + numBasis : (ind - numBasis * (size - 1))
      sitea1pa22(ind) = sitea22(sitea12(ind))

      for (ind, i) in enumerate(lat.sitePositions)   
         if i[3] == 0         
            push!(vertex, [ind, sitea12(ind), sitea1pa22(ind)])
            push!(vertex, [ind, sitea1pa22(ind), sitea22(ind)])
         end
      end
      
   elseif numBasis == 1
      sitea1(ind) = checkbounds(Bool, lat.sitePositions, ind + numBasis * size) ? ind + size  : ind - size * (size - 1)
      sitea2(ind) = (ind % size != 0) ? ind + 1  : (ind - (size - 1))
      sitea1pa2(ind) = sitea2(sitea1(ind))   
      for (ind, i) in enumerate(lat.sitePositions)
            push!(vertex, [ind, sitea1(ind), sitea1pa2(ind)])
            push!(vertex, [ind, sitea1pa2(ind), sitea2(ind)])
      end
   end
   return vertex
end

##Have to organize S(q) stuff out

##Get <S(0) . S(r)> with input lattice and layer and return the values 
function readStructureFactor(layer::Int64, lat::Lattice)
   N = 256
   correlation = getCorrelation(lat, layer + 1) 
   kx = collect(range(-5pi/3,5pi/3,length=N))
   ky = collect(range(-5pi/3,5pi/3,length=N))
   structurefactor = zeros(N,N)
   for i in 1:N
      for j in 1:N
         z = 0.0
         numBasis = length(lat.unitcell.basis)
         # Compute Fourier transformation at momentum (kx, ky). The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
         for k in (layer + 1):numBasis:length(lat)
            z += cos(dot((kx[i],ky[j]),getSitePosition(lat,k)[1:2])) * correlation[k]
         end
         structurefactor[j,i] =  numBasis * z / length(lat)
      end
   end
   
   return structurefactor 
end

##Read a component of <S(r)> HAVE TO FIX THIS STUFF
function readCorrelations(layer::Int64, lat::Lattice)
   N = 256
   kx = collect(range(-5pi/3,5pi/3,length=N))
   ky = collect(range(-5pi/3,5pi/3,length=N))
   structurefactorX = Array{ComplexF64}(undef, (N, N))
   structurefactorY = Array{ComplexF64}(undef, (N, N))
   structurefactorZ = Array{ComplexF64}(undef, (N, N))
   for i in 1:N
      for j in 1:N
         x = 0.0
         y = 0.0
         z = 0.0
         numBasis = length(lat.unitcell.basis)
         # Compute Fourier transformation at momentum (kx, ky). The real-spa  ce position of the i-th spin is obtained via getSitePosition(lattice,i). 
         for k in (layer + 1):numBasis:length(lat)
            x += cos(dot([kx[i], ky[j]], getSitePosition(lat,k)[1:2])) * lat.spins[1,k]
            y += cos(dot([kx[i], ky[j]], getSitePosition(lat,k)[1:2])) * lat.spins[2,k]
            z += cos(dot([kx[i], ky[j]], getSitePosition(lat,k)[1:2])) * lat.spins[3,k]
         end
         structurefactorX[j,i] =  numBasis * x / length(lat)
         structurefactorY[j,i] =  numBasis * y / length(lat)
         structurefactorZ[j,i] =  numBasis * z / length(lat)
      end
   end
   
   return structurefactorX  
end

##function for getting correlation with S(0)
function getCorrelation(lattice::Lattice, spin::Int = 1)
   corr = zeros(length(lattice))
   s0 = getSpin(lattice, spin)
   for i in 1:length(lattice)
       corr[i] = dot(s0, getSpin(lattice, i))
   end
   return corr
end

##Given input HDF5 file, update the lattice spins
function updateSpins!(file, lat::Lattice)
   file = h5open(file)["mc"]
   sites = parse.(Int64,collect(keys(read(file["lattice"]["spins"]))))
   spins = collect(values(read(file["lattice"]["spins"])))

   sorted = sortperm(sites)
   lat.spins = reshape(vcat(spins[sorted]...),(3,lat.length))  

   close(file)
end

##Read the energy of latttice stored in the HDF5 file
function readEnergy(file)
   file = h5open(file)["mc"]
   energy = read(file["observables"]["energyDensity"]["mean"])
   std = read(file["observables"]["energyDensity"]["error"])
   close(file)
   return [energy, std]
end

##Calculate average Sz for a given lattice
function avgSz(lat::Lattice)
   avg = 0
   for i in 1:length(lat.sitePositions)
      avg += lat.spins[3,i]
   end
   return avg/length(lat.sitePositions)
end

##Annealing function for J_Perp vs J_2 phase diagram used for simulations on cluster
function runAnneal(t0,tf,lat,thermSweeps,MeasureSweeps, coolRate, outfile=nothing, init_rewrite=true, extfiled = false)
   ts = [t0 * coolRate^t for t in -500:5000 if t0 >= t0 * coolRate^t >= tf]
   monte = nothing
   for (ind,temp) in enumerate(ts) 
      thermalizationSweeps = thermSweeps
      measurementSweeps = 0
      h = abs.(lat.unitcell.interactionsField[1][3])
        
      if ind == 1
            thermalizationSweeps = 0
            measurementSweeps = MeasureSweeps
            m = MonteCarlo(lat, 1/temp, thermalizationSweeps, measurementSweeps, reportInterval = MeasureSweeps, rewrite = init_rewrite);
            if extfield
                m.lattice.interactionField[1:2:end] = repeat([(0.0, 0.0, h)], length(m.lattice.interactionField[1:2:end]))
                m.lattice.interactionField[2:2:end] = repeat([(0.0, 0.0, -h)], length(m.lattice.interactionField[1:2:end]))
            end

            run_nompi!(m, disableOutput = true)
      else
            if (ind == length(ts)) || (ind == round(length(ts)/2)) 
               thermalizationSweeps = 0
               measurementSweeps = MeasureSweeps
            end
            m = MonteCarlo(monte.lattice, 1/temp, thermalizationSweeps, measurementSweeps, reportInterval = MeasureSweeps, rewrite = false);

            if extfield
                if length(ts)/2 > ind > length(ts)/4
                    m.lattice.interactionField[1:2:end] = repeat([(0.0, 0.0, h/2)], length(m.lattice.interactionField[1:2:end]))
                    m.lattice.interactionField[2:2:end] = repeat([(0.0, 0.0, -h/2)], length(m.lattice.interactionField[2:2:end]))
                elseif 3  * length(ts)/4 > ind >= length(ts)/2
                    m.lattice.interactionField[1:2:end] = repeat([(0.0, 0.0, h/4)], length(m.lattice.interactionField[1:2:end]))
                    m.lattice.interactionField[2:2:end] = repeat([(0.0, 0.0, -h/4)], length(m.lattice.interactionField[2:2:end]))
                elseif ind >= 3 * length(ts)/4
                    m.lattice.interactionField[1:2:end] = repeat([(0.0, 0.0, 0.0)], length(m.lattice.interactionField[1:2:end]))
                    m.lattice.interactionField[2:2:end] = repeat([(0.0, 0.0, 0.0)], length(m.lattice.interactionField[2:2:end]))
                end
            end
            
            if ind != length(ts)
               run_nompi!(m, disableOutput = true)
            else
               run_nompi!(m, outfile = outfile)
            end 

      end
      monte = deepcopy(m);
   end
   return monte
end

##Annealing function for H vs J_Perp phase diagram used for simulations on cluster
function runAnnealTWO(H, t0,tf,lat,thermSweeps,MeasureSweeps, coolRate, outfile=nothing, init_rewrite=true, extfiled = false)
   ts = [t0 * coolRate^t for t in -500:5000 if t0 >= t0 * coolRate^t >= tf]
   monte = nothing
   for (ind,temp) in enumerate(ts) 
      thermalizationSweeps = thermSweeps
      measurementSweeps = 0
      h = abs.(lat.unitcell.interactionsField[1][3])
        
      if ind == 1
            thermalizationSweeps = 0
            measurementSweeps = MeasureSweeps
            m = MonteCarlo(lat, 1/temp, thermalizationSweeps, measurementSweeps, reportInterval = MeasureSweeps, rewrite = init_rewrite);
            if extfield
                m.lattice.interactionField[1:2:end] = repeat([(0.0, 0.0, h + H)], length(m.lattice.interactionField[1:2:end]))
                m.lattice.interactionField[2:2:end] = repeat([(0.0, 0.0, -h + H)], length(m.lattice.interactionField[1:2:end]))
            end

            run_nompi!(m, disableOutput = true)
      else
            if (ind == length(ts)) || (ind == round(length(ts)/2)) 
               thermalizationSweeps = 0
               measurementSweeps = MeasureSweeps
            end
            m = MonteCarlo(monte.lattice, 1/temp, thermalizationSweeps, measurementSweeps, reportInterval = MeasureSweeps, rewrite = false);

            if extfield
                # m.lattice.interactionField[1:2:end] = repeat([(0.0, 0.0, coolRate^(ind) * h)], length(m.lattice.interactionField[1:2:end]))
                # m.lattice.interactionField[2:2:end] = repeat([(0.0, 0.0, -coolRate^(ind) * h)], length(m.lattice.interactionField[2:2:end]))
                if length(ts)/2 > ind > length(ts)/4
                    m.lattice.interactionField[1:2:end] = repeat([(0.0, 0.0, h/2 + H)], length(m.lattice.interactionField[1:2:end]))
                    m.lattice.interactionField[2:2:end] = repeat([(0.0, 0.0, -h/2 + H)], length(m.lattice.interactionField[2:2:end]))
                elseif 3  * length(ts)/4 > ind >= length(ts)/2
                    m.lattice.interactionField[1:2:end] = repeat([(0.0, 0.0, h/4 + H)], length(m.lattice.interactionField[1:2:end]))
                    m.lattice.interactionField[2:2:end] = repeat([(0.0, 0.0, -h/4 + H)], length(m.lattice.interactionField[2:2:end]))
                elseif ind >= 3 * length(ts)/4
                    m.lattice.interactionField[1:2:end] = repeat([(0.0, 0.0, H)], length(m.lattice.interactionField[1:2:end]))
                    m.lattice.interactionField[2:2:end] = repeat([(0.0, 0.0, H)], length(m.lattice.interactionField[2:2:end]))
                end
            end
            
            if ind != length(ts)
               run_nompi!(m, disableOutput = true)
            else
               run_nompi!(m, outfile = outfile)
            end 

      end
      monte = deepcopy(m);
   end
   return monte
end

##Display random lattices to see what spin configuration looks like for them
function randomLattices(N = 20, Hlim = length(Hs), J2lim = length(J2s))
   for i in 1:N
      global H = round(Hs[rand(1:Hlim)], sigdigits=5)
      global J2 = round(J2s[rand(1:J2lim)], sigdigits=5)
      plotMonoReadSpins("C:/Users/tanma/Monolayer_Runs_Take2_36x36/16.01.2024-36x36-Monolayer_H=$H,J2=$J2.h5", tempLattice, vertex)
   end
   return nothing
end