for (jind, jp) in enumerate(Jlls)
    # for k in 1:length(lattice.interactionMatrices)
    #     t = [deepcopy(lattice.interactionMatrices[k])...]
    #     t[1] = SpinMC.InteractionMatrix(jp .* M)
    #    global lattice.interactionMatrices[k] = tuple(t...)
    # end

    global UC.interactions[1] = (1, 2, (0,0,0), jp .* M)

    for (hind, h) in enumerate(Hs)
        global ltemp = nothing;
        m = nothing;
        println("Total done: ", c/(length(JHmat)) * 100)
        global c += 1
        # for k in 1:length(lattice.interactionField)
        #    global lattice.interactionField[k] =  (-h, 0, 0)
        # end

        # setField!(UC, 1, [-h,0,0])
        # setField!(UC, 2, [-h,0,0])

        global UC.interactionsField[1] = [0,0,-h]
        global UC.interactionsField[2] = [0,0,-h]
        
        global ltemp = Lattice(UC, L)

        # MPI.Initialized() || MPI.Init()
        # commSize = MPI.Comm_size(MPI.COMM_WORLD)
        # commRank = MPI.Comm_rank(MPI.COMM_WORLD)
        m = MonteCarlo(ltemp, 1000.0, thermalizationSweeps, measurementSweeps, reportInterval = 10000)
        run!(m)
        # N = 128
        # MPI.Finalize()

        # correlation = mean(m.observables.correlation) # The correlation is measured with respect to the first spin, i.e. the i-th entry is the correlation dot(S_1,S_i). 
        # kx = collect(range(-2pi,2pi,length=N))
        # ky = collect(range(-2pi,2pi,length=N))        
        # structurefactor = zeros(N,N)
        # for a in 1:N
        #     z = 0.0
        #     for k in 1:length(lattice)
        #         z += cos(dot((kx[a],ky[a]),(getSitePosition(lattice,k)[1], getSitePosition(lattice,k)[2]))) * correlation[k]
        #     end
        #     structurefactor[a,a] = z / length(lattice)
        # end
        global JHmat[jind, hind] = ((scal_chiral(1, m) + scal_chiral(0, m))/(size(m.lattice)[1])^2 ≈ 0.0) ? 1.0 : 0.0
        if jind == length(Jlls) && hind == length(Hs)
           global lattice = deepcopy(m.lattice) 
        end
    end
 
end


# Plots.scatter(eachrow(hcat(vals...))..., markersize = 3.5, legend = false, group = round.(vec(JHmat), sigdigits=2))



    ##Magnetic field
    # setField!(UC, i, [-H,0,0])

    # ##onsite anisotropy
    # setInteractionOnsite!(UC, i, A_ion * S_z)


        #2nd NN anti-ferromagnetic interaction
        D1v = D .* [0, -1, 0]
        D2v = D .* [sqrt(3)/2, 1/2, 0] 
        D3v = D .* [-sqrt(3)/2, 1/2, 0] 
        D1ex = [0 -D1v[3] D1v[2]; D1v[3] 0 -D1v[1]; -D1v[2] D1v[1] 0]
        D2ex = [0 -D2v[3] D2v[2]; D2v[3] 0 -D2v[1]; -D2v[2] D2v[1] 0]
        D3ex = [0 -D3v[3] D3v[2]; D3v[3] 0 -D3v[1]; -D3v[2] D3v[1] 0]
            