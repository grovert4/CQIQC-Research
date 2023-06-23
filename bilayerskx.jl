using SpinMC, PyPlot, LinearAlgebra

function isedge(p)
    if (p .+ a1 in lattice.sitePositions) && (p .+ a2 in lattice.sitePositions)
        return false
    else
        return true
    end
end


function scal_chiral(layer)
    chi = 0
    if layer == 0 
        for i in centers
            Si = collect(getSpin(lattice, i[1]))
            Sj = collect(getSpin(lattice, i[2]))
            Sk = collect(getSpin(lattice, i[3]))     
            chi += dot(Si, cross(Sj,Sk))
        end
    else
        for i in centers
            Si = collect(getSpin(lattice, i[1] + 1))
            Sj = collect(getSpin(lattice, i[2] + 1))
            Sk = collect(getSpin(lattice, i[3] + 1))     
            chi += dot(Si, cross(Sj, Sk))
        end    
    end

    return chi
end


a1 = (1/2 , sqrt(3)/2, 0.0)  #/
a2 = (1.0 , 0.0, 0.0)  #-
a3 = (0.0, 0.0, 2.0)  #|
UC = UnitCell(a1 , a2, a3) 

b1 = addBasisSite!(UC, (0.0, 0.0, 0.0)) ##layer A (z = 0)
b2 = addBasisSite!(UC, (0.0, 0.0, 1.0)) ##layer B (z = 1)


J = 1.0
J_ll = 1.0 
A_ion = 0.2
H = 1.0
D = 0.05

M = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
S_z = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]

D1v = D .* [0, -1, 0]
D2v = D .* [sqrt(3)/2, 1/2, 0] 
D3v = D.* [-sqrt(3)/2, 1/2, 0] 
D1ex = [0 -D1v[3] D1v[2]; D1v[3] 0 -D1v[1]; -D1v[2] D1v[1] 0]
D2ex = [0 -D2v[3] D2v[2]; D2v[3] 0 -D2v[1]; -D2v[2] D2v[1] 0]
D3ex = [0 -D3v[3] D3v[2]; D3v[3] 0 -D3v[1]; -D3v[2] D3v[1] 0]

addInteraction!(UC, b1, b2, J_ll .* M, (0,0,0))

##Same layer interactions
for i in 1:length(UC.basis)
    ##onsite anisotropy
    setInteractionOnsite!(UC, i, A_ion * S_z)

    ##Magnetic field
    setField!(UC, i, [-H,0,0])

    ##2nd NN ferromagnetic interaction
    addInteraction!(UC, i, i, -J * M, (2,0,0))
    addInteraction!(UC, i, i, -J * M, (0,2,0))
    addInteraction!(UC, i, i, -J * M, (2,-2,0))
    addInteraction!(UC, i, i, -J * M, (-1,2,0))
    addInteraction!(UC, i, i, -J * M, (2,-1,0))
    addInteraction!(UC, i, i, -J * M, (1,1,0))

    ##DM interaction
    addInteraction!(UC, i, i, (-1)^(i) * D1ex, (0,1,0))
    addInteraction!(UC, i, i, (-1)^(i) * D2ex, (1,-1,0))
    addInteraction!(UC, i, i, (-1)^(i) * D3ex, (-1,0,0))
end

L = (15, 15, 1)
lattice = Lattice(UC, L)

beta_f = 0.001 * J
beta_i = 1
alpha = 0.9993 

Jlls = collect(range(0.025, 0.3, 12))
Hs = collect(range(0, 1.5, 16)) 

thermalizationSweeps = 7500
measurementSweeps = 17500

vals = collect(Iterators.product(Jlls, Hs))
JHmat = zeros(length(Jlls), length(Hs))
c = 1

centers = []
for i in lattice.sitePositions
    if i[3] == 0 && isedge(i) == false
        push!(centers, [findfirst(x -> x== i, lattice.sitePositions), findfirst(x -> x== i .+ a2, lattice.sitePositions), findfirst(x -> x== i .+ a1, lattice.sitePositions)])
    end
    if i[3] == 0 && i .+ a1 .- a2 in lattice.sitePositions
        push!(centers, [findfirst(x -> x== i, lattice.sitePositions), findfirst(x -> x== i .+ a1, lattice.sitePositions), findfirst(x -> x== i .+ a1 .- a2, lattice.sitePositions)])
    end
end



for (jind, jp) in enumerate(Jlls)
    
    for k in 1:length(lattice.interactionMatrices)
        t = [deepcopy(lattice.interactionMatrices[k])...]
        t[1] = SpinMC.InteractionMatrix(jp .* M)
        lattice.interactionMatrices[k] = tuple(t...)
    end

    for (hind, h) in enumerate(Hs)
        println("Total done: ", c/(length(JHmat)) * 100)
        global c += 1
        for k in 1:length(lattice.interactionField)
            lattice.interactionField[k] =  (-h, 0, 0)
        end

        m = MonteCarlo(lattice, (1/beta_f), thermalizationSweeps, measurementSweeps, reportInterval = 50000)
        run!(m)
        N = 128
         
#         # correlation = mean(m.observables.correlation) # The correlation is measured with respect to the first spin, i.e. the i-th entry is the correlation dot(S_1,S_i). 
#         # kx = collect(range(-2pi,2pi,length=N))
#         # ky = collect(range(-2pi,2pi,length=N))        
#         # structurefactor = zeros(N,N)
#         # for a in 1:N
#         #     z = 0.0
#         #     for k in 1:length(lattice)
#         #         z += cos(dot((kx[a],ky[a]),(getSitePosition(lattice,k)[1], getSitePosition(lattice,k)[2]))) * correlation[k]
#         #     end
#         #     structurefactor[a,a] = z / length(lattice)
#         # end
        JHmat[jind, hind] = abs(scal_chiral(1) - scal_chiral(0))
    end
end


# # Plots.scatter(eachrow(hcat(vals...))..., markersize = 0.5, legend = false, group = round.(vec(JHmat), sigdigits=3))


