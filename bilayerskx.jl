using SpinMC, LinearAlgebra, Plots


function getScalarChirality(layer, mc)
    vertex = getVertex()
    chi = 0
    if layer == 0 
        for i in vertex
            Si = collect(getSpin(mc.lattice, i[1]))
            Sj = collect(getSpin(mc.lattice, i[2]))
            Sk = collect(getSpin(mc.lattice, i[3]))            
            chi += dot(Si, cross(Sj,Sk))/((size(mc.lattice)[1])^2  )
        end
    else
        for i in vertex
            Si = collect(getSpin(mc.lattice, i[1] + 1))
            Sj = collect(getSpin(mc.lattice, i[2] + 1))
            Sk = collect(getSpin(mc.lattice, i[3] + 1))     
            chi += dot(Si, cross(Sj, Sk))/((size(mc.lattice)[1])^2)
        end    
    end

    return chi
end

function localChiralities(layer,mc)
    chiralsA = []
    chiralsB = []
    vertex = getVertex()
    chi = 0
    if layer == 0 
        for i in vertex
            Si = collect(getSpin(mc.lattice, i[1]))
            Sj = collect(getSpin(mc.lattice, i[2]))
            Sk = collect(getSpin(mc.lattice, i[3]))
            push!(chiralsA, dot(Si, cross(Sj,Sk)))
        end
        return chiralsA
    else
        for i in vertex
            Si = collect(getSpin(mc.lattice, i[1] + 1))
            Sj = collect(getSpin(mc.lattice, i[2] + 1))
            Sk = collect(getSpin(mc.lattice, i[3] + 1))     
            push!(chiralsB, dot(Si, cross(Sj,Sk)))
        end
        return chiralsB    
    end
end

function getCenters()
    centerpos = []
    for (ind, i) in enumerate(lattice.sitePositions)
        if i[3] == 0 && (i .+ a1 in lattice.sitePositions) && (findmin(x -> norm(i .+ a2 .- x),lattice.sitePositions)[2] == ind + 2)
            push!(centerpos, [i .+ a1 ./ 2 .+ (0,sqrt(3)/6,0)][1][1:2])
        end
        if i[3] == 0 && (findmin(x -> norm(i .+ a2 .- x),lattice.sitePositions)[2] == ind + 2) && (findmin(x -> norm(i .+ a2 .- a1 .- x),lattice.sitePositions)[2] == ind - size(m.lattice)[1] * 2 + 2)
            push!(centerpos, [i .+ (0,sqrt(3)/3,0)][1][1:2])
        end
    end
    return centerpos
end

function getVertex()
    vertex = []
    for (ind, i) in enumerate(lattice.sitePositions)
        if i[3] == 0 && (i .+ a1 in lattice.sitePositions) && (findmin(x -> norm(i .+ a2 .- x),lattice.sitePositions)[2] == ind + 2)
            push!(vertex, [ind, ind + size(lattice)[1], ind + 2])
        end
        if i[3] == 0 && (findmin(x -> norm(i .+ a2 .- x),lattice.sitePositions)[2] == ind + 2) && (findmin(x -> norm(i .+ a2 .- a1 .- x),lattice.sitePositions)[2] == ind - size(m.lattice)[1] * 2 + 2)
            push!(vertex, [ind, ind + 2, size(m.lattice)[1] * 2 + 2])
        end
    end
    return vertex
end

function plotSpins(lat)
    xpos,ypos = [],[]
    for i in 1:2:length(lat.sitePositions)
        push!(xpos, lat.sitePositions[i][1])
        push!(ypos, lat.sitePositions[i][2])
    end
    xspin1,yspin1 = [],[]
    xspin2,yspin2 = [],[]
    for i in 1:2:length(lat.spins[1,:])
        push!(xspin1,lat.spins[1,i])
        push!(yspin1,lat.spins[2,i])
    end
    for i in 2:2:length(lat.spins[1,:])
        push!(xspin2,lat.spins[1,i])
        push!(yspin2,lat.spins[2,i])
    end 
    p1 = quiver(xpos,ypos,quiver = (xspin1,yspin1), legend = false, aspect_ratio = :equal); 
    p2 = quiver(xpos,ypos,quiver = (xspin2,yspin2), legend = false, aspect_ratio = :equal);
    limx = (min(xlims(p1)[1],xlims(p2)[1]), max(xlims(p1)[2],xlims(p2)[2]))
    limy =  (min(ylims(p1)[1],ylims(p2)[1]), max(ylims(p1)[2],ylims(p2)[2]))
    p1 = quiver(xpos,ypos,quiver = (xspin1,yspin1), legend = false, aspect_ratio = :equal, xlims=limx,ylims=limy)
    display(p1)
    p2 = quiver(xpos,ypos,quiver = (xspin2,yspin2), legend = false, aspect_ratio = :equal, xlims=limx,ylims=limy)
    display(p2)
end

a1 = (1.0 , 0.0, 0.0)  #-
a2 = (1/2 , sqrt(3)/2, 0.0)  #/
a3 = (0.0, 0.0, 2.0)  #|
UC = UnitCell(a1 , a2, a3) 

b1 = addBasisSite!(UC, (0.0, 0.0, 0.0)) ##layer A (z = 0)
b2 = addBasisSite!(UC, (0.0, 0.0, 1.0)) ##layer B (z = 1)

J1 = 1.0
J2 = 0.0
J_ll = 0.0
A_ion = 0.0
H = 0.0
D = 0.00

M = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
S_z = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]

D1v = D .* [0, -1, 0]
D2v = D .* [sqrt(3)/2, 1/2, 0] 
D3v = D .* [-sqrt(3)/2, 1/2, 0] 
D1ex = [0 -D1v[3] D1v[2]; D1v[3] 0 -D1v[1]; -D1v[2] D1v[1] 0]
D2ex = [0 -D2v[3] D2v[2]; D2v[3] 0 -D2v[1]; -D2v[2] D2v[1] 0]
D3ex = [0 -D3v[3] D3v[2]; D3v[3] 0 -D3v[1]; -D3v[2] D3v[1] 0]

##inter-layer antiferromagnetic interaction
addInteraction!(UC, b1, b2, J_ll .* M, (0,0,0))
addInteraction!(UC, b2, b1, J_ll .* M, (0,0,1))

##Same layer interactions
for i in 1:length(UC.basis)

    ##1nd NN ferromagnetic interaction
    addInteraction!(UC, i, i, -J1 * M, (1,0,0))
    addInteraction!(UC, i, i, -J1 * M, (0,1,0))
    addInteraction!(UC, i, i, -J1 * M, (-1,1,0))
    
    #2nd NN anti-ferromagnetic interaction
    addInteraction!(UC, i, i, J2 * M, (1,1,0))
    addInteraction!(UC, i, i, J2 * M, (-1,2,0))
    addInteraction!(UC, i, i, J2 * M, (-2,1,0))

    ##DM interaction
    # addInteraction!(UC, i, i, (-1)^(i-1) * D1ex, (1,0,0))
    # addInteraction!(UC, i, i, (-1)^(i-1) * D2ex, (-1,1,0))
    # addInteraction!(UC, i, i, (-1)^(i-1) * D3ex, (0,-1,0))

    ##Magnetic field
    # setField!(UC, i, [-H,0,0])

    # ##onsite anisotropy
    # setInteractionOnsite!(UC, i, A_ion * S_z)

end


L = (24, 24, 1)
global lattice = Lattice(UC, L)

thermalizationSweeps = 25000
measurementSweeps = 75000
m = MonteCarlo(lattice, 1000.0, thermalizationSweeps, measurementSweeps, reportInterval = 50000);
run!(m)
