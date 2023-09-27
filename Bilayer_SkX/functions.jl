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
    p1 = quiver(xpos,ypos,quiver = (xspin1,yspin1), legend = false, aspect_ratio = :equal, xlims=limx,ylims=limy, color=:blue2, )
    display(p1)
    p2 = quiver(xpos,ypos,quiver = (xspin2,yspin2), legend = false, aspect_ratio = :equal, xlims=limx,ylims=limy, color=:orange3)
    display(p2)
end

function addInterLayerInt(UC, strength)
    I = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    addInteraction!(UC, 1, 2, strength .* I, (0,0,0))
    addInteraction!(UC, 2, 1, strength .* I, (0,0,1))
end

function removerInterLayerInt(UC)
    int = findall(x -> x[1]!=x[2], UC.interactions)
    for i in 1:length(int)
        deleteat!(UC.interactions, int[i] - (i-1))
    end
    
end