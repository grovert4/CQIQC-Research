using NPZ
using LinearAlgebra
using JLD2
using Plots
using LaTeXStrings

const t1 = 1.0
const B = 1.0
const beta = 10.0
mu = -4.0
fileName = "./RPA_data/RPA/t1=$(t1)_B=$(B)_beta=$(beta)_mu=$(mu)_suscep.npz"
data = npzread(fileName)

function dress_primitives(data::Dict)
    primitives = data["primitives"]
    primitives = Vector{eltype(data["primitives"])}[eachrow(data["primitives"])...]
    primitives[1] = [primitives[1]; 0.0]
    primitives[2] = [primitives[2]; 0.0]
    push!(primitives, [0.0, 0.0, 1.0])

    return primitives
end

function dress_reciprocal(data::Dict)
    reciprocal = data["reciprocal"]
    reciprocal = Vector{eltype(data["reciprocal"])}[eachrow(data["reciprocal"])...]
    return reciprocal
end

function dress_chis(data::Dict)
    chi = data["chiDD"]

    chis = Matrix{ComplexF64}[]

    for i in 1:size(chi)[1]

        push!(chis, chi[i, :, :])
    end

    return chis
end

function get_reciprocal_ks(data::Dict)
    primitives = dress_primitives(data)
    ks = Vector{eltype(data["ks"])}[eachrow(data["ks"])...]

    k1s = dot.(ks, Ref(primitives[1]))
    k2s = dot.(ks, Ref(primitives[2]))
    k3s = dot.(ks, Ref(primitives[3]))

    Ks = hcat(k1s, k2s, k3s)
    Ks = Vector{eltype(Ks)}[eachrow(Ks)...]

    return Ks
end



function U_NN(U::Float64, k::Vector{Float64} ; primitives = dress_primitives(data))

    ##### Nearest neighbours within unit cell : 21
    NN_00 = [(1, 2), (1, 7), (1, 8), (2, 3), (2, 8), (2, 9), (3, 4), (3, 9), (3, 10),
    (4, 5), (4, 10), (4, 11), (5, 6), (5, 11), (5, 12),
    (6, 12),
    (7, 8), (8, 9), (9, 10), (10, 11), (11, 12)]
    ##### between unit cell at (0, 1) : 11
    NN_01 = [(7, 3), (7, 4), (8, 4), (8, 5), (9, 5), (9, 6), (10, 6), (10, 1), (11, 1), (11, 2), (12, 2)]
    ##### between unit cell at (1, 0) : 3
    NN_10 = [(6, 1), (6, 7), (12, 7)]
    ##### between unit cell at (1, 1) : 1
    NN_11 = [(12, 3)]
    ##### all bonds
    bonds = Dict((0, 0) => NN_00, (0, 1) => NN_01, (1, 0) => NN_10, (1, 1) => NN_11)

    mat = zeros(ComplexF64, 12, 12)
    for (type, neighbours) in bonds
        for (i, j) in neighbours
            mat[i, j] += U * exp( im * dot(k, type[1] * primitives[1] + type[2] * primitives[2]))
            mat[j, i] += U * exp(-im * dot(k, type[1] * primitives[1] + type[2] * primitives[2]))
        end
    end

    return mat
end

function RPAeigs(chis::Vector{Matrix{ComplexF64}}, Us::Vector{Matrix{ComplexF64}})
    mats = inv.(Ref(I) .+ chis .* Us) .* chis
    eigs = eigen.(mats)
    values = getproperty.(eigs, :values)

    vectors = getproperty.(eigs, :vectors)
    return values, vectors
end

function minima(values::Vector{Vector{ComplexF64}}, vectors::Vector{Matrix{ComplexF64}})
    minEigs = getindex.(values, 1)
    value, index = findmin(real.(minEigs))
    return index, value, vectors[index][:, 1]
end

function maxima(values::Vector{Vector{ComplexF64}}, vectors::Vector{Matrix{ComplexF64}})
    minEigs = getindex.(values, 6)
    value, index = findmax(real.(minEigs))
    return index, value, vectors[index][:, 6]
end


function suscep_rpa(data::Dict, U::Float64)
    chis = dress_chis(data)
    ks = Vector{eltype(data["ks"])}[eachrow(data["ks"])...]
    Us = U_NN.(Ref(U), ks)
    rpa = RPAeigs(chis, Us)
    return rpa..., ks
end

function minima(data::Dict, U::Float64)
    values, vectors, ks = suscep_rpa(data, U)
    index, value, vector = minima(values, vectors)
    return Dict("minima index" => index, "lowest eigenvalue" => value,
        "lowest momenta" => ks[index], "lowest eigenvector" => vector)
end

function maxima(data::Dict, U::Float64)
    values, vectors, ks = suscep_rpa(data, U)
    index, value, vector = maxima(values, vectors)
    return Dict("maxima index" => index, "peak eigenvalue" => value,
        "peak momenta" => ks[index], "peak eigenvector" => vector)
end


function find_peak(data::Dict ; steps::Int = 21)

    lower = 0.0
    upper = 2.0
    current = Float64[]
    check = nothing

    for _ in 1:steps
        push!(current, (upper + lower) / 2)
        U = current[end]
        check = minima(data, U)

        if check["lowest eigenvalue"] < 0.0
            upper = current[end]
        else
            lower = current[end]
        end
    end

    if check["lowest eigenvalue"] < 0.0
        U = lower
        check = minima(data, U)
        maxim = maxima(data, U)
        return Dict("U" => lower, check..., maxim...)
    else
        maxim = maxima(data, U)
        return Dict("U" => current[end], check..., maxim...)
    end
end
