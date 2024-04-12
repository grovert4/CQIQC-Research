using NPZ
using LinearAlgebra
using JLD2
using Plots
using LaTeXStrings

fileName = "t1=1.0_t3=0.2_beta=10_suscep.npz"
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



function J_NN(Jmat::Matrix{Float64}, k::Vector{Float64})

    int = zeros(ComplexF64, 6, 6)
    int[1:3, 4:6] = Jmat * (1 + exp(-1.0im * k[1])+ exp(-1.0im * k[2]))
    int[4:6, 1:3] = Jmat' * (1 + exp(1.0im * k[1])+ exp(1.0im * k[2]))

    return int
end

function J_3NN(Jmat::Matrix{Float64}, k::Vector{Float64})

    int = zeros(ComplexF64, 6, 6)
    int[1:3, 4:6] = Jmat * (exp(1.0im * (k[1] - k[2])) + exp(-1.0im * (k[1] - k[2])) + exp(-1.0im * (k[1] + k[2])))
    int[4:6, 1:3] = Jmat' * (exp(-1.0im * (k[1] - k[2])) + exp(1.0im * (k[1] - k[2])) + exp(1.0im * (k[1] + k[2])))

    return int
end

function JMats(J1::Float64, J3::Float64, lambda::Float64, ks::Vector{Vector{Float64}})

    Jmat = diagm([1.0, 1.0, lambda])
    J1s = J_NN.(Ref(Jmat * J1), ks)
    J3s = J_3NN.(Ref(Jmat * J3), ks)

    return J1s+J3s
end

function RPAeigs(chis::Vector{Matrix{ComplexF64}}, Js::Vector{Matrix{ComplexF64}})
    mats = inv.(Ref(I) .+ chis .* Js) .* chis
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


function suscep_rpa(data::Dict, J1::Float64, J3::Float64, lambda::Float64)
    chis = dress_chis(data)
    ks = get_reciprocal_ks(data)
    Js = JMats(J1, J3, lambda, ks)
    rpa = RPAeigs(chis, Js)
    return rpa..., ks
end

function minima(data::Dict, J1::Float64, J3::Float64, lambda::Float64)
    values, vectors, ks = suscep_rpa(data, J1, J3, lambda)
    index, value, vector = minima(values, vectors)
    return Dict("minima index" => index, "lowest eigenvalue" => value,
        "lowest momenta" => ks[index], "eigenvector" => vector)
end

function maxima(data::Dict, J1::Float64, J3::Float64, lambda::Float64)
    values, vectors, ks = suscep_rpa(data, J1, J3, lambda)
    index, value, vector = maxima(values, vectors)
    return Dict("maxima index" => index, "peak eigenvalue" => value,
        "peak momenta" => ks[index], "eigenvector" => vector)
end

function phase_diagram(J1s::Vector{Float64}, J3s::Vector{Float64}, lambda::Float64, data::Dict ;
    title_prefix::String = fileName[1:end-11])

    Qs = Vector{Float64}[]
    minimums = Float64[]
    maximums = Float64[]
    eigenvectors = Vector{ComplexF64}[]

    for (i, J1) in enumerate(J1s)
        for (j, J3) in enumerate(J3s)
            println("J1 = $(J1), J3 = $(J3)")
            values, vectors, ks = suscep_rpa(data, J1, J3, lambda)
            min_index, min_value, min_vector = minima(values, vectors)
            max_index, max_value, max_vector = maxima(values, vectors)

            push!(Qs, ks[max_index])
            push!(minimums, min_value)
            push!(maximums, max_value)
            push!(eigenvectors, max_vector)
        end
    end

    ns = length(J1s)
    @assert length(J1s) == length(J3s)
    Qs = reshape(Qs, ns, ns)
    minimums = reshape(minimums, ns, ns)
    maximums = reshape(maximums, ns, ns)
    eigenvectors = reshape(eigenvectors, ns, ns)

    output = Dict("Qs" => Qs, "minimums" => minimums, "maximums" => maximums, "eigenvectors" => eigenvectors,
    "J1s" => J1s, "J3s" => J3s, "lambda" => lambda)

    save(title_prefix * "_lambda=$(lambda)_suscepRPA_data.jld2", output)
    return output
end

function plot_phase_diagram(output::Dict)
    Qmods = norm.(output["Qs"]) ./ Ref(2*pi)
    ssb = (Ref(1) .+ sign.(output["minimums"])) ./ Ref(2)

    p = plot(framestyle=:box, grid=false, aspect_ratio=:equal)
    heatmap!(output["J1s"], output["J3s"], ssb .* Qmods, c=:inferno, xlabel="J1", ylabel="J3")
    hline!([0.0], c=:white, l=:dash, lw=2, label = "")
    vline!([0.0], c=:white, l=:dash, lw=2, label = "")
    return p
end


function J_polar(J::Float64, theta::Float64)
    J1 = J * cos(theta)
    J3 = J * sin(theta)
    return J1, J3
end

function find_peak(theta::Float64, lambda::Float64, data::Dict ; steps::Int = 21)

    lower = 0.0
    upper = 1.0
    current = Float64[]
    check = nothing

    for _ in 1:steps
        push!(current, (upper + lower) / 2)
        J1, J3 = J_polar(current[end], theta)
        check = minima(data, J1, J3, lambda)

        if check["lowest eigenvalue"] < 0.0
            upper = current[end]
        else
            lower = current[end]
        end
    end

    if check["lowest eigenvalue"] < 0.0
        J1, J3 = J_polar(lower, theta)
        check = minima(data, J1, J3, lambda)
        return Dict("J" => lower, check...)
    else
        return Dict("J" => current[end], check...)
    end
end


const ns = 101
const thetas = collect(LinRange(0.0, 2*pi, ns))
const lambda = 1.0

function find_peak(thetas::Vector{Float64}, lambda::Float64, data::Dict ;
    steps::Int = 21, title_prefix::String = fileName[1:end-11])

    Js = Float64[]
    Qs = Vector{Float64}[]
    minimums = Float64[]
    maximums = Float64[]
    eigenvectors = Vector{ComplexF64}[]

    for theta in thetas
        println("Working on theta = $(theta)...")

        peak = find_peak(theta, lambda, data ; steps = steps)
        push!(Js, peak["J"])
        push!(minimums, peak["lowest eigenvalue"])

        maxim = maxima(data, J_polar(peak["J"], theta)..., lambda)
        push!(maximums, maxim["peak eigenvalue"])
        push!(eigenvectors, maxim["eigenvector"])
        push!(Qs, maxim["peak momenta"])
    end

    output = Dict("Js" => Js, "thetas" => thetas, "lambda" => lambda,
    "Qs" => Qs,
    "minimums" => minimums, "maximums" => maximums,
    "eigenvectors" => eigenvectors)

    save(title_prefix * "_lambda=$(lambda)_suscepRPA_peaks.jld2", output)
    return output
end


function eigenvector_classification(output::Dict)

    ferroX = [1/sqrt(2), 0.0, 0.0, 1/sqrt(2), 0.0, 0.0]
    ferroY = [0.0, 1/sqrt(2), 0.0, 0.0, 1/sqrt(2), 0.0]
    ferroZ = [0.0, 0.0, 1/sqrt(2), 0.0, 0.0, 1/sqrt(2)]
    NeelX = [1/sqrt(2), 0.0, 0.0, -1/sqrt(2), 0.0, 0.0]
    NeelY = [0.0, 1/sqrt(2), 0.0, 0.0, -1/sqrt(2), 0.0]
    NeelZ = [0.0, 0.0, 1/sqrt(2), 0.0, 0.0, -1/sqrt(2)]

    ferroX_weight = dot.(output["eigenvectors"], Ref(ferroX))
    ferroY_weight = dot.(output["eigenvectors"], Ref(ferroY))
    ferroZ_weight = dot.(output["eigenvectors"], Ref(ferroZ))
    NeelX_weight = dot.(output["eigenvectors"], Ref(NeelX))
    NeelY_weight = dot.(output["eigenvectors"], Ref(NeelY))
    NeelZ_weight = dot.(output["eigenvectors"], Ref(NeelZ))

    ferroXY_weight = sqrt.(abs.(ferroX_weight) .^ 2 .+ abs.(ferroY_weight) .^ 2)
    NeelXY_weight = sqrt.(abs.(NeelX_weight) .^ 2 .+ abs.(NeelY_weight) .^ 2)

    return Dict("ferro XY" => ferroXY_weight, "ferro Z" => abs.(ferroZ_weight),
    "Neel XY" => NeelXY_weight, "Neel Z" => abs.(NeelZ_weight))
end


function plot_classification(output::Dict, classification::Dict ; spin_rotation::Bool=false)

    p = plot(framestyle=:box, grid=false, aspect_ratio=:equal)

    if !spin_rotation
        plot!(output["thetas"], classification["ferro XY"], proj=:polar, gridalpha = 0.75, lw=2.0, c=:blue, label="xy FM", marker=:circle)
        plot!(output["thetas"], classification["Neel XY"], proj=:polar, gridalpha = 0.75, lw=2.0, c=:red, label="xy AFM", marker=:square)
        plot!(output["thetas"], classification["ferro Z"], proj=:polar, gridalpha = 0.75, lw=2.0, c=:blue, label="z FM", marker=:circle)
        plot!(output["thetas"], classification["Neel Z"], proj=:polar, gridalpha = 0.75, lw=2.0, c=:red, label="Neel", marker=:square)
    else
        ferros = sqrt.(abs.(classification["ferro XY"]) .^ 2 .+ abs.(classification["ferro Z"]) .^ 2)
        Neels = sqrt.(abs.(classification["Neel XY"]) .^ 2 .+ abs.(classification["Neel Z"]) .^ 2)
        plot!(output["thetas"], ferros, proj=:polar, gridalpha = 0.75, lw=2.0, c=:blue, label="FM", marker=:circle)
        plot!(output["thetas"], Neels, proj=:polar, gridalpha = 0.75, lw=2.0, c=:red, label="AFM", marker=:square)
    end
    return p
end
