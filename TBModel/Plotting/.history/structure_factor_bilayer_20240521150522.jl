using TightBindingToolkit, FixedPointToolkit, MeanFieldToolkit
using Plots, LaTeXStrings, LinearAlgebra
using JLD2


function SSF(values::Vector{Float64}, positions::Vector{Vector{Float64}}, k::Vector{Float64})
    phases = exp.(-im .* dot.(Ref(k), positions))
    return sum((values .- (sum(values) / length(values))) .* phases) / length(values)
end

function SSF(values::Vector{Float64}, positions::Vector{Vector{Float64}}, ks::Matrix{Vector{Float64}})

    return SSF.(Ref(values), Ref(positions), ks)
end

function plot_RS(UC::UnitCell, polarizations::Vector{Float64})

    SpinVec = SpinMats(1 // 2)
    Mz = kron(SpinVec[3], SpinVec[3])
    My = kron(SpinVec[4], SpinVec[2])
    Mx = kron(SpinVec[4], SpinVec[1])

    lookup = Lookup(UC)
    cmp = cgrad(:viridis)

    p = plot(framestyle=:box, aspect_ratio=:equal, xlabel=L"x", ylabel=L"y", grid=false)
    for (index, site) in enumerate(UC.basis)
        scatter!(Tuple(site), label="", c=:black, markersize=6.0, markeralpha=0.5)
        scatter!(Tuple(site), label="", c=:darkorange1, markersize=polarizations[index] * 10.0, markeralpha=0.75)

        mat = lookup[(index, index, [0, 0])]
        field = real.([tr(mat * Mx), tr(mat * My), tr(mat * Mz)])
        field = 0.5 * field
        field[3] = field[3] / norm(field)

        theta = round(acos(field[3]) / pi, digits=3)
        plot!([site[begin] - (field[begin] / 2), site[begin] + (field[begin] / 2)],
            [site[end] - (field[2] / 2), site[end] + (field[2] / 2)],
            arrow=true, color=cmp[(field[3]+1)/2], linewidth=2.0, label="", linealpha=0.75)
    end

    return p
end

function plot_monolayer_Sz(UC::UnitCell, polarizations::Vector{Float64})

    p = plot(framestyle=:box, aspect_ratio=:equal, xlabel=L"x", ylabel=L"y", grid=false)
    for (index, site) in enumerate(UC.basis)
        scatter!(Tuple(site), label="", c=:black, markersize=6.0, markeralpha=0.5)

        if polarizations[index] > 0.0
            scatter!(Tuple(site), label="", c=:red, markersize=polarizations[index] * 10.0, markeralpha=0.6)
        else
            scatter!(Tuple(site), label="", c=:blue, markersize=abs(polarizations[index]) * 10.0, markeralpha=0.6)
        end

    end

    return p
end

# U = 0.0
loc = "/media/andrewhardy/9C33-6BBD/Skyrmion"
outloc = "/home/andrewhardy/Documents/Graduate/Codes/Skyrmion/TBModel/"
##### loading in data files
const layer = "Bilayer"
const directory = loc * "/$(layer)_Data/"
const date = "05.10.2024"
#const date = "04.29.2024"
const date = "05.16.2024"
const date = "05.05.2024"
const date = "05.10.2024"

const n = 27#

const ps = round.(collect(range(1, n, n) / (2 * n)), digits=3)
const Us = collect(LinRange(0.0, 5.0, 11))
const t1 = -1.0

# ps = [0.5]
# Us = [5.0]

const l1 = [1.0, 0]
const l2 = [-0.5, sqrt(3) / 2]

const kSize = 6 * 5 + 3
const kSize_triangle = 6 * 30 + 3

for p in ps
    for U in Us

        filename = directory * "Last_Itr_$(date)_$(layer)_p=$(p)_U=$(U)_t1=$(t1).jld2"
        #filename = directory * "Last_Itr_$(date)_$(layer)_E_p=$(p)_U=$(U)_t1=$(t1).jld2"
        #filename = directory * "Last_Itr_$(date)_$(layer)_4_p=$(p)_U=$(U)_t1=$(t1).jld2"

        data = load(filename)
        println("data loaded" * filename)

        if layer == "Monolayer"
            polarizations = data["Expectations"][2:end]
        elseif layer == "Bilayer"
            polarizations = data["Expectations"]
        end

        UC = data["UC"]

        UC_triangle = UnitCell([l1, l2], 2, 2)
        AddBasisSite!(UC_triangle, [0.0, 0.0])

        bz = BZ([kSize, kSize])
        FillBZ!(bz, UC)

        bz_triangle = BZ([kSize_triangle, kSize_triangle])
        FillBZ!(bz_triangle, UC_triangle)

        kxs = collect(LinRange(-2 * pi, 2 * pi, 401))
        kys = collect(LinRange(-2 * pi, 2 * pi, 401))

        ks = [[kx, ky] for kx in kxs, ky in kys]
        ssf = SSF(polarizations, UC.basis, ks)

        ssf_plot = plot(framestyle=:box, aspect_ratio=:equal, xlabel=L"k_x", ylabel=L"k_y", grid=false)
        heatmap!(kxs, kys, abs.(ssf)', c=:inferno, clim=(0, maximum(abs.(ssf))))
        xlims!(-2 * pi, 2 * pi)
        ylims!(-2 * pi, 2 * pi)
        title!(L"N(k), U = %$(U), \bar{n} = %$(p)")

        symmetry_vectors = Vector{Float64}[]
        for (key, value) in bz_triangle.HighSymPoints
            if key == "M3"
                push!(symmetry_vectors, value - bz_triangle.basis[2])
            elseif key == "-M3"
                push!(symmetry_vectors, value + bz_triangle.basis[2])
            else
                push!(symmetry_vectors, value)
            end
        end

        skyrmion_vectors = [bz.basis[1], bz.basis[2], -bz.basis[1] + bz.basis[2], -bz.basis[1], -bz.basis[2], bz.basis[1] - bz.basis[2]]
        scatter!(getindex.(skyrmion_vectors, 1), getindex.(skyrmion_vectors, 2), label="skyrmion")
        scatter!(getindex.(symmetry_vectors, 1), getindex.(symmetry_vectors, 2), label="lattice")

        savefig(outloc * "/Plotting/Plots/SSF/$(date)_U=$(U)_p=$(p)_NSF_$(layer).png")

        RSPlot = plot_monolayer_Sz(UC, polarizations)
        savefig(RSPlot, outloc * "/Plotting/Plots/SSF/$(date)_U=$(U)_p=$(p)_RS_$(layer).png")
    end
end