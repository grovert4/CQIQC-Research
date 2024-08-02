# package activate skyrmion syntax
using TightBindingToolkit, MeanFieldToolkit, Statistics
using LaTeXStrings, Plots, LinearAlgebra, YAML
using DelimitedFiles, DataFrames, JLD2
cd(@__DIR__)
println(pwd())

############################
function SSF(values::Vector{Float64}, positions::Vector{Vector{Float64}}, k::Vector{Float64})
    phases = exp.(-im .* dot.(Ref(k), positions))
    return sum((values .- (sum(values) / length(values))) .* phases) / length(values)
end

function SSF(values::Vector{Float64}, positions::Vector{Vector{Float64}}, ks::Matrix{Vector{Float64}})

    return SSF.(Ref(values), Ref(positions), ks)
end


function plot_RS(UC::UnitCell, polarizations::Vector{Float64})

    p = plot(framestyle=:box, aspect_ratio=:equal, xlabel=L"x", ylabel=L"y", grid=false)
    for (index, site) in enumerate(UC.basis)
        scatter!(Tuple(site), label="", c=:black, markersize=1.0, markeralpha=0.5)

        if polarizations[index] > 0.0
            scatter!(Tuple(site), label="", c=:purple, markersize=polarizations[index] * 10.0, markeralpha=0.9)
        else
            scatter!(Tuple(site), label="", c=:orange, markersize=abs(polarizations[index]) * 10.0, markeralpha=0.9)
        end

    end
    for (index, site) in enumerate(UC.basis .+ Ref(UC.primitives[1]))
        scatter!(Tuple(site), label="", c=:black, markersize=1.0, markeralpha=0.5)

        if polarizations[index] > 0.0
            scatter!(Tuple(site), label="", c=:purple, markersize=polarizations[index] * 10.0, markeralpha=0.9)
        else
            scatter!(Tuple(site), label="", c=:orange, markersize=abs(polarizations[index]) * 10.0, markeralpha=0.9)
        end

    end
    for (index, site) in enumerate(UC.basis .+ Ref(UC.primitives[2]))
        scatter!(Tuple(site), label="", c=:black, markersize=1.0, markeralpha=0.5)

        if polarizations[index] > 0.0
            scatter!(Tuple(site), label="", c=:purple, markersize=polarizations[index] * 10.0, markeralpha=0.9)
        else
            scatter!(Tuple(site), label="", c=:orange, markersize=abs(polarizations[index]) * 10.0, markeralpha=0.9)
        end

    end

    for (index, site) in enumerate(UC.basis .- Ref(UC.primitives[1]))
        scatter!(Tuple(site), label="", c=:black, markersize=1.0, markeralpha=0.5)

        if polarizations[index] > 0.0
            scatter!(Tuple(site), label="", c=:purple, markersize=polarizations[index] * 10.0, markeralpha=0.9)
        else
            scatter!(Tuple(site), label="", c=:orange, markersize=abs(polarizations[index]) * 10.0, markeralpha=0.9)
        end

    end
    for (index, site) in enumerate(UC.basis .- Ref(UC.primitives[2]))
        scatter!(Tuple(site), label="", c=:black, markersize=1.0, markeralpha=0.5)

        if polarizations[index] > 0.0
            scatter!(Tuple(site), label="", c=:purple, markersize=polarizations[index] * 10.0, markeralpha=0.9)
        else
            scatter!(Tuple(site), label="", c=:orange, markersize=abs(polarizations[index]) * 10.0, markeralpha=0.9)
        end

    end

    for (index, site) in enumerate(UC.basis .+ Ref(UC.primitives[1]) .- Ref(UC.primitives[2]))
        scatter!(Tuple(site), label="", c=:black, markersize=1.0, markeralpha=0.5)

        if polarizations[index] > 0.0
            scatter!(Tuple(site), label="", c=:purple, markersize=polarizations[index] * 10.0, markeralpha=0.9)
        else
            scatter!(Tuple(site), label="", c=:orange, markersize=abs(polarizations[index]) * 10.0, markeralpha=0.9)
        end

    end
    for (index, site) in enumerate(UC.basis .+ Ref(UC.primitives[2]) .- Ref(UC.primitives[1]))
        scatter!(Tuple(site), label="", c=:black, markersize=1.0, markeralpha=0.5)

        if polarizations[index] > 0.0
            scatter!(Tuple(site), label="", c=:purple, markersize=polarizations[index] * 10.0, markeralpha=0.9)
        else
            scatter!(Tuple(site), label="", c=:orange, markersize=abs(polarizations[index]) * 10.0, markeralpha=0.9)
        end

    end

    return p

end


#############################

filename = "04.14-Bloch.2024_Bilayer"
filename = "05.01-0.5.2024_Bilayer"
filename = "05.03-0.375.2024_Bilayer"
filename = "05.03-0.5.2024_Bilayer"
filename = "05.04-0.75.2024_Bilayer"
filename = "05.03-0.33.2024_Bilayer"
filename = "06.17-1.2024_Bilayer"
# filename = "07.04-3.2024_Bilayer"
# filename = "07.04-4.2024_Bilayer"

#filename = "05.04-0.66.2024_Bilayer"

filename = "07.09-25.2024_Bilayer"
filename = "07.31-25.2024_Bilayer"
filename = "07.30-25.2024_Bilayer"

#filename = "07.21-25.2024_Bilayer"

#println(@__DIR__)
params = YAML.load_file("../Input/$(filename).yml")

U_array = collect(LinRange(params["U_min"], params["U_max"], params["U_length"]))
filling_arr = collect(LinRange(params["filling_min"], params["filling_max"], params["filling_length"])) / (params["filling_max"] )
filling_arr = collect(LinRange(params["filling_min"], params["filling_max"], params["filling_length"])) / (params["filling_max"] )
filling_arr = (24 .+ LinRange(params["filling_min"], params["filling_max"], params["filling_length"])) / 48

V_array = collect(LinRange(params["V_min"], params["V_max"], params["V_length"]))
params["V"] = V_array[6]

#J_array = collect(LinRange(params["J_min"], params["J_max"], params["J_length"]))

#params["jh"] = J_array[4]
filling = params["filling"]
println(filling, "filling")
#U_var = U_array[end-1]
#loc = "/Users/ahardy/Library/CloudStorage/GoogleDrive-ahardy@flatironinstitute.org/My Drive/Skyrmion/Bilayer_SkX/TBModel/Monolayer"
loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Bilayer_Data/"
SkXSize = get!(params, "SkXSize", 2)
SkX = get!(params, "SkX", "Neel")
a1 = SkXSize / 2 * [-3.0, sqrt(3)]
a2 = SkXSize / 2 * [3.0, sqrt(3)]
l1 = [1.0, 0]
l2 = [-0.5, sqrt(3) / 2]
dimension = 4
UC = UnitCell([a1, a2], 2, 2)
##Parameters
n = get!(params, "n", 10)
kSize = 6 * n + 3
t = get!(params, "t", 1.0)
jh = get!(params, "jh", -1.0)
U = get!(params, "U", 0.0)
##### Thermodynamic parameters
#filling = get!(params, "filling", 0.5)
T = get!(params, "T", 0.0)
t1 = -t
t1Param = Param(t1, 2)
jhParam = Param(jh, 2)
HoppingParams = [t1Param, jhParam]
su2spin = SpinMats(1 // 2)

##Adding inner-hexagon structure
for j = 0:(SkXSize-1)
    for i = 0:(SkXSize*3-1)
        AddBasisSite!(UC, i .* l1 + j .* l2)
    end
end
Uniform_Status = false
bz = BZ(kSize)
FillBZ!(bz, UC)
if Uniform_Status
    order_parameter = Array{Float64}(undef, (length(U_array), 1))
else
    order_parameter = Array{Float64}(undef, (length(U_array), SkXSize^2 * 3))
end
c_arr = Array{Float64}(undef, (length(U_array), SkXSize^2 * 3 * dimension))
c_fill = Array{Float64}(undef, (length(U_array)))
gap_array = zeros((length(U_array), 2))
ord_array = Array{Float64}(undef, (length(U_array), 2 * SkXSize^2 * 3))
eng_array = Array{Float64}(undef, (length(U_array)))

#for (ind, U_var) in enumerate(U_array[:])
for (ind, V_var) in enumerate(V_array[:])
    U_var = U_array[1]
    params["V"] = V_var
    println(U_var)
    if Uniform_Status == true
        fileName = loc * "Last_Itr_$(filename)_UNIFORM_p=$(round(params["jh"], digits=3))_U=$(round(U_var, digits=2))_t1=$(round(t1, digits=2)).jld2"
    else
        #fileName = loc * "Last_Itr_$(filename)_J=$(round(params["jh"], digits=3))_U=$(round(U_var, digits=2)).jld2"
        fileName = loc * "Last_Itr_$(filename)_V=$(round(params["V"], digits=3))_U=$(round(U_var, digits=2)).jld2"
        #fileName = loc * "Last_Itr_$(filename)_n=$(round(filling, digits=3))_U=$(round(U_var, digits=2)).jld2"


    end
    println(fileName)
    TBResults = load(fileName) #MeanFieldToolkit.MFTResume.ReadMFT(fileName)
    #println(length(TBResults["UC"].basis))
    SkXsize = length(TBResults["UC"].basis)

    gap_array[ind, 1] = U_var
    gap_array[ind, 2] = TBResults["Bands"][1][31]-TBResults["Bands"][1][30] #TBResults["Gap"]
    #println(TBResults["MFT_Energy"])
    eng_array[ind] = TBResults["MFT_Energy"][end]
    #println(TBResults["Gap"])
    c_arr[ind, :] = abs.(TBResults["Chern"])
    c_fill[ind] = abs.(TBResults["Chern Fill"])
    len = length(TBResults["Expectations"][(end-2*SkXsize)+1:(end)])
    order_parameter[ind, :] = TBResults["Expectations"][(end-2*SkXsize)+1:(end-SkXsize)] .- TBResults["Expectations"][(end-SkXsize)+1:end]
    #println(length(ords))
    ord_array[ind, :] = TBResults["Expectations"][(end-2*SkXsize)+1:(end)]#(mean(abs.(ords)))

    #plot = Plot_Band_Data!(TBResults, [L"\Gamma", L"M_2", L"M_3"])
    H = Hamiltonian(TBResults["UC"], bz)
    DiagonalizeHamiltonian!(H)
    Mdl = Model(TBResults["UC"], bz, H; filling=filling)
    SolveModel!(Mdl; get_gap=true)

    bands = Plot_Band_Structure!(Mdl, [bz.HighSymPoints["G"], bz.HighSymPoints["M2"], bz.HighSymPoints["M3"]],[25,26,27,28,29,30,31,32], labels=["G", "M2", "M3"], plot_legend=false)
    plot!(bands, legend=false)
    display(bands)

    p = Plot_Fields!(TBResults["UC"]; use_lookup=true, site_size=1.0,
        field_thickness=1.5, range=1, field_opacity=0.9, scale=0.75,
        cmp=:viridis)
    plot!(p, legend=false)
    #display(p)
    #println(TBResults["UC"].bonds)
    #println(ords)
end
gap_plot = scatter(V_array, gap_array[:, 2], xlabel="U", ylabel="Δ")
display(gap_plot)
savefig(loc * "gap.png")

U_array = V_array
x = scatter(U_array, c_arr)
display(x)
ords = scatter(U_array, ord_array, xlabel="U", ylabel="ΔP", ylims=(0.0, 0.1))
println(ord_array)
display(ords)
energy_plot = scatter(U_array, eng_array, xlabel="U", ylabel="Energy")
display(energy_plot)
ords2 = scatter(U_array, abs.(order_parameter), xlabel="U", ylabel="ΔP")
display(ords2)
println("ords2")

scatter(U_array, [abs.(c_arr[:, 1]), abs.(c_arr[:, 4])], label=["Chern ( first 2 bands)" "Chern ( first 6 bands)"], ylabel="C")
#It's uncertain of what Chern number to use?
xlabel!("U")

C_plot = scatter(U_array, c_fill, xlabel="U", ylabel="σ(0)")
display(C_plot)
savefig(loc * "Chern.png")

kxs = collect(LinRange(-2 * pi, 2 * pi, 201))
kys = collect(LinRange(-2 * pi, 2 * pi, 201))

ks = [[kx, ky] for kx in kxs, ky in kys]
ssf = SSF(ord_array[4, 1:SkXSize^2*3], UC.basis, ks)
ssf_plot = plot(framestyle=:box, aspect_ratio=:equal, xlabel=L"k_x", ylabel=L"k_y", grid=false)
heatmap!(kxs, kys, abs.(ssf)', c=:inferno, clim=(0, maximum(abs.(ssf))))
xlims!(-2 * pi, 2 * pi)
ylims!(-2 * pi, 2 * pi)
title!(L"N(k), U = %$(U), \bar{n} = %$(filling)")

symmetry_vectors = Vector{Float64}[]
for (key, value) in bz.HighSymPoints
    if key == "M3"
        push!(symmetry_vectors, value - bz.basis[2])
    elseif key == "-M3"
        push!(symmetry_vectors, value + bz.basis[2])
    else
        push!(symmetry_vectors, value)
    end
end

skyrmion_vectors = [bz.basis[1], bz.basis[2], -bz.basis[1] + bz.basis[2], -bz.basis[1], -bz.basis[2], bz.basis[1] - bz.basis[2]]
scatter!(getindex.(skyrmion_vectors, 1), getindex.(skyrmion_vectors, 2), label="skyrmion")
scatter!(getindex.(symmetry_vectors, 1), getindex.(symmetry_vectors, 2), label="lattice")
display(ssf_plot)

RSPlot = plot_RS(UC, 20 * ord_array[20, 1:SkXSize^2*3] .- 20 * ord_array[20, SkXSize^2*3+1:SkXSize^2*6].-0.521)
display(RSPlot)
data = 1 * ord_array[15, 1:SkXSize^2*3] .- 1 * ord_array[15, SkXSize^2*3+1:SkXSize^2*6]
# RSPlot = plot_RS(UC, ord_array[1, SkXSize^2*3:SkXSize^2*3*2])
# display(RSPlot)
RSPlot = plot_RS(UC, 10 * ord_array[4, 1:SkXSize^2*3] .- 10 * ord_array[4, SkXSize^2*3+1:SkXSize^2*6])
display(RSPlot)
RSPlot = plot_RS(UC, order_parameter[4,:])
display(RSPlot)
savefig(RSPlot,"file.png")