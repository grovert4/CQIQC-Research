# package activate skyrmion syntax
using TightBindingToolkit, MeanFieldToolkit, Statistics
using LaTeXStrings, Plots, LinearAlgebra, YAML
using DelimitedFiles, DataFrames, JLD2
cd(@__DIR__)
println(pwd())

############################
function Plot_Band_Data!(TBResults, labels, closed::Bool=true, nearest::Bool=true, plot_legend::Bool=true, framestyle::Symbol=:box, guidefontsize::Int64=14, tickfontsize::Int64=12,
    font::String="Helvetica", plot_title::Bool=true)
    bands = TBResults["Bands"]
    label_indices = TBResults["Labels"]
    bzpath = TBResults["BZ_Path"]
    mu = TBResults["mu"]
    plt = plot(grid=false, legend=plot_legend, bg_legend=:transparent,
        framestyle=framestyle, guidefontsize=guidefontsize, tickfontsize=tickfontsize)
    #println(size(bands))
    for i in 1:length(bands[1])
        #println(i)
        plot!(getindex.(bands, i), labels=L"Band : %$i", lw=2.0)
    end
    hline!([mu], linestyle=:dash, label=L"\mu", lw=0.5, linecolor=:black)
    if !closed
        xticks!(label_indices, labels)
        vline!(label_indices, linestyle=:dash, linecolor=:indigo, label="", lw=0.5)
    else
        xticks!(vcat(label_indices, [length(bzpath)]), vcat(labels, [labels[begin]]))
        vline!(label_indices, linestyle=:dash, linecolor=:indigo, label="", lw=0.5)
    end
    ylabel!("Energy", guidefontsize=guidefontsize, guidefont=font)
    if plot_title
        xlabel!("Path", guidefontsize=guidefontsize, guidefont=font)
        title!("Band Structure along path", titlefontsize=12, guidefont=font)
    end

    return plt
end

#############################

filename = "04.18.2024_Monolayer_NN"


#println(@__DIR__)
params = YAML.load_file("../Input/$(filename).yml")

U_array = collect(LinRange(params["U_min"], params["U_max"], params["U_length"]))
filling_arr = collect(LinRange(params["filling_min"], params["filling_max"], params["filling_length"])) / (params["filling_max"] * 2)
filling = filling_arr[12]
println(filling, "filling")
#U_var = U_array[end-1]
#loc = "/Users/ahardy/Library/CloudStorage/GoogleDrive-ahardy@flatironinstitute.org/My Drive/Skyrmion/Bilayer_SkX/TBModel/Monolayer"
loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Monolayer_Data/"
SkXSize = get!(params, "SkXSize", 2)
SkX = get!(params, "SkX", "Neel")
a1 = SkXSize / 2 * [-3.0, sqrt(3)]
a2 = SkXSize / 2 * [3.0, sqrt(3)]
l1 = [1.0, 0]
l2 = [-0.5, sqrt(3) / 2]
dimension = 2
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
    order_parameter = Array{Float64}(undef, (length(U_array), 1 + SkXSize^2 * 3))
end
c_arr = Array{Float64}(undef, (length(U_array), SkXSize^2 * 3 * dimension))
c_fill = Array{Float64}(undef, (length(U_array)))
gap_array = zeros((length(U_array), 2))
ord_array = Array{Float64}(undef, (length(U_array)))
eng_array = Array{Float64}(undef, (length(U_array)))

for (ind, U_var) in enumerate(U_array[:3])
    println(U_var)
    if Uniform_Status == true
        fileName = loc * "Last_Itr_$(filename)_UNIFORM_p=$(round(filling, digits=3))_U=$(round(U_var, digits=2))_t1=$(round(t1, digits=2)).jld2"
    else
        fileName = loc * "Last_Itr_$(filename)_p=$(round(filling, digits=3))_U=$(round(U_var, digits=2))_t1=$(round(t1, digits=2)).jld2"
    end
    println(fileName)
    TBResults = load(fileName) #MeanFieldToolkit.MFTResume.ReadMFT(fileName)
    #println(length(TBResults["UC"].basis))
    gap_array[ind, 1] = U_var
    gap_array[ind, 2] = TBResults["Gap"]
    #println(TBResults["MFT_Energy"])
    eng_array[ind] = TBResults["MFT_Energy"][end]
    #println(TBResults["Gap"])
    c_arr[ind, :] = abs.(TBResults["Chern"])
    c_fill[ind] = abs.(TBResults["Chern Fill"])
    ords = TBResults["Expectations"]
    order_parameter[ind, :] = TBResults["Expectations"]
    println(TBResults["Expectations"])
    println(TBResults["Chern"])
    #println(length(ords))
    ord_array[ind] = (mean(abs.(ords)))

    #plot = Plot_Band_Data!(TBResults, [L"\Gamma", L"M_2", L"M_3"])
    H = Hamiltonian(TBResults["UC"], bz)
    DiagonalizeHamiltonian!(H)
    global Mdl = Model(TBResults["UC"], bz, H; filling=0.5)
    SolveModel!(Mdl; get_gap=true)

    bands = Plot_Band_Structure!(Mdl, [bz.HighSymPoints["G"], bz.HighSymPoints["M1"], bz.HighSymPoints["K1"]], labels=["G", "M1", "K1"], plot_legend=false)
    plot!(bands, legend=false)
    display(bands)

    p = Plot_Fields!(TBResults["UC"]; use_lookup=true, site_size=1.0,
        field_thickness=1.5, range=1, field_opacity=0.9, scale=0.75,
        cmp=:viridis)
    plot!(p, legend=false)
    #display(p)
    #println(TBResults["UC"].bonds)
    display(diag(Mdl.Gr[1, 1]))
    #println(ords)
end
gap_plot = scatter(gap_array[:, 1], gap_array[:, 2], xlabel="U", ylabel=L"\Delta")
display(gap_plot)
savefig(loc * "gap.png")


x = scatter(U_array, c_arr)
display(x)
# ords = scatter(U_array, ord_array, xlabel="U", ylabel="ΔP", ylims=(0.0, 0.1))
# println(ord_array)
# display(ords)
energy_plot = scatter(U_array, eng_array, xlabel="U", ylabel="Energy")
display(energy_plot)
ords = scatter(U_array, abs.(order_parameter[:, 1]), xlabel="U", ylabel="t")
display(ords)
ords2 = scatter(U_array, order_parameter, xlabel="U", ylabel=L"\Delta S_z")
display(ords2)
# scatter(U_array, [abs.(c_arr[:, 1]), abs.(c_arr[:, 4])], label=["Chern ( first 2 bands)" "Chern ( first 6 bands)"], ylabel="C")
# #It's uncertain of what Chern number to use?
# xlabel!("U")

C_plot = scatter(U_array, c_fill, xlabel="U", ylabel=L"\sigma_xy")
display(C_plot)
savefig(loc * "Chern.png")
# plot the bands color code by sign of Chern # 
# or which layer 
#println(size(ords))
#println(ords)
# writedlm(loc * "chern_$(round(U_var, digits=2)).csv", c)
# #f string formatting ? 
# println("Convergence ", U_var)
# c = Array{Float64}(undef, 2 * length(TBModel.uc.basis))
# H = TBModel.Ham
# DiagonalizeHamiltonian!(H)
# for i in 1:2*length(TBModel.uc.basis)
#     c[i] = ChernNumber(H, collect(1:i))
#     #c[i] = ChernNumber(H, :i], check_validity:: = True)
#     # try ChernNumber(H, [i])
#     # except ChernNumber(H,[i-1,i+1])
#     # or just always do ChernNumber[4] first 
#     # first pass, make another file, uniform: 
#     println(round(c[i]), "Chern")
# end
#TBModel = TBResults["MFT"].model
#display(TBModel.gap

#display(TBResults["Convergence"]) # This is currently missing from the sort I just did, but will add it back in later.

# for i in 2:25
#     order_parameter[i-1] = TBResults["MFT"].HoppingOrders[i].value[end]
# end
# order_parameter[13] = 0
# for i in 1:12
#     norm = sqrt(order_parameter[i]^2 + order_parameter[12+i]^2)
#     println(norm)
# end
# #PlotTB.Plot_Fields!(UC, 
#display(Plot_Fields!(TBModel.uc; use_lookup=true))
#fileName = loc * "Bilayer_Uniform_11.06.2023=$(round(filling, digits=3))_U=$(round(0, digits=2))_t1=$(round(t1, digits=2)).jld2"
#TBResults = MeanFieldToolkit.MFTResume.ReadMFT(fileName)
#TBModel = TBResults["MFT"].model
#
# for (ind, U_var) in enumerate(U_array)
#     c = readdlm(loc * "chern_$(round(U_var, digits=2)).csv")
#     c_arr[ind, :] = c
# end