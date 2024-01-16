using TightBindingToolkit, MeanFieldToolkit
using LaTeXStrings, Plots, LinearAlgebra, YAML
using DelimitedFiles, DataFrames, JLD2
filling = 0.232
t1 = -1.0
filename = "11.27.2023_Bilayer"
cd(@__DIR__)
println(pwd())
println(@__DIR__)
params = YAML.load_file("../Input/$(filename).yml")
filename = "Bilayer_11.09.2023"

U_array = collect(LinRange(params["U_min"], params["U_max"], params["U_length"]))
U_arr = U_array#append!(U_array_1, U_array_2)
#U_var = U_array[end-1]
#loc = "/Users/ahardy/Library/CloudStorage/GoogleDrive-ahardy@flatironinstitute.org/My Drive/Skyrmion/Bilayer_SkX/TBModel/Monolayer"
loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Bilayer_Data/"
#date = "11.09.2023"
type = "_Uniform"
#############################
function Plot_Band_Data!(TBResults, labels, closed::Bool=true, nearest::Bool=true, plot_legend::Bool=true, framestyle::Symbol=:box, guidefontsize::Int64=14, tickfontsize::Int64=12,
    font::String="Helvetica", plot_title::Bool=true)
    bands = TBResults["Bands"]
    label_indices = TBResults["Labels"]
    bzpath = TBResults["BZ_Path"]
    mu = TBResults["mu"]
    plt = plot(grid=false, legend=plot_legend, bg_legend=:transparent,
        framestyle=framestyle, guidefontsize=guidefontsize, tickfontsize=tickfontsize)
    println(size(bands))
    for (i, band) in enumerate(bands[1:48])
        println(i)
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
gap_array = zeros((length(U_array), 2))
const a1 = [-3.0, sqrt(3)]
const a2 = [3.0, sqrt(3)]

const l1 = [1.0, 0]
const l2 = [-0.5, sqrt(3) / 2]

UC = UnitCell([a1, a2], 4)
order_parameter = Array{Float64}(undef, 24)
c_arr = Array{Float64}(undef, (length(U_array), 24))
c_fill = Array{Float64}(undef, (length(U_array)))


for (ind, U_var) in enumerate(U_arr)

    fileName = loc * "Last_Itr/Last_Itr_$(filename)=$(round(filling, digits=3))_U=$(round(U_var, digits=2))_t1=$(round(t1, digits=2)).jld2"
    println(fileName)
    TBResults = load(fileName) #MeanFieldToolkit.MFTResume.ReadMFT(fileName)
    println(keys(TBResults))

    gap_array[ind, 1] = U_var
    gap_array[ind, 2] = TBResults["Gap"]
    print(TBResults["Gap"])
    c_arr[ind, :] = abs.(TBResults["Chern"])
    c_fill[ind] = abs.(TBResults["Chern Fill"])
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

    plot = Plot_Band_Data!(TBResults, [L"\Gamma", L"M_2", L"M_3"])
    #plot = Plot_Band_Structure!(TBModel, [TBModel.bz.HighSymPoints["G"], TBModel.bz.HighSymPoints["M2"], TBModel.bz.HighSymPoints["M3"]]; labels=[L"\Gamma", L"M_2", L"M_3"])
    display(plot)
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
end
gap_plot = scatter(gap_array[:, 1], gap_array[:, 2], xlabel="U", ylabel="Δ")
display(gap_plot)
savefig(loc * "gap.png")

#fileName = loc * "Bilayer_Uniform_11.06.2023=$(round(filling, digits=3))_U=$(round(0, digits=2))_t1=$(round(t1, digits=2)).jld2"
#TBResults = MeanFieldToolkit.MFTResume.ReadMFT(fileName)
#TBModel = TBResults["MFT"].model
#
# for (ind, U_var) in enumerate(U_array)
#     c = readdlm(loc * "chern_$(round(U_var, digits=2)).csv")
#     c_arr[ind, :] = c
# end
x = scatter(U_arr, c_arr)
display(x)
scatter(U_arr, [abs.(c_arr[:, 1]), abs.(c_arr[:, 4])], label=["Chern ( first 2 bands)" "Chern ( first 6 bands)"])
#It's uncertain of what Chern number to use?
xlabel!("U")

C_plot = scatter(U_arr, c_fill, xlabel="U", ylabel="σ")
display(C_plot)
savefig(loc * "Chern.png")
# plot the bands color code by sign of Chern # 
# or which layer 