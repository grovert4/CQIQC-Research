using TightBindingToolkit
using MeanFieldToolkit, LaTeXStrings
using JLD2, Plots, LinearAlgebra, YAML
using DelimitedFiles, DataFrames
filling = 0.5
t1 = -1.0
filename = "11.27.2023_Bilayer"
cd(@__DIR__)
println(pwd())
println(@__DIR__)
params = YAML.load_file("../Input/$(filename).yml")
filename = "Bilayer_11.09.2023"

U_array = collect(LinRange(params["U_min"], params["U_max"], params["U_length"]))
U_array = [8.0]
U_arr = U_array#append!(U_array_1, U_array_2)
#U_var = U_array[end-1]
#loc = "/Users/ahardy/Library/CloudStorage/GoogleDrive-ahardy@flatironinstitute.org/My Drive/Skyrmion/Bilayer_SkX/TBModel/Monolayer"
loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Bilayer_Data/"
#date = "11.09.2023"
type = "_Uniform"
gap_array = zeros((length(U_array), 2))


const a1 = [-3.0, sqrt(3)]
const a2 = [3.0, sqrt(3)]

const l1 = [1.0, 0]
const l2 = [-0.5, sqrt(3) / 2]

UC = UnitCell([a1, a2], 4)
order_parameter = Array{Float64}(undef, 24)
for (ind, U_var) in enumerate(U_arr)
    println(U_var)
    fileName = loc * "$(filename)=$(round(filling, digits=3))_U=$(round(U_var, digits=2))_t1=$(round(t1, digits=2)).jld2"
    TBResults = MeanFieldToolkit.MFTResume.ReadMFT(fileName)
    TBModel = TBResults["MFT"].model
    plot = Plot_Band_Structure!(TBModel, [TBModel.bz.HighSymPoints["G"], TBModel.bz.HighSymPoints["M2"], TBModel.bz.HighSymPoints["M3"]]; labels=[L"\Gamma", L"M_2", L"M_3"])
    display(plot)
    #display(TBModel.gap
    gap_array[ind, 1] = U_var
    gap_array[ind, 2] = TBModel.gap
    c = Array{Float64}(undef, 2 * length(TBModel.uc.basis))
    H = TBResults["MFT"].model.Ham
    DiagonalizeHamiltonian!(H)
    for i in 1:2*length(TBModel.uc.basis)
        c[i] = ChernNumber(H, collect(1:i))
        #c[i] = ChernNumber(H, :i], check_validity:: = True)
        # try ChernNumber(H, [i])
        # except ChernNumber(H,[i-1,i+1])
        # or just always do ChernNumber[4] first 
        # first pass, make another file, uniform: 
        println(round(c[i]), "Chern")
    end
    writedlm(loc * "chern_$(round(U_var, digits=2)).csv", c)
    #f string formatting ? 
    println("Convergence ", U_var)
    display(TBResults["Convergence"])

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
fileName = loc * "Bilayer_11.09.2023=$(round(filling, digits=3))_U=$(round(8, digits=2))_t1=$(round(t1, digits=2)).jld2"
TBResults = MeanFieldToolkit.MFTResume.ReadMFT(fileName)
ResumeMFT!(fileName; max_iter=200, tol=1e-6)
TBModel = TBResults["MFT"].model

SolveModel!(TBModel; get_gap=true)
c_arr = Array{Float64}(undef, (length(U_array), 2 * length(TBModel.uc.basis)))
for (ind, U_var) in enumerate(U_array)
    c = readdlm(loc * "chern_$(round(U_var, digits=2)).csv")
    c_arr[ind, :] = c
end
x = scatter(U_arr, c_arr)
display(x)
scatter(U_arr, [abs.(c_arr[:, 2]), abs.(c_arr[:, 6])], label=["Chern ( first 2 bands)" "Chern ( first 6 bands)"])
xlabel!("U")
# plot the bands color code by sign of Chern # 
# or which layer 