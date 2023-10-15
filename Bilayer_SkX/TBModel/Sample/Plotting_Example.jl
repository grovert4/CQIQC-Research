using TightBindingToolkit
using MeanFieldToolkit, LaTeXStrings
using JLD2
fillings = collect(LinRange(0.1, 0.5, 20))

H =  load("./SquareHubbard/filling="*string(fillings[end])*"_U=4.0_t1=-1.0.jld2")
TBModel = H["function args"][1].TightBindingModel

plot = Plot_Band_Structure!(TBModel , [TBModel.bz.HighSymPoints["G"], TBModel.bz.HighSymPoints["M2"], TBModel.bz.HighSymPoints["M3"]] ; labels = [L"\Gamma", L"M_2", L"M_3"] )
display(plot)
display(TBModel.gap)