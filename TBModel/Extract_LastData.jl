using FixedPointToolkit, JLD2, MeanFieldToolkit, TightBindingToolkit, LinearAlgebra
##Script that reads in a jld2 MFT data file and extracts the last of iteration of the data

function extract_data!(folderpath::String, substring::String=".jld2")
    file_list = sort(readdir(folderpath))
    for file in file_list
        if occursin(substring, string(file))
            if isfile(folderpath * "/Last_Itr/Last_Itr_" * string(file))
                println("FILE EXISTS : " * folderpath * "/Last_Itr/Last_Itr_" * string(file))
                continue
            else
                try
                    println("TRYING TO LOAD " * folderpath * "/" * string(file))
                    data_entry = load(folderpath * "/" * string(file))
                    println("SUCCESFULLY LOADED " * folderpath * "/" * string(file))
                    dict = Dict()
                    dict["Iterations"] = data_entry["Self-consistency params"][:iter]
                    dict["MFT Energy"] = data_entry["function args"][1].MFTEnergy
                    dict["Hopping Block"] = data_entry["function args"][1].HoppingOrders
                    TBModel = data_entry["function args"][1].model
                    path = [TBModel.bz.HighSymPoints["G"], TBModel.bz.HighSymPoints["M2"], TBModel.bz.HighSymPoints["M3"]]
                    bzpath = CombinedBZPath(TBModel.bz, path; nearest=true, closed=true)
                    path_index = GetQIndex.(bzpath, Ref(TBModel.bz); nearest=true)
                    bands_from_index = getindex.(Ref(TBModel.Ham.bands), CartesianIndex.(Tuple.(path_index)))
                    label_indices = getindex.(findmin.([norm.(Ref(ReduceQ(x, TBModel.bz)) .- bzpath) for x in path]), 2)
                    #Can you save these data arrays? 
                    c = Array{Float32}(undef, 2 * length(TBModel.uc.basis))
                    for i in 1:2*length(TBModel.uc.basis)
                        c[i] = PartialChernNumber(TBModel.Ham, i, TBModel.mu)
                        println(round(c[i]), "Chern")
                        c_fill = FilledChernNumber(TBModel.Ham, TBModel.mu)
                    end
                    dict["Bands"] = bands_from_index
                    println(bands_from_index)
                    dict["Labels"] = label_indices
                    dict["BZ_Path"] = bzpath
                    dict["UC"] = TBModel.uc
                    dict["Gap"] = TBModel.gap # I seem to have done this wrong? 
                    dict["mu"] = TBModel.mu
                    dict["Outputs"] = data_entry["outputs"][end]
                    #plot = Plot_Band_Structure!(TBModel, [TBModel.bz.HighSymPoints["G"], TBModel.bz.HighSymPoints["M2"], TBModel.bz.HighSymPoints["M3"]]; labels=[L"\Gamma", L"M_2", L"M_3"])

                    dict["Chern"] = c
                    save(folderpath * "/Last_Itr/Last_Itr_" * string(file), dict)
                catch e
                    println("Error Loading $file")
                end
            end
        end
    end
end
pwd()
extract_data!("/scratch/a/aparamek/andykh/Data/Bilayer_Data")
#dict["Gr"] = data_entry["function args"][1].model.Gr
#dict["Convergence"] = [maximum(norm.(data_entry["outputs"][i] - data_entry["inputs"][i])) for i in 1:length(data_entry["inputs"])]
#dict["Pairing Block"] = Lookup(data_entry["function args"][1].PairingOrders)
