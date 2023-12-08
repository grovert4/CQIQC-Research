using FixedPointToolkit,JLD2,MeanFieldToolkit
##Script that reads in a jld2 MFT data file and extracts the last of iteration of the data

function extract_data!(folderpath::String, substring::String=".jld2")
    file_list = sort(readdir(folderpath))
    for file in file_list
        if occursin(substring,string(file))
        try
            data_entry = load(folderpath * "/" * string(file))
            jldopen(folderpath * "/LastItr" *  "/" * "Last_Itr_" * file , "w" ) do dict 
                dict["Iterations"] = data_entry["Self-consistency params"][:iter]
                dict["Convergence"] = [maximum(norm.(data_entry["outputs"][i] - data_entry["inputs"][i])) for i in 1:length(data_entry["inputs"])]
                dict["MFT Energy"] = data_entry["function args"][1].MFTEnergy
                dict["Pairing Block"] = Lookup(data_entry["function args"][1].PairingOrders)
                dict["Hopping Block"] = Lookup(data_entry["function args"][1].HoppingOrders)
                dict["Hopping UC"] = data_entry["function args"][1].model.uc_hop
                dict["Gr"] = data_entry["function args"][1].model.Gr
                dict["outputs"] = data_entry["outputs"][end]
            end
        catch e
        println("Error Loading $file")
        end
        end
    end
end
pwd()
extract_data!(pwd())
