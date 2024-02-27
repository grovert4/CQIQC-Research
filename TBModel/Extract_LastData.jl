using FixedPointToolkit, JLD2, MeanFieldToolkit, TightBindingToolkit, LinearAlgebra
##Script that reads in a jld2 MFT data file and extracts the last of iteration of the data
function KuboChern(Ham::Hamiltonian, bz::BZ, mu::Float64)

    Vx = conj.(permutedims.(Ham.states)) .* Ham.velocity[1] .* Ham.states
    Vy = conj.(permutedims.(Ham.states)) .* Ham.velocity[2] .* Ham.states

    chern = 0.0 + im * 0.0
    for k in eachindex(Ham.bands)
        Es = Ham.bands[k]
        vx = Vx[k]
        vy = Vy[k]

        ind = searchsortedfirst(Es, mu)
        if ind == 1 || ind == length(Es)
            continue
        else
            for i in 1:ind-1
                for j in ind:length(Es)
                    chern += (vx[i, j] * vy[j, i] - vx[j, i] * vy[i, j]) / ((Es[j] - Es[i])^2)
                end
            end
        end

    end

    b1 = [bz.basis[1]; 0.0]
    b2 = [bz.basis[2]; 0.0]
    bzUnitArea = cross(b1, b2)[3] / (4 * pi^2)

    return imag(chern) * bzUnitArea * 2 * pi / length(Ham.bands)

end
# how to make this MPI compatible ? 
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
                    dict["MFT_Energy"] = data_entry["function args"][1].MFTEnergy
                    dict["Hopping_Block"] = data_entry["function args"][1].HoppingOrders
                    dict["Expectations"] = data_entry["outputs"][end]
                    order_param = last.(getproperty.(dict["Hopping_Block"], :value))
                    dict["Order_Parameter"] = order_param
                    TBModel = data_entry["function args"][1].model
                    SolveModel!(TBModel; get_gap=true)
                    path = [TBModel.bz.HighSymPoints["G"], TBModel.bz.HighSymPoints["M2"], TBModel.bz.HighSymPoints["M3"]]
                    bzpath = CombinedBZPath(TBModel.bz, path; nearest=true, closed=true)
                    path_index = GetQIndex.(bzpath, Ref(TBModel.bz); nearest=true)
                    bands_from_index = getindex.(Ref(TBModel.Ham.bands), CartesianIndex.(Tuple.(path_index)))
                    label_indices = getindex.(findmin.([norm.(Ref(ReduceQ(x, TBModel.bz)) .- bzpath) for x in path]), 2)
                    n = 10
                    kSize = 6 * n + 3
                    TBModel.bz = BZ(kSize)
                    FillBZ!(TBModel.bz, TBModel.uc)
                    TBModel.Ham = Hamiltonian(TBModel.uc, TBModel.bz)
                    DiagonalizeHamiltonian!(TBModel.Ham)
                    c = Array{Float32}(undef, 2 * length(TBModel.uc.basis))
                    for i in 1:2*length(TBModel.uc.basis)
                        c[i] = PartialChernNumber(TBModel.Ham, i, TBModel.mu)
                        println(round(c[i]), "Chern")
                    end
                    GetVelocity!(TBModel.Ham, TBModel.bz)
                    c_fill = KuboChern(TBModel.Ham, TBModel.bz, TBModel.mu)
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
                    dict["Chern Fill"] = c_fill
                    dict["Convergence"] = norm(data_entry["inputs"][end] - data_entry["outputs"][end]) #[maximum(norm.(data_entry["outputs"][i] - data_entry["inputs"][i])) for i in 1:length(data_entry["inputs"])]#
                    # save convergences
                    save(folderpath * "/Last_Itr/Last_Itr_" * string(file), dict)
                catch e
                    println("Error Loading $file")
                    println(e)
                    rethrow(e)
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
