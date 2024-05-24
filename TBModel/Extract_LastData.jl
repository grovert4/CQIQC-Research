using FixedPointToolkit, JLD2, MeanFieldToolkit, TightBindingToolkit, LinearAlgebra
using MPI
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
function SSF(values::Vector{Float64}, positions::Vector{Vector{Float64}}, k::Vector{Float64})
    phases = exp.(-im .* dot.(Ref(k), positions))
    return sum((values .- (sum(values) / length(values))) .* phases) / length(values)
end

function SSF(values::Vector{Float64}, positions::Vector{Vector{Float64}}, ks::Matrix{Vector{Float64}})

    return SSF.(Ref(values), Ref(positions), ks)
end
function N(max)

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
end


# how to make this MPI compatible ? 
function extract_data!(folderpath::String, date, layer="Bilayer", substring::String=".jld2")
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    file_list = sort(readdir(folderpath))
    file_list = sort(readdir(folderpath))
    files_per_process = length(file_list) ÷ size
    start_index = rank * files_per_process + 1
    end_index = (rank == size - 1) ? length(file_list) : start_index + files_per_process - 1

    for i in start_index:end_index
        file = file_list[i]
        if occursin(substring, string(file)) && occursin(date, string(file))
            if isfile(folderpath * "/Last_Itr/Last_Itr_" * string(file))
                println("FILE EXISTS : " * string(file))
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
                    order_param = last.(getproperty.(dict["Hopping_Block"], :value))
                    dict["Order_Parameter"] = order_param
                    TBModel = data_entry["function args"][1].model
                    n = 5
                    kSize = 6 * n + 3
                    TBModel.bz = BZ(kSize)
                    FillBZ!(TBModel.bz, TBModel.uc)
                    TBModel.Ham = Hamiltonian(TBModel.uc, TBModel.bz)
                    DiagonalizeHamiltonian!(TBModel.Ham)
                    SolveModel!(TBModel; get_gap=true)
                    idx = IsBandGapped(TBModel.Ham)
                    band_list = 1:length(idx[1, :])
                    c = Array{Float32}(undef, length(idx[1, :]))
                    for i in 1:length(idx[1, :])
                        println(.!idx[i, :])
                        c[i] = ChernNumber(TBModel.Ham, band_list[.!idx[i, :]])#, TBModel.mu)
                        #c[i] = PartialChernNumber(TBModel.Ham, i)#, TBModel.mu)
                        #println(round(c[i]), "Chern")
                        # This has an error for some reason ? 
                    end

                    GetVelocity!(TBModel.Ham, TBModel.bz)
                    c_fill = KuboChern(TBModel.Ham, TBModel.bz, TBModel.mu)
                    path = [TBModel.bz.HighSymPoints["G"], TBModel.bz.HighSymPoints["M2"], TBModel.bz.HighSymPoints["M3"]]
                    bzpath = CombinedBZPath(TBModel.bz, path; nearest=true, closed=true)
                    path_index = GetQIndex.(bzpath, Ref(TBModel.bz); nearest=true)
                    bands_from_index = getindex.(Ref(TBModel.Ham.bands), CartesianIndex.(Tuple.(path_index)))
                    label_indices = getindex.(findmin.([norm.(Ref(ReduceQ(x, TBModel.bz)) .- bzpath) for x in path]), 2)
                    dict["Bands"] = bands_from_index
                    dict["Labels"] = label_indices
                    dict["BZ_Path"] = bzpath
                    dict["UC"] = TBModel.uc
                    dict["Gap"] = TBModel.gap # I seem to have done this wrong? 
                    dict["mu"] = TBModel.mu
                    dict["Expectations"] = data_entry["outputs"][end]
                    dict["Chern"] = c
                    dict["Chern Fill"] = c_fill
                    kxs = collect(LinRange(-2 * pi, 2 * pi, 101))
                    kys = collect(LinRange(-2 * pi, 2 * pi, 101))
                    ks = [[kx, ky] for kx in kxs, ky in kys]
                    if layer == "Monolayer"
                        polarizations = dict["Expectations"][3:end]
                    elseif layer == "Bilayer"
                        len = length(dict["Expectations"][3:end])
                        polarization_up = dict["Expectations"][3:div(len, 2)+2]
                        polarization_dn = dict["Expectations"][3+div(len, 2), end]

                    end
                    dict["ssf_up"] = SSF(polarization_up, TBModel.uc.basis, ks)
                    dict["ssf_dn"] = SSF(polarization_dn, TBModel.uc.basis, ks)

                    dict["Convergence"] = norm(data_entry["inputs"][end] - data_entry["outputs"][end]) #[maximum(norm.(data_entry["outputs"][i] - data_entry["inputs"][i])) for i in 1:length(data_entry["inputs"])]#
                    # save convergences
                    save(folderpath * "/Last_Itr/Last_Itr_" * string(file), dict)
                    println("COMPLETED : " * string(file))
                catch e
                    println("Error Loading $file")
                    println(e)
                    rethrow(e)
                end
            end
        end
    end
    MPI.Finalize()
end
pwd()
extract_data!("$(ARGS[1])", "$(ARGS[2])")
