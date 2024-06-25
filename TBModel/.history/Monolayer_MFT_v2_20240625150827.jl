using Plots, LinearAlgebra, ColorSchemes
using MeanFieldToolkit, TightBindingToolkit, FixedPointToolkit
using YAML, LazyGrids, JLD2

function MFT(params, filename)
    loc = get!(params, "loc", "/scratch/a/aparamek/andykh/Data/Monolayer_Data")
    ##Triangular Lattice
    SkXSize = get!(params, "SkXSize", 2)
    SkX = get!(params, "SkX", "Neel")
    a1 = SkXSize / 2 * [-3.0, sqrt(3)]
    a2 = SkXSize / 2 * [3.0, sqrt(3)]
    l1 = [1.0, 0]
    l2 = [-0.5, sqrt(3) / 2]
    UC = UnitCell([a1, a2], 2, 2)
    ##Parameters
    n = get!(params, "n", 10)
    kSize = 6 * n + 3
    t = get!(params, "t", 1.0)
    jh = get!(params, "jh", -1.0)
    U = get!(params, "U", 0.0)
    ##### Thermodynamic parameters
    filling = get!(params, "filling", 0.5)
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
    AddIsotropicBonds!(t1Param, UC, 1.0, su2spin[4], "t1", checkOffsetRange=1)
    ##Functions that will be useful for adding anisotropic bonds
    weiss_neel(v) = [sin(pi * (norm(v) / (SkXSize))) * v[1] / norm(v), sin(pi * (norm(v) / (SkXSize))) * v[2] / norm(v), cos(pi * (norm(v) / (SkXSize)))]
    weiss_bloch(v) = [sin(pi * (norm(v) / (SkXSize))) * v[2] / norm(v), sin(pi * (norm(v) / (SkXSize))) * -v[1] / norm(v), cos(pi * (norm(v) / (SkXSize)))]
    weiss = Dict("Neel" => weiss_neel, "Bloch" => weiss_bloch)
    sigmav(i, j) = 2 .* [su2spin[1][i, j], su2spin[2][i, j], su2spin[3][i, j]]
    s11 = sigmav(1, 1)
    s12 = sigmav(1, 2)
    s21 = sigmav(2, 1)
    s22 = sigmav(2, 2)

    intermat(s) = [dot(s, s11) dot(s, s12); dot(s, s21) dot(s, s22)]

    ##Adding anisotropic bonds and normalizing if needed

    bz = BZ(kSize)
    FillBZ!(bz, UC)
    ##Adding anisotropic bonds and normalizing if needed
    # Adding MFT Parameters
    n_up = [1.0 0.0; 0.0 0.0]
    n_down = [0.0 0.0; 0.0 1.0]
    n_tot = [1.0 0.0; 0.0 1.0]
    Hubbard = DensityToPartonCoupling(n_tot, n_tot)
    UParam = Param(1.0, 4)
    Nexp = []
    tParam = Param(1.0, 2)
    UParam.value = [U]
    AddIsotropicBonds!(UParam, UC, 1.0, Hubbard, "Hubbard Interaction", checkOffsetRange=1) # Do I need to add this to all sites?
    for (ind, bas) in enumerate(UC.basis)
        closest = [bas, bas - a1, bas - a2, bas - a1 - a2, bas + a1, bas + a2, bas + a1 + a2, bas + a1 - a2, bas - a1 + a2]
        minimal = findmin(x -> norm(x), closest)[2]
        if (SkXSize - 1) < norm(closest[minimal]) < SkXSize
            mat = intermat(normalize(weiss[SkX](closest[minimal]) + weiss[SkX](-closest[minimal])))
        else
            spn = weiss[SkX](closest[minimal])
            replace!(spn, NaN => 0.0)
            mat = intermat(normalize(spn))
        end
        AddAnisotropicBond!(jhParam, UC, ind, ind, [0, 0], mat, 0.0, "Hunds")
    end
    CreateUnitCell!(UC, HoppingParams)
    hopping = []
    for (ind,bond) in enumerate(t1Param.unitBonds)
        AddAnisotropicBonds!(hopping[], UC, bond.base, bond.target,bond.offset, su2spin[4],bond.dist, "texp-" * string(ind))


    end
    for (ind, bas) in enumerate(UC.basis)
        push!(Nexp, Param(1.0, 2))
        #AddAnisotropicBond!(Sz[ind], UC, ind, ind, [0, 0], su2spin[3], 0.0, "Sz-" * string(ind))
        AddAnisotropicBond!(Nexp[ind], UC, ind, ind, [0, 0], n_tot, 0.0, "Ntotal-" * string(ind))
    end

    ChiParams = vcat(tParam, Nexp)
    ChiParams = Vector{Param{2,Float64}}(ChiParams)
    ##Creating BZ and Hamiltonian Model
    bz = BZ(kSize)
    FillBZ!(bz, UC)
    H = Hamiltonian(UC, bz)
    DiagonalizeHamiltonian!(H)
    Mdl = Model(UC, bz, H; filling=filling, T=T) # Does T matter, don't I want 0 T, or is that technically impossible?
    SolveModel!(Mdl)
    mft = TightBindingMFT(Mdl, ChiParams, [UParam], InterQuarticToHopping)
    # add filename to input
    #fileName = loc * "/$(filename)_J=$(round(jh, digits=3))_U=$(round(U, digits=2)).jld2"
    fileName = loc * "/$(filename)_n=$(round(filling, digits=3))_U=$(round(U , digits=2)).jld2"

    GC.gc()
    rand_noise = rand(SkXSize^2 * 3) .- 0.5
    rand_noise = 0.05 .* filling .* (rand_noise .- sum(rand_noise) / (SkXSize^2 * 3))
    init_guess = vcat([1.0], fill(filling, SkXSize^2 * 3) .+ rand_noise)   # some random # which
    if isfile(fileName)
        println("TRYING TO LOAD " * fileName)
        try
            data = load(fileName)
            println("SUCCESFULLY LOADED " * fileName)
            #ResumeMFT!(fileName; max_iter=params["max_iter"], tol=params["tol"])#, Update=BroydenMixing)
        catch e
            println("Error Loading $fileName")
            if haskey(params, "U_prev")
                oldfile = loc * "/$(filename)_J=$(round(jh, digits=3))_U=$(round(params["U_prev"], digits=2)).jld2"
                oldfile = loc * "/$(filename)_n=$(round(filling, digits=3))_U=$(round(params["U_prev"] , digits=2)).jld2"

                init_guess = load(oldfile)["outputs"][end] .+ vcat([0.0001], rand_noise)
                SolveMFT!(mft, init_guess, fileName; max_iter=params["max_iter"], tol=params["tol"])
            else
                SolveMFT!(mft, init_guess, fileName; max_iter=params["max_iter"], tol=params["tol"])
            end
        end
    else
        if haskey(params, "U_prev")
            oldfile = loc * "/$(filename)_J=$(round(jh, digits=3))_U=$(round(params["U_prev"], digits=2)).jld2"
            oldfile = loc * "/$(filename)_n=$(round(filling, digits=3))_U=$(round(params["U_prev"] , digits=2)).jld2"

            init_guess = load(oldfile)["outputs"][end] .+ vcat([0.001], rand_noise)
            SolveMFT!(mft, init_guess, fileName; max_iter=params["max_iter"], tol=params["tol"])
        else
            SolveMFT!(mft, init_guess, fileName; max_iter=params["max_iter"], tol=params["tol"])
        end
    end
    print("finished run!")
end
# filename = "01.01.2024_Monolayer_NN_test"
# params = YAML.load_file("../Input/$(filename).yml")
# data = MFT(params, filename)
#data = load("/media/andrewhardy/9C33-6BBD/Skyrmion/Monolayer_Data/01.01.2024_Monolayer_NN_test_p=0.5_U=1.0_t1=-1.0.jld2")