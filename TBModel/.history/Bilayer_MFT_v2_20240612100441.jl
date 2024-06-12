using Plots, LinearAlgebra, ColorSchemes
using MeanFieldToolkit, TightBindingToolkit, FixedPointToolkit
using YAML, LazyGrids, JLD2

loc = "/scratch/a/aparamek/andykh/Data/Bilayer_Data"
#loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Bilayer_Data"
function MFT(params, filename)
    ##Triangular Lattice 
    SkXSize = get!(params, "SkXSize", 2)
    SkX = get!(params, "SkX", "Neel")
    a1 = SkXSize / 2 * [-3.0, sqrt(3)]
    a2 = SkXSize / 2 * [3.0, sqrt(3)]
    l1 = [1.0, 0]
    l2 = [-0.5, sqrt(3) / 2]

    UC = UnitCell([a1, a2], 4)
    ##Parameters

    n = get!(params, "n", 10)
    kSize = 6 * n + 3
    t = get!(params, "t", 1.0)
    t_inter = get!(params, "t_inter", 0.0)
    jh = get!(params, "jh", 1.0)
    U = get!(params, "U", 0.0)
    t_density = get!(params, "t_density", 0.0)
    ##### Thermodynamic parameters
    filling = get!(params, "filling", 0.5)
    T = get!(params, "T", 0.0)
    t1 = -t
    tinter_param = Param(t_inter, 2)
    t1Param = Param(t1, 2)
    jhParam = Param(jh, 2)
    tdParam = Param(t_density, 2)
    tiParam = Param(t_inter, 2)
    HoppingParams = [t1Param, tdParam, tiParam, jhParam]
    su2spin = SpinMats(1 // 2)
    su4spin = SpinMats(3 // 2)
    ##Adding inner-hexagon structure  
    for j = 0:(SkXSize-1)
        for i = 0:(SkXSize*3-1)
            AddBasisSite!(UC, i .* l1 + j .* l2)
        end
    end

    ##Istrotropic bonds
    AddIsotropicBonds!(t1Param, UC, 1.0, su4spin[4], "t1", checkOffsetRange=1)
    AddIsotropicBonds!(tiParam, UC, 0.0, 2 * kron(su2spin[1], su2spin[4]), "interlayer")
    AddIsotropicBonds!(tdParam, UC, 0.0, 2 * kron(su2spin[3], su2spin[4]), "imbalance")

    ##Functions that will be useful for adding anisotropic bonds
    weiss_neel(v) = [sin(pi * (norm(v) / (SkXSize))) * v[1] / norm(v), sin(pi * (norm(v) / (SkXSize))) * v[2] / norm(v), cos(pi * (norm(v) / (SkXSize)))]
    weiss_bloch(v) = [sin(pi * (norm(v) / (SkXSize))) * v[2] / norm(v), sin(pi * (norm(v) / (SkXSize))) * -v[1] / norm(v), cos(pi * (norm(v) / (SkXSize)))]
    weiss = Dict("Neel" => weiss_neel, "Bloch" => weiss_bloch)
    sigmav(i, j) = 2 .* [su2spin[1][i, j], su2spin[2][i, j], su2spin[3][i, j]]

    s11 = sigmav(1, 1)
    s12 = sigmav(1, 2)
    s21 = sigmav(2, 1)
    s22 = sigmav(2, 2)
    intermat(s1, s2) = [dot(s1, s11) dot(s1, s12) 0 0; dot(s1, s21) dot(s1, s22) 0 0; 0 0 dot(s2, s11) dot(s2, s12); 0 0 dot(s2, s21) dot(s2, s22)]
    ##Creating BZ and Hamiltonian Model
    bz = BZ(kSize)
    FillBZ!(bz, UC)
    n_top = real.(kron([1.0 0.0; 0.0 0.0], su2spin[4]))
    n_bottom = real.(kron([0.0 0.0; 0.0 1.0], su2spin[4]))

    Hubbard = DensityToPartonCoupling(n_top, n_bottom)
    Hubbard_V_up = DensityToPartonCoupling(n_top, n_top)
    Hubbard_V_dn = DensityToPartonCoupling(n_bottom, n_bottom)

    UParam = Param(1.0, 4)
    VParam = Param(1.0, 4)
    F1Param = Param(1.0, 2)
    F2Param = Param(1.0, 2)

    AddIsotropicBonds!(UParam, UC, 0.0, Hubbard, "Hubbard Interaction")
    AddIsotropicBonds!(VParam, UC, 0.0, Hubbard_V_up + Hubbard_V_dn, "Hubbard_V Interaction", checkOffsetRange=1)

    for (ind, bas) in enumerate(UC.basis)
        closest = [bas, bas - a1, bas - a2, bas - a1 - a2, bas + a1, bas + a2, bas + a1 + a2, bas + a1 - a2, bas - a1 + a2]
        minimal = findmin(x -> norm(x), closest)[2]
        if (SkXSize - 1) < norm(closest[minimal]) < SkXSize
            mat = intermat(normalize(weiss[SkX](closest[minimal]) + weiss[SkX](-closest[minimal])), normalize(weiss[SkX](closest[minimal]) .* [1, 1, -1] + weiss[SkX](-closest[minimal]) .* [1, 1, -1]))
        else
            spn = weiss[SkX](closest[minimal])
            replace!(spn, NaN => 0.0)
            mat = intermat(normalize(spn), normalize(spn .* [1, 1, -1]))
        end
        AddAnisotropicBond!(jhParam, UC, ind, ind, [0, 0], mat, 0.0, "Hunds")
    end
    CreateUnitCell!(UC, HoppingParams)
    Density_up = []
    Density_dn = []
    UParam.value = [U]
    VParam.value = [params["V"]]
    V = params["V"]
    for (ind, bas) in enumerate(UC.basis)
        push!(Density_up, Param(1.0, 2))
        push!(Density_dn, Param(1.0, 2))
        AddAnisotropicBond!(Density_up[ind], UC, ind, ind, [0, 0], n_top, 0.0, "Dens_up-" * string(ind))
        AddAnisotropicBond!(Density_dn[ind], UC, ind, ind, [0, 0], n_bottom, 0.0, "Dens_dn-" * string(ind))

        #AddAnisotropicBond!(Sz[ind], UC, ind, ind, [0, 0], 2 * kron(su2spin[4], su2spin[3]), 0.0, "Sz-" * string(ind))

    end
    AddIsotropicBonds!(F1Param, UC, 1.0, su4spin[4], "texp")
    AddIsotropicBonds!(F2Param, UC, 1.0, 2 * kron(su2spin[1], su2spin[4]), "tdexp")

    ChiParams = vcat(F1Param, F2Param, Density_up, Density_dn)
    ChiParams = Vector{Param{2,Float64}}(ChiParams)
    H = Hamiltonian(UC, bz)
    DiagonalizeHamiltonian!(H)
    Mdl = Model(UC, bz, H; filling=filling, T=T) # Does T matter, don't I want 0 T, or is that technically impossible? 
    SolveModel!(Mdl; get_gap=true)
    mft = TightBindingMFT(Mdl, ChiParams, [UParam, VParam], [IntraQuarticToHopping, InterQuarticToHopping])
    # add filename to input 
    #fileName = loc * "/$(filename)_p=$(round(filling, digits=3))_U=$(round(U, digits=2))_t1=$(round(t1, digits=2)).jld2"
    fileName = loc * "/$(filename)_J=$(round(jh, digits=3))_U=$(round(U, digits=2)).jld2"
    #fileName = loc * "/$(filename)_V=$(round(V, digits=3))_U=$(round(U, digits=2)).jld2"
    GC.gc()
    rand_noise = rand(SkXSize^2 * 3) .- 0.5
    rand_noise = 0.01 .* (rand_noise .- sum(rand_noise) / (SkXSize^2 * 3))

    init_up = fill(filling, SkXSize^2 * 3) .+ rand_noise .- 0.02
    init_dn = fill(filling, SkXSize^2 * 3) .- rand_noise .+ 0.02
    init_guess = vcat([1.0], [0.0001], init_up, init_dn)
    if isfile(fileName)
        println("TRYING TO LOAD " * fileName)
        try
            println("SUCCESFULLY LOADED " * fileName)
            ResumeMFT!(fileName; max_iter=params["max_iter"], tol=params["tol"])#, Update=BroydenMixing)
        catch e
            println("Error Loading $fileName")
            if haskey(params, "U_prev")
                #oldfile = loc * "/$(filename)_p=$(round(filling, digits=3))_U=$(round(params["U_prev"], digits=2))_t1=$(round(t1, digits=2)).jld2"
                oldfile = loc * "/$(filename)_J=$(round(jh, digits=3))_U=$(round(params["U_prev"], digits=2)).jld2"
                #oldfile = loc * "/$(filename)_V=$(round(V, digits=3))_U=$(round(params["U_prev"], digits=2)).jld2"

                init_guess = load(oldfile)["outputs"][end] .+ vcat([0.0], [0.0], rand_noise, -1 .* rand_noise)
                SolveMFT!(mft, init_guess, fileName; max_iter=params["max_iter"], tol=params["tol"])
            else
                SolveMFT!(mft, init_guess, fileName; max_iter=params["max_iter"], tol=params["tol"])
            end
        end
    else
        if haskey(params, "U_prev")
            #oldfile = loc * "/$(filename)_p=$(round(filling, digits=3))_U=$(round(params["U_prev"], digits=2))_t1=$(round(t1, digits=2)).jld2"
            oldfile = loc * "/$(filename)_J=$(round(jh, digits=3))_U=$(round(params["U_prev"], digits=2)).jld2"
            #oldfile = loc * "/$(filename)_V=$(round(V, digits=3))_U=$(round(params["U_prev"], digits=2)).jld2"

            init_guess = load(oldfile)["outputs"][end] .+ vcat([0.0], [0.0], rand_noise, -1 .* rand_noise)
            SolveMFT!(mft, init_guess, fileName; max_iter=params["max_iter"], tol=params["tol"])
        else
            SolveMFT!(mft, init_guess, fileName; max_iter=params["max_iter"], tol=params["tol"])
        end
    end
end