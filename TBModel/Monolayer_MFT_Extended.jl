using Plots, LinearAlgebra, ColorSchemes
using MeanFieldToolkit, TightBindingToolkit, FixedPointToolkit
loc = "/scratch/a/aparamek/andykh/Data/Bilayer_Data"
#loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Bilayer_Data"
function MFT(params, filename)
    ##Triangular Lattice 
    a1 = [-3 / 2, 3 * sqrt(3) / 2]
    a2 = [9.0, 0]
    l1 = [1.0, 0]
    l2 = [-0.5, sqrt(3) / 2]
    UC = UnitCell([a1, a2], 2, 2)
    ##Parameters
    n = get!(params, "n", 10)
    kSize = 6 * n + 3
    t = get!(params, "t", 1.0)
    jh = get!(params, "jh", 0.0)
    U = get!(params, "U", 0.0)
    ##### Thermodynamic parameters
    filling = get!(params, "filling", 0.5)
    T = get!(params, "T", 0.0)
    t1 = -t
    t1Param = Param(t1, 2)
    jhParam = Param(0.0, 2)
    HoppingParams = [t1Param, jhParam]
    su2spin = SpinMats(1 // 2)

    ##Adding inner-hexagon structure  
    for j = 0:2
        for i = 0:8
            AddBasisSite!(UC, i .* l1 + j .* l2)
        end
    end
    AddIsotropicBonds!(t1Param, UC, 1.0, su2spin[4], "t1", checkOffsetRange=1)
    ##Functions that will be useful for adding anisotropic bonds
    weiss1(v) = [sin(pi * (1 - norm(v) / 3)) * v[1] / norm(v), sin(pi * (1 - norm(v) / 3)) * v[2] / norm(v), cos(pi * (1 - norm(v) / 3))]
    sigmav(i, j) = 2 .* [su2spin[1][i, j], su2spin[2][i, j], su2spin[3][i, j]]
    s11 = sigmav(1, 1)
    s12 = sigmav(1, 2)
    s21 = sigmav(2, 1)
    s22 = sigmav(2, 2)

    intermat(s) = [dot(s, s11) dot(s, s12); dot(s, s21) dot(s, s22)]

    bz = BZ(kSize)
    FillBZ!(bz, UC)
    ##Adding anisotropic bonds and normalizing if needed
    # Adding MFT Parameters
    HoppingParams = [t1Param]
    n_up = [1.0 0.0; 0.0 0.0]
    n_down = [0.0 0.0; 0.0 1.0]
    Hubbard = DensityToPartonCoupling(n_up, n_down)
    UParam = Param(1.0, 4)
    Sz = []
    tParam = Param(1.0, 2)
    UParam.value = [U]
    AddIsotropicBonds!(UParam, UC, 0.0, Hubbard, "Hubbard Interaction") # Do I need to add this to all sites?
    for (ind, bas) in enumerate(UC.basis)
        if 1 < norm(bas) < 3
            mat = intermat(normalize(weiss1(bas) + weiss1(-bas)))
        else
            closest = [bas, bas - a1, bas - a2]
            spn = weiss1(closest[findmin(x -> norm(x), closest)[2]])
            replace!(spn, NaN => 0.0)
            mat = intermat(spn)
        end
        AddAnisotropicBond!(jhParam, UC, ind, ind, [0, 0], mat, 0.0, "Hunds")
    end
    CreateUnitCell!(UC, HoppingParams)
    AddIsotropicBonds!(tParam, UC, 1.0, su2spin[4], "s_H") # Am I not double counting the hopping ?? 
    for (ind, bas) in enumerate(UC.basis)
        push!(Sz, Param(1.0, 2))
        AddAnisotropicBond!(Sz[ind], UC, ind, ind, [0, 0], su2spin[3], 0.0, "Sz-" * string(ind))
        #AddAnisotropicBond!(Nd[ind], UC, ind, ind, [0, 0], n_down, 0.0, "Ndown-" * string(ind))
    end

    ChiParams = vcat(tParam, Sz)
    ChiParams = Vector{Param{2,Float64}}(ChiParams)
    ##Creating BZ and Hamiltonian Model
    bz = BZ(kSize)
    FillBZ!(bz, UC)
    path = CombinedBZPath(bz, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]]; nearest=true)
    H = Hamiltonian(UC, bz)
    DiagonalizeHamiltonian!(H)
    Mdl = Model(UC, bz, H; filling=filling, T=T) # Does T matter, don't I want 0 T, or is that technically impossible? 
    mft = TightBindingMFT(Mdl, ChiParams, [UParam], IntraQuarticToHopping)
    # add filename to input 
    fileName = loc * "/$(filename)_p=$(round(filling, digits=3))_U=$(round(U, digits=2))_t1=$(round(t1, digits=2)).jld2"
    GC.gc()
    init_guess = fill(0.01, 1 + 9 * 3)
    if isfile(fileName)
        println("TRYING TO LOAD " * fileName)
        try
            println("SUCCESFULLY LOADED " * fileName)
            ResumeMFT!(fileName; max_iter=params["max_iter"], tol=params["tol"])#, Update=BroydenMixing)
        catch e
            println("Error Loading $fileName")
            SolveMFT!(mft, init_guess, fileName; max_iter=params["max_iter"], tol=params["tol"])
        end
    else
        SolveMFT!(mft, init_guess, fileName; max_iter=params["max_iter"], tol=params["tol"])
    end
    for i in 1:2*length(UC.basis)
        c = ChernNumber(H, [i])
        println(round(c))

    end
end