using Plots, LinearAlgebra, ColorSchemes
using MeanFieldToolkit, TightBindingToolkit, FixedPointToolkit
loc = "/scratch/a/aparamek/andykh/Data/Bilayer_Data"
#loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Bilayer_Data"
function MFT(params, filename)
    ##Triangular Lattice 

    a1 = [-3.0, sqrt(3)]
    a2 = [3.0, sqrt(3)]

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
    SpinVec = SpinMats(1 // 2)
    ##### Thermodynamic parameters
    filling = get!(params, "filling", 0.5)
    T = get!(params, "T", 0.0)
    tinter_param = Param(t_inter, 2)
    t1 = -t
    t1Param = Param(t1, 2)
    jhParam = Param(jh, 2)
    tdParam = Param(t_density, 2)
    tiParam = Param(t_inter, 2)
    HoppingParams = [t1Param, tdParam, tiParam, jhParam]

    ##Adding inner-hexagon structure  
    for j = 1:2
        for i = -1:4
            AddBasisSite!(UC, i .* l1 + j .* l2)
        end
    end

    ##Istrotropic bonds
    t1 = -t
    t1Param = Param(t1, 2)
    jhParam = Param(-jh, 2)
    HoppingParams = [t1Param]

    AddIsotropicBonds!(t1Param, UC, 1.0, SpinVec[4], "t1")

    ##Functions that will be useful for adding anisotropic bonds
    weiss1(v) = [sin(pi * (1 - norm(v) / 2)) * v[1] / norm(v), sin(pi * (1 - norm(v) / 2)) * v[2] / norm(v), cos(pi * (1 - norm(v) / 2))]
    weiss2(v) = [0, 0, 1]
    #weiss(v) = ep .* weiss1(v) .+ (1 - ep) .* weiss2(v)
    weiss(v) = weiss1(v)
    sigmav(i, j) = 2 .* [SpinVec[1][i, j], SpinVec[2][i, j], SpinVec[3][i, j]]
    s11 = sigmav(1, 1)
    s12 = sigmav(1, 2)
    s21 = sigmav(2, 1)
    s22 = sigmav(2, 2)

    intermat(s) = [dot(s, s11) dot(s, s12); dot(s, s21) dot(s, s22)]

    bz = BZ(kSize)
    FillBZ!(bz, UC)
    path = CombinedBZPath(bz, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]]; nearest=true)

    ##Adding anisotropic bonds and normalizing if needed
    CreateUnitCell!(UC, HoppingParams)


    # Adding MFT Parameters
    HoppingParams = [t1Param]
    n_up = [1.0 0.0; 0.0 0.0]
    n_down = [0.0 0.0; 0.0 1.0]
    Hubbard = DensityToPartonCoupling(n_up, n_down)
    UParam = Param(1.0, 4)
    AddIsotropicBonds!(UParam, UC, 0.0, Hubbard, "Hubbard Interaction") # Do I need to add this to all sites?
    t_s = Param(1.0, 2)
    AddIsotropicBonds!(t_s, UC, 1.0, SpinVec[4], "s Hopping") # Am I not double counting the hopping ?? 
    Dx = []
    Dy = []
    #Dz = []
    UParam.value = [U]

    for (ind, bas) in enumerate(UC.basis)
        if 1 < norm(bas) < 2
            mat = intermat(normalize(weiss(bas) + weiss(-bas)))
        else
            closest = [bas, bas - a1, bas - a2]
            spn = weiss(closest[findmin(x -> norm(x), closest)[2]])
            replace!(spn, NaN => 0.0)
            mat = intermat(spn)
        end
        push!(Dx, Param(1.0, 2))
        push!(Dy, Param(1.0, 2))
        #push!(Dz, Param(1.0, 2))
        AddAnisotropicBond!(jhParam, UC, ind, ind, [0, 0], mat, 0.0, "interaction")
        AddAnisotropicBond!(Dx[ind], UC, ind, ind, [0, 0], SpinVec[1], 0.0, "Sx-" * string(ind))
        if ind > 1
            AddAnisotropicBond!(Dy[ind], UC, ind, ind, [0, 0], SpinVec[2], 0.0, "Sy-" * string(ind))
            #AddAnisotropicBond!(Dz[ind], UC, ind, ind, [0, 0], SpinVec[3], 0.0, "Sz-" * string(ind))
            # on one site only Sx 
        end

    end
    ChiParams = vcat(t_s, Dx, Dy)
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
    init_guess = fill(0.25, 12)
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