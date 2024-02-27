using MPI
MPI.Initialized() || MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)

using SpinMC_more_more, LinearAlgebra, JSON, LazyGrids
include("functions.jl")

inputFile = JSON.parsefile("./Input_Files/"*ARGS[1]*".json")
J1 = inputFile["J_1"]
D = inputFile["D"]
A_ion = inputFile["A_ion"]
t0 = inputFile["t_max"]
tf = inputFile["t_min"]
thermSweeps = inputFile["thermalizationSweeps"]
measureSweeps = inputFile["measurementSweeps"]
extfield = inputFile["extfield"]

#Unit Cell Construction
a1 = (1.0 , 0.0, 0.0)  #-
a2 = (-1/2 , sqrt(3)/2, 0.0)  #/
a3 = (0.0, 0.0, 2.0)  #|
UCglobal = UnitCell(a1, a2, a3)

b1 = addBasisSite!(UCglobal, (0.0, 0.0, 0.0)) ##layer A (z = 0)
b2 = addBasisSite!(UCglobal, (0.0, 0.0, 1.0)) ##layer A  (z = 1)


#Helpful Matrices
I = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
Sz = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]

D1v = D .* [1, 0, 0]
D2v = D .* [-1/2,sqrt(3)/2 , 0] 
D3v = D .* [-1/2, -sqrt(3)/2, 0] 

ExchangeD(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]


##Same layer interactions
for i in 1:length(UCglobal.basis)
    #1nd NN ferromagnetic interaction
    addInteraction!(UCglobal, i, i, -J1 * I, (1,0,0))
    addInteraction!(UCglobal, i, i, -J1 * I, (0,1,0))
    addInteraction!(UCglobal, i, i, -J1 * I, (-1,-1,0))

    ##onsite anisotropy
    setInteractionOnsite!(UCglobal, i, A_ion * Sz)

    ##DM interaction
    addInteraction!(UCglobal, i, i, (-1)^(i+1) * ExchangeD(D1v), (1,0,0))
    addInteraction!(UCglobal, i, i, (-1)^(i+1) * ExchangeD(D2v), (0,1,0)) #(0,1)
    addInteraction!(UCglobal, i, i, (-1)^(i+1) * ExchangeD(D3v), (-1,-1,0)) #
end


(Jperparr,J2arr) = ndgrid(range(inputFile["Jperp_min"],inputFile["Jperp_max"],inputFile["Jperp_length"]),range(inputFile["J2_min"],inputFile["J2_max"],inputFile["J2_length"]) )
Jperps = collect(Iterators.flatten(Jperparr))
J2s = collect(Iterators.flatten(J2arr))


L = (inputFile["System_Size"], inputFile["System_Size"], 1)
gridsize =inputFile["Jperp_length"]*inputFile["J2_length"]

elements_per_process = div(gridsize, commSize)
remainder = rem(gridsize, commSize)
 
start_index = commRank * elements_per_process + min(commRank, remainder) + 1
end_index = start_index + elements_per_process - 1 + (commRank < remainder ? 1 : 0)

# path = "/home/grovert4/scratch/Data/Monolayer_Runs_Take2/13.10.2024-Monolayer_"
# Hs = collect(range(0,1.0, 31))
# J2s2 = collect(range(0,-0.5,31))

#println(commSize, " commSize?")
for i in start_index:end_index
   jperp = round(Jperps[i],sigdigits=5)
   j2 = round(J2s[i],sigdigits=5)

   #println("Rank " , commRank , " working on h = " , h, " working on j2 = ", j2) 
   filename = "/scratch/grovert4/Data/Bilayer_Runs2/noguess/"*ARGS[1]*"_Jperp=$jperp,J2=$j2.h5"
   #println(filename)
   if isfile(filename) 
      println("Already Completed "*filename)
   else
      UClocal = deepcopy(UCglobal)
      for i in 1:length(UCglobal.basis)
         #Add J2 2NN AF interaction 
         addInteraction!(UClocal, i, i, -j2 * I, (-1,1,0))
         addInteraction!(UClocal, i, i, -j2 * I, (1,2,0))
         addInteraction!(UClocal, i, i, -j2 * I, (2,1,0))

         # setField!(UClocal, i, [0,0,(-1)^(i + 1) * jperp/4])
      end
      addInteraction!(UClocal, b1, b2, -jperp * Sz , (0,0,0))
      latticeLocal = Lattice(UClocal, L)

      # temph = findmin(abs.(jperp .+ Hs))[2] 
      # tempj = findmin(abs.(j2 .- J2s))[2]
      # file = h5open(path * "H=$(Hs[temph]),J2=$(J2s2[tempj]).h5")["mc"]
      # sites = parse.(Int64,collect(keys(read(file["lattice"]["spins"]))))
      # spins = collect(values(read(file["lattice"]["spins"])))
      # sorted = sortperm(sites)
      # sortedspins = reshape(vcat(spins[sorted]...),(3,latticeLocal.length./2))
      # latticeLocal[1:2:end] = sortedspins
      # sortedspins[3,:] = -1 .* sortedspins[3,:]
      # latticeLocal[2:2:end] = sortedspins
      mc = runAnneal(t0,tf,latticeLocal,thermSweeps,measureSweeps,inputFile["coolRate"],filename,true, extfield);
      # mc = runAnneal(t0,tf,latticeLocal,thermSweeps,measureSweeps,inputFile["coolRate"],filename,false, extfield);
   end
end
