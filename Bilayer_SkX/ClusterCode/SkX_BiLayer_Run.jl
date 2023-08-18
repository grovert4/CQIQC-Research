using Pkg
Pkg.add(url="https://github.com/grovert4/SpinMC_more_more.jl")
using SpinMC_more_more, LinearAlgebra, Plots
include("functions.jl")

using MPI
MPI.Initialized() || MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)

J1 = inputFile["J_1"]
D = inputFile["D"]
A_ion = inputFile["A_ion"]
t0 = inputFile["t_max"]
tf = inputFile["t_min"]
thermSweeps = inputFile["thermalizationSweeps"]
measureSweeps = inputFile["measurementSweeps"]
J_ll = inputFile["J_perp"]
cores = commSize

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

D1ex = ExchangeD(D1v)
D2ex = ExchangeD(D2v)
D3ex = ExchangeD(D3v)

addInteraction!(UCglobal, b1, b2, -J_ll * Sz , (0,0,0))

##Same layer interactions
for i in 1:length(UCglobal.basis)
    #1nd NN ferromagnetic interaction
    addInteraction!(UCglobal, i, i, -J1 * I, (1,0,0))
    addInteraction!(UCglobal, i, i, -J1 * I, (0,1,0))
    addInteraction!(UCglobal, i, i, -J1 * I, (-1,-1,0))


    ##onsite anisotropy
    setInteractionOnsite!(UCglobal, i, A_ion * Sz)


    ##DM interaction
    addInteraction!(UCglobal, i, i, (-1)^(i+1) * D1ex, (1,0,0))
    addInteraction!(UCglobal, i, i, (-1)^(i+1) * D2ex, (0,1,0)) #(0,1)
    addInteraction!(UCglobal, i, i, (-1)^(i+1) * D3ex, (-1,-1,0)) #
end


L = (inputFile["System_Size"], inputFile["System_Size"])

(Harr,J2arr) = ndgrid(range(inputFile["H_min"],inputFile["H_max"],inputFile["H_length"]),range(inputFile["J2_min"],inputFile["J2_max"],inputFile["J2_length"]) )
Hs = collect(Iterators.flatten(Harr))
J2s = collect(Iterators.flatten(J2arr))

gridsize =inputFile["H_length"]*inputFile["J2_length"]

elements_per_process = div(gridsize, commSize)
remainder = rem(gridsize, commSize)

SkXnumberPhase = zeros(length(Hs),length(J2s))

 
start_index = commRank * elements_per_process + min(commRank, remainder) + 1
end_index = start_index + elements_per_process - 1 + (commRank < remainder ? 1 : 0)

#println(commSize, " commSize?")
for (j2idx, j2) in enumerate(J2s[start_index:end_index])
   println(j2idx, "index")
   h = round(Hs[j2idx],sigdigits=3)
   println(h, Hs[j2idx], "confused")
   j2 = round(j2,sigdigits=3)
   println("Rank " , commRank , " working on h = " , h, " working on j2 = ", j2) 
   filename = "/scratch/andykh/02_Data/Bilayer_Runs/"*ARGS[1]*"_H=$h,J2=$j2.h5"
   #println(filename)
   if isfile(filename) 
        println("Already Completed "*filename)
   else
      UClocal = deepcopy(UCglobal)

      #Add J2 2NN AF interaction 
      addInteraction!(UClocal, 1, 1, -j2 * I, (-1,1))
      addInteraction!(UClocal, 1, 1, -j2 * I, (1,2))
      addInteraction!(UClocal, 1, 1, -j2 * I, (2,1))
      
      #Local Magnetic field
      setField!(UClocal, 1, [0,0,-h])

      latticeLocal = Lattice(UClocal, L)

      mc = runAnneal(t0,tf,latticeLocal,thermSweeps,measureSweeps,inputFile["coolRate"], h, j2,filename);

   end
end


