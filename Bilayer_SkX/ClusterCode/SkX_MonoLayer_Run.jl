#using Pkg
#Pkg.add(url="https://github.com/grovert4/SpinMC_more_more.jl")
using MPI
MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)

using SpinMC_more_more, LinearAlgebra, LazyGrids, JSON
include("functions.jl")

inputFile = JSON.parsefile("./Input_Files/"*ARGS[1]*".json")

#Parameters

J1 = inputFile["J_1"]
D = inputFile["D"]
A_ion = inputFile["A_ion"]
t0 = inputFile["t_max"]
tf = inputFile["t_min"]
thermSweeps = inputFile["thermalizationSweeps"]
measureSweeps = inputFile["measurementSweeps"]
cores = commSize

#Unit Cell Construction
a1 = (1.0 , 0.0)  #-
a2 = (-1/2 , sqrt(3)/2)  #/
UCglobal = UnitCell(a1 , a2)

b1 = addBasisSite!(UCglobal, (0.0, 0.0)) ##layer A (z = 0)

#Helpful Matrices
I = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
Sz = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]

#DM Vector
D1v =  D .* [1, 0, 0]
D2v =  D .* [-1/2,sqrt(3)/2 , 0]
D3v =  D .* [-1/2, -sqrt(3)/2, 0]
Dex(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]

##Global ayer interactions

#1nd NN ferromagnetic interaction
addInteraction!(UCglobal, 1, 1, -J1 * I, (1,0))
addInteraction!(UCglobal, 1, 1, -J1 * I, (0,1))
addInteraction!(UCglobal, 1, 1, -J1 * I, (-1,-1))

##onsite anisotropy
setInteractionOnsite!(UCglobal, 1, A_ion * Sz)

##DM interaction
addInteraction!(UCglobal, 1, 1, Dex(D1v), (1,0))
addInteraction!(UCglobal, 1, 1, Dex(D2v), (0,1)) #(0,1)
addInteraction!(UCglobal, 1, 1, Dex(D3v), (-1,-1)) #

L = (inputFile["System_Size"], inputFile["System_Size"])

(Harr,J2arr) = ndgrid(range(inputFile["H_min"],inputFile["H_max"],inputFile["H_length"]),range(inputFile["J2_min"],inputFile["J2_max"],inputFile["J2_length"]) )
Hs = collect(Iterators.flatten(Harr))
J2s = collect(Iterators.flatten(J2arr))

gridsize =inputFile["H_length"]*inputFile["J2_length"]

elements_per_process = div(gridsize, commSize)
remainder = rem(gridsize, commSize)

start_index = commRank * elements_per_process + min(commRank, remainder) + 1
end_index = start_index + elements_per_process - 1 + (commRank < remainder ? 1 : 0)

for i in start_index:end_index
   h = round(Hs[i],sigdigits=5)
   j2 = round(J2s[i],sigdigits=5)
   filename = "/scratch/grovert4/Data/Monolayer_Runs_Take2_36x36/"*ARGS[1]*"_H=$h,J2=$j2.h5"

   if isfile(filename) 
      println("Already Completed "*filename)
   else
      println("Rank " , commRank , " working on h = " , h, " working on j2 = ", j2) 

      UClocal = deepcopy(UCglobal)

      #Add J2 2NN AF interaction 
      addInteraction!(UClocal, 1, 1, -j2 * I, (-1,1))
      addInteraction!(UClocal, 1, 1, -j2 * I, (1,2))
      addInteraction!(UClocal, 1, 1, -j2 * I, (2,1))
      
      #Local Magnetic field
      setField!(UClocal, 1, [0,0,-h])

      latticeLocal = Lattice(UClocal, L)

      mc = runAnneal(t0,tf,latticeLocal,thermSweeps,measureSweeps,inputFile["coolRate"],filename);

   end
end



