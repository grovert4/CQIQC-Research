using SpinMC_more_more, LinearAlgebra, LazyGrids, JSON
include("functions.jl")
using MPI
MPI.Initialized() || MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)
inputFile = JSON.parsefile("./Input_Files/"*ARGS[1]*".json")

#Parameters

J1 = inputFile["J_1"]
A_ion = inputFile["A_ion"]
T = inputFile["temperature"]
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

##Global ayer interactions

#1nd NN ferromagnetic interaction
addInteraction!(UCglobal, 1, 1, -J1 * I, (1,0))
addInteraction!(UCglobal, 1, 1, -J1 * I, (0,1))
addInteraction!(UCglobal, 1, 1, -J1 * I, (-1,-1))

##onsite anisotropy
setInteractionOnsite!(UCglobal, 1, A_ion * Sz)

L = (inputFile["System_Size"], inputFile["System_Size"])

(Harr,J3arr) = ndgrid(range(inputFile["H_min"],inputFile["H_max"],inputFile["H_length"]),range(inputFile["J3_min"],inputFile["J3_max"],inputFile["J3_length"]) )
Hs = collect(Iterators.flatten(Harr))
J3s = collect(Iterators.flatten(J3arr))

gridsize =inputFile["H_length"]*inputFile["J3_length"]

elements_per_process = div(gridsize, commSize)
remainder = rem(gridsize, commSize)

start_index = commRank * elements_per_process + min(commRank, remainder) + 1
end_index = start_index + elements_per_process - 1 + (commRank < remainder ? 1 : 0)

for i in start_index:end_index
   h = round(Hs[i],sigdigits=5)
   j3 = round(J3s[i],sigdigits=5)
   filename = "/scratch/andykh/02_Data/Monolayer_Runs/"*ARGS[1]*"_H=$h,J3=$j3.h5"

   if isfile(filename) 
      println("Already Completed "*filename)
   else
      println("Rank " , commRank , " working on h = " , h, " working on j3 = ", j3) 

      UClocal = deepcopy(UCglobal)

      #Add J3 3NN AF interaction 
      addInteraction!(UClocal, 1, 1, -j3 * I, (2,0))
      addInteraction!(UClocal, 1, 1, -j3 * I, (0,2))
      addInteraction!(UClocal, 1, 1, -j3 * I, (-2,-2))
      
      #Local Magnetic field
      setField!(UClocal, 1, [0,0,-h])

      latticeLocal = Lattice(UClocal, L)

      mc = MonteCarlo(latticeLocal, 1/T, thermSweeps, measureSweeps, reportInterval = 250000, rewrite = true);
      run_nompi!(mc, outfile = filename)

   end
end
