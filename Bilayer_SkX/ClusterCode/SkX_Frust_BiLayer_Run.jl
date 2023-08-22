using SpinMC_more_more, LinearAlgebra, JSON, LazyGrids
include("functions.jl")
using MPI
MPI.Initialized() || MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)
inputFile = JSON.parsefile("./Input_Files/"*ARGS[1]*".json")
J1 = inputFile["J_1"]
A_ion = inputFile["A_ion"]
T = inputFile["temperature"]
thermSweeps = inputFile["thermalizationSweeps"]
measureSweeps = inputFile["measurementSweeps"]
J_ll = inputFile["J_perp"]

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

addInteraction!(UCglobal, b1, b2, -J_ll * Sz , (0,0,0))

##Same layer interactions
for i in 1:length(UCglobal.basis)
    #1nd NN ferromagnetic interaction
    addInteraction!(UCglobal, i, i, -J1 * I, (1,0,0))
    addInteraction!(UCglobal, i, i, -J1 * I, (0,1,0))
    addInteraction!(UCglobal, i, i, -J1 * I, (-1,-1,0))


    ##onsite anisotropy
    setInteractionOnsite!(UCglobal, i, A_ion * Sz)
end


(Harr,J3arr) = ndgrid(range(inputFile["H_min"],inputFile["H_max"],inputFile["H_length"]),range(inputFile["J3_min"],inputFile["J3_max"],inputFile["J3_length"]) )
Hs = collect(Iterators.flatten(Harr))
J3s = collect(Iterators.flatten(J3arr))


L = (inputFile["System_Size"], inputFile["System_Size"], 1)
gridsize =inputFile["H_length"]*inputFile["J3_length"]

elements_per_process = div(gridsize, commSize)
remainder = rem(gridsize, commSize)
 
start_index = commRank * elements_per_process + min(commRank, remainder) + 1
end_index = start_index + elements_per_process - 1 + (commRank < remainder ? 1 : 0)

#println(commSize, " commSize?")
for i in start_index:end_index
   h = round(Hs[i],sigdigits=5)
   j3 = round(J3s[i],sigdigits=5)

   println("Rank " , commRank , " working on h = " , h, " working on j3 = ", j3) 
   filename = "/scratch/andykh/02_Data/Bilayer_Runs/"*ARGS[1]*"_H=$h,J3=$j3.h5"
   #println(filename)
   if isfile(filename) 
      println("Already Completed "*filename)
   else
      UClocal = deepcopy(UCglobal)
      for i in 1:length(UCglobal.basis)
         #Add J3 2NN AF interaction 
         addInteraction!(UClocal, i, i, -j3 * I, (2,0,0))
         addInteraction!(UClocal, i, i, -j3 * I, (0,2,0))
         addInteraction!(UClocal, i, i, -j3 * I, (-2,-2,0))
         
         #Local Magnetic field
         setField!(UClocal, i, [0,0,-h])
      end
      latticeLocal = Lattice(UClocal, L)
      mc = MonteCarlo(latticeLocal, 1/T, thermSweeps, measureSweeps, reportInterval = 250000, rewrite = true);
      run_nompi!(mc, outfile = filename)
   end
end