using MPI
MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)

using SpinMC_more_more, LinearAlgebra, LazyGrids, JSON
include("functions.jl")

inputFile = JSON.parsefile("./Input_Files/"*ARGS[1]*".json")

#Parameters 

J1 = inputFile["J_1"]
# D = inputFile["D"]
# A_ion = inputFile["A_ion"]
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
D1v =  [1, 0, 0]
D2v =  [-1/2,sqrt(3)/2 , 0]
D3v =  [-1/2, -sqrt(3)/2, 0]
Dex(v, D) = D .* [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]

##Global ayer interactions

#1nd NN ferromagnetic interaction
addInteraction!(UCglobal, 1, 1, -J1 * I, (1,0))
addInteraction!(UCglobal, 1, 1, -J1 * I, (0,1))
addInteraction!(UCglobal, 1, 1, -J1 * I, (-1,-1))

L = (inputFile["System_Size"], inputFile["System_Size"])

Hs = union(collect(range(inputFile["H_min"], inputFile["H_max"],length = inputFile["H_length"])))
J2s = collect(range(inputFile["J2_max"],inputFile["J2_min"],length = inputFile["J2_length"]))
Ds = union(collect(range(inputFile["D_min"], inputFile["D_max"],length = inputFile["D_length"])))
A_ions = union(collect(range(inputFile["A_ion_min"], inputFile["A_ion_max"],length = inputFile["A_ion_length"])))


# (Harr,J2arr, Darr, Aarr) = ndgrid(range(0,0.8,5),range(0,-0.4,5), range(0, 0.8, 5), range(0, 0.4, 5))
(Harr,J2arr, Darr, Aarr) = ndgrid(range(inputFile["H_min"],inputFile["H_max"],inputFile["H_length"]),range(inputFile["J2_min"],inputFile["J2_max"],inputFile["J2_length"]),range(inputFile["D_min"],inputFile["D_max"],inputFile["D_length"]),range(inputFile["A_ion_min"],inputFile["A_ion_max"],inputFile["A_ion_length"]))

Hs = collect(Iterators.flatten(Harr))
J2s = collect(Iterators.flatten(J2arr))
Ds = collect(Iterators.flatten(Darr))
As = collect(Iterators.flatten(Aarr))

gridsize = inputFile["H_length"]*inputFile["J2_length"]*inputFile["D_length"]*inputFile["A_ion_length"]

elements_per_process = div(gridsize, commSize)
remainder = rem(gridsize, commSize)

start_index = commRank * elements_per_process + min(commRank, remainder) + 1
end_index = start_index + elements_per_process - 1 + (commRank < remainder ? 1 : 0)

for i in start_index:end_index
   h = round(Hs[i],sigdigits=5)
   j2 = round(J2s[i],sigdigits=5)
   d = round(Ds[i],sigdigits=5)
   aion = round(As[i],sigdigits=5)

   filename = "/scratch/grovert4/Data/Parameter_Search_Runs/"*ARGS[1]*"_H=$h,J2=$j2,D=$d,A=$aion.h5"

   if isfile(filename) 
      println("Already Completed "*filename)
   else
      println("Rank " , commRank , " working on h = " , h, " working on j2 = ", j2, " working on D = " , d, " working on A_ion = ", aion) 

      UClocal = deepcopy(UCglobal)

      #Add J2 2NN AF interaction 
      addInteraction!(UClocal, 1, 1, -j2 * I, (-1,1))
      addInteraction!(UClocal, 1, 1, -j2 * I, (1,2))
      addInteraction!(UClocal, 1, 1, -j2 * I, (2,1))
      
      #Local Magnetic field
      setField!(UClocal, 1, [0,0,-h])
      
      ##onsite anisotropy
      setInteractionOnsite!(UClocal, 1, aion * Sz)

      ##DM interaction
      addInteraction!(UClocal, 1, 1, Dex(D1v, D), (1,0))
      addInteraction!(UClocal, 1, 1, Dex(D2v, D), (0,1)) #(0,1)
      addInteraction!(UClocal, 1, 1, Dex(D3v, D), (-1,-1)) #
      
      latticeLocal = Lattice(UClocal, L)

      mc = runAnneal(t0,tf,latticeLocal,thermSweeps,measureSweeps,inputFile["coolRate"],filename);

   end
end



