#using Pkg
#Pkg.add(url="https://github.com/grovert4/SpinMC_more_more.jl")
using SpinMC_more_more, LinearAlgebra, LazyGrids
include("functions.jl")
using MPI
MPI.Initialized() || MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)

#Unit Cell Construction
a1 = (1.0 , 0.0)  #-
a2 = (-1/2 , sqrt(3)/2)  #/
UCglobal = UnitCell(a1 , a2)

b1 = addBasisSite!(UCglobal, (0.0, 0.0)) ##layer A (z = 0)

#Parameters
J1 = 1.0
A_ion = 0.2
D = 0.25
J_ll = 0.0
Hlength = 30
Jlength = 30
size = 36
cores = commSize
t0 = 1/1E-1
tf = 1/1E-3
thermSweeps = 2500
measureSweeps = 250000


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

L = (size, size)

num = Hlength * Jlength
(Harr, J2arr) = ndgrid(range(0.0,1.0,Hlength),range(0.0,-0.5,Jlength) )
Hs = collect(Iterators.flatten(Harr))
J2s = collect(Iterators.flatten(J2arr))

gridsize = Hlength*Jlength

elements_per_process = div(gridsize, commSize)
remainder = rem(gridsize, commSize)

SkXnumberPhase = zeros(length(Hs),length(J2s))

 
start_index = commRank * elements_per_process + min(commRank, remainder) + 1
end_index = start_index + elements_per_process - 1 + (commRank < remainder ? 1 : 0)

println(commSize, "commsize?")
for (j2idx, j2) in enumerate(J2s[start_index:end_index])
   for (hidx,h) in enumerate(Hs[start_index:end_index])
      h = round(h,sigdigits=3)
      j2 = round(j2,sigdigits=3)
      println("Rank " , commRank , " working on h = " , h, " working on j2 = ", j2) 
      filename = "/scratch/andykh/02_Data/Monolayer_Runs/H=$h,J2=$j2.hdf"
      #println(filename)
      if isfile(filename) 
           println("Already Completed "*filename)
      else
         println("Inside Here")
         UClocal = deepcopy(UCglobal)

         #Add J2 2NN AF interaction 
         addInteraction!(UClocal, 1, 1, -j2 * I, (-1,1))
         addInteraction!(UClocal, 1, 1, -j2 * I, (1,2))
         addInteraction!(UClocal, 1, 1, -j2 * I, (2,1))
         
         #Local Magnetic field
         setField!(UClocal, 1, [0,0,-h])

         latticeLocal = Lattice(UClocal, L)
         println("Running Simulation")

         mc = runAnneal(t0,tf,latticeLocal,thermSweeps,measureSweeps,0.99, h, j2,"/scratch/andykh/02_Data/Monolayer_Runs/H=$h,J2=$j2.hdf");
         println("Finished Simulation")

      end
   end
end



