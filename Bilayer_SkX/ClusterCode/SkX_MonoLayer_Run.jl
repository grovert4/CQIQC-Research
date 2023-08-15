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
cores = 40
t0 = 1E-1
tf = 1E-3
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


#UClocal0 = deepcopy(UCglobal)
#Lattice0 = Lattice(UClocal0, L)
#vertex=getVertex(Lattice0)

SkXnumberPhase = zeros(length(Hs),length(J2s))


# need to combine looped for loop 
# How to use multiple nodes? 
println(commSize)
for (j2idx, j2) in enumerate(J2s[1+ceil(Int64, num/cores)*commRank:1+ceil(Int64,num/cores)*(commRank+1)])
   for (hidx,h) in enumerate(Hs[1+ceil(Int64, num/cores)*commRank:1+ceil(Int64,num/cores)*(commRank+1)])
      println("Rank " , commRank , " working on h = " , h, "working on j2 = ", j2) 

      filename = "/scratch/andykh/02_Data/Monolayer_Runs/H=$h,J2=$j2.hdf"
      println(filename)

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

         mc = runAnneal(t0,tf,latticeLocal,thermSweeps,measureSweeps,0.99, h, j2,"/scratch/andykh/02_Data/Monolayer_Runs/H=$h,J2=$j2.hdf");
            # DetailedMonoPlot(mc,mc.lattice,vertex)
            # SkXnumberPhase[hidx, j2idx] = round(getSkyrmionNumber(0,mc.lattice,vertex),digits=1)
      end
   end
end



