using Pkg
Pkg.add(url="https://github.com/grovert4/SpinMC_more_more.jl")
using SpinMC_more_more, LinearAlgebra, Plots
include("functions.jl")

using MPI
MPI.Initialized() || MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)

#Unit Cell Construction
a1 = (1.0 , 0.0, 0.0)  #-
a2 = (-1/2 , sqrt(3)/2, 0.0)  #/
a3 = (0.0, 0.0, 2.0)  #|
UCglobal = UnitCell(a1, a2, a3)

b1 = addBasisSite!(UCglobal, (0.0, 0.0, 0.0)) ##layer A (z = 0)
b2 = addBasisSite!(UCglobal, (0.0, 0.0, 1.0)) ##layer A  (z = 1)

#Parameters
J2 = 0.0
J1 = 1.0
A_ion = 0.2
H = 0.0
J_ll = -0.2
D = 0.25

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


Hs = union(collect(range(0.0,1.1,length = 36)))
J2s = collect(range(0.0,-0.5,length = 36))

# addInteraction!(UC, b1, b2, -J_ll * Sz, (0,0,-1))

L = (24, 24, 1)

UClocal0 = deepcopy(UCglobal)
Lattice0 = Lattice(UClocal0, L)
vertex=getVertex(Lattice0)

for (j2idx, j2) in enumerate(J2s)
   global J2 = j2
   for (hidx,h) in enumerate(Hs)
      
      global H = h
      UClocal = deepcopy(UCglobal)

      for i in 1:length(UCglobal.basis)
         #Add J2 2NN AF interaction 
         addInteraction!(UClocal, i, i, -J2 * I, (-1,1,0))
         addInteraction!(UClocal, i, i, -J2 * I, (1,2,0))
         addInteraction!(UClocal, i, i, -J2 * I, (2,1,0))
         
         #Local Magnetic field
         setField!(UClocal, i, [0,0,-H])
      end

      latticeLocal = Lattice(UClocal, L)

      mc = runAnneal(69,680,latticeLocal,2500,250000,0.99, H, J2,"Official-Cluster-Run-1-BiLayer/H=$h,J2=$j2");
      # DetailedMonoPlot(mc,mc.lattice,vertex)
      # SkXnumberPhase[hidx, j2idx] = round(getSkyrmionNumber(0,mc.lattice,vertex),digits=1)
   end
end


