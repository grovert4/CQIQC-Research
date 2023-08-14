using Pkg
Pkg.add(url="https://github.com/grovert4/SpinMC_more_more.jl")
using SpinMC_more_more, LinearAlgebra, Plots
include("functions.jl")

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

L = (24, 24)


Hs = union(collect(range(0.0,1.1,length = 36)))
J2s = collect(range(0.0,-0.5,length = 36))

UClocal0 = deepcopy(UCglobal)
Lattice0 = Lattice(UClocal0, L)
vertex=getVertex(Lattice0)

SkXnumberPhase = zeros(length(Hs),length(J2s))

for (j2idx, j2) in enumerate(J2s)
   global J2 = j2
   for (hidx,h) in enumerate(Hs)

      global H = h
      UClocal = deepcopy(UCglobal)

      #Add J2 2NN AF interaction 
      addInteraction!(UClocal, 1, 1, -J2 * I, (-1,1))
      addInteraction!(UClocal, 1, 1, -J2 * I, (1,2))
      addInteraction!(UClocal, 1, 1, -J2 * I, (2,1))
      
      #Local Magnetic field
      setField!(UClocal, 1, [0,0,-H])

      latticeLocal = Lattice(UClocal, L)

      mc = runAnneal(69,680,latticeLocal,2500,250000,0.99, H, J2,"Official-Cluster-Run-1-MonoLayer/H=$h,J2=$j2");
      # DetailedMonoPlot(mc,mc.lattice,vertex)
      # SkXnumberPhase[hidx, j2idx] = round(getSkyrmionNumber(0,mc.lattice,vertex),digits=1)
   end
end



