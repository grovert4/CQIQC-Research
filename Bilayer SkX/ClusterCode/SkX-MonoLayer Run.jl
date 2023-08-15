using SpinMC_more_more, LinearAlgebra, Plots, JSON
include("functions.jl")


inputFile = JSON.parsefile("inputParametersMonoLayer.json")

#Unit Cell Construction
a1 = (1.0 , 0.0)  #-
a2 = (-1/2 , sqrt(3)/2)  #/
UCglobal = UnitCell(a1 , a2)

b1 = addBasisSite!(UCglobal, (0.0, 0.0)) ##layer A (z = 0)

J1 = inputFile["J_1"]
D = inputFile["D"]
A_ion = inputFile["A_ion"]

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


Hs = union(collect(range(inputFile["H_min"], inputFile["H_max"],length = inputFile["H_length"])))
J2s = collect(range(inputFile["J2_max"],inputFile["J2_min"],length = inputFile["J2_length"]))

UClocal0 = deepcopy(UCglobal)
Lattice0 = Lattice(UClocal0, L)
vertex=getVertex(Lattice0)

SkXnumberPhase = zeros(length(Hs),length(J2s))

t0 = inputFile["t_max"]
tf = inputFile["t_min"]
coolRate = inputFile["coolRate"]
thermalizationSweeps = inputFile["thermalizationSweeps"]
measurementSweeps = inputFile["measurementSweeps"]

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

      mc = runAnneal(t0,tf,latticeLocal,thermalizationSweeps,measurementSweeps,coolRate, H, J2,"Official-Cluster-Run-1-MonoLayer/H=$h,J2=$j2");
      # DetailedMonoPlot(mc,mc.lattice,vertex)
      # SkXnumberPhase[hidx, j2idx] = round(getSkyrmionNumber(0,mc.lattice,vertex),digits=1)
   end
end



