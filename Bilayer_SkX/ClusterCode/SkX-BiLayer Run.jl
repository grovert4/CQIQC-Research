using SpinMC_more_more, LinearAlgebra, Plots, JSON
include("functions.jl")

inputFile = JSON.parsefile("inputParametersBiLayer.json")

#Unit Cell Construction
a1 = (1.0 , 0.0, 0.0)  #-
a2 = (-1/2 , sqrt(3)/2, 0.0)  #/
a3 = (0.0, 0.0, 2.0)  #|
UCglobal = UnitCell(a1, a2, a3)

b1 = addBasisSite!(UCglobal, (0.0, 0.0, 0.0)) ##layer A (z = 0)
b2 = addBasisSite!(UCglobal, (0.0, 0.0, 1.0)) ##layer A  (z = 1)

#Parameters
J1 = inputFile["J_1"]
D = inputFile["D"]
A_ion = inputFile["A_ion"]
J_ll = inputFile["J_perp"]

#Helpful Matrices
I = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
Sz = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]

D1v = D .* [1, 0, 0]
D2v = D .* [-1/2,sqrt(3)/2 , 0] 
D3v = D .* [-1/2, -sqrt(3)/2, 0] 

ExchangeD(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]

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
    addInteraction!(UCglobal, i, i, (-1)^(i+1) * ExchangeD(D1v), (1,0,0))
    addInteraction!(UCglobal, i, i, (-1)^(i+1) * ExchangeD(D2v), (0,1,0)) #(0,1)
    addInteraction!(UCglobal, i, i, (-1)^(i+1) * ExchangeD(D3v), (-1,-1,0)) #
end


Hs = union(collect(range(inputFile["H_min"], inputFile["H_max"],length = inputFile["H_length"])))
J2s = collect(range(inputFile["J2_max"],inputFile["J2_min"],length = inputFile["J2_length"]))

# addInteraction!(UC, b1, b2, -J_ll * Sz, (0,0,-1))

L = (inputFile["System_Size"], inputFile["System_Size"], 1)

UClocal0 = deepcopy(UCglobal)
Lattice0 = Lattice(UClocal0, L)
vertex=getVertex(Lattice0)

t0 = inputFile["t_max"]
tf = inputFile["t_min"]
coolRate = inputFile["coolRate"]
thermalizationSweeps = inputFile["thermalizationSweeps"]
measurementSweeps = inputFile["measurementSweeps"]


for (j2idx, j2) in enumerate(J2s)
   for (hidx,h) in enumerate(Hs)
      
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

      mc = runAnneal(t0,tf,latticeLocal,thermalizationSweeps, measurementSweeps, coolRate, H, J2,"Official-Cluster-Run-1-BiLayer/H=$h,J2=$j2");
      # DetailedMonoPlot(mc,mc.lattice,vertex)
      # SkXnumberPhase[hidx, j2idx] = round(getSkyrmionNumber(0,mc.lattice,vertex),digits=1)
   end
end
