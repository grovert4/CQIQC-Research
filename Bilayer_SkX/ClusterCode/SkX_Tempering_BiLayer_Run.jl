using MPI
MPI.Initialized() || MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)

using SpinMC_more_more, LinearAlgebra, JSON, LazyGrids
include("functions.jl")

inputFile = JSON.parsefile("./Input_Files/Tempering_Input_File/"*ARGS[1]*".json")
J1 = inputFile["J_1"]
D = inputFile["D"]
A_ion = inputFile["A_ion"]
t0 = inputFile["t_max"]
tf = inputFile["t_min"]
thermSweeps = inputFile["thermalizationSweeps"]
measureSweeps = inputFile["measurementSweeps"]

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

D1v = D .* [1, 0, 0]
D2v = D .* [-1/2,sqrt(3)/2 , 0] 
D3v = D .* [-1/2, -sqrt(3)/2, 0] 

ExchangeD(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]


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


(Jperparr,J2arr) = ndgrid(range(inputFile["Jperp_min"],inputFile["Jperp_max"],inputFile["Jperp_length"]),range(inputFile["J2_min"],inputFile["J2_max"],inputFile["J2_length"]) )
Jperps = collect(Iterators.flatten(Jperparr))
J2s = collect(Iterators.flatten(J2arr))


L = (inputFile["System_Size"], inputFile["System_Size"], 1)
gridsize =inputFile["Jperp_length"]*inputFile["J2_length"]

j2 = inputFile["J2"]

for j in 1:length(Jperps)
      jperp = round(Jperps[j],sigdigits=5)

      filename = "/scratch/grovert4/Data/Tempering/longertry/"*ARGS[1]*"_Jperp=$(jperp),J2=$(j2).h5"
      if isfile(filename) 
         println("Already Completed "*filename)
      else
         UClocal = deepcopy(UCglobal)
         for i in 1:length(UClocal.basis)
            #Add J2 2NN AF interaction 
            addInteraction!(UClocal, i, i, -j2 * I, (-1,1,0))
            addInteraction!(UClocal, i, i, -j2 * I, (1,2,0))
            addInteraction!(UClocal, i, i, -j2 * I, (2,1,0))
            end
         addInteraction!(UClocal, b1, b2, -jperp * Sz , (0,0,0))
         latticeLocal = Lattice(UClocal, L)
   
         mc = MPIrunAnneal(inputFile["tmax"],inputFile["tmin"],inputFile["exchangeRate"],t0,tf,latticeLocal,thermSweeps,measureSweeps,inputFile["coolRate"],filename,true);
      end
end
