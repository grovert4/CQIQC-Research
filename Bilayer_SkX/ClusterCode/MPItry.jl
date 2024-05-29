using MPI, LinearAlgebra, Plots, SpinMC_more_more
MPI.Initialized() || MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)

include("functions.jl")

#Unit Cell Construction
a1 = (1.0 , 0.0, 0.0)  #-
a2 = (-1/2 , sqrt(3)/2, 0.0)  #/
a3 = (0.0, 0.0, 2.0)  #|
UC = UnitCell(a1 , a2, a3) 


const b1 = addBasisSite!(UC, (0.0, 0.0, 0.0)) ##layer A (z = 0)
const b2 = addBasisSite!(UC, (0.0, 0.0, 1.0)) ##layer A  (z = 1)

#Parameters
J2 = -0.32
const J1 = 1.0
const A_ion = 0.2
J_ll = -0.85
const H = 0.0
D = 0.25

const I = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
# I2 = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0]
const Sz = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]

# D1v =  [0, -1, 0]
# D2v =  [sqrt(3)/2, 1/2, 0] 
# D3v =  [-sqrt(3)/2, 1/2, 0] 
const D1v = D .* [1, 0, 0]
const D2v = D .* [-1/2,sqrt(3)/2 , 0] 
const D3v = D .* [-1/2, -sqrt(3)/2, 0] 

ExchangeD(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]

const D1ex = ExchangeD(D1v)
const D2ex = ExchangeD(D2v)
const D3ex = ExchangeD(D3v)

#Same layer interactions
for i in 1:length(UC.basis)
    #1nd NN ferromagnetic interaction
    addInteraction!(UC, i, i, -J1 * I, (1,0,0))
    addInteraction!(UC, i, i, -J1 * I, (0,1,0))
    addInteraction!(UC, i, i, -J1 * I, (-1,-1,0))

    addInteraction!(UC, i, i, -J2 * I, (-1,1,0))
    addInteraction!(UC, i, i, -J2 * I, (1,2,0))
    addInteraction!(UC, i, i, -J2 * I, (2,1,0))

    ##onsite anisotropy
    setInteractionOnsite!(UC, i, A_ion * Sz)

    #Local Magnetic field
   #  setField!(UC, i, [0,0,0.0])

    ##DM interaction
    addInteraction!(UC, i, i, (-1)^(i+1) * D1ex, (1,0,0))
    addInteraction!(UC, i, i, (-1)^(i+1) * D2ex, (0,1,0)) #(0,1)
    addInteraction!(UC, i, i, (-1)^(i+1) * D3ex, (-1,-1,0)) #
end

addInteraction!(UC, b1, b2, -J_ll * Sz , (0,0,0))

L = (24, 24, 1)
global lattice = Lattice(UC, L);
   
const tmax = 10.0
const tmin = 1.0
coolRate = 0.99
# const thermSweeps = 25000
# const measureSweeps = 50000
const extfield = false


# m = MonteCarlo(lattice, beta, thermSweeps, measureSweeps,  replicaExchangeRate=10, reportInterval = thermSweeps)
# run!(m, outfile="MPITry/simulation.h5") 

ts = []
for i in 1:625
   push!(ts, i)
end

beta = (commSize == 1) ? 1.0/tmin : 1.0 / (reverse([ tmax * (tmin / tmax)^(n/(commSize-1)) for n in 0:commSize-1 ])[commRank+1])

for (ind,i) in enumerate(ts)
   thermalizationSweeps = 3000
   measurementSweeps = 0
   m = nothing   
   # h = lattice.unitcell.interactionsField[1][3]
   if ind == 1
       m = MonteCarlo(lattice, beta*(1/coolRate)^ind, thermalizationSweeps, measurementSweeps,  replicaExchangeRate=5, reportInterval = thermalizationSweeps)
       # if extfield
       #     m.lattice.interactionField[1:2:end] = repeat([(0.0, 0.0, h)], length(m.lattice.interactionField[1:2:end]))
       #     m.lattice.interactionField[2:2:end] = repeat([(0.0, 0.0, -h)], length(m.lattice.interactionField[1:2:end]))
       # end
   else
       if ind == length(ts)
           thermalizationSweeps = 0
           measurementSweeps = 25000
       end
       m = MonteCarlo(monte.lattice, beta*(1/coolRate)^ind, thermalizationSweeps, measurementSweeps,  replicaExchangeRate=10, reportInterval = 10000, rewrite = false);
   end

   if ind == length(ts)
      run!(m, outfile="/scratch/grovert4/Data/MPITry/simulation.h5")
   else
      run!(m, disableOutput=true)
   end
   global monte = deepcopy(m)
end

