using MPI, MPIPreferences
# Initialize MPI
MPI.Init()
commRank = MPI.Comm_rank(MPI.COMM_WORLD)
commSize = MPI.Comm_size(MPI.COMM_WORLD)
println("Hello, I am $(commRank) of $(commSize)")