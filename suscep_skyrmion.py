from triqs_tprf.tight_binding import TBLattice
import numpy as np
import matplotlib.pyplot as plt
from triqs_tprf.lattice_utils import k_space_path
import inspect
import string
from triqs_tprf.lattice import lindhard_chi00
from triqs.gf import MeshImFreq, Idx
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

###### Pauli matrices ###########
s1 = np.matrix([[0,1],[1,0]])
s2 = np.matrix([[0,-1j],[1j,0]])
s3 = np.matrix([[1,0],[0,-1]])
sup = np.matrix([[1,0],[0,0]])
sdo = np.matrix([[0,0],[0,1]])
paulis = [s1,s2,s3]

##### return the Hunds coupling matrix
def hunds(B: list) -> np.matrix:
    return B[0]*s1 + B[1]*s2 + B[2]*s3

##Triangular Lattice primitives
l1 = np.array([1.0, 0.0, 0.0])
l2 = np.array([-0.5, np.sqrt(3)/2, 0.0])

##### returns the positions of the atoms in the triangular lattice for a n1xn2 supercell
def positions(n1: int, n2:int, orbitals:int) -> list:
    pos = []
    names = []
    orbital_names = list(string.ascii_lowercase[0:orbitals])
    for i in range(n1):
        for j in range(n2):
            for orbital in range(orbitals):
                pos.append((i*l1[0] + j*l2[0], i*l1[1] + j*l2[1], 0.0))
                names.append(str(n2*i+j+1) + orbital_names[orbital])
    return pos, names
    
######### returns the correct hunds coupled matrix ##############
def skyrmion():
    B_center = [0.0, 0.0, 1.0]
    B_edge = [0.0, 0.0, -1.0]
    B_mid1 = [-1.0, 0.0, 0.0]
    B_mid2 = [-0.5, np.sqrt(3)/2, 0.0]
    B_mid3 = [0.5, np.sqrt(3)/2, 0.0]
    B_mid4 = [1.0, 0.0, 0.0]
    B_mid5 = [0.5, -np.sqrt(3)/2, 0.0]
    B_mid6 = [-0.5, -np.sqrt(3)/2, 0.0]
    
    Bs = {"center" : B_center, "edge" : B_edge, 
                "mid1" : B_mid1, "mid2" : B_mid2, "mid3" : B_mid3, 
                "mid4" : B_mid4, "mid5" : B_mid5, "mid6" : B_mid6}
    
    positions = {12 : "center", 
                        1 : "edge", 8 : "edge", 9 : "edge", 10 : "edge", 4 : "edge",
                        7 : "mid1", 6 : "mid2", 5 : "mid3", 11: "mid4", 2 : "mid5", 3 : "mid6"}
    
    mat = np.zeros((24, 24), dtype=np.complex128)
    
    for i in range(12):
        mat[2*i:2*i+2, 2*i:2*i+2] += hunds(Bs[positions[i+1]])
        
    return mat
    
####### returns the nearest neighbour hopping matrix ##########
def hopping(orbital_space: np.matrix, neighbours: list) -> np.matrix:
    mat = np.zeros((24, 24), dtype=np.complex128)
    for neighbour in neighbours:
        i,j = neighbour
        mat[2*(i-1):2*(i-1)+2, 2*(j-1):2*(j-1)+2] += orbital_space
    return mat

##### Nearest neighbours within unit cell : 21
NN_00 = [(1, 2), (1, 7), (1, 8), (2, 3), (2, 8), (2, 9), (3, 4), (3, 9), (3, 10), 
        (4, 5), (4, 10), (4, 11), (5, 6), (5, 11), (5, 12), 
        (6, 12), 
        (7, 8), (8, 9), (9, 10), (10, 11), (11, 12)] 
##### between unit cell at (0, 1) : 11
NN_01 = [(7, 3), (7, 4), (8, 4), (8, 5), (9, 5), (9, 6), (10, 6), (10, 1), (11, 1), (11, 2), (12, 2)]
##### between unit cell at (1, 0) : 3
NN_10 = [(6, 1), (6, 7), (12, 7)]
##### between unit cell at (1, 1) : 1
NN_11 = [(12, 3)]
##### all bonds
bonds = {(0, 0) : NN_00, (0, 1) : NN_01, (1, 0) : NN_10, (1, 1) : NN_11}

t1 = 1.0
B = 1.0 * t1

pos, names = positions(6, 2, 2)
########## tight binding model #################
model = TBLattice(
    units = [(-3.0, np.sqrt(3), 0.0), [3.0, np.sqrt(3), 0.0]],##### Primitive vectors
    orbital_positions = pos,##### 2 sublattice x 2 spins
    orbital_names = names,
    hoppings = {
        # Intra-diamond hopping
        ( 0, 0): -t1 * hopping(np.eye(2), bonds[(0, 0)]) - t1 * np.conj(np.transpose(hopping(np.eye(2), bonds[(0, 0)]))) - B * skyrmion(),##### within the unit cell
        (+1, 0): -t1 * hopping(np.eye(2), bonds[(1, 0)]),
        (-1, 0): -t1 * np.conj(np.transpose(hopping(np.eye(2), bonds[(1, 0)]))),
        ( 0,+1): -t1 * hopping(np.eye(2), bonds[(0, 1)]),
        ( 0,-1): -t1 * np.conj(np.transpose(hopping(np.eye(2), bonds[(0, 1)]))),
        ( 1, 1): -t1 * hopping(np.eye(2), bonds[(1, 1)]),
        (-1,-1): -t1 * np.conj(np.transpose(hopping(np.eye(2), bonds[(1, 1)]))),

        },
    )

kmesh = model.get_kmesh(n_k=(8, 8, 1))
e_k = model.fourier(kmesh)

G = [0.0, 0.0, 0.0]
K1 = [1/3, 1/3, 0.0]
K2 = [2/3, 2/3, 0.0]
M1 = [0.5, 0.0, 0.0]
M2 = [0.0, 0.5, 0.0]
M3 = [0.5, 0.5, 0.0]

paths = [(G, K1), (K1, M2), (M2, G)]
k_vecs, k_plot, k_ticks = k_space_path(paths, bz=model.bz)

def energies(k, e_k):
    return np.linalg.eigvalsh(e_k(k))
##### Band structure
# bands_plot = plt.plot(k_plot, [energies(k, e_k) for k in k_vecs])
# plt.xticks(k_ticks, [r'$\Gamma$',r'$K_1$',r'$M_2$',r'$\Gamma$'])
# plt.ylabel(r'$\epsilon(\mathbf{k})$')
# plt.grid(True)
# plt.show()

##### fermi distribution function
def fermi(e, beta, mu):
    return 1.0 / (np.exp(beta * (e-mu)) + 1.0)
##### returns band energies at each k-value
def bands(e_k, kmesh, N):
    band = [np.linalg.eigvalsh(e_k(kmesh[i].value)) for i in range(N**2)]
    band = np.concatenate(band, axis = 0)
    
    return band

##### filling at fixed temperature and chemical potential
def filling(band, beta, mu):
    return np.sum(fermi(band, beta, mu))/len(band)

##### Temperature
beta = 10.0
mus = np.linspace(-7.0, 4.0, 111)
fillings = np.zeros(len(mus))
print("calculating chemical potential vs filling...")
for (i, mu) in enumerate(mus):
    fillings[i] = filling(bands(e_k, kmesh, 32), beta, mu)
    
# mu_vs_filling = plt.plot(fillings, mus)
# plt.ylabel(r'$\mu$')
# plt.xlabel(r'$n$')
# plt.hlines(-0.4, 0.0, 1.0, ls='--', color='orange')
# plt.vlines(0.3, -2.5, 2.5, ls='--', color='orange')
# plt.show(mu_vs_filling)
print("calculating bare bubble...")
ind = 40
mu = mus[ind]
filling = fillings[ind]
wmesh = MeshImFreq(beta=10.0, S='Fermion', n_max=30)
g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=1)