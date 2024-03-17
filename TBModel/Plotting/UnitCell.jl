using Plots, LinearAlgebra, ColorSchemes
using MeanFieldToolkit, TightBindingToolkit, FixedPointToolkit
# a1 = [-3 / 2, 3 * sqrt(3) / 2]
# a2 = [9.0, 0]
a1 = [-1.0, sqrt(3)]
a2 = [6.0, 0.0]
l1 = [1.0, 0]
l2 = [-0.5, sqrt(3) / 2]
UC = UnitCell([a1, a2], 2, 2)
t1Param = Param(1.0, 2)
su2spin = SpinMats(1 // 2)
for j = 0:1
    for i = 0:5
        AddBasisSite!(UC, i .* l1 + j .* l2)
    end
end
AddIsotropicBonds!(t1Param, UC, 1.0, su2spin[4], "t1", checkOffsetRange=1)
CreateUnitCell!(UC, t1Param)
Plot_UnitCell!(UC, range=0, bond_cmp=:viridis)
