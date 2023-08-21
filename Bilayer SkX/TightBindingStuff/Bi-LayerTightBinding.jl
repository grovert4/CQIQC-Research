using Plots, TightBindingToolkit, LinearAlgebra, ColorSchemes

#Skyrmion Lattice
# a1 = [3.0 , sqrt(3)]
# a2 = [3.0 , -sqrt(3)]
a1 = [-3.0, sqrt(3)]
a2 = [3.0, sqrt(3)]
UC = UnitCell([a1 , a2] , 4) 

##Triangular Lattice 
l1 = [1.0, 0]
l2 = [-0.5, sqrt(3)/2]

##Parameters
t = 1.0
jh = -5.0
su2spin = SpinMats(1//2)
su4spin = SpinMats(3//2)

##Adding inner-hexagon structure  
for j = 1:2
    for i = -1:4
        AddBasisSite!(UC, i .* l1 + j .* l2)
    end
end

##Adding boundary points of outer-hexagon 
# AddBasisSite!(UC, 2 * l2)
# AddBasisSite!(UC, 2 * l2 + l1)
# AddBasisSite!(UC, 2 * l2 + 2 * l1)
# AddBasisSite!(UC, l2 + 2 * l1)
# AddBasisSite!(UC, -l1 + l2)

##Istrotropic bonds
AddIsotropicBonds!(UC, 1.0, -t * su4spin[4], "iso", checkOffsetRange = 1)

##Functions that will be useful for adding anisotropic bonds
tau0(v) = [sin(pi * (1 - norm(v)/2)) * v[1]/norm(v), sin(pi * (1 - norm(v)/2 )) * v[2]/norm(v), cos(pi * (1 - norm(v)/2 ))]
tau1(v) = [sin(pi * (1 - norm(v)/2)) * v[1]/norm(v), sin(pi * (1 - norm(v)/2 )) * v[2]/norm(v), -cos(pi * (1 - norm(v)/2 ))]
sigmav(i,j) = 2 .* [su2spin[1][i,j], su2spin[2][i,j], su2spin[3][i,j]]

s11 = sigmav(1,1)
s12 = sigmav(1,2)
s21 = sigmav(2,1)
s22 = sigmav(2,2)
intermat(s1, s2) = [dot(s1, s11) dot(s1, s12) 0 0;dot(s1, s21) dot(s1, s22) 0 0; 0 0 dot(s2, s11) dot(s2, s12);0 0 dot(s2, s21) dot(s2, s22)]

##Adding anisotropic bonds and normalizing if needed
for (ind, bas) in enumerate(UC.basis)
    if 1 < norm(bas) < 2
        mat = intermat(normalize(tau0(bas) + tau0(-bas)), normalize(tau1(bas) + tau1(-bas))) 
    else 
        closest = [bas, bas-a1, bas-a2]
        clv = closest[findmin(x -> norm(x), closest)[2]] 
        mat = intermat(tau0(clv), tau1(clv))
    end
    replace!(mat, NaN + NaN * im => 0.0 + 0.0 * im)
    AddAnisotropicBond!(UC, ind, ind, [0,0], -jh * mat, 0.0, "interaction")
end

##Plotting the unit cell
plot_UC = Plot_UnitCell!(UC);
# display(plot_UC)

##Creating BZ and Hamiltonian Model
kSize = 6 * 15 + 3  
bz = BZ(kSize, 2)
FillBZ!(bz, UC)
path = CombinedBZPath(bz, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]] ; nearest=true)
H = Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)
Mdl = Model(UC, bz, H)
SolveModel!(Mdl)

##Plotting the band structure
bands = Plot_Band_Structure!(Mdl, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]] , labels = ["G", "K1", "M2"], plot_legend=false);
# plot(bands, legend = false)
# display(bands)

##Calculating Chern Numbers for bands
# totC = 0
# for i in 1:2*length(UC.basis)
#     c = ChernNumber(H, [i])
#     println(round(c))
# end
println(ChernNumber(H, collect(1:2*length(UC.basis))))


function spinsSkx()
    xspin, yspin, zspin =[],[],[]
    xpos,ypos = [],[]
    for (ind, bas) in enumerate(UC.basis)
    end

    for j = -3:3
        for i = -3:3
            bas = i .* l1 + j .* l2
            if norm(bas) <= 2
                # if 1 < norm(bas) < 2
                #     push!(xspin, normalize(tau(bas) + tau(-bas))[1]) 
                #     push!(yspin, normalize(tau(bas) + tau(-bas))[2]) 
                #     push!(zspin, normalize(tau(bas) + tau(-bas))[3]) 
                #     push!(xpos, bas[1])
                #     push!(ypos, bas[2])
                # else 
                push!(xpos, bas[1])
                push!(ypos, bas[2])
                push!(xspin, tau(bas)[1])
                push!(yspin, tau(bas)[2])
                push!(zspin, tau(bas)[3])
            end
        end
    end
    replace!(xspin, NaN=> 0.0)
    replace!(yspin, NaN=> 0.0)
    replace!(zspin, NaN=> 0.0)


    display(Plots.quiver(xpos,ypos,quiver=(xspin,yspin), line_z=repeat([[zspin zspin]'...], inner=2), c=:bwr, legend=false, aspect_ratio=:equal, background_color = :black))
end


function BandColor(bs)
    for i in 1:2 * length(UC.basis)
        c = round(ChernNumber(H, [i]))

        bs.series_list[i].plotattributes[:linecolor] = ColorSchemes.oslo10[findfirst(x -> x==c, cs)]
    end
end