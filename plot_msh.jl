
using ApproxOperator, GLMakie

import Gmsh: gmsh

ndiv = 8
gmsh.initialize()
# gmsh.open("./msh/plate_with_hole_tri3_"*string(ndiv)*".msh")
# gmsh.open("./msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate/SquarePlate_quad_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate/SquarePlate_tri6_"*string(ndiv)*".msh")
gmsh.open("./msh/circular/circular_"*string(ndiv)*".msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes,entities["Ω"])
# elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"])
elements["Γᵇ"] = getElements(nodes,entities["Γᵇ"])
elements["Γˡ"] = getElements(nodes,entities["Γˡ"])
# elements["Γʳ"] = getElements(nodes,entities["Γʳ"])
# elements["∂Ω"] = elements["Γᵇ"]∪elements["Γᵗ"]∪elements["Γˡ"]∪elements["Γʳ"]

elements["Γᵉ"] = getElements(nodes,entities["Γᵉ"])
elements["∂Ω"] = elements["Γᵇ"]∪elements["Γˡ"]∪elements["Γᵉ"]

# gmsh.finalize()

f = Figure()

# axis
ax = Axis3(f[1, 1], perspectiveness = 0.8, aspect = :data, azimuth = -0.5*pi, elevation = 0.5*pi, xlabel = " ", ylabel = " ", zlabel = " ", xticksvisible = false,xticklabelsvisible=false, yticksvisible = false, yticklabelsvisible=false, zticksvisible = false, zticklabelsvisible=false, protrusions = (0.,0.,0.,0.))
hidespines!(ax)
hidedecorations!(ax)

x =  nodes.x
y = nodes.y
z = 0
ps = Point3f.(x,y,z)
scatter!(ps, 
    marker=:circle,
    markersize = 20,
    color = :black
)

# elements
for elm in elements["Ω"]
    x = [x.x for x in elm.𝓒[[1,2,3,1]]]
    y = [x.y for x in elm.𝓒[[1,2,3,1]]]

    lines!(x,y,linestyle = :dash, linewidth = 0.5, color = :black)
end

# # boundaries
for elm in elements["∂Ω"]
    ξ¹ = [x.x for x in elm.𝓒]
    ξ² = [x.y for x in elm.𝓒]
    x =  [x.x for x in elm.𝓒]
    y =  [x.y for x in elm.𝓒]
    lines!(x,y,linewidth = 1.5, color = :black)
end

save("./png/Circular_"*string(ndiv)*"_msh.png",f)
# save("./png/plate_with_hole_tri3_"*string(ndiv)*"_msh.png",f)
# save("./png/SquarePlate_quad_"*string(ndiv)*"_msh.png",f)
# save("./png/MorleysAcuteSkewPlate_"*string(ndiv)*"_msh.png",f)
f