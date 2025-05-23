using ApproxOperator,CairoMakie,Tensors, BenchmarkExample, Statistics
import Gmsh: gmsh
lwb = 2;lwm =1.5;mso =20;msx =1;ppu = 2.5;α = 0.7;
# filename1 = "./msh/cantilever_nonuniform_16.msh"
filename1 = "./msh/square_tri3_8.msh"
# filename2 = "./msh/cantilever_bubble_4165.msh"
# savename = "./png/cantilever_nonuniform_16.png"
savename = "./png/square_tri3_8.png"
# filename1 = "./msh/plate_with_hole_45.msh"
# filename2 = "./msh/plate_with_hole_2.msh"
# savename = "./png/plate_with_hole_2_45.png"

gmsh.initialize()
gmsh.open(filename1)

entities = getPhysicalGroups()
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
integrationOrder = 2
elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)


# gmsh.open(filename1)
nodes_p = get𝑿ᵢ()
xᵖ = nodes_p.x
yᵖ = nodes_p.y
zᵖ = nodes_p.z
elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"])


if occursin("quad",filename1)
    index = [1,2,3,4,1]
else
    index = [1,2,3,1]
end

f = Figure(backgroundcolor = :transparent)
ax = Axis(f[1,1],aspect = DataAspect(),backgroundcolor = :transparent)
hidespines!(ax)
hidedecorations!(ax)
L = 1.
b = 1.

for elm in elements["Ω"]
    id = [node.𝐼 for node in elm.𝓒]
    # lines!(x[id[index]],y[id[index]], linewidth = lwm, color = :grey)
end
lines!([0.0,L,L,0.0,0.0],[-b/2,-b/2,b/2,b/2,-b/2], linewidth = lwb, color = :black)
scatter!(x,y,marker = :circle, markersize = mso, color = :black)

# for elm in elements["Ωᵖ"]
#     id = [node.𝐼 for node in elm.𝓒]
#     # lines!(xᵖ[id[[1,2,3,1]]],yᵖ[id[[1,2,3,1]]], linewidth = lwm, color = :blue)
# end
# scatter!(xᵖ,yᵖ,marker = :xcross, markersize = msx, color = (:blue, α))

save(savename,f,px_per_unit = ppu)
f