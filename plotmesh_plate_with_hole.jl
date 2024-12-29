using ApproxOperator,CairoMakie
import Gmsh: gmsh
n=2
# lwb = 2.0;lwm = 0.5;mso = 7;msx = 7;ppu = 2.5;α = 0.7;
lwb = 1;lwm = 2;mso = 18;msx = 18;ppu = 2.5;α = 0.7;
filename1 = "./msh/plate_with_hole_tri3_"*string(n)*".msh"
filename2 = "./msh/plate_with_hole_tri3_"*string(n)*".msh"
savename = "./png/2.png"

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

if occursin("quad",filename1)
    index = [1,2,3,4,1]
else
    index = [1,2,3,1]
end

f = Figure(backgroundcolor = :transparent)
ax = Axis(f[1,1],aspect = DataAspect(),backgroundcolor = :transparent)
hidespines!(ax)
hidedecorations!(ax)
# L = 48.
# b = 12.
# lines!([0.0,L,L,0.0,0.0],[-b/2,-b/2,b/2,b/2,-b/2], linewidth = lwb, color = :black)

for elm in elements["Ω"]
    id = [node.𝐼 for node in elm.𝓒]
    # lines!(x[id[index]],y[id[index]], linewidth = lwm, color = :black)
end
scatter!(x,y,marker = :circle, markersize = mso, color = :black)
# lines!([1.0,5.0,5.0,0.0],[0.0,0.0,5.0,5.0], linewidth = lwb, color = :black)
# Circle!([0.0,0.0,1.0],[1.0,0.0,0.0], linewidth = lwb, color = :black)
# for elm in elms_p["Ω"]
#     id = [i for i in elm.i]
#     lines!(xᵖ[id[[1,2,3,1]]],yᵖ[id[[1,2,3,1]]], linewidth = lwm, color = :blue)
# end
# scatter!(xᵖ,yᵖ,marker = :xcross, markersize = msx, color = (:blue, α))
save(savename,f,px_per_unit = ppu)
f