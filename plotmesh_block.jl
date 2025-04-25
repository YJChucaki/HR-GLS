using ApproxOperator,CairoMakie,Tensors, BenchmarkExample, Statistics
import Gmsh: gmsh
lwb = 1.5;lwm =1.0;mso =8;msx =1;ppu = 2.5;α = 0.7;
ndiv = 2
ploy= "nonuniform"
# ploy= "tet4"
filename = "./msh/block_"*ploy*"_"*string(ndiv)*".msh"
savename = "./png/block_"*ploy*"_"*string(ndiv)*".png"

gmsh.initialize()
gmsh.open(filename)

entities = getPhysicalGroups()
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes, entities["Ω"])
elements["Ω"] = getElements(nodes,entities["Ω"])
elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"])
elements["Γᵍ"] = getElements(nodes,entities["Γᵍ"])
elements["Γʳ"] = getElements(nodes,entities["Γʳ"])
elements["∂Ω"] = elements["Γʳ"]∪elements["Γᵍ"]∪elements["Γᵗ"]
# elements["∂Ω"] = elements["Γᵍ"]

# gmsh.open(filename1)


tetra_indices = [
    [1,2,3,1],    # 底面
    [1,4], [2,4], [3,4]  # 侧棱
]

if occursin("quad",filename)
    index = [1,2,3,4,1]
else
    index = [1,2,3,1]
end
f = Figure(backgroundcolor = :transparent)


ax = Axis3( # 3 4
    f[1,1],
    xlabel = " ",
    ylabel = " ",
    xticksvisible = false,
    yticksvisible = false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
    backgroundcolor = :transparent,
    aspect = :data,
    # aspect = (1,1,0.3), # 0.3(3) 0.2(4)
    # azimuth = 1.5π, # 3
    # elevation = 0.1π, # 3
)
hidespines!(ax)
hidedecorations!(ax)

ps = Point3f.(x,y,z)


# 绘制三维边界框
L = 1.0
b = 1.0
h = 1.0  


cube_lines = [
    # 底面
    [[0,0,0], [L,0,0], [L,b,0], [0,b,0], [0,0,0]],
    # 顶面
    [[0,0,h], [L,0,h], [L,b,h], [0,b,h], [0,0,h]],
    # 侧棱
    [[0,0,0], [0,0,h]],
    [[L,0,0], [L,0,h]],
    [[L,b,0], [L,b,h]],
    [[0,b,0], [0,b,h]]
]

# cube_vertices = [
#     [0.0, 0.0, 0.0],  # 1
#     [L, 0.0, 0.0],    # 2
#     [L, b, 0.0],      # 3
#     [0.0, b, 0.0],    # 4
#     [0.0, 0.0, h],    # 5
#     [L, 0.0, h],      # 6
#     [L, b, h],        # 7
#     [0.0, b, h]       # 8
# ]

cube_vertices = [
    Point3f(0.0, 0.0, 0.0),  # 1
    Point3f(L, 0.0, 0.0),   # 2
    Point3f(L, b, 0.0),  # 3
    Point3f(0.0, b, 0.0), # 4
    Point3f(0.0, 0.0, h),   # 5
    Point3f(L, 0.0, h),    # 6
    Point3f(L, b, h),   # 7
    Point3f(0.0, b, h)   # 8
]
cube_faces = [
    [1, 2, 3, 4],  # 底面 (z=0)
    [5, 6, 7, 8],   # 顶面 (z=h)
    [1, 2, 6, 5],   # 前面 (y=0)
    [3, 4, 8, 7],   # 后面 (y=b)
    [2, 3, 7, 6],   # 右面 (x=L)
    [1, 4, 8, 5]    # 左面 (x=0)
]

# # 绘制所有四边形面
# for face in quad_faces
#     points = cube_vertices[face]
#     # 使用poly!绘制四边形
#     poly!(ax, points, color=:blue, strokecolor=:black, strokewidth=2,
#           transparency=true, alpha=0.5)
# end

# cube_faces = [
#     # [1, 2, 3],  # 底面（正确顺序）
#     # [3, 4, 1],
#     [5, 6, 7],  # 顶面
#     [7, 8, 5],  
#     [1, 2, 6],  # 前面（1→2→6→5）
#     [6, 5, 1],
#     # [2, 3, 7],  # 右面（2→3→7→6）
#     # [7, 6, 2],
#     # [3, 4, 8],  # 后面（3→4→8→7）
#     # [8, 7, 3],
#     [4, 1, 5],   # 左面（4→1→5→8）
#     [5, 8, 4],   # 左面（4→1→5→8）
# ]

# 绘制立方体的外表面
# for (i, face) in enumerate(cube_faces)
#     # 提取当前面的顶点坐标
#     x_face = [cube_vertices[v][1] for v in face]
#     y_face = [cube_vertices[v][2] for v in face]
#     z_face = [cube_vertices[v][3] for v in face]
    
#     # 绘制面
#     poly!(ax, [Point3f.(x_face, y_face, z_face)], color=:gray, transparency=true, alpha=0.2)

#     # mesh!(ax, [Point3f(x_face[j], y_face[j], z_face[j]) for j in 1:4], color = :gray,depth_shift = -1f-4,transparency = true, alpha = 0.2 , shading = true,overdraw = false )
# end


for lines in cube_lines
    for i in 1:length(lines)-1
        seg = lines[i:i+1]
        x = [p[1] for p in seg]
        y = [p[2] for p in seg]
        z = [p[3] for p in seg]
        lines!(ax, x, y, z, linewidth=1, color=:black)
    end
end

# elements


# 绘制三维网格元素
for elm in elements["Ω"]
    id = [node.𝐼 for node in elm.𝓒]
    # 根据元素类型选择索引（示例为六面体）
    if length(id) == 8
        indices = hexa_indices
    else
        indices = tetra_indices
    end
    for seg in indices
        xs = [x[id[i]] for i in seg]
        ys = [y[id[i]] for i in seg]
        zs = [z[id[i]] for i in seg]
        # lines!(ax, xs, ys, zs, linewidth=1, color=:gray,transparency = true, alpha = 0.5)
    end
end


# # 绘制点
# scatter!(ps, 
#     marker=:circle,
#     markersize = 4,
#     overdraw = true,
#     depth_shift = 1f-1,
#     color = :gray,
#     transparency = true, 
#     alpha = 0.3
# )

# for elm in elements["∂Ω"]
#     id = [node.𝐼 for node in elm.𝓒]
#     # lines!(x[id[index]],y[id[index]],z[id[index]], linewidth = lwm, color = :black)
#     scatter!(x[id[index]],y[id[index]],z[id[index]],markersize = 4,color = :black,overdraw = true,depth_shift = 1f-1)
# end


for elm in elements["∂Ω"]
    id = [node.𝐼 for node in elm.𝓒]
    # lines!(x[id[index]],y[id[index]],z[id[index]], linewidth = lwm, color = :black)
    scatter!(x[id[index]],y[id[index]],z[id[index]],markersize = 4,color = :black,overdraw = true,depth_shift = 1f-1)
end


# for elm in elements["Γᵗ"]
#     id = [node.𝐼 for node in elm.𝓒]
#     # lines!(x[id[index]],y[id[index]],z[id[index]], linewidth = lwm, color = :black)
#     scatter!(x[id[index]],y[id[index]],z[id[index]],markersize = 4,color = :black,overdraw = true,depth_shift = 1f-1)
# end


# for elm in elements["Γʳ"]
#     id = [node.𝐼 for node in elm.𝓒]
#     # lines!(x[id[index]],y[id[index]],z[id[index]], linewidth = lwm, color = :black)
#     scatter!(x[id[index]],y[id[index]],z[id[index]],markersize = 4,color = :black,overdraw = true,depth_shift = 1f-1)
# end

save(savename,f,px_per_unit = ppu)

f