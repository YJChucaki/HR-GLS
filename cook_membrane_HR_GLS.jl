
using TimerOutputs 
using SparseArrays, Pardiso, Printf
using CairoMakie, XLSX, WriteVTK
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric,∫∫σᵢⱼσₖₗdxdy_PlaneStrian,∫∫σᵢⱼσₖₗdxdy,∫σᵢⱼnⱼuᵢds,∫∫∇σᵢⱼuᵢdxdy,∫σᵢⱼnⱼgᵢds,∫∫τ∇σᵢⱼ∇σᵢₖdxdy,∫∫τ∇σᵢⱼ∇σᵢₖdxdy_new,∫∫τ∇σᵢⱼ∇σᵢₖdxdy_Taylor,∫∫σᵢⱼσₖₗdxdy_Taylor

include("import_cook.jl")
include("wirteVTK.jl")

const to = TimerOutput()
ps = MKLPardisoSolver()
nn = [ 4 8 12 16 ]
for i in 1:4
ndiv = nn[i]
ndiv2 = nn[i]
n = nn[i]
# ndiv = 4
# ndiv2 = 4
# n = 4
poly = "tri3"
# poly = "tri6"
# poly = "quad"
@timeit to "import data" begin
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_c_"*string(nₚ)*".msh")
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly*"_"*string(ndiv2)*".msh")
# elements, nodes, sp, type, Ω, nodes_c= import_HR_GLS("./msh/cook_"*poly*"_"*string(ndiv)*".msh","./msh/cook_"*poly*"_"*string(ndiv2)*".msh",n)
# elements, nodes,  sp, type = import_HR_GLS_reduced("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly*"_"*string(ndiv2)*".msh")
# elements, nodes, sp, type = import_HR_GLS("./msh/cantilever_nonuniform_"*string(ndiv)*".msh","./msh/cantilever_nonuniform_"*string(ndiv2)*".msh")
elements, nodes, sp, type, Ω, nodes_c= import_HR_GLS("./msh/cook_membrane_"*poly*"_"*string(ndiv)*".msh","./msh/cook_membrane_"*poly*"_"*string(ndiv2)*".msh",n)
end

nₑ = length(elements["Ωᵘ"])
nₛ = 6
nᵤ = length(nodes)
# nₚ = length(nodes_p)
# nₚ = length(nodes)

L = 48.0
D = 44
P = 6.25
ℎ = D/ndiv

# Ē = 3e6
# ν̄  = 0.3
E = 70.0
# ν = 0.3
ν = 0.5-1e-8
Ē = E/(1.0-ν^2)
ν̄ = ν/(1.0-ν)
I = D^3/12
EI = Ē*I
Cᵢᵢᵢᵢ = E/(1+ν)/(1-2*ν)*(1-ν)
Cᵢᵢⱼⱼ = E/(1+ν)/(1-2*ν)*ν
Cᵢⱼᵢⱼ = E/(1+ν)/2
𝐺 = E/(1+ν)/2
K=E/3/(1-2ν )



β =0.1*ℎ^2/2/𝐺
prescribe!(elements["Ωˢ"],:τ=>(x,y,z)->β)
prescribe!(elements["Ωˢ"],:ℎ=>(x,y,z)->ℎ) 
prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν)

prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν)

prescribe!(elements["Ωˢ"],:b₁=>(x,y,z)->0.0)
prescribe!(elements["Ωˢ"],:b₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P) 
# prescribe!(elements["Γʳ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
# prescribe!(elements["Γʳ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂) 
prescribe!(elements["Γᵍᵘ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍᵘ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₁₂=>(x,y,z)->0.0)
# prescribe!(elements["Ωᵍᵘ"],:u=>(x,y,z)->u(x,y))
# prescribe!(elements["Ωᵍᵘ"],:v=>(x,y,z)->v(x,y))
# prescribe!(elements["Ωᵍᵘ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
# prescribe!(elements["Ωᵍᵘ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
# prescribe!(elements["Ωᵍᵘ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
# prescribe!(elements["Ωᵍᵘ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))


𝑎 =∫∫σᵢⱼσₖₗdxdy_PlaneStrian=>elements["Ωˢ"]
# 𝑎 =∫∫σᵢⱼσₖₗdxdy_Taylor=>elements["Ωˢ"]
𝑏 = [
    ∫σᵢⱼnⱼuᵢds=>(elements["∂Ωˢ"],elements["∂Ωᵘ"]),
    ∫∫∇σᵢⱼuᵢdxdy=>(elements["Ωˢ"],elements["Ωᵘ"]),
    # ∫σᵢⱼnⱼuᵢds_Taylor=>(elements["∂Ωˢ"],elements["∂Ωᵘ"]),
    # ∫∫∇σᵢⱼuᵢdxdy_Taylor=>(elements["Ωˢ"],elements["Ωᵘ"]),
    ]

𝑏ᵅ = ∫σᵢⱼnⱼgᵢds=>(elements["Γᵍˢ"],elements["Γᵍᵘ"])
# 𝑏ᵅ = ∫σᵢⱼnⱼgᵢds_Taylor=>(elements["Γᵍˢ"],elements["Γᵍᵘ"])

𝑏ᵝ = ∫∫τ∇σᵢⱼ∇σᵢₖdxdy=>elements["Ωˢ"]
# 𝑏ᵝ = ∫∫τ∇σᵢⱼ∇σᵢₖdxdy_new=>elements["Ωˢ"]
# 𝑏ᵝ = ∫∫τ∇σᵢⱼ∇σᵢₖdxdy_Taylor=>elements["Ωˢ"]
𝑓 = ∫vᵢtᵢds=>elements["Γᵗ"]
# 𝑓 = [
#     ∫vᵢtᵢds=>elements["Γᵗ"]∪elements["Γʳ"],
#     ∫∫vᵢbᵢdxdy=>elements["Ωᵘ"]
# ]

kˢˢ = zeros(3*nₛ*nₑ,3*nₛ*nₑ)
kˢᵘ = zeros(3*nₛ*nₑ,2*nᵤ)
kˢᵘⁿ  = zeros(3*nₛ*nₑ,2*nᵤ)
fˢ = zeros(3*nₛ*nₑ)
fᵘ = zeros(2*nᵤ)



@timeit to "assembly" begin

    𝑎(kˢˢ)
    𝑏(kˢᵘ)
    𝑏ᵅ(kˢᵘ,fˢ)
    𝑏ᵝ(kˢˢ,fˢ)
    𝑓(fᵘ)
    end
    
    # k = sparse([kᵖᵖ kᵖᵘ;kᵖᵘ' zeros(2*nₚ,2*nₚ)])
    # set_matrixtype!(ps,-2)
    # k = get_matrix(ps,k,:N)
    # f = [fᵖ;-fᵘ]
    # @timeit to "solve" pardiso(ps,d,k,f)
    d = [kˢˢ kˢᵘ;kˢᵘ' zeros(2*nᵤ,2*nᵤ)]\[fˢ;-fᵘ]
    d₁ = d[3*nₛ*nₑ+1:2:end]
    d₂ = d[3*nₛ*nₑ+2:2:end]
    push!(nodes,:d₁=>d₁,:d₂=>d₂)
    
    # 𝐿₂ = L₂(elements["Ωᵍ"])
    # 𝐻ₑ, 𝐿₂ = Hₑ_PlaneStress(elements["Ωᵍᵘ"])
    # println(log10(𝐿₂))
    # println(log10(𝐻ₑ))
# println(log10(Hₑ_dev))
# println(log10(L₂_𝑝))
# eval(VTK_Guass_point)
# eval(displacement_stress)


α = 1.0
nc = length(nodes_c)
# vertices = [[node.x+α*node.d₁ for node in nodes] [node.y+α*node.d₂ for node in nodes]]
colors = zeros(nc)
x = zeros(nc)
y = zeros(nc)
𝗠 = zeros(21)
for (i,node_c) in enumerate(nodes_c)
    xs = node_c.x
    ys = node_c.y
    indices = sp(xs,ys,0.0)
    ni = length(indices)
    𝓒 = [nodes[i] for i in indices]
    # data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
    data = Dict([:x=>(2,[xs]),:y=>(2,[ys]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:∂𝝭∂x=>(4,zeros(ni)),:∂𝝭∂y=>(4,zeros(ni)),:𝗠=>(0,zeros(21)),:∂𝗠∂x=>(0,zeros(21)),:∂𝗠∂y=>(0,zeros(21))])
    ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
    𝓖 = [ξ]
    a = type(𝓒,𝓖)
    set∇𝝭!(a)
    p = 0.0
    d₁ = 0.0
    d₂ = 0.0
    u₁ = 0.0
    u₂ = 0.0
    𝝭 = ξ[:𝝭]
    B₁ = ξ[:∂𝝭∂x]
    B₂ = ξ[:∂𝝭∂y]
    ε₁₁ = 0.0
    ε₂₂ = 0.0
    ε₁₂ = 0.0
    N = ξ[:𝝭]
    for (k,xₖ) in enumerate(𝓒)
        ε₁₁ += B₁[k]*xₖ.d₁
        ε₂₂ += B₂[k]*xₖ.d₂
        ε₁₂ += B₁[k]*xₖ.d₂ + B₂[k]*xₖ.d₁
        u₁ += 𝝭[k]*xₖ.d₁
        u₂ += 𝝭[k]*xₖ.d₂
    end
    p=K*(ε₁₁+ε₂₂)
    x[i] = xs+α*u₁
    y[i] = ys+α*u₂
    colors[i] = p
end

# fig = Figure(figure_padding = 1,size = (400,600))
# ind = 100
# ax = Axis(fig[1,1], 
#     aspect = DataAspect(), 
#     xticksvisible = false,
#     xticklabelsvisible=false, 
#     yticksvisible = false, 
#     yticklabelsvisible=false,
#     backgroundcolor = :transparent,
# )
# hidespines!(ax)
# hidedecorations!(ax)

# x = [node.x+α*node.d₁ for node in nodes]
# y = [node.y+α*node.d₂ for node in nodes]
# # contourf!(ax,x,y,colors,levels=collect(-60:5:20), colormap=Reverse(:deep))
# contourf!(ax,x,y,colors,levels=5, colormap=Reverse(:deep))

# save("./png/cook_mix_"*poly*"_"*string(ndiv)*"_"*string(ndiv)*".png",fig, px_per_unit = 10.0)



# points = [[node.x+α*node.d₁ for node in nodes]';[node.y+α*node.d₂ for node in nodes]';zeros(1,nᵤ)]
points = [x';y';zeros(1,nc)]
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
vtk_grid("./vtk/cook_GLS_"*poly*"_"*string(ndiv)*"_"*string(ndiv)*"_"*string(β),points,cells) do vtk
    vtk["𝑝"] = colors
end

println(y[3] - nodes_c[3].y)
# @timeit to "plot figure" begin
# fig = Figure()
# ind = 100
# ax = Axis(fig[1,1], 
#     aspect = DataAspect(), 
#     xticksvisible = false,
#     xticklabelsvisible=false, 
#     yticksvisible = false, 
#     yticklabelsvisible=false,
# )
# hidespines!(ax)
# hidedecorations!(ax)
# xs = LinRange(0, 48, 4*ind)
# ys = LinRange(-6, 6, ind)
# zs = zeros(4*ind,ind)
# 𝗠 = zeros(6)
# for (i,x) in enumerate(xs)
#     for (j,y) in enumerate(ys)
#         indices = sp(x,y,0.0)
#         ni = length(indices)
#         𝓒 = [nodes[i] for i in indices]
#         data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:∂𝝭∂x=>(4,zeros(ni)),:∂𝝭∂y=>(4,zeros(ni)),:𝗠=>(0,𝗠),:∂𝗠∂x=>(0,𝗠),:∂𝗠∂y=>(0,𝗠)])
#         ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
#         𝓖 = [ξ]
#         a = type(𝓒,𝓖)
#         # set𝝭!(a)
#         set∇𝝭!(a)
#         d₁ = 0.0
#         d₂ = 0.0
#         B₁ = ξ[:∂𝝭∂x]
#         B₂ = ξ[:∂𝝭∂y]
#         ε₁₁ = 0.0
#         ε₂₂ = 0.0
#         ε₁₂ = 0.0
#         for (k,xₖ) in enumerate(𝓒)
#             ε₁₁ += B₁[k]*xₖ.d₁
#             ε₂₂ += B₂[k]*xₖ.d₂
#             ε₁₂ += B₁[k]*xₖ.d₂ + B₂[k]*xₖ.d₁
#         end
#         p=K*(ε₁₁+ε₂₂)
#         zs[i,j] = p
#     end
# end

# surface!(xs,ys,zeros(4*ind,ind),color=zs,shading=NoShading,colormap=:lightrainbow)
# contour!(xs,ys,zs,levels=-1e3:200:1e3,color=:azure)
# # Colorbar(fig[1,2], limits=(-900,900), colormap=:lightrainbow)
# save("./png/cantilever_mix_tri3_"*string(ndiv)*"_"*string(ndiv)*"_ls.png",fig, px_per_unit = 10.0)
# # end
# fig
show(to)

end
