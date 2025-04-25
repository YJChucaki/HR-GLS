
using TimerOutputs 
using SparseArrays, Pardiso, Printf
using CairoMakie, XLSX, WriteVTK
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric,∫∫σᵢⱼσₖₗdxdy_PlaneStrian,∫∫σᵢⱼσₖₗdxdy,∫σᵢⱼnⱼuᵢds,∫∫∇σᵢⱼuᵢdxdy,∫σᵢⱼnⱼgᵢds,∫∫τ∇σᵢⱼ∇σᵢₖdxdy,∫∫τ∇σᵢⱼ∇σᵢₖdxdy_new,∫∫τ∇σᵢⱼ∇σᵢₖdxdy_Taylor,∫∫σᵢⱼσₖₗdxdy_Taylor

include("import_cook.jl")
include("wirteVTK.jl")

const to = TimerOutput()
ps = MKLPardisoSolver()
nn = [ 2 4 8 16 ]
for i in 4:4
ndiv = nn[i]
ndiv2 = nn[i]
n = nn[i]

# ndiv = 18
# ndiv2 = 18
# n = 18
# poly = "tri3"
test = "cook"
poly = "nouniform"
# poly = "tri6"
# poly = "quad"
@timeit to "import data" begin
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_c_"*string(nₚ)*".msh")
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly*"_"*string(ndiv2)*".msh")
# elements, nodes, sp, type, Ω, nodes_c= import_HR_GLS("./msh/cook_"*poly*"_"*string(ndiv)*".msh","./msh/cook_"*poly*"_"*string(ndiv2)*".msh",n)
# elements, nodes,  sp, type = import_HR_GLS_reduced("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly*"_"*string(ndiv2)*".msh")
# elements, nodes, sp, type, Ω, nodes_c = import_HR_GLS("./msh/cook_membrane_nonuniform_"*string(ndiv)*".msh","./msh/cook_membrane_nonuniform_"*string(ndiv2)*".msh",n)
elements, nodes, sp, type, Ω, nodes_c= import_HR_GLS("./msh/cook_membrane_"*poly*"_"*string(ndiv)*".msh","./msh/cook_membrane_"*poly*"_"*string(ndiv2)*".msh",n)
end
nc = length(nodes_c)
nₑ = length(elements["Ωᵘ"])
nₑₛ = length(elements["Ωˢ"])
nₛ = 3
nᵤ = length(nodes)
# nₚ = length(nodes_p)
# nₚ = length(nodes)
ni = 6
L = 48.0
D = 44
P = 6.25
ℎ = D/ndiv
# ℎ = 1.0
# Ē = 3e6
# ν̄  = 0.3
E = 70.0
# ν = 0.3 
# ν = 0.5-1e-6
ν = 0.49999
Ē = E/(1.0-ν^2)
ν̄ = ν/(1.0-ν)
I = D^3/12
EI = Ē*I
Cᵢᵢᵢᵢ = E/(1+ν)/(1-2*ν)*(1-ν)
Cᵢᵢⱼⱼ = E/(1+ν)/(1-2*ν)*ν
Cᵢⱼᵢⱼ = E/(1+ν)/2
𝐺 = E/(1+ν)/2
K=E/3/(1-2ν )

ℎ = zeros(nₑₛ*nₛ)
β = zeros(nₑₛ*nₛ)
for (i,elm) in enumerate(elements["Ω"])
    𝓒 = elm.𝓒
    𝓒ₚ = elements["Ωˢ"][i].𝓒
    x1 = 𝓒[1].x
    x2 = 𝓒[2].x
    x3 = 𝓒[3].x
    y1 = 𝓒[1].y
    y2 = 𝓒[2].y
    y3 = 𝓒[3].y
    A = 0.5*sqrt((x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))^2)
    # A = 0.5*abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))
    hₑ= sqrt(A)
     for (j,xⱼ) in enumerate(𝓒ₚ)
        J = xⱼ.𝐼
        ℎ[J] =hₑ
        β[J] =0.1*hₑ^2/2/𝐺
     end
end
for elm in elements["Ωˢ"]
    𝓒ₚ = elm.𝓒
        push!(𝓒ₚ,:β=>β)
        push!(𝓒ₚ,:ℎ=>ℎ)
end


# β =0.1*ℎ^2/2/𝐺
# prescribe!(elements["Ωˢ"],:τ=>(x,y,z)->β)
# prescribe!(elements["Ωˢ"],:ℎ=>(x,y,z)->ℎ) 
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
    # 𝑏ᵝ(kˢˢ,fˢ)
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
    dₛ₁₁ = d[1:3:3*nₛ*nₑₛ]
    dₛ₂₂ = d[2:3:3*nₛ*nₑₛ]
    dₛ₁₂ = d[3:3:3*nₛ*nₑₛ]
    push!(nodes,:d₁=>d₁,:d₂=>d₂)
    for elm in elements["Ωˢ"]
        𝓒ₚ = elm.𝓒
        𝓖 = elm.𝓖
            push!(𝓒ₚ,:dₛ₁₁=>dₛ₁₁,:dₛ₂₂=>dₛ₂₂,:dₛ₁₂=>dₛ₁₂)
    end

    pₑ = zeros(nₑ)
    
    for (i,elm) in enumerate(elements["Ωˢ"])
        𝓒ₚ = elm.𝓒
        𝓖 = elm.𝓖
        𝓒 = elements["Ω"][i].𝓒
        a = length(𝓒)
        x = 0.0
        y = 0.0
        for j in 𝓒 
            x += j.x
            y += j.y
        end
        xc = x/a
        yc = y/a
        if nₛ==3
        σ₁₁ = 𝓒ₚ[1].dₛ₁₁+𝓒ₚ[2].dₛ₁₁*xc+𝓒ₚ[3].dₛ₁₁*yc
        σ₂₂ = 𝓒ₚ[1].dₛ₂₂+𝓒ₚ[2].dₛ₂₂*xc+𝓒ₚ[3].dₛ₂₂*yc
        elseif nₛ==6
            σ₁₁ = 𝓒ₚ[1].dₛ₁₁+𝓒ₚ[2].dₛ₁₁*xc+𝓒ₚ[3].dₛ₁₁*yc+𝓒ₚ[4].dₛ₁₁*xc^2+𝓒ₚ[6].dₛ₁₁*yc^2+𝓒ₚ[5].dₛ₁₁*xc*yc
            σ₂₂ = 𝓒ₚ[1].dₛ₂₂+𝓒ₚ[2].dₛ₂₂*xc+𝓒ₚ[3].dₛ₂₂*yc+𝓒ₚ[4].dₛ₂₂*xc^2+𝓒ₚ[6].dₛ₂₂*yc^2+𝓒ₚ[5].dₛ₂₂*xc*yc
        end
        σ₃₃ = ν*(σ₁₁ + σ₂₂)
        # pₑ[i] = σ₁₁
        pₑ[i] = (σ₁₁ + σ₂₂ + σ₃₃)/3 
    end
    
    p_node = zeros(nc)
    
    w = zeros(nc)
    for (i,elm) in enumerate(elements["Ω"])
        𝓒 = elm.𝓒
         for (j,xⱼ) in enumerate(𝓒)
            J = xⱼ.𝐼
            p_node[J] +=pₑ[i]
            w[J] +=1 
         end
    end
    
pc = zeros(nc)
for (i,elm) in enumerate(elements["Ω"])
    𝓒 = elm.𝓒
     for (j,xⱼ) in enumerate(𝓒)
        J = xⱼ.𝐼
        ξ¹ = xⱼ.x
        ξ² = xⱼ.y
        ∂ū₁∂x = -P/EI*(L-ξ¹)*ξ²
        ∂ū₁∂y = -P/6/EI*((6*L-3*ξ¹)*ξ¹ + (2+ν )*(3*ξ²^2-D^2/4))
        ∂ū₂∂x = P/6/EI*((6*L-3*ξ¹)*ξ¹ - 3*ν *ξ²^2 + (4+5*ν )*D^2/4)
        ∂ū₂∂y = P/EI*(L-ξ¹)*ξ²*ν 
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        pc[J] = K*(ε̄₁₁+ε̄₂₂)
        
     end
end

for elm in elements["Ω"]
    𝓒ₚ = elm.𝓒
    push!(𝓒ₚ,:pc=>pc)
end
    # eval(VTK_HR_displacement_pressure)
    # eval(VTK_HR_displacement_pressure_smoothing)
    # eval(VTK_HR_displacement_pressure_sigma11)


α = 1.0

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
    data = Dict([:x=>(2,[xs]),:y=>(2,[ys]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:∂𝝭∂x=>(4,zeros(ni)),:∂𝝭∂y=>(4,zeros(ni)),:𝗠=>(0,zeros(60)),:∂𝗠∂x=>(0,zeros(60)),:∂𝗠∂y=>(0,zeros(60))])
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
    σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ +Cᵢᵢⱼⱼ*ε₂₂
    σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁ +Cᵢᵢᵢᵢ*ε₂₂
    σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
    σ₃₃ = ν*(σ₁₁ + σ₂₂)
    p = (σ₁₁ + σ₂₂ + σ₃₃)/3 
    # p=K*(ε₁₁+ε₂₂)
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
vtk_grid("./vtk/cook_GLS_"*poly*"_"*string(ndiv)*"_"*string(ndiv),points,cells) do vtk
    vtk["𝑝"] = colors
    vtk["p_element"] = pₑ
    vtk["p_node"] = pc
end

println(y[3] - nodes_c[3].y)

show(to)

end
