
using TimerOutputs 
using SparseArrays, Pardiso, Printf, LinearAlgebra
using CairoMakie, WriteVTK
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric,∫∫σᵢⱼσₖₗdxdy_PlaneStrian,∫∫σᵢⱼσₖₗdxdy,∫σᵢⱼnⱼuᵢds,∫∫∇σᵢⱼuᵢdxdy,∫σᵢⱼnⱼgᵢds,∫∫τ∇σᵢⱼ∇σᵢₖdxdy,∫∫τ∇σᵢⱼ∇σᵢₖdxdy_new,∫∫τ∇σᵢⱼ∇σᵢₖdxdy_Taylor,∫∫σᵢⱼσₖₗdxdy_Taylor, Hₑ_PlaneStrain_Dil, 𝐿₂_PlaneStrain_Pressure, Hₑ_PlaneStrain_Deviatoric

include("import_cantilever.jl")
include("wirteVTK.jl")

const to = TimerOutput()
ps = MKLPardisoSolver()
# n = [ 2 4 8 16 ]
# for i in 1:4
# ndiv = n[i]
# ndiv2 = n[i]
ndiv = 8
ndiv2 = 8
poly = "tri3"
test = "cantilever"
# poly = "tri6"
# poly = "quad"
# poly = "nonuniform"
@timeit to "import data" begin
elements, nodes, sp, type, Ω, nodes_c = import_HR_GLS("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly*"_"*string(ndiv2)*".msh")
# elements, nodes, sp, type, Ω, nodes_c = import_HR_GLS("./msh/cantilever_nonuniform_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_nonuniform_"*poly*"_"*string(ndiv2)*".msh")
end
nₑ = length(elements["Ωᵘ"])
nₑₛ = length(elements["Ωˢ"])
nₛ = 3
ni = 3
nᵤ = length(nodes)
nc = length(nodes_c)

L = 48.0
D = 12.0
P = 1000
# ℎ = D/ndiv
ℎ = 1.0

Ē = 3e6
# ν̄  = 0.3
ν̄  = 0.5-1e-5
# E = 3e6
# ν = 0.3
# ν = 0.5-1e-4
E =Ē/(1.0-ν̄ ^2)
ν = ν̄ /(1.0-ν̄ )
I = D^3/12
EI = E*I
Cᵢᵢᵢᵢ = Ē/(1+ν̄ )/(1-2*ν̄ )*(1-ν̄ )
Cᵢᵢⱼⱼ = Ē/(1+ν̄ )/(1-2*ν̄ )*ν̄ 
Cᵢⱼᵢⱼ = Ē/(1+ν̄ )/2
𝐺 = Ē/(1+ν̄ )/2
K=Ē/3/(1-2ν̄  )


u(x,y) = -P*y/6/EI*((6*L-3*x)*x + (2+ν)*(y^2-D^2/4))
v(x,y) = P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2)
∂u∂x(x,y) = -P/EI*(L-x)*y
∂u∂y(x,y) = -P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4))
∂v∂x(x,y) = P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4)
∂v∂y(x,y) = P/EI*(L-x)*y*ν

ε₁₁(x,y) = ∂u∂x(x,y)
ε₂₂(x,y) = ∂v∂y(x,y)
ε₁₂(x,y) = ∂u∂y(x,y) + ∂v∂x(x,y)
σ₁₁(x,y) = -P*(L-x)*y/I
σ₂₂(x,y) = 0.0
σ₃₃(x,y) = Cᵢᵢⱼⱼ*ε₁₁(x,y) + Cᵢᵢⱼⱼ*ε₂₂(x,y)
σ₁₂(x,y) = P/2/I*(D^2/4-y^2)
p(x,y) = (σ₁₁(x,y)+σ₂₂(x,y)+σ₃₃(x,y))/3

β =1*ℎ^2/2/𝐺
# β =1/2/𝐺
prescribe!(elements["Ωˢ"],:τ=>(x,y,z)->β)
prescribe!(elements["Ωˢ"],:ℎ=>(x,y,z)->ℎ) 
prescribe!(elements["Ωˢ"],:Ē=>(x,y,z)->Ē)
prescribe!(elements["Ωˢ"],:ν̄ =>(x,y,z)->ν̄ )

prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν)

prescribe!(elements["Ωˢ"],:b₁=>(x,y,z)->0.0)
prescribe!(elements["Ωˢ"],:b₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂) 
# prescribe!(elements["Γʳ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
# prescribe!(elements["Γʳ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂) 
prescribe!(elements["Γᵍᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍᵘ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γᵍᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ωᵍᵘ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍᵘ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
prescribe!(elements["Ωˢ"],:p=>(x,y,z)->p(x,y))
prescribe!(elements["Ωˢ"],:p=>(x,y,z)->p(x,y))

𝑎 =∫∫σᵢⱼσₖₗdxdy=>elements["Ωˢ"]
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

kˢˢ = zeros(3*nₛ*nₑₛ,3*nₛ*nₑₛ)
kˢᵘ = zeros(3*nₛ*nₑₛ,2*nᵤ)
kˢᵘⁿ  = zeros(3*nₛ*nₑₛ,2*nᵤ)
fˢ = zeros(3*nₛ*nₑₛ)
fᵘ = zeros(2*nᵤ)



@timeit to "assembly" begin

    𝑎(kˢˢ)
    𝑏(kˢᵘ)
    𝑏ᵅ(kˢᵘ,fˢ)
    𝑏ᵝ(kˢˢ,fˢ)
    𝑓(fᵘ)
end
    
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

    # for elm in elements["Ω"]
    #     𝓒ₚ = elm.𝓒
    #     𝓖 = elm.𝓖
    #         push!(𝓒ₚ,:dₛ₁₁=>dₛ₁₁,:dₛ₂₂=>dₛ₂₂,:dₛ₁₂=>dₛ₁₂)
    # end
    
                  
𝐻ₑ, 𝐿₂ = Hₑ_PlaneStress(elements["Ωᵍᵘ"])
Hₑ_dev = Hₑ_PlaneStrain_Deviatoric(elements["Ωᵍᵘ"])
Hₑ_dil = Hₑ_PlaneStrain_Dil(elements["Ωᵍᵘ"])
L₂_𝑝 = 𝐿₂_PlaneStrain_Pressure(elements["Ωᵍᵘ"])
# L₂𝑝_stress = L₂𝑝_stress(elements["Ωˢ"])
println(log10(𝐿₂))
println(log10(𝐻ₑ))
println(log10(Hₑ_dev))
println(log10(Hₑ_dil))
println(log10(L₂_𝑝))
# println(log10(L₂_stress))


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
    pₑ[i]= (σ₁₁ + σ₂₂ + σ₃₃)/3 
    
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

eval(VTK_HR_displacement_pressure)
eval(VTK_HR_displacement_pressure_smoothing)
# eval(VTK_exact_pressure)


α = 1.0
nc = length(nodes_c)
# vertices = [[node.x+α*node.d₁ for node in nodes] [node.y+α*node.d₂ for node in nodes]]
colors = zeros(nc)
x = zeros(nc)
y = zeros(nc)
for (i,node) in enumerate(nodes)
    xs = node.x
    ys = node.y
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
    # p=K*(ε₁₁+ε₂₂)
    
    σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
    σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂
    σ₃₃ = ν̄ *(σ₁₁ + σ₂₂)
    p = (σ₁₁ + σ₂₂ + σ₃₃ )/3
    x[i] = xs+α*u₁
    y[i] = ys+α*u₂
    colors[i] = p
end


# points = [[node.x+α*node.d₁ for node in nodes]';[node.y+α*node.d₂ for node in nodes]';zeros(1,nᵤ)]
points = [x';y';zeros(1,nc)]
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
vtk_grid("./vtk/cantilever_GLS_"*poly*"_"*string(ndiv)*"_"*string(ndiv),points,cells) do vtk
    vtk["𝑝"] = colors
end
show(to)
# fig
# end
