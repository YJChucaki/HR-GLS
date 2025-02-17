using SparseArrays, Pardiso
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric,∫∫τ∇σᵢⱼ∇σᵢₖdxdy,∫∫σᵢⱼσₖₗdxdy, ∫∫∇σᵢⱼuᵢdxdy, ∫σᵢⱼnⱼuᵢds, ∫σᵢⱼnⱼgᵢds,  ∫∫τ∇σᵢⱼ∇σᵢₖdxdy_Real

include("import_patchtest.jl")

ndiv = 8
ndiv2 = 8

# nₚ = 60
# elements, nodes, nodes_p = import_patchtest_elasticity_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_c_"*string(nₚ)*".msh")
# elements, nodes, nodes_p = import_patchtest_elasticity_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndiv)*".msh")
elements, nodes = import_patchtest_mix("msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndiv2)*".msh");
nₑ = length(elements["Ωᵘ"])
nₛ = 3*nₑ
nᵤ = length(nodes)

E = 1.0
ν = 0.3
# ν = 0.4999999
ℎ = 1.0/ndiv
𝐺 = E/(1+ν)/2

set𝝭!(elements["Ωᵘ"])
set𝝭!(elements["∂Ωᵘ"])
set∇𝝭!(elements["Ωᵍᵘ"])
set𝝭!(elements["Γᵘ"])
set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["∂Ωˢ"])

n = 2
# u(x,y) = x^5
# v(x,y) = - 5*x^4*y
# ∂u∂x(x,y) = 5*x^4
# ∂u∂y(x,y) = 0.0
# ∂v∂x(x,y) = -20*x^3*y
# ∂v∂y(x,y) = -5*x^4
# ∂²u∂x²(x,y)  = 20*x^3
# ∂²u∂x∂y(x,y) = 0.0
# ∂²u∂y²(x,y)  = 0.0
# ∂²v∂x²(x,y)  = - 60*x^2*y
# ∂²v∂x∂y(x,y) = - 20*x^3
# ∂²v∂y²(x,y)  = 0.0
u(x,y) = (1+2*x+3*y)^n
v(x,y) = (4+5*x+6*y)^n
∂u∂x(x,y) = 2*n*(1+2*x+3*y)^abs(n-1)
∂u∂y(x,y) = 3*n*(1+2*x+3*y)^abs(n-1)
∂v∂x(x,y) = 5*n*(4+5*x+6*y)^abs(n-1)
∂v∂y(x,y) = 6*n*(4+5*x+6*y)^abs(n-1)
∂²u∂x²(x,y)  = 4*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²u∂x∂y(x,y) = 6*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²u∂y²(x,y)  = 9*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²v∂x²(x,y)  = 25*n*(n-1)*(4+5*x+6*y)^abs(n-2)
∂²v∂x∂y(x,y) = 30*n*(n-1)*(4+5*x+6*y)^abs(n-2)
∂²v∂y²(x,y)  = 36*n*(n-1)*(4+5*x+6*y)^abs(n-2)

ε₁₁(x,y) = ∂u∂x(x,y)
ε₂₂(x,y) = ∂v∂y(x,y)
ε₁₂(x,y) = 0.5*(∂u∂y(x,y) + ∂v∂x(x,y))
σ₁₁(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*ε₁₁(x,y) + ν*ε₂₂(x,y))
σ₂₂(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + (1-ν)*ε₂₂(x,y))
σ₃₃(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + ν*ε₂₂(x,y))
σ₁₂(x,y) = E/(1+ν)*ε₁₂(x,y)
∂ε₁₁∂x(x,y) = ∂²u∂x²(x,y)
∂ε₁₁∂y(x,y) = ∂²u∂x∂y(x,y)
∂ε₂₂∂x(x,y) = ∂²v∂x∂y(x,y)
∂ε₂₂∂y(x,y) = ∂²v∂y²(x,y)
∂ε₁₂∂x(x,y) = 0.5*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
∂ε₁₂∂y(x,y) = 0.5*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))

∂σ₁₁∂x(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂x(x,y) + ν*∂ε₂₂∂x(x,y))
∂σ₁₁∂y(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂y(x,y) + ν*∂ε₂₂∂y(x,y))
∂σ₂₂∂x(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y) + (1-ν)*∂ε₂₂∂x(x,y))
∂σ₂₂∂y(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y) + (1-ν)*∂ε₂₂∂y(x,y))
∂σ₁₂∂x(x,y) = E/(1+ν)*∂ε₁₂∂x(x,y)
∂σ₁₂∂y(x,y) = E/(1+ν)*∂ε₁₂∂y(x,y)
b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)


prescribe!(elements["Ωˢ"],:τ=>(x,y,z)->-1*ℎ^2/2/𝐺, index=:𝑔)

prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν, index=:𝑔)

prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωᵘ"],:b₁=>(x,y,z)->b₁(x,y))
prescribe!(elements["Ωᵘ"],:b₂=>(x,y,z)->b₂(x,y))

prescribe!(elements["Ωˢ"],:b₁=>(x,y,z)->b₁(x,y))
prescribe!(elements["Ωˢ"],:b₂=>(x,y,z)->b₂(x,y))
prescribe!(elements["Γ¹ᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ¹ᵘ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ²ᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ²ᵘ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ³ᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ³ᵘ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ⁴ᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ⁴ᵘ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ¹ᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ¹ᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ¹ᵘ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ²ᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ²ᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ²ᵘ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ³ᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ³ᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ³ᵘ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ⁴ᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴ᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴ᵘ"],:n₁₂=>(x,y,z)->0.0)

prescribe!(elements["Ωᵍᵘ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍᵘ"],:v=>(x,y,z)->v(x,y))

prescribe!(elements["Ωᵍᵘ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))






𝑎 = ∫∫σᵢⱼσₖₗdxdy=>elements["Ωˢ"]
# 𝑎 = ∫∫σᵢⱼσₖₗdxdy_PlaneStrian=>elements["Ωˢ"]
𝑏 = [
    ∫σᵢⱼnⱼuᵢds=>(elements["∂Ωˢ"],elements["∂Ωᵘ"]),
    ∫∫∇σᵢⱼuᵢdxdy=>(elements["Ωˢ"],elements["Ωᵘ"]),
]
𝑏ᵅ = ∫σᵢⱼnⱼgᵢds=>(elements["Γˢ"],elements["Γᵘ"])

𝑏ᵝ = ∫∫τ∇σᵢⱼ∇σᵢₖdxdy=>elements["Ωˢ"]
𝑓 = ∫∫vᵢbᵢdxdy=>elements["Ωᵘ"]

kᵖᵖ = zeros(3*nₛ,3*nₛ)
fᵖ = zeros(3*nₛ)
kᵖᵘ = zeros(3*nₛ,2*nᵤ)
fᵘ = zeros(2*nᵤ)

𝑎(kᵖᵖ)
𝑏(kᵖᵘ)
𝑏ᵅ(kᵖᵘ,fᵖ)
# 𝑏ᵝ(kᵖᵖ,fᵖ)
𝑓(fᵘ)


k = [kᵖᵖ kᵖᵘ;kᵖᵘ' zeros(2*nᵤ,2*nᵤ)]
f = [fᵖ;-fᵘ]

d = k\f
d₁ = d[3*nₛ+1:2:end]
d₂ = d[3*nₛ+2:2:end]
push!(nodes,:d₁=>d₁,:d₂=>d₂)

# 𝐿₂ = L₂(elements["Ωᵍᵘ"])
𝐿₂, 𝐻ₑ = Hₑ_PlaneStress(elements["Ωᵍᵘ"])
println(𝐿₂)
println(𝐻ₑ)