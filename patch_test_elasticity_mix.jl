
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫εᵢⱼσᵢⱼdxdy_new, ∫∫σᵢⱼσₖₗdxdy, ∫∫∇σᵢⱼuᵢdxdy, ∫σᵢⱼnⱼuᵢds, ∫σᵢⱼnⱼgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, Hₑ_PlaneStress,∫∫τ∇σᵢⱼ∇σᵢₖdxdy, ∫∫τ∇σᵢⱼ∇σᵢₖdxdy_Taylor,  ∫∫τ∇σᵢⱼ∇σᵢₖdxdy_new, ∫∫σᵢⱼσₖₗdxdy_Taylor

include("import_patchtest.jl")

# nₚ = 49
ndivs = 8
ndiv = 8
# elements, nodes = import_patchtest_mix("msh/patchtest_u_"*string(nₚ)*".msh","./msh/patchtest_"*string(ndiv)*".msh");
elements, nodes = import_patchtest_mix("msh/patchtest_"*string(ndivs)*".msh","./msh/patchtest_"*string(ndiv)*".msh");
# elements, nodes = import_patchtest_mix("msh/patchtest_uniform_"*string(ndivs)*".msh","./msh/patchtest_uniform_"*string(ndiv)*".msh");
nₛ = 3
nₚ = length(nodes)
nₑ = length(elements["Ω"])

Ē = 1.0
ν̄  = 0.3
# ν̄  = 0.4999999
E = Ē/(1.0-ν̄ ^2)
ν = ν̄ /(1.0-ν̄ )
ℎ = 1.0/ndiv
𝐺 = E/(1+ν)/2

set∇𝝭!(elements["Ω"])
set𝝭!(elements["∂Ω"])
set∇𝝭!(elements["Ωᵍ"])
set𝝭!(elements["Γ"])
set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["∂Ωˢ"])

n = 1
# u(x,y) = (1+2*x+3*y)^n
# v(x,y) = (4+5*x+6*y)^n
# ∂u∂x(x,y) = 2*n*(1+2*x+3*y)^abs(n-1)
# ∂u∂y(x,y) = 3*n*(1+2*x+3*y)^abs(n-1)
# ∂v∂x(x,y) = 5*n*(4+5*x+6*y)^abs(n-1)
# ∂v∂y(x,y) = 6*n*(4+5*x+6*y)^abs(n-1)
# ∂²u∂x²(x,y)  = 4*n*(n-1)*(1+2*x+3*y)^abs(n-2)
# ∂²u∂x∂y(x,y) = 6*n*(n-1)*(1+2*x+3*y)^abs(n-2)
# ∂²u∂y²(x,y)  = 9*n*(n-1)*(1+2*x+3*y)^abs(n-2)
# ∂²v∂x²(x,y)  = 25*n*(n-1)*(4+5*x+6*y)^abs(n-2)
# ∂²v∂x∂y(x,y) = 30*n*(n-1)*(4+5*x+6*y)^abs(n-2)
# ∂²v∂y²(x,y)  = 36*n*(n-1)*(4+5*x+6*y)^abs(n-2)

u(x,y) = (x+y)^n
v(x,y) = -(x+y)^n
∂u∂x(x,y) = n*(x+y)^abs(n-1)
∂u∂y(x,y) = n*(x+y)^abs(n-1)
∂v∂x(x,y) = -n*(x+y)^abs(n-1)
∂v∂y(x,y) =- n*(x+y)^abs(n-1)
∂²u∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
∂²u∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²v∂x²(x,y)  = -n*(n-1)*(x+y)^abs(n-2)
∂²v∂x∂y(x,y) = -n*(n-1)*(x+y)^abs(n-2)
∂²v∂y²(x,y)  = -n*(n-1)*(x+y)^abs(n-2)


ε₁₁(x,y) = ∂u∂x(x,y)
ε₂₂(x,y) = ∂v∂y(x,y)
ε₁₂(x,y) = 0.5*(∂u∂y(x,y) + ∂v∂x(x,y))
σ₁₁(x,y) = E/(1-ν^2)*(ε₁₁(x,y) + ν*ε₂₂(x,y))
σ₂₂(x,y) = E/(1-ν^2)*(ν*ε₁₁(x,y) + ε₂₂(x,y))
σ₁₂(x,y) = E/(1+ν)*ε₁₂(x,y)
∂ε₁₁∂x(x,y) = ∂²u∂x²(x,y)
∂ε₁₁∂y(x,y) = ∂²u∂x∂y(x,y)
∂ε₂₂∂x(x,y) = ∂²v∂x∂y(x,y)
∂ε₂₂∂y(x,y) = ∂²v∂y²(x,y)
∂ε₁₂∂x(x,y) = 0.5*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
∂ε₁₂∂y(x,y) = 0.5*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))

∂σ₁₁∂x(x,y) = E/(1-ν^2)*(∂ε₁₁∂x(x,y) + ν*∂ε₂₂∂x(x,y))
∂σ₁₁∂y(x,y) = E/(1-ν^2)*(∂ε₁₁∂y(x,y) + ν*∂ε₂₂∂y(x,y))
∂σ₂₂∂x(x,y) = E/(1-ν^2)*(ν*∂ε₁₁∂x(x,y) + ∂ε₂₂∂x(x,y))
∂σ₂₂∂y(x,y) = E/(1-ν^2)*(ν*∂ε₁₁∂y(x,y) + ∂ε₂₂∂y(x,y))
∂σ₁₂∂x(x,y) = E/(1+ν)*∂ε₁₂∂x(x,y)
∂σ₁₂∂y(x,y) = E/(1+ν)*∂ε₁₂∂y(x,y)
b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)

prescribe!(elements["Ωˢ"],:τ=>(x,y,z)->0.00001*ℎ^2/2/𝐺, index=:𝑔)
prescribe!(elements["Ωˢ"],:ℎ=>(x,y,z)->ℎ, index=:𝑔) 
prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωˢ"],:b₁=>(x,y,z)->b₁(x,y))
prescribe!(elements["Ωˢ"],:b₂=>(x,y,z)->b₂(x,y))
prescribe!(elements["Ωᵍ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵍ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ω"],:b₁=>(x,y,z)->b₁(x,y))
prescribe!(elements["Ω"],:b₂=>(x,y,z)->b₂(x,y))
prescribe!(elements["Γ¹"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ¹"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ²"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ²"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ³"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ³"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ⁴"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ⁴"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ¹"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ¹"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ¹"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ²"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ²"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ²"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ³"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ³"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ³"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ⁴"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))

𝑎 = ∫∫σᵢⱼσₖₗdxdy=>elements["Ωˢ"]
𝑎ᵝ = ∫∫εᵢⱼσᵢⱼdxdy_new=>(elements["Ωˢ"],elements["Ω"])
𝑏 = [
    ∫σᵢⱼnⱼuᵢds=>(elements["∂Ωˢ"],elements["∂Ω"]),
    ∫∫∇σᵢⱼuᵢdxdy=>(elements["Ωˢ"],elements["Ω"]),
]
𝑏ᵅ = ∫σᵢⱼnⱼgᵢds=>(elements["Γˢ"],elements["Γ"])
𝑏ᵝ = ∫∫τ∇σᵢⱼ∇σᵢₖdxdy=>elements["Ωˢ"]
# 𝑏ᵝ = ∫∫τ∇σᵢⱼ∇σᵢₖdxdy_Taylor=>elements["Ωˢ"]
# 𝑏ᵝ = ∫∫τ∇σᵢⱼ∇σᵢₖdxdy_new=>elements["Ωˢ"]
𝑓 = ∫∫vᵢbᵢdxdy=>elements["Ω"]

kᵖᵖ = zeros(3*nₛ*nₑ,3*nₛ*nₑ)
fᵖ = zeros(3*nₛ*nₑ)
kᵖᵘ = zeros(3*nₛ*nₑ,2*nₚ)
kᵖᵘᵅ = zeros(3*nₛ*nₑ,2*nₚ)
fᵖᵅ = zeros(3*nₛ*nₑ)
fᵘ = zeros(2*nₚ)
kᵖᵖˢ = zeros(3*nₛ*nₑ,3*nₛ*nₑ)
fᵖˢ = zeros(3*nₛ*nₑ)

# 𝑎(kᵖᵖ)
𝑎ᵝ(kᵖᵘ)
𝑏(kᵖᵘ)
𝑏ᵅ(kᵖᵘᵅ,fᵖᵅ)
# 𝑏ᵝ(kᵖᵖˢ,fᵖˢ)
# 𝑏ᵝ(kᵖᵖ,fᵖ)
𝑓(fᵘ)

d = [(kᵖᵖ + kᵖᵖˢ) (kᵖᵘ+kᵖᵘᵅ);(kᵖᵘ+kᵖᵘᵅ)' zeros(2*nₚ,2*nₚ)]\[(fᵖ+fᵖˢ+fᵖᵅ);-fᵘ]
# d = [kᵖᵖ  kᵖᵘ;kᵖᵘ' zeros(2*nₚ,2*nₚ)]\[fᵖ;-fᵘ]

d₁ = d[3*nₛ*nₑ+1:2:end]
d₂ = d[3*nₛ*nₑ+2:2:end]
push!(nodes,:d₁=>d₁,:d₂=>d₂)

# 𝐿₂ = L₂(elements["Ωᵍ"])
𝐻ₑ,𝐿₂= Hₑ_PlaneStress(elements["Ωᵍ"])
println(log10(𝐿₂))
println(log10(𝐻ₑ))

dₛ = zeros(3*nₛ*nₑ)


dᵤ = zeros(2*nₚ)
for (i,node) in enumerate(nodes)
    x = node.x
    y = node.y
    dᵤ[2*i-1] = u(x,y)
    dᵤ[2*i]   = v(x,y)
end

dₛ = zeros(3*nₛ*nₑ)
for i in 1:nₑ
    dₛ[3*i-2] = Ē/(1+ν̄ )/(1-2*ν̄ )*((1-ν̄ )*2 + ν̄ *6)
    dₛ[3*i-1] = Ē/(1+ν̄ )/(1-2*ν̄ )*(ν̄ *2 + (1-ν̄ )*6)
    dₛ[3*i]   = Ē/(1+ν̄ )*4
end

# for (i,elm) in enumerate(elements["Ωˢ"])
#     for j in 1:nₛ
#        𝓒ₑ = elm.𝓒
#        dₛ[3*i*j-2] = Ē/(1+ν̄ )/(1-2*ν̄ )*((1-ν̄ )*2* + ν̄ *6)
#        dₛ[3*i*j-1] = Ē/(1+ν̄ )/(1-2*ν̄ )*(ν̄ *2 + (1-ν̄ )*6)
#        dₛ[3*i*j]   = Ē/(1+ν̄ )*4
#      end 
# end

err1 = kᵖᵖ*dₛ + kᵖᵘ*dᵤ - fᵖ
err2 = kᵖᵘ'*dₛ + fᵘ
# err3 = kᵖᵖˢ*dₛ  - fᵖˢ
err4 = kᵖᵘᵅ*dᵤ - fᵖᵅ