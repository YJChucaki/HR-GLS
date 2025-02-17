
using TimerOutputs 
using SparseArrays, Pardiso
using CairoMakie
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫εᵢⱼσᵢⱼdxdy, ∫∫εᵢⱼσᵢⱼdxdy_PlaneStrian,∫σᵢⱼnⱼgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, Hₑ_PlaneStress,  ∫vᵢgᵢds

include("import_cantilever.jl")


const to = TimerOutput()
# ps = MKLPardisoSolver()
n = [ 2 4 8 16 32 ]
for i in 1:4
ndiv = n[i]
ndiv2 = n[i]
# ndiv = 4


poly = "tri3"
# poly = "nonuniform"
# poly = "tri6"
# poly = "quad"
@timeit to "import data" begin

# elements, nodes = import_MF_Gauss("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh")
elements, nodes = import_MF_Gauss("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh")

end

nₑ = length(elements["Ω"])

nₚ = length(nodes)
# nₚ = length(nodes_p)
# nₚ = length(nodes)

L = 48.0
D = 12.0
P = 1000
ℎ = D/ndiv

# Ē = 3e6
# ν̄  = 0.3
E = 3e6
# ν = 0.3
ν = 0.5-1e-7
Ē = E/(1.0-ν^2)
ν̄ = ν/(1.0-ν)
I = D^3/12
EI = Ē*I


u(x,y) = -P*y/6/EI*((6*L-3*x)*x + (2+ν̄)*(y^2-D^2/4))
v(x,y) = P/6/EI*(3*ν̄*y^2*(L-x) + (4+5*ν̄)*D^2*x/4 + (3*L-x)*x^2)
∂u∂x(x,y) = -P/EI*(L-x)*y
∂u∂y(x,y) = -P/6/EI*((6*L-3*x)*x + (2+ν̄)*(3*y^2-D^2/4))
∂v∂x(x,y) = P/6/EI*((6*L-3*x)*x - 3*ν̄*y^2 + (4+5*ν̄)*D^2/4)
∂v∂y(x,y) = P/EI*(L-x)*y*ν̄

σ₁₁(x,y) = -P*(L-x)*y/I
σ₂₂(x,y) = 0.0
σ₁₂(x,y) = P/2/I*(D^2/4-y^2)

prescribe!(elements["Ω"],:E=>(x,y,z)->Ē,index=:𝑔)
prescribe!(elements["Ω"],:ν=>(x,y,z)->ν̄,index=:𝑔)
prescribe!(elements["Ωᵍ"],:E=>(x,y,z)->E,index=:𝑔)
prescribe!(elements["Ωᵍ"],:ν=>(x,y,z)->ν,index=:𝑔)
prescribe!(elements["Γᵍ"],:E=>(x,y,z)->E)
prescribe!(elements["Γᵍ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->v(x,y))

prescribe!(elements["Γᵍ"],:α=>(x,y,z)->1e3*E,index=:𝑔)
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))

# 𝑎 = ∫∫εᵢⱼσᵢⱼdxdy_PlaneStrian=>elements["Ω"]
𝑎 = ∫∫εᵢⱼσᵢⱼdxdy=>elements["Ω"]
𝑎ᵅ = ∫vᵢgᵢds=>elements["Γᵍ"]
# 𝑎ᵅ = ∫σᵢⱼnⱼgᵢds=>elements["Γᵍ"]
𝑓 = ∫vᵢtᵢds=>elements["Γᵗ"]


k = zeros(2*nₚ,2*nₚ)
kᵍ = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

@timeit to "assembly matrix" begin

𝑎(k)
𝑎ᵅ(kᵍ,f)
𝑓(f)
end


d = (k+kᵍ)\f
d₁ = d[1:2:end]
d₂ = d[2:2:end]
push!(nodes,:d₁=>d₁,:d₂=>d₂)


𝐻ₑ, 𝐿₂ = Hₑ_PlaneStress(elements["Ωᵍ"])

println(log10(𝐿₂))
println(log10(𝐻ₑ))

show(to)
end
