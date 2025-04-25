
using TimerOutputs 
using Pardiso
using SparseArrays
using SharedArrays, Distributed
using LinearAlgebra
using WriteVTK
using ApproxOperator
using ApproxOperator.Elasticity: ∫vᵢbᵢdΩ, ∫vᵢtᵢdΓ,∫∫σᵢⱼσₖₗdΩ,∫σᵢⱼnⱼuᵢdΓ,∫∫∇σᵢⱼuᵢdΩ,∫σᵢⱼnⱼgᵢdΓ, Hₑ, ∫∫τ∇σᵢⱼ∇σᵢₖdΩ
         

# addprocs(3)
# println(nprocs())
println(Threads.nthreads())

include("import_patchtest3d.jl")

const to = TimerOutput()
ps = MKLPardisoSolver()

ndiv = 1
ndiv2 = 1
poly = "tet4"
# poly = "hex8"
test = "PatchTest3D"
# test = "block"
@timeit to "import data" begin
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/block_"*string(ndiv)*".msh","./msh/block_"*string(ndiv2)*".msh",ndiv2)
elements, nodes, sp, type, Ω, nodes_c = import_HR_GLS("./msh/"*test*"_"*poly*"_"*string(ndiv)*".msh","./msh/"*test*"_"*poly*"_"*string(ndiv2)*".msh",ndiv)
end

nᵤ = length(nodes)
nc = length(nodes_c)
nₑ = length(elements["Ωᵘ"])
nₑₛ = length(elements["Ωˢ"])
nₛ = 4
ℎ = 1.0/ndiv
E = 240.56839
# ν = 0.5-1e-8
ν = 0.3
P = 80.0
𝐺 = E/(1+ν)/2
β =0.001*ℎ^2/2/𝐺
n₁₁(n₁,n₂,n₃) = n₃ ≈ 1.0 || n₁ ≈ -1.0 ? 1.0 : 0.0
n₂₂(n₁,n₂,n₃) = n₃ ≈ 1.0 || n₂ ≈ -1.0 ? 1.0 : 0.0
n₃₃(n₁,n₂,n₃) = n₃ ≈ -1.0 ? 1.0 : 0.0
prescribe!(elements["Ωˢ"],:τ=>(x,y,z)->β)
prescribe!(elements["Ωˢ"],:ℎ=>(x,y,z)->ℎ) 
prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωˢ"],:b₁=>(x,y,z)->0.0)
prescribe!(elements["Ωˢ"],:b₂=>(x,y,z)->0.0)


n = 2
u(x,y,z) = (x+y+z)^n
v(x,y,z) = (x+y+z)^n
w(x,y,z) = (x+y+z)^n
∂u∂x(x,y,z) = n*(x+y+z)^abs(n-1)
∂u∂y(x,y,z) = n*(x+y+z)^abs(n-1)
∂u∂z(x,y,z) = n*(x+y+z)^abs(n-1)
∂v∂x(x,y,z) = n*(x+y+z)^abs(n-1)
∂v∂y(x,y,z) = n*(x+y+z)^abs(n-1)
∂v∂z(x,y,z) = n*(x+y+z)^abs(n-1)
∂w∂x(x,y,z) = n*(x+y+z)^abs(n-1)
∂w∂y(x,y,z) = n*(x+y+z)^abs(n-1)
∂w∂z(x,y,z) = n*(x+y+z)^abs(n-1)
∂²u∂x²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²u∂y²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²u∂z²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²u∂x∂y(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²u∂x∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²u∂y∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²v∂x²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²v∂y²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²v∂z²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²v∂x∂y(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²v∂x∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²v∂y∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²w∂x²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²w∂y²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²w∂z²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²w∂x∂y(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²w∂x∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²w∂y∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)

ε₁₁(x,y,z) = ∂u∂x(x,y,z)
ε₂₂(x,y,z) = ∂v∂y(x,y,z)
ε₃₃(x,y,z) = ∂w∂z(x,y,z)
ε₁₂(x,y,z) = 0.5*(∂u∂y(x,y,z) + ∂v∂x(x,y,z))
ε₁₃(x,y,z) = 0.5*(∂u∂z(x,y,z) + ∂w∂x(x,y,z))
ε₂₃(x,y,z) = 0.5*(∂v∂z(x,y,z) + ∂w∂y(x,y,z))
∂ε₁₁∂x(x,y,z) = ∂²u∂x²(x,y,z)
∂ε₁₁∂y(x,y,z) = ∂²u∂x∂y(x,y,z)
∂ε₁₁∂z(x,y,z) = ∂²u∂x∂z(x,y,z)
∂ε₂₂∂x(x,y,z) = ∂²v∂x∂y(x,y,z)
∂ε₂₂∂y(x,y,z) = ∂²v∂y²(x,y,z)
∂ε₂₂∂z(x,y,z) = ∂²v∂y∂z(x,y,z)
∂ε₃₃∂x(x,y,z) = ∂²w∂x∂z(x,y,z)
∂ε₃₃∂y(x,y,z) = ∂²w∂y∂z(x,y,z)
∂ε₃₃∂z(x,y,z) = ∂²w∂z²(x,y,z)
∂ε₁₂∂x(x,y,z) = 0.5*(∂²u∂x∂y(x,y,z) + ∂²v∂x²(x,y,z))
∂ε₁₂∂y(x,y,z) = 0.5*(∂²u∂y²(x,y,z) + ∂²v∂x∂y(x,y,z))
∂ε₁₂∂z(x,y,z) = 0.5*(∂²u∂y∂z(x,y,z) + ∂²v∂x∂z(x,y,z))
∂ε₁₃∂x(x,y,z) = 0.5*(∂²u∂x∂z(x,y,z) + ∂²w∂x²(x,y,z))
∂ε₁₃∂y(x,y,z) = 0.5*(∂²u∂y∂z(x,y,z) + ∂²w∂x∂y(x,y,z))
∂ε₁₃∂z(x,y,z) = 0.5*(∂²u∂z²(x,y,z) + ∂²w∂x∂z(x,y,z))
∂ε₂₃∂x(x,y,z) = 0.5*(∂²v∂x∂z(x,y,z) + ∂²w∂x∂y(x,y,z))
∂ε₂₃∂y(x,y,z) = 0.5*(∂²v∂y∂z(x,y,z) + ∂²w∂y²(x,y,z))
∂ε₂₃∂z(x,y,z) = 0.5*(∂²v∂z²(x,y,z) + ∂²w∂y∂z(x,y,z))
σ₁₁(x,y,z) = E/(1+ν)/(1-2*ν)*((1-ν)*ε₁₁(x,y,z) + ν*ε₂₂(x,y,z) + ν*ε₃₃(x,y,z))
σ₂₂(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y,z) + (1-ν)*ε₂₂(x,y,z) + ν*ε₃₃(x,y,z))
σ₃₃(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y,z) + ν*ε₂₂(x,y,z) + (1-ν)*ε₃₃(x,y,z))
σ₁₂(x,y,z) = E/(1+ν)*ε₁₂(x,y,z)
σ₁₃(x,y,z) = E/(1+ν)*ε₁₃(x,y,z)
σ₂₃(x,y,z) = E/(1+ν)*ε₂₃(x,y,z)
𝑝(x,y,z) = (σ₁₁(x,y,z)+σ₂₂(x,y,z)+σ₃₃(x,y,z))/3
∂σ₁₁∂x(x,y,z) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂x(x,y,z) + ν*∂ε₂₂∂x(x,y,z) + ν*∂ε₃₃∂x(x,y,z))
∂σ₁₁∂y(x,y,z) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂y(x,y,z) + ν*∂ε₂₂∂y(x,y,z) + ν*∂ε₃₃∂y(x,y,z))
∂σ₁₁∂z(x,y,z) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂z(x,y,z) + ν*∂ε₂₂∂z(x,y,z) + ν*∂ε₃₃∂z(x,y,z))
∂σ₂₂∂x(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y,z) + (1-ν)*∂ε₂₂∂x(x,y,z) + ν*∂ε₃₃∂x(x,y,z))
∂σ₂₂∂y(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y,z) + (1-ν)*∂ε₂₂∂y(x,y,z) + ν*∂ε₃₃∂y(x,y,z))
∂σ₂₂∂z(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂z(x,y,z) + (1-ν)*∂ε₂₂∂z(x,y,z) + ν*∂ε₃₃∂z(x,y,z))
∂σ₃₃∂x(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y,z) + ν*∂ε₂₂∂x(x,y,z) + (1-ν)*∂ε₃₃∂x(x,y,z))
∂σ₃₃∂y(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y,z) + ν*∂ε₂₂∂y(x,y,z) + (1-ν)*∂ε₃₃∂y(x,y,z))
∂σ₃₃∂z(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂z(x,y,z) + ν*∂ε₂₂∂z(x,y,z) + (1-ν)*∂ε₃₃∂z(x,y,z))
∂σ₁₂∂x(x,y,z) = E/(1+ν)*∂ε₁₂∂x(x,y,z)
∂σ₁₂∂y(x,y,z) = E/(1+ν)*∂ε₁₂∂y(x,y,z)
∂σ₁₂∂z(x,y,z) = E/(1+ν)*∂ε₁₂∂z(x,y,z)
∂σ₁₃∂x(x,y,z) = E/(1+ν)*∂ε₁₃∂x(x,y,z)
∂σ₁₃∂y(x,y,z) = E/(1+ν)*∂ε₁₃∂y(x,y,z)
∂σ₁₃∂z(x,y,z) = E/(1+ν)*∂ε₁₃∂z(x,y,z)
∂σ₂₃∂x(x,y,z) = E/(1+ν)*∂ε₂₃∂x(x,y,z)
∂σ₂₃∂y(x,y,z) = E/(1+ν)*∂ε₂₃∂y(x,y,z)
∂σ₂₃∂z(x,y,z) = E/(1+ν)*∂ε₂₃∂z(x,y,z)
b₁(x,y,z) = - ∂σ₁₁∂x(x,y,z) - ∂σ₁₂∂y(x,y,z) - ∂σ₁₃∂z(x,y,z)
b₂(x,y,z) = - ∂σ₁₂∂x(x,y,z) - ∂σ₂₂∂y(x,y,z) - ∂σ₂₃∂z(x,y,z)
b₃(x,y,z) = - ∂σ₁₃∂x(x,y,z) - ∂σ₂₃∂y(x,y,z) - ∂σ₃₃∂z(x,y,z)

prescribe!(elements["Ωᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵘ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωˢ"],:b₁=>b₁)
prescribe!(elements["Ωˢ"],:b₂=>b₂)
prescribe!(elements["Ωˢ"],:b₃=>b₃)

prescribe!(elements["Ωᵘ"],:b₁=>b₁)
prescribe!(elements["Ωᵘ"],:b₂=>b₂)
prescribe!(elements["Ωᵘ"],:b₃=>b₃)

prescribe!(elements["Γᵍᵘ"],:α=>(x,y,z)->1e12*E)
prescribe!(elements["Γᵍᵘ"],:g₁=>u)
prescribe!(elements["Γᵍᵘ"],:g₂=>v)
prescribe!(elements["Γᵍᵘ"],:g₃=>w)
prescribe!(elements["Γᵍᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₃₃=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍᵘ"],:n₁₃=>(x,y,z)->0.0)
prescribe!(elements["Γᵍᵘ"],:n₂₃=>(x,y,z)->0.0)

prescribe!(elements["Γᵍᵘ"],:t₁=>(x,y,z,n₁,n₂,n₃)->σ₁₁(x,y,z)*n₁+σ₁₂(x,y,z)*n₂+σ₁₃(x,y,z)*n₃)
prescribe!(elements["Γᵍᵘ"],:t₂=>(x,y,z,n₁,n₂,n₃)->σ₁₂(x,y,z)*n₁+σ₂₂(x,y,z)*n₂+σ₂₃(x,y,z)*n₃)
prescribe!(elements["Γᵍᵘ"],:t₃=>(x,y,z,n₁,n₂,n₃)->σ₁₃(x,y,z)*n₁+σ₂₃(x,y,z)*n₂+σ₃₃(x,y,z)*n₃)

prescribe!(elements["Ωᵍᵘ"],:u₁=>u)
prescribe!(elements["Ωᵍᵘ"],:u₂=>v)
prescribe!(elements["Ωᵍᵘ"],:u₃=>w)
prescribe!(elements["Ωᵍᵘ"],:∂u₁∂x=>∂u∂x)
prescribe!(elements["Ωᵍᵘ"],:∂u₁∂y=>∂u∂y)
prescribe!(elements["Ωᵍᵘ"],:∂u₁∂z=>∂u∂z)
prescribe!(elements["Ωᵍᵘ"],:∂u₂∂x=>∂v∂x)
prescribe!(elements["Ωᵍᵘ"],:∂u₂∂y=>∂v∂y)
prescribe!(elements["Ωᵍᵘ"],:∂u₂∂z=>∂v∂z)
prescribe!(elements["Ωᵍᵘ"],:∂u₃∂x=>∂w∂x)
prescribe!(elements["Ωᵍᵘ"],:∂u₃∂y=>∂w∂y)
prescribe!(elements["Ωᵍᵘ"],:∂u₃∂z=>∂w∂z)


𝑎 =∫∫σᵢⱼσₖₗdΩ=>elements["Ωˢ"]

𝑏 = [
    ∫σᵢⱼnⱼuᵢdΓ=>(elements["∂Ωˢ"],elements["∂Ωᵘ"]),
    ∫∫∇σᵢⱼuᵢdΩ=>(elements["Ωˢ"],elements["Ωᵘ"]),
   
    ]

𝑏ᵅ = ∫σᵢⱼnⱼgᵢdΓ=>(elements["Γᵍˢ"],elements["Γᵍᵘ"])


𝑏ᵝ = ∫∫τ∇σᵢⱼ∇σᵢₖdΩ=>elements["Ωˢ"]


𝑓 = [
    # ∫vᵢtᵢdΓ=>elements["Γᵗ"]∪elements["Γʳ"]∪elements["Γᵍᵘ"],
    # ∫vᵢtᵢdΓ=>elements["Γᵗ"]∪elements["Γʳ"],
    # ∫vᵢtᵢdΓ=>elements["Γᵗ"],
    ∫vᵢbᵢdΩ=>elements["Ωᵘ"],
]

kˢˢ = zeros(6*nₛ*nₑ,6*nₛ*nₑ)
kˢᵘ = zeros(6*nₛ*nₑ,3*nᵤ)
kˢᵘⁿ  = zeros(6*nₛ*nₑ,3*nᵤ)
fˢ = zeros(6*nₛ*nₑ)
fˢⁿ = zeros(6*nₛ*nₑ)
fᵘ = zeros(3*nᵤ)


# kᵘᵘ = SharedMatrix{Float64}(3*nᵤ,3*nᵤ)
# kᵖᵖ = SharedMatrix{Float64}(nₚ,nₚ)
# kᵖᵘ = SharedMatrix{Float64}(nₚ,3*nᵤ)
# fᵖ  = SharedVector{Float64}(nₚ)
# fᵘ  = SharedVector{Float64}(3*nᵤ)

@timeit to "assembly" begin
    𝑎(kˢˢ)
    𝑏(kˢᵘ)
    𝑏ᵅ(kˢᵘⁿ,fˢⁿ)
    # 𝑏ᵝ(kˢˢ,fˢ)
    𝑓(fᵘ)
end



d = [kˢˢ (kˢᵘ+kˢᵘⁿ);(kˢᵘ+kˢᵘⁿ)' zeros(3*nᵤ,3*nᵤ)]\[(-fˢ+fˢⁿ);-fᵘ]
d₁ = d[6*nₛ*nₑ+1:3:end]
d₂ = d[6*nₛ*nₑ+2:3:end]
d₃ = d[6*nₛ*nₑ+3:3:end]
dₛ₁₁ = d[1:6:6*nₛ*nₑₛ]
dₛ₂₂ = d[2:6:6*nₛ*nₑₛ]
dₛ₃₃ = d[3:6:6*nₛ*nₑₛ]
dₛ₁₂ = d[4:6:6*nₛ*nₑₛ]
dₛ₂₃ = d[5:6:6*nₛ*nₑₛ]
dₛ₁₃ = d[6:6:6*nₛ*nₑₛ]


push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
for elm in elements["Ωˢ"]
    𝓒ₚ = elm.𝓒
    𝓖 = elm.𝓖
        push!(𝓒ₚ,:dₛ₁₁=>dₛ₁₁,:dₛ₂₂=>dₛ₂₂,:dₛ₃₃=>dₛ₃₃,:dₛ₁₂=>dₛ₁₂,:dₛ₂₃=>dₛ₂₃,:dₛ₁₃=>dₛ₁₃)
end
Hₑ_𝒖, L₂_𝒖 = Hₑ(elements["Ωᵍᵘ"])


println(log10(L₂_𝒖))
println(log10(Hₑ_𝒖))

# colors = zeros(nᵤ)
# 𝗠 = zeros(10)
# for (i,node) in enumerate(nodes)
#     x = node.x
#     y = node.y
#     z = node.z
#     indices = sp(x,y,z)
#     ni = length(indices)
#     𝓒 = [nodes_p[i] for i in indices]
#     data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[z]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
#     ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
#     𝓖 = [ξ]
#     a = type(𝓒,𝓖)
#     set𝝭!(a)
#     p = 0.0
#     N = ξ[:𝝭]
#     for (k,xₖ) in enumerate(𝓒)
#         p += N[k]*xₖ.p
#     end
#     colors[i] = p
# end
# α = 1.0
# points = [[node.x+α*node.u₁ for node in nodes]';[node.y+α*node.u₂ for node in nodes]';[node.z+α*node.u₃ for node in nodes]']
# # cells = [MeshCell(VTKCellTypes.VTK_TETRA,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# vtk_grid("./vtk/block_"*poly*"_"*string(ndiv)*"_"*string(nₚ),points,cells) do vtk
#     vtk["u"] = (𝑢₁,𝑢₂,𝑢₃, VTKPointData())
#     vtk["𝑝",] = (colors, VTKCelldData())
# end

# println(nodes[5])



dₛ = zeros(6*nₛ*nₑ)


dᵤ = zeros(3*nᵤ)
for (i,node) in enumerate(nodes)
    x = node.x
    y = node.y
    z = node.z
    dᵤ[3*i-2] = u(x,y,z)
    dᵤ[3*i-1] = v(x,y,z)
    dᵤ[3*i]   = w(x,y,z)
end

dₛ = zeros(6*nₛ*nₑ)
for i in 1:nₑ
    dₛ[6*i-5] = E/(1+ν)/(1-2*ν)*((1-ν ) + ν + ν)
    dₛ[6*i-4] = E/(1+ν)/(1-2*ν)*(ν + (1-ν )+ν)
    dₛ[6*i-3] = E/(1+ν)/(1-2*ν)*(ν + ν+ (1-ν ))
    dₛ[6*i-2] = E/(1+ν)
    dₛ[6*i-1] = E/(1+ν)
    dₛ[6*i]   = E/(1+ν)
   

end
err1 = kˢˢ*dₛ + (kˢᵘ+kˢᵘⁿ)*dᵤ - fˢ
err2 = (kˢᵘ+kˢᵘⁿ)' *dₛ + fᵘ

err4 = kˢᵘⁿ*dᵤ - fˢⁿ
show(to)