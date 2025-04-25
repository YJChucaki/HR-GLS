
using TimerOutputs 
using SparseArrays, Pardiso
using CairoMakie, LinearAlgebra, WriteVTK
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫εᵢⱼσᵢⱼdxdy, ∫∫εᵢⱼσᵢⱼdxdy_PlaneStrian,∫σᵢⱼnⱼgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, Hₑ_PlaneStress,  ∫vᵢgᵢds

include("import_cantilever.jl")


const to = TimerOutput()
ps = MKLPardisoSolver()
n = [4 8 16 32 ]
for i in 1:4
ndiv = n[i]
ndiv2 = n[i]
# ndiv = 8                                                                                


poly = "tri3"
# poly = "nonuniform"
# poly = "tri6"
# poly = "quad"
@timeit to "import data" begin

# elements, nodes = import_MF_Gauss("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh")
elements, nodes, sp, type , nodes_c= import_MF_Gauss("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh")

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
Ē = 1.0
ν̄  = 0.3
# ν̄  = 0.5-1e-8
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

σ₁₁(x,y) = -P*(L-x)*y/I
σ₂₂(x,y) = 0.0
σ₁₂(x,y) = P/2/I*(D^2/4-y^2)

prescribe!(elements["Ω"],:E=>(x,y,z)->E,index=:𝑔)
prescribe!(elements["Ω"],:ν=>(x,y,z)->ν,index=:𝑔)
prescribe!(elements["Ωᵍ"],:E=>(x,y,z)->E,index=:𝑔)
prescribe!(elements["Ωᵍ"],:ν=>(x,y,z)->ν,index=:𝑔)
prescribe!(elements["Γᵍ"],:E=>(x,y,z)->E)
prescribe!(elements["Γᵍ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->v(x,y))

prescribe!(elements["Γᵍ"],:α=>(x,y,z)->1e6*Ē,index=:𝑔)
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



α = 0.0
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


points = [[node.x+α*node.d₁ for node in nodes]';[node.y+α*node.d₂ for node in nodes]';zeros(1,nₚ)]
# points = [x';y';zeros(1,nc)]
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["ΩC"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
vtk_grid("./vtk/cantilever_Gauss_"*poly*"_"*string(ndiv)*"_"*string(ndiv),points,cells) do vtk
    vtk["𝑝"] = colors
 
end
show(to)
# show(to)
end
