
using TimerOutputs 
using SparseArrays, Pardiso
using CairoMakie, WriteVTK
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫qpdxdy, ∫∫sᵢⱼsᵢⱼdxdy, ∫∫p∇udxdy, ∫∫sᵢⱼεᵢⱼdxdy, ∫pnᵢgᵢds, ∫sᵢⱼnⱼgᵢds, ∫∫τ∇q∇pdxdy, ∫∫τ∇sᵢⱼ∇sᵢₖdxdy, ∫∫τ∇sᵢⱼ∇pdxdy, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric, ∫∫∇sᵢⱼuᵢdxdy, ∫∫∇puᵢdxdy ,∫pnᵢuᵢds,∫sᵢⱼnⱼuᵢds 
using  Printf
include("import_cantilever.jl")
include("wirteVTK.jl")
const to = TimerOutput()
ps = MKLPardisoSolver()

ndiv = 4
ndiv2 = 4
# nₚ =  40
poly1 = "tri3"
# poly1 = "tri6"
# poly1 = "quad8"
# poly1 = "quad"
@timeit to "import data" begin
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_c_"*string(nₚ)*".msh")
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly*"_"*string(ndiv2)*".msh")
# elements, nodes,nodes_p,  sp, type, Ω = import_HR_GLS("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly*"_"*string(ndiv2)*".msh")
elements, nodes,  sp, type,Ω, nodes_c = import_HR_GLS_MPP("./msh/cantilever_"*poly1*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly1*"_"*string(ndiv2)*".msh")
end

nₑ = length(elements["Ωᵘ"])
nₛ = 3*nₑ
nᵤ = length(nodes)
nₚ = 1*nₑ
nₑₚ = length(Ω)
# nₚ = length(nodes_p)
# nₚ = length(nodes)

L = 48.0
D = 12.0
P = 1000
ℎ = D/ndiv
E = 3e6
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

u(x,y) = -P*y/6/EI*((6*L-3*x)*x + (2+ν̄)*(y^2-D^2/4))
v(x,y) = P/6/EI*(3*ν̄*y^2*(L-x) + (4+5*ν̄)*D^2*x/4 + (3*L-x)*x^2)
∂u∂x(x,y) = -P/EI*(L-x)*y
∂u∂y(x,y) = -P/6/EI*((6*L-3*x)*x + (2+ν̄)*(3*y^2-D^2/4))
∂v∂x(x,y) = P/6/EI*((6*L-3*x)*x - 3*ν̄*y^2 + (4+5*ν̄)*D^2/4)
∂v∂y(x,y) = P/EI*(L-x)*y*ν̄

ε₁₁(x,y) = ∂u∂x(x,y)
ε₂₂(x,y) = ∂v∂y(x,y)
ε₁₂(x,y) = ∂u∂y(x,y) + ∂v∂x(x,y)
σ₁₁(x,y) = -P*(L-x)*y/I
σ₂₂(x,y) = 0.0
σ₃₃(x,y) = Cᵢᵢⱼⱼ*ε₁₁(x,y) + Cᵢᵢⱼⱼ*ε₂₂(x,y)
σ₁₂(x,y) = P/2/I*(D^2/4-y^2)

β = 0.01*ℎ^2/2/𝐺
prescribe!(elements["Ωˢ"],:τ=>(x,y,z)->β, index=:𝑔)
prescribe!(elements["Ωᵖ"],:τ=>(x,y,z)->β, index=:𝑔)
prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωᵖ"],:b₁=>(x,y,z)->0.0)
prescribe!(elements["Ωᵖ"],:b₂=>(x,y,z)->0.0)
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

prescribe!(elements["Γᵍᵘᵖ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍᵘᵖ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γᵍᵘᵖ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘᵖ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘᵖ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ωᵍᵘᵖ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍᵘᵖ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍᵘᵖ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍᵘᵖ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍᵘᵖ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍᵘᵖ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
prescribe!(elements["Ωᵍᵖ"],:p=>(x,y,z)->(σ₁₁(x,y)+σ₂₂(x,y)+σ₃₃(x,y))/3)


𝑎ˢ = ∫∫sᵢⱼsᵢⱼdxdy=>elements["Ωˢ"]
𝑎ᵖ = ∫∫qpdxdy=>elements["Ωᵖ"]
𝑏ˢ = ∫∫∇sᵢⱼuᵢdxdy=>(elements["Ωˢ"],elements["Ωᵘ"])
𝑏ᵖ = ∫∫∇puᵢdxdy=>(elements["Ωᵖ"],elements["Ωᵘ"])
# 𝑏ˢ = ∫∫sᵢⱼεᵢⱼdxdy=>(elements["Ωˢ"],elements["Ωᵘ"])
# 𝑏ᵖ = ∫∫p∇udxdy=>(elements["Ωᵖ"],elements["Ωᵘ"])
𝑏ˢᵅ = ∫sᵢⱼnⱼgᵢds=>(elements["Γᵍˢ"],elements["Γᵍᵘ"])
𝑏ˢᵇ = ∫sᵢⱼnⱼuᵢds=>(elements["∂Ωˢ"],elements["∂Ωᵘ"])
𝑏ᵖᵅ = ∫pnᵢgᵢds=>(elements["Γᵍᵖ"],elements["Γᵍᵘ"])
𝑏ᵖᵇ = ∫pnᵢuᵢds=>(elements["∂Ωᵖ"],elements["∂Ωᵘ"])
𝑏ˢᵝ = ∫∫τ∇sᵢⱼ∇sᵢₖdxdy=>elements["Ωˢ"]
𝑏ᵖᵝ = ∫∫τ∇q∇pdxdy=>elements["Ωᵖ"]
𝑏ˢᵖᵝ = ∫∫τ∇sᵢⱼ∇pdxdy=>(elements["Ωˢ"],elements["Ωᵖ"])
𝑓 = ∫vᵢtᵢds=>elements["Γᵗ"]
# 𝑓 = [
#     ∫vᵢtᵢds=>elements["Γᵗ"]∪elements["Γʳ"],
#     ∫∫vᵢbᵢdxdy=>elements["Ωᵘ"]
# ]

kˢˢ = zeros(4*nₛ,4*nₛ)
kᵖᵖ = zeros(nₚ,nₚ)
kˢᵖ = zeros(4*nₛ,nₚ)
kˢᵘ = zeros(4*nₛ,2*nᵤ)
kᵖᵘ = zeros(nₚ,2*nᵤ)
fˢ = zeros(4*nₛ)
fᵖ = zeros(nₚ)
fᵘ = zeros(2*nᵤ)


@timeit to "assembly" begin
𝑎ˢ(kˢˢ)
𝑎ᵖ(kᵖᵖ)
𝑏ˢ(kˢᵘ)
𝑏ᵖ(kᵖᵘ)
𝑏ˢᵅ(kˢᵘ,fˢ)
𝑏ˢᵇ(kˢᵘ)
𝑏ᵖᵅ(kᵖᵘ,fᵖ)
𝑏ᵖᵇ(kᵖᵘ)

𝑏ˢᵝ(kˢˢ,fˢ)
𝑏ᵖᵝ(kᵖᵖ,fᵖ)
𝑏ˢᵖᵝ(kˢᵖ)
𝑓(fᵘ)
end

k = [zeros(2*nᵤ,2*nᵤ) kᵖᵘ' kˢᵘ';kᵖᵘ kᵖᵖ kˢᵖ';kˢᵘ kˢᵖ kˢˢ]
f = [-fᵘ;fᵖ;fˢ]
d = zeros(2*nᵤ+nₚ+4*nₛ)
d = k\f


𝑢₁ = d[1:2:2*nᵤ]
𝑢₂ = d[2:2:2*nᵤ]
𝑝 = d[2*nᵤ+1:2*nᵤ+nₚ]
push!(nodes,:d₁=>𝑢₁)
push!(nodes,:d₂=>𝑢₂)

for elm in elements["Ωᵍᵖ"]
    𝓒ₚ = elm.𝓒
    push!(𝓒ₚ,:p=>𝑝)
    end

# eval(VTK_HR_MPP)
@timeit to "compute error" begin
Hₑ_𝒖, L₂_𝒖 = Hₑ_PlaneStress(elements["Ωᵍᵘ"])
Hₑ_dev = Hₑ_PlaneStrain_Deviatoric(elements["Ωᵍᵘ"])
L₂_𝑝 = L₂𝑝(elements["Ωᵍᵖ"])
end

println(log10(L₂_𝒖))
println(log10(Hₑ_𝒖))
println(log10(Hₑ_dev))
println(log10(L₂_𝑝))

α = 1.0

colors = zeros(nᵤ)
x = zeros(nᵤ)
y = zeros(nᵤ)
# type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
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
    p=K*(ε₁₁+ε₂₂)
    x[i] = xs+α*u₁
    y[i] = ys+α*u₂
    colors[i] = p
end


# points = [[node.x+α*node.d₁ for node in nodes]';[node.y+α*node.d₂ for node in nodes]';zeros(1,nᵤ)]
points = [x';y';zeros(1,nᵤ)]
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
vtk_grid("./vtk/cantilever_GLSMPP_"*poly1*"_"*string(ndiv)*"_"*string(ndiv),points,cells) do vtk
    vtk["𝑝"] = colors
end


show(to)
# fig


# dᵤ = zeros(2*nᵤ)
# for (i,node) in enumerate(nodes)
#     x = node.x
#     y = node.y
#     dᵤ[2*i-1] = u(x,y)
#     dᵤ[2*i]   = v(x,y)
# end
# dₚ = zeros(nₚ)


# dₛ = zeros(4*nₛ)
# for i in 1:3*nₑ
#     dₛ[4*i-3] = E/(1+ν)/(1-2*ν)*((1-ν)*2 + ν*6)
#     dₛ[4*i-2] = E/(1+ν)/(1-2*ν)*(ν*2 + (1-ν)*6)
#     dₛ[4*i-1] = E/(1+ν)/(1-2*ν)*(ν*2 + ν*6)
#     dₛ[4*i]   = E/(1+ν)*4
# end
# dₛₚ = zeros(4*nₛ)
# dₚ = zeros(nₚ)
# for i in 1:3*nₑ
#     dₛₚ[4*i-3] = (2*dₛ[4*i-3]-dₛ[4*i-2]-dₛ[4*i-1])/3
#     dₛₚ[4*i-2] = (-dₛ[4*i-3]+2*dₛ[4*i-2]-dₛ[4*i-1])/3
#     dₛₚ[4*i-1] = (-dₛ[4*i-3]-dₛ[4*i-2]+2*dₛ[4*i-1])/3
#     dₛₚ[4*i]   = dₛ[4*i]
# end

# for i in 1:3*nₑ
#     dₚ[i]   =(dₛ[4*i-3]+dₛ[4*i-2]+dₛ[4*i-1])/3
# end

# err1 = kᵖᵘ'*dₚ + kˢᵘ'*dₛₚ + fᵘ
# err2 = kᵖᵘ*dᵤ + kᵖᵖ*dₚ - fᵖ
# err3 = kˢᵘ*dᵤ + kˢˢ*dₛₚ - fˢ