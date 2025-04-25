
using TimerOutputs 
using SparseArrays, Pardiso
using CairoMakie, WriteVTK
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫εᵢⱼσᵢⱼdxdy, ∫∫εᵢⱼσᵢⱼdxdy_PlaneStrian,∫σᵢⱼnⱼgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, Hₑ_PlaneStress,  ∫vᵢgᵢds

include("import_cook.jl")


const to = TimerOutput()
# ps = MKLPardisoSolver()
nn = [ 4 8 12 16 20 24]
for i in 5:6
ndiv = nn[i]
ndiv2 = nn[i]
n = nn[i]

# ndiv = 4
# ndiv2 = 4
# n = 4
poly = "tri3"
# poly = "nonuniform"
# poly = "tri6"
# poly = "quad"
@timeit to "import data" begin

# elements, nodes = import_MF_Gauss("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh")
elements, nodes, sp, type = import_MF_Gauss("./msh/cook_membrane_"*poly*"_"*string(ndiv)*".msh",n)

end

nₑ = length(elements["Ω"])

nₚ = length(nodes)
# nₚ = length(nodes_p)
# nₚ = length(nodes)

L = 48.0
D = 44.0
P = 6.25
ℎ = D/ndiv

# Ē = 3e6
# ν̄  = 0.3
E = 70
# ν = 0.3
ν = 0.5-1e-5
Ē = E/(1.0-ν^2)
ν̄ = ν/(1.0-ν)
I = D^3/12
EI = Ē*I
Cᵢᵢᵢᵢ = E/(1+ν)/(1-2*ν)*(1-ν)
Cᵢᵢⱼⱼ = E/(1+ν)/(1-2*ν)*ν
Cᵢⱼᵢⱼ = E/(1+ν)/2
𝐺 = E/(1+ν)/2
K=E/3/(1-2ν )


prescribe!(elements["Ω"],:E=>(x,y,z)->Ē,index=:𝑔)
prescribe!(elements["Ω"],:ν=>(x,y,z)->ν̄,index=:𝑔)
# prescribe!(elements["Ωᵍ"],:E=>(x,y,z)->E,index=:𝑔)
# prescribe!(elements["Ωᵍ"],:ν=>(x,y,z)->ν,index=:𝑔)
# prescribe!(elements["Γᵍ"],:E=>(x,y,z)->E)
# prescribe!(elements["Γᵍ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->P)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P)
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)

prescribe!(elements["Γᵍ"],:α=>(x,y,z)->1e3*E,index=:𝑔)
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)


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





α = 1.0

# vertices = [[node.x+α*node.d₁ for node in nodes] [node.y+α*node.d₂ for node in nodes]]
colors = zeros(nₚ)
x = zeros(nₚ)
y = zeros(nₚ)

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
points = [x';y';zeros(1,nₚ)]
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["ΩC"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
vtk_grid("./vtk/cook_GI13_"*poly*"_"*string(ndiv)*"_"*string(ndiv),points,cells) do vtk
    vtk["𝑝"] = colors
end

println(y[3] - nodes[3].y)




show(to)




end
