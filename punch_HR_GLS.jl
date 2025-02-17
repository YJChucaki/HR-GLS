
using TimerOutputs 
using SparseArrays, Pardiso, Printf
using CairoMakie, WriteVTK
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric,∫∫σᵢⱼσₖₗdxdy_PlaneStrian,∫∫σᵢⱼσₖₗdxdy,∫σᵢⱼnⱼuᵢds,∫∫∇σᵢⱼuᵢdxdy,∫σᵢⱼnⱼgᵢds,∫∫τ∇σᵢⱼ∇σᵢₖdxdy,∫∫τ∇σᵢⱼ∇σᵢₖdxdy_new,∫∫τ∇σᵢⱼ∇σᵢₖdxdy_Taylor,∫∫σᵢⱼσₖₗdxdy_Taylor, Hₑ_PlaneStrain_Dil, 𝐿₂_PlaneStrain_Pressure, Hₑ_PlaneStrain_Deviatoric

include("import_punch.jl")
include("wirteVTK.jl")

const to = TimerOutput()
ps = MKLPardisoSolver()
# n = [ 2 4 8 16 ]
# for i in 1:4
# ndiv = n[i]
# ndiv2 = n[i]
ndiv = 2
ndiv2 = 2
poly = "tri3"
test = "punch"
# poly = "tri6"
# poly = "quad"
@timeit to "import data" begin
elements, nodes, sp, type, Ω, nodes_c = import_HR_GLS("./msh/punch_"*poly*"_"*string(ndiv)*".msh","./msh/punch_"*poly*"_"*string(ndiv2)*".msh")
 end
nₑ = length(elements["Ωᵘ"])
nₑₛ = length(elements["Ωˢ"])
nₛ = 3
ni = 6
nᵤ = length(nodes)


L = 48.0
D = 12.0
P = 20
ℎ = D/ndiv

# Ē = 3e6
# ν̄  = 0.3
E = 70.0
ν = 0.3
# ν = 0.5-1e-5
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
prescribe!(elements["Ωˢ"],:τ=>(x,y,z)->β, index=:𝑔)
prescribe!(elements["Ωˢ"],:ℎ=>(x,y,z)->ℎ, index=:𝑔) 
prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν, index=:𝑔)

prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν, index=:𝑔)

prescribe!(elements["Ωˢ"],:b₁=>(x,y,z)->0.0)
prescribe!(elements["Ωˢ"],:b₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->-20.0) 
# prescribe!(elements["Γʳ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
# prescribe!(elements["Γʳ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂) 
prescribe!(elements["Γᵍᵘ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍᵘ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₁₂=>(x,y,z)->0.0)


𝑎 =∫∫σᵢⱼσₖₗdxdy_PlaneStrian=>elements["Ωˢ"]
# 𝑎 =∫∫σᵢⱼσₖₗdxdy_Taylor=>elements["Ωˢ"]
𝑏 = [
    ∫σᵢⱼnⱼuᵢds=>(elements["∂Ωˢ"],elements["∂Ωᵘ"]),
    ∫∫∇σᵢⱼuᵢdxdy=>(elements["Ωˢ"],elements["Ωᵘ"]),
    # ∫σᵢⱼnⱼuᵢds_Taylor=>(elements["∂Ωˢ"],elements["∂Ωᵘ"]),
    # ∫∫∇σᵢⱼuᵢdxdy_Taylor=>(elements["Ωˢ"],elements["Ωᵘ"]),
    ]

𝑏ᵅ = ∫σᵢⱼnⱼgᵢds=>(elements["Γᵍˢ"]∪elements["Γᵗˢ"],elements["Γᵍᵘ"]∪elements["Γᵗ"])
# 𝑏ᵅ = ∫σᵢⱼnⱼgᵢds_Taylor=>(elements["Γᵍˢ"],elements["Γᵍᵘ"])

𝑏ᵝ = ∫∫τ∇σᵢⱼ∇σᵢₖdxdy=>elements["Ωˢ"]
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
    # 𝑏ᵝ(kˢˢ,fˢ)
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



pₑ = zeros(nₑ)
for (i,elm) in enumerate(elements["Ωˢ"])
    𝓒ₚ = elm.𝓒
    𝓖 = elm.𝓖
    𝓒 = elements["Ω"][i].𝓒
    x = (𝓒[1].x+𝓒[2].x+𝓒[3].x)/3
    x = (𝓒[1].x+𝓒[2].x+𝓒[3].x)/3
    σ₁₁ = 𝓒ₚ[1].dₛ₁₁+𝓒ₚ[2].dₛ₁₁*𝓒[1].x+𝓒ₚ[3].dₛ₁₁*𝓒[1].y
    σ₂₂ = 𝓒ₚ[1].dₛ₂₂+𝓒ₚ[2].dₛ₂₂*𝓒[1].x+𝓒ₚ[3].dₛ₂₂*𝓒[1].y
    σ₃₃ = ν*(σ₁₁ + σ₂₂)
    pₑ[i]= (σ₁₁ + σ₂₂ + σ₃₃)/3 
end

p_node = zeros(nᵤ)
w = zeros(nᵤ)
for (i,elm) in enumerate(elements["Ω"])
    𝓒 = elm.𝓒
     for (j,xⱼ) in enumerate(𝓒)
        J = xⱼ.𝐼
        p_node[J] +=pₑ[i]
        w[J] +=1 
     end
end

eval(VTK_HR_displacement_pressure)
eval(VTK_HR_displacement_pressure_smoothing)



# α = 1.0
# nc = length(nodes_c)
# # vertices = [[node.x+α*node.d₁ for node in nodes] [node.y+α*node.d₂ for node in nodes]]
# colors = zeros(nc)
# x = zeros(nc)
# y = zeros(nc)
# for (i,node) in enumerate(nodes)
#     xs = node.x
#     ys = node.y
#     indices = sp(xs,ys,0.0)
#     ni = length(indices)
#     𝓒 = [nodes[i] for i in indices]
#     # data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
#     data = Dict([:x=>(2,[xs]),:y=>(2,[ys]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:∂𝝭∂x=>(4,zeros(ni)),:∂𝝭∂y=>(4,zeros(ni)),:𝗠=>(0,zeros(21)),:∂𝗠∂x=>(0,zeros(21)),:∂𝗠∂y=>(0,zeros(21))])
#     ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
#     𝓖 = [ξ]
#     a = type(𝓒,𝓖)
#     set∇𝝭!(a) 
#     p = 0.0
#     d₁ = 0.0
#     d₂ = 0.0
#     u₁ = 0.0
#     u₂ = 0.0
#     𝝭 = ξ[:𝝭]
#     B₁ = ξ[:∂𝝭∂x]
#     B₂ = ξ[:∂𝝭∂y]
#     ε₁₁ = 0.0
#     ε₂₂ = 0.0
#     ε₁₂ = 0.0
#     N = ξ[:𝝭]
#     for (k,xₖ) in enumerate(𝓒)
#         ε₁₁ += B₁[k]*xₖ.d₁
#         ε₂₂ += B₂[k]*xₖ.d₂
#         ε₁₂ += B₁[k]*xₖ.d₂ + B₂[k]*xₖ.d₁
#         u₁ += 𝝭[k]*xₖ.d₁
#         u₂ += 𝝭[k]*xₖ.d₂
#     end
#     p=K*(ε₁₁+ε₂₂)
#     x[i] = xs+α*u₁
#     y[i] = ys+α*u₂
#     colors[i] = p
# end


# # points = [[node.x+α*node.d₁ for node in nodes]';[node.y+α*node.d₂ for node in nodes]';zeros(1,nᵤ)]
# points = [x';y';zeros(1,nc)]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# # # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # # cells = [MeshCell(VTKCellTypes.VTK_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# vtk_grid("./vtk/cantilever_GLS_"*poly*"_"*string(ndiv)*"_"*string(ndiv),points,cells) do vtk
#     vtk["𝑝"] = colors
# end
show(to)
# fig
# end
