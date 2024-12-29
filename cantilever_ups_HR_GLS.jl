
using TimerOutputs 
using SparseArrays, Pardiso
using CairoMakie
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫qpdxdy, ∫∫sᵢⱼsᵢⱼdxdy, ∫∫p∇udxdy, ∫∫sᵢⱼεᵢⱼdxdy, ∫pnᵢgᵢds, ∫sᵢⱼnⱼgᵢds, ∫∫τ∇q∇pdxdy, ∫∫τ∇sᵢⱼ∇sᵢₖdxdy, ∫∫τ∇sᵢⱼ∇pdxdy, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric, ∫∫∇sᵢⱼuᵢdxdy, ∫∫∇puᵢdxdy ,∫pnᵢuᵢds,∫sᵢⱼnⱼuᵢds 
using  Printf
include("import_cantilever.jl")
include("wirteVTK.jl")
const to = TimerOutput()
ps = MKLPardisoSolver()

ndiv = 2
ndiv2 = 2
nₚ =  40
poly = "tri3"
# poly = "quad"
@timeit to "import data" begin
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_c_"*string(nₚ)*".msh")
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly*"_"*string(ndiv2)*".msh")
elements, nodes,nodes_p,  sp, type, Ω = import_HR_GLS_MMP("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly*"_"*string(ndiv2)*".msh")
# elements, nodes,  sp, type = import_HR_GLS_reduced("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly*"_"*string(ndiv2)*".msh")
end

nₑ = length(elements["Ωᵘ"])
nₛ = 3*nₑ
nᵤ = length(nodes)
# nₚ = 3*nₑ
nₑₚ = length(Ω)
# nₚ = length(nodes_p)
nₚ = length(nodes)

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

β = 1*ℎ^2/2/𝐺
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

# kˢˢ = spzeros(4*nₛ*nₑ,4*nₛ*nₑ)
# kᵖᵖ = spzeros(nₚ,nₚ)
# kˢᵘ = spzeros(4*nₛ*nₑ,2*nᵤ)
# kᵖᵘ = spzeros(nₚ,2*nᵤ)
# fˢ = spzeros(4*nₛ*nₑ)
# fᵖ = spzeros(nₚ)
# fᵘ = spzeros(2*nᵤ)

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
# d = zeros(2*nᵤ+nₚ+4*nₛ)
d = k\f



# @timeit to "solve" d = k\f

# @timeit to "solve" d = solve(ps, k, f)

𝑢₁ = d[1:2:2*nᵤ]
𝑢₂ = d[2:2:2*nᵤ]
𝑝 = d[2*nᵤ+1:2*nᵤ+nₚ]
push!(nodes,:d₁=>𝑢₁)
push!(nodes,:d₂=>𝑢₂)
push!(nodes_p,:p=>𝑝)

# eval(VTK_mix_pressure_HR)
@timeit to "compute error" begin
Hₑ_𝒖, L₂_𝒖 = Hₑ_PlaneStress(elements["Ωᵍᵘ"])
Hₑ_dev = Hₑ_PlaneStrain_Deviatoric(elements["Ωᵍᵘ"])
L₂_𝑝 = L₂𝑝(elements["Ωᵍᵖ"])
end

println(log10(L₂_𝒖))
println(log10(Hₑ_𝒖))
println(log10(Hₑ_dev))
println(log10(L₂_𝑝))

# @timeit to "plot figure" begin
# fig = Figure()
# ind = 100
# ax = Axis(fig[1,1], 
#     aspect = DataAspect(), 
#     xticksvisible = false,
#     xticklabelsvisible=false, 
#     yticksvisible = false, 
#     yticklabelsvisible=false,
# )
# hidespines!(ax)
# hidedecorations!(ax)
# xs = LinRange(0, 48, 4*ind)
# ys = LinRange(-6, 6, ind)
# zs = zeros(4*ind,ind)
# 𝗠 = zeros(21)
# for (i,x) in enumerate(xs)
#     for (j,y) in enumerate(ys)
#         indices = sp(x,y,0.0)
#         ni = length(indices)
#         𝓒 = [nodes_p[i] for i in indices]
#         data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
#         ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
#         𝓖 = [ξ]
#         a = type(𝓒,𝓖)
#         set𝝭!(a)
#         p = 0.0
#         N = ξ[:𝝭]
#         for (k,xₖ) in enumerate(𝓒)
#             p += N[k]*xₖ.p
#         end
#         zs[i,j] = p
#     end
# end
# surface!(xs,ys,zeros(4*ind,ind),color=zs,shading=NoShading,colormap=:lightrainbow)
# contour!(xs,ys,zs,levels=-1e3:200:1e3,color=:azure)
# # Colorbar(fig[1,2], limits=(-900,900), colormap=:lightrainbow)
# save("./png/cantilever_mix_tri3_"*string(ndiv)*"_"*string(nₚ)*"_ls.png",fig, px_per_unit = 10.0)
# end




show(to)
# fig