# using Revise
using TimerOutputs 
using SparseArrays, Pardiso
using CairoMakie
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫qpdxdy, ∫∫sᵢⱼsᵢⱼdxdy, ∫∫p∇udxdy, ∫∫sᵢⱼεᵢⱼdxdy, ∫pnᵢgᵢds, ∫sᵢⱼnⱼgᵢds, ∫∫τ∇q∇pdxdy, ∫∫τ∇sᵢⱼ∇sᵢₖdxdy, ∫∫τ∇sᵢⱼ∇pdxdy, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric, ∫∫∇sᵢⱼuᵢdxdy, ∫∫∇puᵢdxdy ,∫pnᵢuᵢds,∫sᵢⱼnⱼuᵢds 
using  Printf
include("import_plate_with_hole.jl")

const to = TimerOutput()
ps = MKLPardisoSolver()

ndiv = 8
# nₚ = 243
# poly = "tri6"
poly = "tri3"
# poly = "quad"
@timeit to "import data" begin
n = 8
# poly = "tri6"
elements, nodes, nodes_p = import_elasticity_linear_mix("./msh/plate_with_hole_"*poly*"_"*string(ndiv)*".msh","./msh/plate_with_hole_"*poly*"_"*string(ndiv)*".msh",n)
# elements, nodes, nodes_p = import_elasticity_quadratic_mix("./msh/plate_with_hole_tri6_"*string(ndiv)*".msh","./msh/plate_with_hole_tri3_"*string(n)*".msh",n)
# nx = 131;ny = 32
# elements, nodes, nodes_p = import_linear_mix("./msh/plate_with_hole_"*poly*"_"*string(ndiv)*".msh","./msh/plate_with_hole_"*poly*"_"*string(ndiv)*".msh",n)
nₚ = length(nodes_p)
end

nₑ = length(elements["Ωᵘ"])
nₛ = 3*nₑ
nᵤ = length(nodes)

# T = 1.0e3
# E = 3.0e6
T = 1.0
E = 1.0e2
ν = 0.4995
# ν = 0.5-1e-8
Ē = E/(1.0-ν^2)
ν̄ = ν/(1.0-ν)
Cᵢᵢᵢᵢ = E/(1+ν)/(1-2*ν)*(1-ν)
Cᵢᵢⱼⱼ = E/(1+ν)/(1-2*ν)*ν
Cᵢⱼᵢⱼ = E/(1+ν)/2
𝐺 = E/(1+ν)/2

a = 1
b = 5
D = 5.0
ℎ = D/ndiv
# n = 1
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

r(x,y) = (x^2+y^2)^0.5
θ(x,y) = atan(y/x)
u(x,y) = T*a*(1+ν̄)/2/Ē*(r(x,y)/a*2/(1+ν̄)*cos(θ(x,y)) + a/r(x,y)*(4/(1+ν̄)*cos(θ(x,y))+cos(3*θ(x,y))) - a^3/r(x,y)^3*cos(3*θ(x,y)))
v(x,y) = T*a*(1+ν̄)/2/Ē*( -r(x,y)/a*2*ν̄/(1+ν̄)*sin(θ(x,y)) - a/r(x,y)*(2*(1-ν̄)/(1+ν̄)*sin(θ(x,y))-sin(3*θ(x,y))) - a^3/r(x,y)^3*sin(3*θ(x,y)) )
∂u∂x(x,y) = T/Ē*(1 + a^2/2/r(x,y)^2*((ν̄-3)*cos(2*θ(x,y))-2*(1+ν̄)*cos(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν̄)*cos(4*θ(x,y)))
∂u∂y(x,y) = T/Ē*(-a^2/r(x,y)^2*((ν̄+5)/2*sin(2*θ(x,y))+(1+ν̄)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν̄)*sin(4*θ(x,y)))
∂v∂x(x,y) = T/Ē*(-a^2/r(x,y)^2*((ν̄-3)/2*sin(2*θ(x,y))+(1+ν̄)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν̄)*sin(4*θ(x,y)))
∂v∂y(x,y) = T/Ē*(-ν̄ - a^2/2/r(x,y)^2*((1-3*ν̄)*cos(2*θ(x,y))-2*(1+ν̄)*cos(4*θ(x,y))) - 3*a^4/2/r(x,y)^4*(1+ν̄)*cos(4*θ(x,y)))

ε₁₁(x,y) = ∂u∂x(x,y)
ε₂₂(x,y) = ∂v∂y(x,y)
ε₁₂(x,y) = 0.5*(∂u∂y(x,y) + ∂v∂x(x,y))
σ₁₁(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*ε₁₁(x,y) + ν*ε₂₂(x,y))
σ₂₂(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + (1-ν)*ε₂₂(x,y))
σ₃₃(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + ν*ε₂₂(x,y))
σ₁₂(x,y) = E/(1+ν)*ε₁₂(x,y)

# ∂ε₁₁∂x(x,y) = ∂²u∂x²(x,y)
# ∂ε₁₁∂y(x,y) = ∂²u∂x∂y(x,y)
# ∂ε₂₂∂x(x,y) = ∂²v∂x∂y(x,y)
# ∂ε₂₂∂y(x,y) = ∂²v∂y²(x,y)
# ∂ε₁₂∂x(x,y) = 0.5*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
# ∂ε₁₂∂y(x,y) = 0.5*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))

# ∂σ₁₁∂x(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂x(x,y) + ν*∂ε₂₂∂x(x,y))
# ∂σ₁₁∂y(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂y(x,y) + ν*∂ε₂₂∂y(x,y))
# ∂σ₂₂∂x(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y) + (1-ν)*∂ε₂₂∂x(x,y))
# ∂σ₂₂∂y(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y) + (1-ν)*∂ε₂₂∂y(x,y))
# ∂σ₁₂∂x(x,y) = E/(1+ν)*∂ε₁₂∂x(x,y)
# ∂σ₁₂∂y(x,y) = E/(1+ν)*∂ε₁₂∂y(x,y)
# b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
# b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)
p(x,y) = (σ₁₁(x,y)+σ₂₂(x,y)+σ₃₃(x,y))/3


β = 10*ℎ^2/2/𝐺
prescribe!(elements["Ωˢ"],:τ=>(x,y,z)->β, index=:𝑔)
prescribe!(elements["Ωᵖ"],:τ=>(x,y,z)->β, index=:𝑔)

prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν)

# prescribe!(elements["Ωᵘ"],:b₁=>(x,y,z)->b₁(x,y))
# prescribe!(elements["Ωᵘ"],:b₂=>(x,y,z)->b₂(x,y))

prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂) 
prescribe!(elements["Γᵍᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍᵘ"],:g₂=>(x,y,z)->v(x,y))
# prescribe!(elements["Γᵍᵘ"],:n₁₁=>(x,y,z,n₁,n₂)->1.0)
# prescribe!(elements["Γᵍᵘ"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₁₁=>(x,y,z,n₁,n₂)->abs(n₁))
prescribe!(elements["Γᵍᵘ"],:n₂₂=>(x,y,z,n₁,n₂)->abs(n₂))
prescribe!(elements["Γᵍᵘ"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Ωᵍᵘ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍᵘ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
prescribe!(elements["Ωᵍᵖ"],:p=>(x,y,z)->p(x,y))

prescribe!(elements["Ωᵖ"],:b₁=>(x,y,z)->0.0)
prescribe!(elements["Ωᵖ"],:b₂=>(x,y,z)->0.0)
prescribe!(elements["Ωˢ"],:b₁=>(x,y,z)->0.0)
prescribe!(elements["Ωˢ"],:b₂=>(x,y,z)->0.0)

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
    
# 𝑏ˢᵝ(kˢˢ,fˢ)
# 𝑏ᵖᵝ(kᵖᵖ,fᵖ)
# 𝑏ˢᵖᵝ(kˢᵖ)
𝑓(fᵘ)
end
# k = [zeros(2*nᵤ,2*nᵤ) kᵖᵘ' kˢᵘ';kᵖᵘ kᵖᵖ zeros(nₚ,4*nₛ*nₑ);kˢᵘ zeros(4*nₛ*nₑ,nₚ) kˢˢ]
k = sparse([zeros(2*nᵤ,2*nᵤ) kᵖᵘ' kˢᵘ';kᵖᵘ kᵖᵖ zeros(nₚ,4*nₛ);kˢᵘ zeros(4*nₛ,nₚ) kˢˢ])
f = [-fᵘ;fᵖ;fˢ]
d = zeros(2*nᵤ+nₚ+4*nₛ)
# d = k\f

set_matrixtype!(ps, -2)
k = get_matrix(ps,k,:N)
@timeit to "solve" pardiso(ps,d,k,f)

𝑢₁ = d[1:2:2*nᵤ]
𝑢₂ = d[2:2:2*nᵤ]
𝑝 = d[2*nᵤ+1:2*nᵤ+nₚ]
push!(nodes,:d₁=>𝑢₁)
push!(nodes,:d₂=>𝑢₂)
push!(nodes_p,:p=>𝑝)

@timeit to "compute error" begin
Hₑ_𝒖, L₂_𝒖 = Hₑ_PlaneStress(elements["Ωᵍᵘ"])
Hₑ_dev = Hₑ_PlaneStrain_Deviatoric(elements["Ωᵍᵘ"])
L₂_𝑝 = L₂𝑝(elements["Ωᵍᵖ"])
end

println(log10(L₂_𝒖))
println(log10(Hₑ_𝒖))
println(log10(Hₑ_dev))
println(log10(L₂_𝑝))


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
# Colorbar(fig[1,2], limits=(-900,900), colormap=:lightrainbow)
# save("./png/cantilever_mix_"*poly*"_"*string(ndiv)*"_"*string(nₚ)*".png",fig, px_per_unit = 10.0)
# end

show(to)
# fig