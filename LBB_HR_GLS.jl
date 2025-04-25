
using TimerOutputs 
using SparseArrays, Pardiso
using CairoMakie, LinearAlgebra
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric,∫∫σᵢⱼσₖₗdxdy_PlaneStrian,∫∫σᵢⱼσₖₗdxdy,∫σᵢⱼnⱼuᵢds,∫∫∇σᵢⱼuᵢdxdy,∫σᵢⱼnⱼgᵢds,∫∫τ∇σᵢⱼ∇σᵢₖdxdy,∫∫uᵢⱼuₖₗdxdy

include("import_cantilever.jl")


const to = TimerOutput()
ps = MKLPardisoSolver()
n = [ 2 4 8 16 32 ]
for i in 5:5
ndiv = n[i]
ndiv2 = n[i]
# ndiv = 4
# ndiv2 = 4
poly = "tri3"
test = "square" 
# test = "cantilever" 
# poly = "tri6"
# poly = "quad"
@timeit to "import data" begin
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_c_"*string(nₚ)*".msh")
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly*"_"*string(ndiv2)*".msh")
elements, nodes, sp, type = import_HR_GLS("./msh/"*test*"_"*poly*"_"*string(ndiv)*".msh","./msh/"*test*"_"*poly*"_"*string(ndiv2)*".msh")
# elements, nodes,  sp, type = import_HR_GLS_reduced("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*poly*"_"*string(ndiv2)*".msh")
# elements, nodes, sp, type = import_HR_GLS("./msh/cantilever_nonuniform_"*string(ndiv)*".msh","./msh/cantilever_nonuniform_"*string(ndiv2)*".msh")
end

nₑ = length(elements["Ωᵘ"])
nₛ = 3
nᵤ = length(nodes)
# nₚ = length(nodes_p)
# nₚ = length(nodes)

L = 1.0
D = 1.0
P = 0
ℎ = D/ndiv

# Ē = 3e6
Ē = 1.0
# ν̄  = 0.3
ν̄  = 0.5-1e-7
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


β =1*ℎ^2/2/𝐺
prescribe!(elements["Ωˢ"],:τ=>(x,y,z)->β)
prescribe!(elements["Ωˢ"],:ℎ=>(x,y,z)->ℎ) 
prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν)

prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν)

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


𝑎 =∫∫σᵢⱼσₖₗdxdy=>elements["Ωˢ"]

𝑏 = [
    ∫σᵢⱼnⱼuᵢds=>(elements["∂Ωˢ"],elements["∂Ωᵘ"]),
    ∫∫∇σᵢⱼuᵢdxdy=>(elements["Ωˢ"],elements["Ωᵘ"]),
   
    ]

𝑏ᵅ = ∫σᵢⱼnⱼgᵢds=>(elements["Γᵍˢ"],elements["Γᵍᵘ"])


𝑏ᵝ = ∫∫τ∇σᵢⱼ∇σᵢₖdxdy=>elements["Ωˢ"]
c= ∫∫uᵢⱼuₖₗdxdy=>elements["Ωᵘ"]

𝑓 = ∫vᵢtᵢds=>elements["Γᵗ"]
# 𝑓 = [
#     ∫vᵢtᵢds=>elements["Γᵗ"]∪elements["Γʳ"],
#     ∫∫vᵢbᵢdxdy=>elements["Ωᵘ"]
# ]

kᵘᵘ = zeros(2*nᵤ,2*nᵤ)
kˢˢ = zeros(3*nₛ*nₑ,3*nₛ*nₑ)

kˢᵘ = zeros(3*nₛ*nₑ,2*nᵤ)
kˢᵘⁿ  = zeros(3*nₛ*nₑ,2*nᵤ)
fˢ = zeros(3*nₛ*nₑ)
fᵘ = zeros(2*nᵤ)



@timeit to "assembly" begin

    𝑎(kˢˢ)
    𝑏(kˢᵘ)
    𝑏ᵅ(kˢᵘ,fˢ)
    # 𝑏ᵝ(kˢˢ,fˢ)
    c(kᵘᵘ)
    𝑓(fᵘ)
end

# svdB = svd(kˢᵘ')
# S = svdB.S
# V = svdB.V

# # 设置一个容差
# tol = 1e-12
# # 找出奇异值小于 tol 的列，这些列对应零空间
# # null_inds = findall(x -> x < tol, S)
# null_inds = count(S .> tol)
# if isempty(null_inds)
#     # 如果没有近似为零的奇异值，则认为 B 的秩为 m，零空间基为 V 的后 n-m 列
#     Z = V[:, 2*nᵤ+1:end]
# else
#     Z = V[:, null_inds]
# end

qr_fact = qr(kˢᵘ)
# qr_fact = qr(kˢᵘ, full=true)
Q = qr_fact.Q

kˢˢ_p = Q' * kˢˢ *Q
kˢˢ_new=zeros(3*nₛ*nₑ-2*nᵤ,3*nₛ*nₑ-2*nᵤ) 
kˢˢ_new= kˢˢ_p[2*nᵤ+1:end,2*nᵤ+1:end]
println("size(kˢˢ_p) = ", size(kˢˢ_p))
println("size(kˢˢ_new) = ", size(kˢˢ_new))
println("size(Q) = ", size(Q))

# val = eigvals(kˢᵘ'*(kˢˢ\kˢᵘ),kᵘᵘ)
# val = eigvals(kˢᵘ'*((inv(kˢˢ)\kˢᵘ),kᵘᵘ)

val = eigvals(kˢᵘ*(kᵘᵘ\kˢᵘ'),kˢˢ)
# val = eigvals(kˢᵘ*inv(kᵘᵘ)*kˢᵘ',kˢˢ)

λ_min = minimum(eigvals(Symmetric(kˢˢ_new)))

# println(minimum(infsup_val))
# println(infsup_val[[3*nₛ*nₑ-2*nᵤ+1]])



val_sign = zeros(3*nₛ*nₑ)
for (ii,v) in enumerate(val)
    if v isa Real
        val_sign[ii] = sign(v)
    else
        val_sign[ii] = sign(v.re) < -1e-5 ? -1.0 : 1.0
    end
end


val_real = val_sign .* abs.(val)
val_abs = abs.(val)
val_sort = sort(val_abs)
n_eig_positive = count(x-> isa(x,Real) ? x > 1e-5 : x.re > 1e-6,val)
n_eig_nonzeros = count(x-> x > 1e-5,val_sort)
min_eig_positive = isa(val[3*nₛ*nₑ - n_eig_positive + 1],Real) ? val[3*nₛ*nₑ - n_eig_positive + 1] : val[3*nₛ*nₑ - n_eig_positive + 1].re
min_eig_nonzeros = val_sort[3*nₛ*nₑ - n_eig_nonzeros + 1]
min_eig_real = min(val_real[abs.(val_real).>1e-5]...)




# val_sign = zeros(2*nᵤ)
# for (ii,v) in enumerate(val)
#     if v isa Real
#         val_sign[ii] = sign(v)
#     else
#         val_sign[ii] = sign(v.re) < -1e-10 ? -1.0 : 1.0
#     end
# end

# val_real = val_sign .* abs.(val)
# val_abs = abs.(val)
# val_sort = sort(val_abs)
# n_eig_positive = count(x-> isa(x,Real) ? x > 1e-10 : x.re > 1e-10,val)
# n_eig_nonzeros = count(x-> x > 1e-10,val_sort)
# min_eig_positive = isa(val[2*nᵤ - n_eig_positive + 1],Real) ? val[2*nᵤ - n_eig_positive + 1] : val[2*nᵤ - n_eig_positive + 1].re
# min_eig_nonzeros = val_sort[2*nᵤ - n_eig_nonzeros + 1]
# min_eig_real = min(val_real[abs.(val_real).>1e-10]...)



println(λ_min)
println(min_eig_nonzeros)
println(min_eig_positive)
show(to)
# fig
end
