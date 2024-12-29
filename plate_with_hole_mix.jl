
using ApproxOperator, XLSX, TimerOutputs
using SparseArrays, Pardiso
using ApproxOperator.Elasticity: ∫∫σᵢⱼσₖₗdxdy,∫∫σᵢⱼσₖₗdxdy_PlaneStrian, ∫∫∇σᵢⱼuᵢdxdy, ∫σᵢⱼnⱼuᵢds, ∫σᵢⱼnⱼgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, Hₑ_PlaneStress,∫∫τ∇σᵢⱼ∇σᵢₖdxdy,∫∫τ∇σᵢⱼ∇σᵢₖdxdy_new,∫∫τ∇σᵢⱼ∇σᵢₖdxdy_Taylor 

include("import_plate_with_hole.jl")


const to = TimerOutput()
ps = MKLPardisoSolver()
# n = [2 4 8 16]
# for i in 1:4
# ndivs = n[i]
# ndiv = n[i]
ndivs = 8
ndiv = 8
# elements, nodes = import_patchtest_mix("msh/patchtest_u_"*string(nₚ)*".msh","./msh/patchtest_"*string(ndiv)*".msh");
elements, nodes, ds₂, ds₁ = import_plate_with_hole_mix("msh/PlateWithHole_"*string(ndivs)*".msh","./msh/PlateWithHole_"*string(ndiv)*".msh",2*ndiv,0.91);

nₛ = 6
nₚ = length(nodes)
nₑ = length(elements["Ω"])
@timeit to "shape function" begin 
set𝝭!(elements["Ω"])
set𝝭!(elements["∂Ω"])
set∇𝝭!(elements["Ωᵍ"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])
set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["∂Ωˢ"])
end
T = 1000.0
E = 3e6
ν = 0.5-1e-5
# ν = 0.3 

Ē = E/(1.0-ν^2)
ν̄ = ν/(1.0-ν)
Cᵢᵢᵢᵢ = E/(1+ν)/(1-2*ν)*(1-ν)
Cᵢᵢⱼⱼ = E/(1+ν)/(1-2*ν)*ν
Cᵢⱼᵢⱼ = E/(1+ν)/2
𝐺 = E/(1+ν)/2
a = 1.0
# r(x,y) = (x^2+y^2)^0.5
# θ(x,y) = atan(y/x)
# u(x,y) = T*a*(1+ν)/2/E*( r(x,y)/a*2/(1+ν)*cos(θ(x,y)) + a/r(x,y)*(4/(1+ν)*cos(θ(x,y))+cos(3*θ(x,y))) - a^3/r(x,y)^3*cos(3*θ(x,y)) )
# v(x,y) = T*a*(1+ν)/2/E*( -r(x,y)/a*2*ν/(1+ν)*sin(θ(x,y)) - a/r(x,y)*(2*(1-ν)/(1+ν)*sin(θ(x,y))-sin(3*θ(x,y))) - a^3/r(x,y)^3*sin(3*θ(x,y)) )
# ∂u∂x(x,y) = T/E*(1 + a^2/2/r(x,y)^2*((ν-3)*cos(2*θ(x,y))-2*(1+ν)*cos(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*cos(4*θ(x,y)))
# ∂u∂y(x,y) = T/E*(-a^2/r(x,y)^2*((ν+5)/2*sin(2*θ(x,y))+(1+ν)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*sin(4*θ(x,y)))
# ∂v∂x(x,y) = T/E*(-a^2/r(x,y)^2*((ν-3)/2*sin(2*θ(x,y))+(1+ν)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*sin(4*θ(x,y)))
# ∂v∂y(x,y) = T/E*(-ν - a^2/2/r(x,y)^2*((1-3*ν)*cos(2*θ(x,y))-2*(1+ν)*cos(4*θ(x,y))) - 3*a^4/2/r(x,y)^4*(1+ν)*cos(4*θ(x,y)))
# σ₁₁(x,y) = T - T*a^2/r(x,y)^2*(3/2*cos(2*θ(x,y))+cos(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
# σ₂₂(x,y) = -T*a^2/r(x,y)^2*(1/2*cos(2*θ(x,y))-cos(4*θ(x,y))) - T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
# σ₁₂(x,y) = -T*a^2/r(x,y)^2*(1/2*sin(2*θ(x,y))+sin(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*sin(4*θ(x,y))




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

# β(x,y) = 0.1*((r(x,y)/2^0.5)*ℎ)^2/2/𝐺
ℎ(x,y) = (ds₁ + (r(x,y)-1)/4/2^0.5*(ds₂-ds₁))
β(x,y) = 0.1*(ds₁ + (r(x,y)-1)/4/2^0.5*(ds₂-ds₁))^2/2/𝐺
prescribe!(elements["Ωˢ"],:τ=>(x,y,z)->β(x,y), index=:𝑔)
prescribe!(elements["Ωˢ"],:ℎ=>(x,y,z)->ℎ(x,y), index=:𝑔) 
prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵍ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵍ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z,n₁,n₂)->(1-abs(n₂))*abs(n₁))
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z,n₁,n₂)->(1-abs(n₁))*abs(n₂))
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))

prescribe!(elements["Ωˢ"],:b₁=>(x,y,z)->0.0)
prescribe!(elements["Ωˢ"],:b₂=>(x,y,z)->0.0)

# 𝑎 = ∫∫σᵢⱼσₖₗdxdy=>elements["Ωˢ"]
𝑎 = ∫∫σᵢⱼσₖₗdxdy_PlaneStrian=>elements["Ωˢ"]
𝑏 = [
    ∫σᵢⱼnⱼuᵢds=>(elements["∂Ωˢ"],elements["∂Ω"]),
    ∫∫∇σᵢⱼuᵢdxdy=>(elements["Ωˢ"],elements["Ω"]),
]
𝑏ᵅ = ∫σᵢⱼnⱼgᵢds=>(elements["Γᵍˢ"],elements["Γᵍ"])

𝑏ᵝ = ∫∫τ∇σᵢⱼ∇σᵢₖdxdy=>elements["Ωˢ"]
# 𝑏ᵝ = ∫∫τ∇σᵢⱼ∇σᵢₖdxdy_Taylor=>elements["Ωˢ"]
𝑓 =  ∫vᵢtᵢds=>elements["Γᵗ"]

@timeit to "assembly matrix" begin

kᵖᵖ = spzeros(3*nₛ*nₑ,3*nₛ*nₑ)
fᵖ = zeros(3*nₛ*nₑ)
kᵖᵘ = spzeros(3*nₛ*nₑ,2*nₚ)
fᵘ = zeros(2*nₚ)
# d = zeros(3*nₛ*nₑ+2*nₚ)

𝑎(kᵖᵖ)
𝑏(kᵖᵘ)
𝑏ᵅ(kᵖᵘ,fᵖ)
𝑏ᵝ(kᵖᵖ,fᵖ)
𝑓(fᵘ)
end

# k = sparse([kᵖᵖ kᵖᵘ;kᵖᵘ' zeros(2*nₚ,2*nₚ)])
# set_matrixtype!(ps,-2)
# k = get_matrix(ps,k,:N)
# f = [fᵖ;-fᵘ]
# @timeit to "solve" pardiso(ps,d,k,f)
d = [kᵖᵖ kᵖᵘ;kᵖᵘ' zeros(2*nₚ,2*nₚ)]\[fᵖ;-fᵘ]
d₁ = d[3*nₛ*nₑ+1:2:end]
d₂ = d[3*nₛ*nₑ+2:2:end]
push!(nodes,:d₁=>d₁,:d₂=>d₂)

# 𝐿₂ = L₂(elements["Ωᵍ"])
𝐻ₑ, 𝐿₂ = Hₑ_PlaneStress(elements["Ωᵍ"])
println(log10(𝐿₂))
println(log10(𝐻ₑ))

# XLSX.openxlsx("./xlsx/platewithhole.xlsx", mode="rw") do xf
# index = 3,4,5,6,7,8,10,289,81,1089,32,16,4225,19,628,2272,78,12,1081,4255,25,23
#     Sheet = xf[3]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = nₑ
#     Sheet["B"*string(ind)] = log10(𝐿₂)
#     Sheet["C"*string(ind)] = log10(𝐻ₑ)
# end
show(to)

# end