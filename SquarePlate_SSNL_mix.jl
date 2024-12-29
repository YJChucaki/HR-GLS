using ApproxOperator, JLD, XLSX, Printf
using CairoMakie
using SparseArrays, Pardiso
import BenchmarkExample: BenchmarkExample
include("import_SquarePlate.jl")
include("wirteVTK.jl")
ndiv  = 16
ndivs = 16
ndivs2 = 16
# elements, nodes, nodes_s, Ω, sp, type = import_SquarePlate_mix("msh/SquarePlate/SquarePlate_"*string(ndiv)*".msh","msh/SquarePlate/SquarePlate_"*string(ndivs)*".msh");
elements, nodes, nodes_s, Ω, sp, type = import_SquarePlate_mix("msh/SquarePlate/SquarePlate_quad8_"*string(ndiv)*".msh","msh/SquarePlate/SquarePlate_quad8_"*string(ndivs)*".msh");
# elements, nodes, nodes_s, Ω, sp, type = import_SquarePlate_mix("msh/SquarePlate/SquarePlate_tri6_"*string(ndiv)*".msh","msh/SquarePlate/SquarePlate_tri6_"*string(ndivs)*".msh");
# elements, nodes, nodes_s, Ω, sp, type = import_SquarePlate_mix("msh/SquarePlate/SquarePlate_"*string(ndiv)*".msh","msh/SquarePlate/SquarePlate_"*string(ndivs)*"_"*string(ndivs2)*".msh");
# elements, nodes, nodes_s, Ω , sp, type= import_SquarePlate_mix("msh/SquarePlate/SquarePlate_quad_"*string(ndiv)*".msh","msh/SquarePlate/SquarePlate_quad_"*string(ndivs)*".msh");
# elements, nodes, nodes_s, Ω = import_SquarePlate_mix("msh/SquarePlate/SquarePlate_quad_"*string(ndiv)*".msh","msh/SquarePlate/SquarePlate_quad_"*string(ndivs)*"_"*string(ndivs2)*".msh");

nᵇ = length(nodes)
nˢ = length(nodes_s)
nₑ = length(elements["Ω"])
nₑₛ = length(Ω)
E = 10.92e6
ν = 0.3
h = 0.001
L = 1.0
# ps = MKLPardisoSolver()

Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))
w(x,y) = 1/3*x^3*(x-1)^3*y^3*(y-1)^3-2*h^2/(5*(1-ν))*(y^3*(y-1)^3*x*(x-1)*(5*x^2-5*x+1)+x^3*(x-1)^3*y*(y-1)*(5*y^2-5*y+1))
θ₁(x,y) = y^3*(y-1)^3*x^2*(x-1)^2*(2*x-1)
θ₂(x,y) = x^3*(x-1)^3*y^2*(y-1)^2*(2*y-1)
F(x,y) = E*h^3/(12*(1-ν^2))*(12*y*(y-1)*(5*x^2-5*x+1)*(2*y^2*(y-1)^2+x*(x-1)*(5*y^2-5*y+1))+12*x*(x-1)*(5*y^2-5*y+1)*(2*x^2*(x-1)^2+y*(y-1)*(5*x^2-5*x+1)))

w₁(x,y) = (x-1)^2*x^2*(2*x-1)*(y-1)^3*y^3-2*h^2/(5*(1-ν))*((20*x^3-30*x^2+12*x-1)*(y-1)^3*y^3+3*(x-1)^2*x^2*(2*x-1)*(y-1)*y*(5*y^2-5*y+1))
w₂(x,y) = (x-1)^3*x^3*(y-1)^2*y^2*(2*y-1)-2*h^2/(5*(1-ν))*(3*(x-1)*x*(5*x^2-5*x+1)*(y-1)^2*y^2*(2*y-1)+x^3*(x-1)^3*(20*y^3-30*y^2+12*y-1))
# θ₁₁(x,y) = 2*(x-1)*x*(5*x^2-5*x+1)*(y-1)^3*y^3
# θ₁₂(x,y) = 3*(x-1)^2*x^2*(2*x-1)*(y-1)^2*y^2*(2*y-1)
# θ₂₂(x,y) = 2*(x-1)^3*x^3*(y-1)*y*(5*y^2-5*y+1)
# M₁₁(x,y)= -Dᵇ*(θ₁₁(x,y)+ν*θ₂₂(x,y))
# M₁₂(x,y)= -Dᵇ*(1-ν)*θ₁₂(x,y)
# M₂₂(x,y)= -Dᵇ*(ν*θ₁₁(x,y)+θ₂₂(x,y))
Q₁(x,y) = Dˢ*(w₁(x,y)-θ₁(x,y))
Q₂(x,y) = Dˢ*(w₂(x,y)-θ₂(x,y))
eval(prescribeForSSNonUniformLoading)
# eval(prescribeForSimpleSupported)
# eval(prescribeForCantilever)
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωˢ"])
set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["Γᵇ"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γˡ"])
set𝝭!(elements["Γʳ"])

ops = [
    Operator{:∫κMdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫wQdΩ}(),
    Operator{:∫QQdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫wqdΩ}(),
    Operator{:∫vwdΓ}(:α=>1e13*E),
    Operator{:∫vθ₁dΓ}(:α=>1e13*E),
    Operator{:∫vθ₂dΓ}(:α=>1e13*E),
    Operator{:L₂_ThickPlate}(:E=>E,:ν=>ν),
    Operator{:L₂_ThickPlate_Q}(:E=>E,:ν=>ν),
    Operator{:∫θM₁dΓ}(),
    Operator{:∫θM₂dΓ}(),
    Operator{:∫wVdΓ}(),
]
kᵇ = spzeros(3*nᵇ,3*nᵇ)
kʷˢ = spzeros(3*nᵇ,2*nˢ)
kˢˢ = spzeros(2*nˢ,2*nˢ)
f = zeros(3*nᵇ)
# d = zeros(3*nᵇ+2*nˢ)

ops[1](elements["Ω"],kᵇ)
ops[2](elements["Ω"],elements["Ωˢ"],kʷˢ)
ops[3](elements["Ωˢ"],kˢˢ)
ops[4](elements["Ω"],f)
ops[5](elements["Γᵇ"],kᵇ,f)
ops[5](elements["Γᵗ"],kᵇ,f)
ops[5](elements["Γˡ"],kᵇ,f)
ops[5](elements["Γʳ"],kᵇ,f)
ops[6](elements["Γᵇ"],kᵇ,f)
ops[6](elements["Γᵗ"],kᵇ,f)
ops[6](elements["Γˡ"],kᵇ,f)
ops[6](elements["Γʳ"],kᵇ,f)
ops[7](elements["Γᵇ"],kᵇ,f)
ops[7](elements["Γᵗ"],kᵇ,f)
ops[7](elements["Γˡ"],kᵇ,f)
ops[7](elements["Γʳ"],kᵇ,f)
# ops[9](elements["Γᵇ"],f)
# ops[9](elements["Γᵗ"],f)
# ops[9](elements["Γˡ"],f)
# ops[9](elements["Γʳ"],f)
# ops[10](elements["Γᵇ"],f)
# ops[10](elements["Γᵗ"],f)
# ops[10](elements["Γʳ"],f)
# ops[11](elements["Γᵇ"],f) 
# ops[11](elements["Γᵗ"],f)
# ops[11](elements["Γʳ"],f)

k = [kᵇ kʷˢ;kʷˢ' kˢˢ]
# k = sparse([kᵇ kʷˢ;kʷˢ' kˢˢ])
f = [f;zeros(2*nˢ)]

# k = kʷˢ*inv(kˢˢ)*kʷˢ'
# k = -kʷˢ*(kˢˢ\kʷˢ')
# a = eigvals(k)
# println(log10(a[3*nᵇ-2nˢ+1]))
# println(a[3*nᵇ-2nˢ+1])

d = k\f
# pardiso(ps,d,k,f)
d₁ = d[1:3:3*nᵇ]
d₂ = d[2:3:3*nᵇ] 
d₃ = d[3:3:3*nᵇ]
s₁ = d[3*nᵇ+1:2:3*nᵇ+2*nˢ]
s₂ = d[3*nᵇ+2:2:3*nᵇ+2*nˢ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
push!(nodes_s,:q₁=>s₁,:q₂=>s₂)
# eval(VTK_mix_pressure)

set𝝭!(elements["Ωᵍ"])
set∇𝝭!(elements["Ωᵍ"])
set𝝭!(elements["Ωᵍˢ"])
set∇𝝭!(elements["Ωᵍˢ"])

prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->w(x,y))
prescribe!(elements["Ωᵍ"],:θ₁=>(x,y,z)->θ₁(x,y))
prescribe!(elements["Ωᵍ"],:θ₂=>(x,y,z)->θ₂(x,y))
prescribe!(elements["Ωᵍˢ"],:Q₁=>(x,y,z)->Q₁(x,y))
prescribe!(elements["Ωᵍˢ"],:Q₂=>(x,y,z)->Q₂(x,y))
L₂_u = ops[8](elements["Ωᵍ"])
L₂_q = ops[9](elements["Ωᵍˢ"])
a = log10(L₂_u)
b = log10(L₂_q)
println(a)
println(b)
# index = 1:20
# XLSX.openxlsx("./xlsx/SquarePlate.xlsx", mode="rw") do xf
#     Sheet = xf[2]
#     ind = findfirst(n->n==ndivs,index)+2
#     Sheet["A"*string(ind)] = nˢ
#     Sheet["B"*string(ind)] = a
#     Sheet["C"*string(ind)] = b
# end

# println(wᶜ)
# e = abs(wᶜ[1]-𝑣)

# index = [8,16,32,64]
# XLSX.openxlsx("./xlsx/SquarePlate.xlsx", mode="rw") do xf
#     Sheet = xf[3]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["B"*string(ind)] = log10(1/ndiv)
#     Sheet["C"*string(ind)] = a
# end

fig = Figure()
ind = 100
ax = Axis(fig[1,1], 
    aspect = DataAspect(), 
    xticksvisible = false,
    xticklabelsvisible=false, 
    yticksvisible = false, 
    yticklabelsvisible=false,
)
hidespines!(ax)
hidedecorations!(ax)
xs = LinRange(0, 1, ind)
ys = LinRange(0, 1, ind)
zs = zeros(ind,ind)
𝗠 = zeros(21)
for (i,x) in enumerate(xs)
    for (j,y) in enumerate(ys)
        indices = sp(x,y,0.0)
        ni = length(indices)
        𝓒 = [nodes_s[i] for i in indices]
        data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
        ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
        𝓖 = [ξ]
        a = type(𝓒,𝓖)
        set𝝭!(a)
        q = 0.0
        N = ξ[:𝝭]
        for (k,xₖ) in enumerate(𝓒)
            q += N[k]*xₖ.q₁
            # q += N[k]*xₖ.q₂
        end
        zs[i,j] = q
    end
end
surface!(xs,ys,zeros(ind,ind),color=zs,colorrange=(-0.000025,0.000025),colormap=:summer)
contour!(xs,ys,zs,levels=-0.000025:0.00000715:0.000025,color=:azure)
# Colorbar(fig[1,2], limits=(-900,900), colormap=:lightrainbow)
# save("./png/SquarePlate_mix_tri3_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)
# save("./png/SquarePlate_mix_tri3_q2_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)
save("./png/SquarePlate_mix_quad8_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)
# save("./png/SquarePlate_mix_quad8_q2_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)

fig