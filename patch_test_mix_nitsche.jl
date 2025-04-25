
using ApproxOperator, Tensors, JLD, LinearAlgebra

import BenchmarkExample: BenchmarkExample
include("import_prescrible_ops.jl")
include("import_patch_test.jl")
elements, nodes = import_roof_mix_nitsche("msh/patchtest.msh")

nₘ = 55
# nₘ = 21
𝑅 = BenchmarkExample.PatchTestThinShell.𝑅
𝐿 = BenchmarkExample.PatchTestThinShell.𝐿
E = BenchmarkExample.PatchTestThinShell.𝐸
ν = BenchmarkExample.PatchTestThinShell.𝜈
h = BenchmarkExample.PatchTestThinShell.ℎ

# cs = BenchmarkExample.sphericalCoordinate(𝑅)
# cs = BenchmarkExample.cylindricalCoordinate(𝑅)
# cs = BenchmarkExample.cartesianCoordinate()

nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
s = 2.5*0.1*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

n = 1
u(x) = Vec{3}(((x[1]+2*x[2])^n,(3*x[1]+4*x[2])^n,(5*x[1]+6*x[2])^n))
vs = BenchmarkExample.PatchTestThinShell.variables(cs,u)

eval(prescribeForMix)
eval(prescribeVariables)
eval(prescribeForNitsche)

set𝝭!(elements["Ω"])
set∇²𝝭!(elements["Ωₚ"])
set∇𝝭!(elements["Γₚ"])
set∇𝝭!(elements["Ωₘ"])
set∇𝝭!(elements["Γₘ"])
set∇̂³𝝭!(elements["Γ¹"])
set∇̂³𝝭!(elements["Γ²"])
set∇̂³𝝭!(elements["Γ³"])
set∇̂³𝝭!(elements["Γ⁴"])
set∇𝝭!(elements["Γ¹"])
set∇𝝭!(elements["Γ²"])
set∇𝝭!(elements["Γ³"])
set∇𝝭!(elements["Γ⁴"])

eval(opsMix)

opForce = Operator{:∫vᵢbᵢdΩ}()
f = zeros(3*nₚ)
opForce(elements["Ωₘ"],f)

kᵋᵋ = zeros(6*nᵥ,6*nᵥ)
kᴺᵋ = zeros(6*nᵥ,6*nᵥ)
kᴺᵛ = zeros(6*nᵥ,3*nₚ)
kᵏᵏ = zeros(9*nᵥ,9*nᵥ)
kᴹᵏ = zeros(9*nᵥ,9*nᵥ)
kᴹᵛ = zeros(9*nᵥ,3*nₚ)

ops[1](elements["Ω"],kᵋᵋ)
ops[2](elements["Ω"],kᴺᵋ)
ops[3](elements["Γₚ"],elements["Γₘ"],kᴺᵛ)
ops[4](elements["Ωₚ"],elements["Ωₘ"],kᴺᵛ)

ops[5](elements["Ω"],kᵏᵏ)
ops[6](elements["Ω"],kᴹᵏ)
ops[7](elements["Γₚ"],elements["Γₘ"],kᴹᵛ)
ops[8](elements["Γₚ"],elements["Γₘ"],kᴹᵛ)
ops[9](elements["Γₚ"],elements["Γₘ"],kᴹᵛ)
ops[10](elements["Ωₚ"],elements["Ωₘ"],kᴹᵛ)

αᵥ = 1e5
αᵣ = 1e3
eval(opsPenalty)
eval(opsNitsche)
kᵛ = zeros(3*nₚ,3*nₚ)
fᵛ = zeros(3*nₚ)
opsv[1](elements["Γ¹"],kᵛ,fᵛ)
opsv[1](elements["Γ²"],kᵛ,fᵛ)
opsv[1](elements["Γ³"],kᵛ,fᵛ)
opsv[1](elements["Γ⁴"],kᵛ,fᵛ)
opsv[2](elements["Γ¹"],kᵛ,fᵛ)
opsv[2](elements["Γ²"],kᵛ,fᵛ)
opsv[2](elements["Γ³"],kᵛ,fᵛ)
opsv[2](elements["Γ⁴"],kᵛ,fᵛ)
opsv[3](elements["Γ¹"],kᵛ,fᵛ)
opsv[3](elements["Γ²"],kᵛ,fᵛ)
opsv[3](elements["Γ³"],kᵛ,fᵛ)
opsv[3](elements["Γ⁴"],kᵛ,fᵛ)
opsv[4](elements["Γ¹"],kᵛ,fᵛ)
opsv[4](elements["Γ²"],kᵛ,fᵛ)
opsv[4](elements["Γ³"],kᵛ,fᵛ)
opsv[4](elements["Γ⁴"],kᵛ,fᵛ)

kᵅ = zeros(3*nₚ,3*nₚ)
fᵅ = zeros(3*nₚ)
opsα[1](elements["Γ¹"],kᵅ,fᵅ)
opsα[1](elements["Γ²"],kᵅ,fᵅ)
opsα[1](elements["Γ³"],kᵅ,fᵅ)
opsα[1](elements["Γ⁴"],kᵅ,fᵅ)
opsα[2](elements["Γ¹"],kᵅ,fᵅ)
opsα[2](elements["Γ²"],kᵅ,fᵅ)
opsα[2](elements["Γ³"],kᵅ,fᵅ)
opsα[2](elements["Γ⁴"],kᵅ,fᵅ)

kᵋᵛ = kᴺᵋ\kᴺᵛ
kᵏᵛ = kᴹᵏ\kᴹᵛ
d = (kᵋᵛ'*kᵋᵋ*kᵋᵛ + kᵏᵛ'*kᵏᵏ*kᵏᵛ - kᵛ - kᵅ)\(-f - fᵛ - fᵅ)

d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
@save compress=true "jld/patchtest_mix_nitsche.jld" d₁ d₂ d₃

set∇²𝝭!(elements["Ωᵍ"])
opE = Operator{:Hₑ_ThinShell}(:E=>E,:ν=>ν,:h=>h)
Hₑ, L₂ = opE(elements["Ωᵍ"])

# dᵗ = zeros(3*nₚ)
# for (i,node) in enumerate(nodes)
#     x = Vec{3}((node.x,node.y,node.z))
#     dᵗ[3*i-2] = 0.0
#     dᵗ[3*i-1] = 0.0
#     dᵗ[3*i]   = (node.x+node.y)^2
# end

# dᵛ = zeros(3*nₚ)
# for (i,node) in enumerate(nodes)
#     x = Vec{3}((node.x,node.y,node.z))
#     u_ = u(x)
#     dᵛ[3*i-2] = u_[1]
#     dᵛ[3*i-1] = u_[2]
#     dᵛ[3*i]   = u_[3]
# end

# d₁ = zeros(nₚ)
# d₂ = zeros(nₚ)
# d₃ = zeros(nₚ)
# for (i,node) in enumerate(nodes)
#     x = Vec{3}((node.x,node.y,node.z))
#     u_ = u(x)
#     d₁[i] = u_[1]
#     d₂[i] = u_[2]
#     d₃[i] = u_[3]
# end

# push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)

# opsv[1](elements["Γ¹"],kᵛ,fᵛ)
# opsv[1](elements["Γ²"],kᵛ,fᵛ)
# opsv[1](elements["Γ³"],kᵛ,fᵛ)
# opsv[1](elements["Γ⁴"],kᵛ,fᵛ)

# kᵛ1 = zeros(3*nₚ,3*nₚ)
# fᵛ1 = zeros(3*nₚ)
# kᵛ2 = zeros(3*nₚ,3*nₚ)
# fᵛ2 = zeros(3*nₚ)
# kᵛ3 = zeros(3*nₚ,3*nₚ)
# fᵛ3 = zeros(3*nₚ)
# opsv[1](elements["Γ¹"],kᵛ1,fᵛ1)
# opsv[1](elements["Γ²"],kᵛ1,fᵛ1)
# opsv[1](elements["Γ³"],kᵛ1,fᵛ1)
# opsv[1](elements["Γ⁴"],kᵛ1,fᵛ1)
# opsv[2](elements["Γ¹"],kᵛ2,fᵛ2)
# opsv[2](elements["Γ²"],kᵛ2,fᵛ2)
# opsv[2](elements["Γ³"],kᵛ2,fᵛ2)
# opsv[2](elements["Γ⁴"],kᵛ2,fᵛ2)
# opsv[3](elements["Γ¹"],kᵛ2,fᵛ2)
# opsv[3](elements["Γ²"],kᵛ2,fᵛ2)
# opsv[3](elements["Γ³"],kᵛ2,fᵛ2)
# opsv[3](elements["Γ⁴"],kᵛ2,fᵛ2)
# opsv[4](elements["Γ¹"],kᵛ2,fᵛ2)
# opsv[4](elements["Γ²"],kᵛ2,fᵛ2)
# opsv[4](elements["Γ³"],kᵛ2,fᵛ2)
# opsv[4](elements["Γ⁴"],kᵛ2,fᵛ2)
# opsv[4](elements["Γ¹"],kᵛ3,fᵛ3)
# opsv[4](elements["Γ²"],kᵛ3,fᵛ3)
# opsv[4](elements["Γ³"],kᵛ3,fᵛ3)
# opsv[4](elements["Γ⁴"],kᵛ3,fᵛ3)
# fex = (kᵋᵛ'*kᵋᵋ*kᵋᵛ + kᵏᵛ'*kᵏᵏ*kᵏᵛ - kᵛ)*dᵛ - (-f - fᵛ)
# fex1 = (kᵋᵛ'*kᵋᵋ*kᵋᵛ- kᵛ1)*dᵛ - (-f - fᵛ1)
# fex2 = (kᵏᵛ'*kᵏᵏ*kᵏᵛ-kᵛ2)*dᵛ - (- fᵛ2)
# fex2 = (kᵏᵛ'*kᵏᵏ*kᵏᵛ)*dᵛ
# fex1 = (kᵋᵛ'*kᵋᵋ*kᵋᵛ)*dᵛ - (-f)
# fex2 = kᵛ1*dᵛ - fᵛ1
# fex1 = kᵛ2*dᵛ - fᵛ2
# dᵗ'*fᵛ2/2
# fex1 = kᵛ*dᵛ
# fex2 = fᵛ
# fex1 = kᵏᵛ'*kᵏᵏ*kᵏᵛ*dᵛ
# fex2 = kᵛ2*dᵛ - fᵛ2
# fex3 = kᵛ3*dᵛ - fᵛ3
# fex3 = fex1-fex2
# fex1 = kᵛ2*dᵛ
# println(norm(fex))
# println(norm(fex1))
# println(norm(fex2))

# err = ones(3*nₚ)'*(kᵛ*dᵛ - fᵛ)

# dᵋ = zeros(6*nᵥ)
# for i in 1:length(elements["Ω"])
#     dᵋ[18*i-17] = 1.0
#     dᵋ[18*i-14] = 2.0
# end

# dᴺ = zeros(6*nᵥ)
# for i in 1:length(elements["Ω"])
#     dᴺ[18*i-17] = E*h
#     dᴺ[18*i-16] = E*h
#     dᴺ[18*i-14] = E*h
# end

# dᵏ = zeros(9*nᵥ)
# for i in 1:length(elements["Ω"])
#     dᵏ[27*i-24] = -2.
#     dᵏ[27*i-21] = -2.
#     dᵏ[27*i-18] = -4.
#     dᵏ[27*i-24] = 1.
# end

# dᴹ = zeros(9*nᵥ)
# for i in 1:length(elements["Ω"])
#     dᴹ[27*i-24] = -E*h^3/12*2
#     dᴹ[27*i-21] = -E*h^3/12*2
#     dᴹ[27*i-18] = -E*h^3/12*2
#     dᴹ[27*i-24] = E*h^3/12
# end
# println(norm(kᴹᵛ'*dᴹ))
# println(norm(kᵏᵏ*dᵏ+kᴹᵏ'*dᴹ))
# println(norm(kᴹᵛ*dᵛ+kᴹᵏ*dᵏ-fᴹ))
# f1 = kᴹᵛ*dᵛ
# f2 = kᴹᵏ*dᵏ

# println(norm(kᴺᵛ'*dᴺ+kᴹᵛ'*dᴹ))
# println(norm(kᵋᵋ*dᵋ+kᴺᵋ'*dᴺ))
# println(norm(kᵏᵏ*dᵏ+kᴹᵏ'*dᴹ))
# println(norm(kᴺᵛ*dᵛ+kᴺᵋ*dᵋ-fᴺ))
# println(norm(kᴹᵛ*dᵛ+kᴹᵏ*dᵏ-fᴹ))
