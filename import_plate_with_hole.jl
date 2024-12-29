
using Gmsh, Statistics

function import_fem(filename::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    elements["Ω"] = getElements(nodes,entities["Ω"])
    elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"],normal=true)
    elements["Γᵍ"] = getElements(nodes,entities["Γᵍ"],normal=true)

    push!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
    push!(elements["Γᵗ"],:𝝭)
    push!(elements["Γᵍ"],:𝝭)

    gmsh.finalize()

    set∇𝝭!(elements["Ω"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γᵍ"])

    return elements, nodes
end

function import_linear_mix(filename1::String,filename2::String,n::Int)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_u = get𝑿ᵢ()
    xᵘ = nodes_u.x
    yᵘ = nodes_u.y
    zᵘ = nodes_u.z
    s = zeros(length(nodes_u))
    
    for (i,node) in enumerate(nodes_u) 
        xᵢ = node.x
        yᵢ = node.y
        r = (xᵢ^2+yᵢ^2)^0.5
        θ = atan(yᵢ/xᵢ)
        s₀ = 0.25π*r/n
        s[i] = s₀ + 2.0*s₀*(cos(θ)+sin(θ)-1.0)
    end
    s .*= 1.5
    push!(nodes_u,:s₁=>s,:s₂=>s,:s₃=>s)

    integrationOrder_Ω = 2
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 2

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ωᵖ"] = getElements(nodes,entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍᵖ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["∂Ωᵖ"] = getElements(nodes, entities["Γ"],   integrationOrder_Γ, normal = true)
    elements["Γᵍᵖ"] = getElements(nodes,entities["Γᵍ"], integrationOrder_Γ, normal = true)

    push!(elements["Ωᵖ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍᵖ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["∂Ωᵖ"],:𝝭)
    push!(elements["Γᵍᵖ"],:𝝭)

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵘ,yᵘ,zᵘ,n = 3,γ = 5)
    elements["Ωᵘ"] = getElements(nodes_u, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["∂Ωᵘ"] = getElements(nodes_u, entities["Γ"], type, integrationOrder_Γ, sp, normal = true)
    elements["Ωᵍᵘ"] = getElements(nodes_u, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γᵗ"] = getElements(nodes_u, entities["Γᵗ"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γᵍᵘ"] = getElements(nodes_u, entities["Γᵍ"],type,  integrationOrder_Γ, sp, normal = true)

    nₘ = 6
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵘ"], :𝝭)
    push!(elements["∂Ωᵘ"], :𝝭)
    push!(elements["Γᵗ"], :𝝭)
    push!(elements["Γᵍᵘ"], :𝝭)
    push!(elements["Ωᵘ"],  :𝗠=>𝗠)
    push!(elements["∂Ωᵘ"], :𝗠=>𝗠)
    push!(elements["Γᵗ"], :𝗠=>𝗠)
    push!(elements["Γᵍᵘ"], :𝗠=>𝗠)
    push!(elements["Ωᵍᵘ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y, :∂𝝭∂z)
    push!(elements["Ωᵍᵘ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    set∇𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["∂Ωᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵍᵖ"])
    set𝝭!(elements["Ωᵘ"])
    set𝝭!(elements["∂Ωᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γᵍᵘ"])

    gmsh.finalize()

    return elements, nodes, nodes_u
end

function import_elasticity_linear_mix(filename1::String,filename2::String,n::Int)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    s = zeros(length(nodes_p))
    
    for (i,node) in enumerate(nodes_p) 
        xᵢ = node.x
        yᵢ = node.y
        r = (xᵢ^2+yᵢ^2)^0.5
        θ = atan(yᵢ/xᵢ)
        s₀ = 0.25π*r/n
        s[i] = s₀ + 2.0*s₀*(cos(θ)+sin(θ)-1.0)
    end
    s .*= 3.0
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    integrationOrder_Ω = 4
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 4

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ωᵘ"] = getElements(nodes,entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γᵍᵘ"] = getElements(nodes,entities["Γᵍ"], integrationOrder_Γ, normal = true)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],  integrationOrder_Γ, normal = true)

    push!(elements["Ωᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Γᵍᵘ"],:𝝭)
    push!(elements["Γᵗ"], :𝝭)

    set∇𝝭!(elements["Ωᵘ"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γᵍᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"],type,  integrationOrder_Γ, sp, normal = true)

    
    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γᵍᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Γᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    set∇𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["Γᵍᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])

    elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"],   integrationOrder_Γ, normal = true)
    push!(elements["∂Ωᵘ"],:𝝭)
    set𝝭!(elements["∂Ωᵘ"])
    elements["∂Ωᵖ"] = getElements(nodes_p, entities["Γ"], type, integrationOrder_Γ, sp)
    push!(elements["∂Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    set∇𝝭!(elements["∂Ωᵖ"])

    # type = PiecewisePolynomial{:Constant}
    type = PiecewisePolynomial{:Linear2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], type, integrationOrder_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], type, integrationOrder_Γ)
    elements["Γᵍˢ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωˢ"])
    
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)

    set∇𝝭!(elements["Ωˢ"])
    set𝝭!(elements["∂Ωˢ"])

    gmsh.finalize()

    return elements, nodes, nodes_p
end


function import_HR_mix(filename1::String,filename2::String,n::Int)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    nodes_p = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    Ω = getElements(nodes, entities["Ω"])
    sᵤ = zeros(length(nodes))
    sₚ = zeros(length(nodes))
    for (i,node) in enumerate(nodes) 
        xᵢ = node.x
        yᵢ = node.y
        r = (xᵢ^2+yᵢ^2)^0.5
        θ = atan(yᵢ/xᵢ)
        s₀ = 0.25π*r/n
        sᵤ[i] = s₀ + 0.9*s₀*(cos(θ)+sin(θ)-1.0)
        sₚ[i] = s₀ + 0.9*s₀*(cos(θ)+sin(θ)-1.0)
    end
    sᵤ .*= 1.5
    sₚ .*= 1.5
    push!(nodes,:s₁=>sᵤ,:s₂=>sᵤ,:s₃=>sᵤ)
    push!(nodes_p,:s₁=>sₚ,:s₂=>sₚ,:s₃=>sₚ)

    integrationOrder_Ω = 6
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 6

   
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    elements["Ωᵘ"] = getElements(nodes,entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["Γᵍᵘ"] = getElements(nodes,entities["Γᵍ"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"], type,  integrationOrder_Γ, sp, normal = true)

    # elements["Ωᵘ"] = getElements(nodes,entities["Ω"], integrationOrder_Ω)
    # elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    # elements["Γᵍᵘ"] = getElements(nodes,entities["Γᵍ"],  integrationOrder_Γ,  normal = true)
    # elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],   integrationOrder_Γ,  normal = true)
    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Γᵗ"],:𝝭)
    push!(elements["Γᵍᵘ"],:𝝭)

    push!(elements["Ωᵘ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵘ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γᵗ"],:𝗠=>𝗠)
    push!(elements["Γᵍᵘ"],:𝗠=>𝗠)

    set∇𝝭!(elements["Ωᵘ"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γᵍᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])

    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"],type,  integrationOrder_Γ, sp, normal = true)

    
    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γᵍᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Γᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    set∇𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["Γᵍᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])

    elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"], type,  integrationOrder_Γ,sp, normal = true)
    # elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"],  integrationOrder_Γ, normal = true)
    push!(elements["∂Ωᵘ"],:𝝭)
    push!(elements["∂Ωᵘ"],:𝗠=>𝗠)
    set𝝭!(elements["∂Ωᵘ"])
    elements["∂Ωᵖ"] = getElements(nodes_p, entities["Γ"], type, integrationOrder_Γ, sp)
    push!(elements["∂Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    set∇𝝭!(elements["∂Ωᵖ"])

    types = PiecewisePolynomial{:Constant}
    # types = PiecewisePolynomial{:Linear2D}
    # types = PiecewisePolynomial{:Quadratic2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], types, integrationOrder_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], types, integrationOrder_Γ)
    elements["Γᵍˢ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωˢ"])
    
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)

    set∇𝝭!(elements["Ωˢ"])
    set𝝭!(elements["∂Ωˢ"])

    gmsh.finalize()

    return elements, nodes, nodes_p,Ω
end
function import_plate_with_hole_mix(filename1::String,filename2::String,n,c)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    w = 0.0
    for i in 0:n-1
        w += c^i
    end
    ds₂ = 4*2^0.5/w
    ds₁ = ds₂*c^(n-1)
    s = zeros(length(nodes))
    for (i,node) in enumerate(nodes) 
        xᵢ = node.x
        yᵢ = node.y
        r = (xᵢ^2+yᵢ^2)^0.5
        s[i] = ds₁ + (r-1)/4/2^0.5*(ds₂-ds₁)
    end
    s .*= 2.5
    push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

    integration_Ω = 4
    integrationOrder_Ωᵍ = 8
    integration_Γ = 4

    gmsh.open(filename2)
    entities = getPhysicalGroups()

    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    elements["Ω"] = getElements(nodes, entities["Ω"], type, integration_Ω, sp)
    elements["∂Ω"] = getElements(nodes, entities["Γ"], type, integration_Γ, sp, normal = true)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"],type, integration_Γ, sp, normal = true)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],type, integration_Γ, sp, normal = true)

    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ω"], :𝝭)
    push!(elements["∂Ω"], :𝝭)
    push!(elements["Γᵍ"], :𝝭)
    push!(elements["Γᵗ"], :𝝭)
    push!(elements["Ω"],  :𝗠=>𝗠)
    push!(elements["∂Ω"], :𝗠=>𝗠)
    push!(elements["Γᵍ"], :𝗠=>𝗠)
    push!(elements["Γᵗ"], :𝗠=>𝗠)
    push!(elements["Ωᵍ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)



    # type = PiecewisePolynomial{:Linear2D}
    type = PiecewisePolynomial{:Quadratic2D}
    println(entities)
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], type, integration_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], type, integration_Γ)
    elements["Γᵍˢ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γᵗˢ"] = getElements(entities["Γᵗ"],entities["Γ"], elements["∂Ωˢ"])
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)

    # gmsh.finalize()

    return elements, nodes, ds₂, ds₁
end
function import_HR_mix_xyt(filename1::String,filename2::String,n,c)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    w = 0.0
    for i in 0:n-1
        w += c^i
    end
    ds₂ = 4*2^0.5/w
    ds₁ = ds₂*c^(n-1)
    s = zeros(length(nodes))
    for (i,node) in enumerate(nodes) 
        xᵢ = node.x
        yᵢ = node.y
        r = (xᵢ^2+yᵢ^2)^0.5
        s[i] = ds₁ + (r-1)/4/2^0.5*(ds₂-ds₁)
    end
    s .*= 2.5
    push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

    integration_Ω = 2
    integrationOrder_Ωᵍ = 8
    integration_Γ = 2

    gmsh.open(filename2)
    entities = getPhysicalGroups()

    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    elements["Ωᵘ"] = getElements(nodes, entities["Ω"], type, integration_Ω, sp)
    elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"], type, integration_Γ, sp, normal = true)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["Γᵍᵘ"] = getElements(nodes, entities["Γᵍ"],type, integration_Γ, sp, normal = true)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],type, integration_Γ, sp, normal = true)

    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵘ"], :𝝭)
    push!(elements["∂Ωᵘ"], :𝝭)
    push!(elements["Γᵍᵘ"], :𝝭)
    push!(elements["Γᵗ"], :𝝭)
    push!(elements["Ωᵘ"],  :𝗠=>𝗠)
    push!(elements["∂Ωᵘ"], :𝗠=>𝗠)
    push!(elements["Γᵍᵘ"], :𝗠=>𝗠)
    push!(elements["Γᵗ"], :𝗠=>𝗠)
    push!(elements["Ωᵍᵘ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍᵘ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)



    type = PiecewisePolynomial{:Linear2D}
    # type = PiecewisePolynomial{:Quadratic2D}
    println(entities)
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], type, integration_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], type, integration_Γ)
    elements["Γᵍˢ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γᵗˢ"] = getElements(entities["Γᵗ"],entities["Γ"], elements["∂Ωˢ"])
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)

    # gmsh.finalize()

    return elements, nodes
end

function import_HR_mix_old(filename1::String,filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    nodes_p = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    Ω = getElements(nodes, entities["Ω"])
    
    s, var𝐴 = cal_area_support(Ω)
    sᵤ = 1.5*s*ones(length(nodes))
    sₚ= 1.0*s*ones(length(nodes))
    push!(nodes,:s₁=>sᵤ,:s₂=>sᵤ,:s₃=>sᵤ)
    push!(nodes_p,:s₁=>sₚ,:s₂=>sₚ,:s₃=>sₚ)

    integrationOrder_Ω = 4
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 4

   
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    elements["Ωᵘ"] = getElements(nodes,entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["Γ¹ᵍᵘ"] = getElements(nodes,entities["Γᵍ₁"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ²ᵍᵘ"] = getElements(nodes,entities["Γᵍ₂"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ³ᵍᵘ"] = getElements(nodes,entities["Γᵍ₃"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ¹ᵗᵘ"] = getElements(nodes, entities["Γᵗ₁"], type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ²ᵗᵘ"] = getElements(nodes, entities["Γᵗ₂"], type,  integrationOrder_Γ, sp, normal = true)
    elements["Γᵍᵘ"] = elements["Γ¹ᵍᵘ"]∪elements["Γ²ᵍᵘ"]∪elements["Γ³ᵍᵘ"]
    elements["Γᵗ"] = elements["Γ¹ᵗᵘ"]∪elements["Γ²ᵗᵘ"]

    # elements["Ωᵘ"] = getElements(nodes,entities["Ω"], integrationOrder_Ω)
    # elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    # elements["Γᵍᵘ"] = getElements(nodes,entities["Γᵍ"],  integrationOrder_Γ,  normal = true)
    # elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],   integrationOrder_Γ,  normal = true)
    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    
    # push!(elements["Γ¹ᵗᵘ"], :𝝭=>:𝑠)
    # push!(elements["Γ²ᵗᵘ"], :𝝭=>:𝑠)
    # push!(elements["Γ¹ᵍᵘ"], :𝝭=>:𝑠)
    # push!(elements["Γ²ᵍᵘ"], :𝝭=>:𝑠)
    # push!(elements["Γ³ᵍᵘ"], :𝝭=>:𝑠)
    
    # push!(elements["Γ¹ᵗᵘ"], :𝗠=>𝗠)
    # push!(elements["Γ²ᵗᵘ"], :𝗠=>𝗠)
    # push!(elements["Γ¹ᵍᵘ"], :𝗠=>𝗠)
    # push!(elements["Γ²ᵍᵘ"], :𝗠=>𝗠)
    # push!(elements["Γ³ᵍᵘ"], :𝗠=>𝗠)
   
    push!(elements["Ωᵘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵘ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵘ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γ¹ᵍᵖ"] = getElements(nodes_p, entities["Γᵍ₁"], type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ²ᵍᵖ"] = getElements(nodes_p, entities["Γᵍ₂"], type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ³ᵍᵖ"] = getElements(nodes_p, entities["Γᵍ₃"], type,  integrationOrder_Γ, sp, normal = true)
    elements["Γᵍᵖ"] = elements["Γ¹ᵍᵖ"]∪elements["Γ²ᵍᵖ"]∪elements["Γ³ᵍᵖ"]

    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γ¹ᵍᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Γ¹ᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γ²ᵍᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Γ²ᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γ³ᵍᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Γ³ᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    set∇𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["Γᵍᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])

    # elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"], type,  integrationOrder_Γ,sp, normal = true)
    elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"],  integrationOrder_Γ, normal = true)
    push!(elements["∂Ωᵘ"],:𝝭)
    push!(elements["∂Ωᵘ"],:𝗠=>𝗠)
    set𝝭!(elements["∂Ωᵘ"])
    elements["∂Ωᵖ"] = getElements(nodes_p, entities["Γ"], type, integrationOrder_Γ, sp)
    push!(elements["∂Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    set∇𝝭!(elements["∂Ωᵖ"])

    # types = PiecewisePolynomial{:Constant}
    # types = PiecewisePolynomial{:Linear2D}
    types = PiecewisePolynomial{:Quadratic2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], types, integrationOrder_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], types, integrationOrder_Γ)
    elements["Γ¹ᵍˢ"] = getElements(entities["Γᵍ₁"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γ²ᵍˢ"] = getElements(entities["Γᵍ₂"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γ³ᵍˢ"] = getElements(entities["Γᵍ₃"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γᵍˢ"] =  elements["Γ¹ᵍˢ"]∪elements["Γ²ᵍˢ"]∪elements["Γ³ᵍˢ"]

    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)

    set∇𝝭!(elements["Ωˢ"])
    set𝝭!(elements["∂Ωˢ"])

    gmsh.finalize()

    return elements, nodes, nodes_p
end

function import_elasticity_quadratic_mix(filename1::String,filename2::String,n::Int)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename2)
    # entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    s = zeros(length(nodes_p))
    
    for (i,node) in enumerate(nodes_p) 
        xᵢ = node.x
        yᵢ = node.y
        r = (xᵢ^2+yᵢ^2)^0.5
        θ = atan(yᵢ/xᵢ)
        s₀ = 0.25π*r/n
        s[i] = s₀ + 2.0*s₀*(cos(θ)+sin(θ)-1.0)
    end
    s .*= 2.5
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    integrationOrder_Ω = 4
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 4

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ωᵘ"] = getElements(nodes,entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γᵍᵘ"] = getElements(nodes,entities["Γᵍ"], integrationOrder_Γ, normal = true)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],  integrationOrder_Γ, normal = true)

    push!(elements["Ωᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Γᵍᵘ"],:𝝭)
    push!(elements["Γᵗ"], :𝝭)
    set∇𝝭!(elements["Ωᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])
    set𝝭!(elements["Γᵍᵘ"])
    set𝝭!(elements["Γᵗ"])

    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"],type,  integrationOrder_Γ, sp, normal = true)
    # elements["∂Ωᵖ"] = getElements(nodes_p, entities["Γ"], type, integrationOrder_Γ, sp)
    
    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωᵖ"], :𝝭)
    push!(elements["Γᵍᵖ"], :𝝭)
    push!(elements["Ωᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["∂Ωᵖ"], :𝗠=>𝗠)
    push!(elements["Γᵍᵖ"], :𝗠=>𝗠)

    push!(elements["Ωᵍᵖ"], :𝝭)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠)

    set∇𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵍᵖ"])

    filename1s = split(filename1,"_")
    if filename1s[2] == "quad8"
        filename3 = replace(filename1,"quad8"=>"quad")
        gmsh.open(filename3)
        entities = getPhysicalGroups()
    end

    elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"],   integrationOrder_Γ, normal = true)
    push!(elements["∂Ωᵘ"],:𝝭)
    set𝝭!(elements["∂Ωᵘ"])
    elements["∂Ωᵖ"] = getElements(nodes_p, entities["Γ"], type, integrationOrder_Γ, sp)
    push!(elements["∂Ωᵖ"], :𝝭)
    push!(elements["∂Ωᵖ"], :𝗠=>𝗠)
    set𝝭!(elements["∂Ωᵖ"])

    type = PiecewisePolynomial{:Linear2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], type, integrationOrder_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], type, integrationOrder_Γ)
    elements["Γᵍˢ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωˢ"])
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)

    set∇𝝭!(elements["Ωˢ"])
    set𝝭!(elements["∂Ωˢ"])

    gmsh.finalize()

    return elements, nodes, nodes_p
end

function import_plate_with_hole_gauss(filename::String,n,c)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename)
    entities = getPhysicalGroups()
    
    nodes = get𝑿ᵢ()
    Ω = getElements(nodes, entities["Ω"])
    x = nodes.x
    y = nodes.y
    z = nodes.z 
    w = 0.0
    for i in 0:n-1
        w += c^i
    end
    ds₂ = 4*2^0.5/w
    ds₁ = ds₂*c^(n-1)
    s = zeros(length(nodes))
    for (i,node) in enumerate(nodes) 
        xᵢ = node.x
        yᵢ = node.y
        r = (xᵢ^2+yᵢ^2)^0.5
        s[i] = ds₁ + (r-1)/4/2^0.5*(ds₂-ds₁)
    end
    s .*= 2.5
    push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

    integration_Ω = 7
    integrationOrder_Ωᵍ = 8
    integration_Γ = 7

    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    elements["Ω"] = getElements(nodes, entities["Ω"], type, integration_Ω, sp)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"],type, integration_Γ, sp, normal = true)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],type, integration_Γ, sp, normal = true)

    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    ∂𝗠∂z = zeros(nₘ)
    push!(elements["Ω"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Γᵍ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Γᵗ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ω"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γᵍ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γᵗ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y, :∂𝝭∂z)
    push!(elements["Ωᵍ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂𝗠∂z=>∂𝗠∂z)

    set∇𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    set∇𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])
    # gmsh.finalize()

    return elements, nodes, Ω
end

function cal_area_support(elms::Vector{ApproxOperator.AbstractElement})
    𝐴s = zeros(length(elms))
    for (i,elm) in enumerate(elms)
        x₁ = elm.𝓒[1].x
        y₁ = elm.𝓒[1].y
        x₂ = elm.𝓒[2].x
        y₂ = elm.𝓒[2].y
        x₃ = elm.𝓒[3].x
        y₃ = elm.𝓒[3].y
        𝐴s[i] = 0.5*(x₁*y₂ + x₂*y₃ + x₃*y₁ - x₂*y₁ - x₃*y₂ - x₁*y₃)
    end
    avg𝐴 = mean(𝐴s)
    var𝐴 = var(𝐴s)
    s = (4/3^0.5*avg𝐴)^0.5
    return s, var𝐴
end
