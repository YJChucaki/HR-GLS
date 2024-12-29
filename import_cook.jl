
using Gmsh, Statistics

function import_fem(filename::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    elements["Ω"] = getElements(nodes,entities["Ω"])
    elements["Ωᵍ"] = getElements(nodes,entities["Ω"],8)
    elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"],normal=true)
    elements["Γᵍ"] = getElements(nodes,entities["Γᵍ"],normal=true)
    elements["Γʳ"] = getElements(nodes,entities["Γʳ"],normal=true)

    push!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Γᵗ"],:𝝭)
    push!(elements["Γᵍ"],:𝝭)
    push!(elements["Γʳ"],:𝝭)

    gmsh.finalize()

    set∇𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γʳ"])

    return elements, nodes
end

function import_linear_mix(filename1::String,filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    Ω = getElements(nodes_p, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 2.5*s*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    integrationOrder_Ω = 3
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 2

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ωᵘ"] = getElements(nodes,entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"],   integrationOrder_Γ, normal = true)
    elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"], integrationOrder_Γ, normal = true)
    elements["Γʳ"] = getElements(nodes,entities["Γʳ"], integrationOrder_Γ, normal = true)
    elements["Γᵍᵘ"] = getElements(nodes,entities["Γᵍ"], integrationOrder_Γ, normal = true)

    push!(elements["Ωᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["∂Ωᵘ"],:𝝭)
    push!(elements["Γᵗ"],:𝝭)
    push!(elements["Γʳ"],:𝝭)
    push!(elements["Γᵍᵘ"],:𝝭)

    set∇𝝭!(elements["Ωᵘ"])
    set𝝭!(elements["∂Ωᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γʳ"])
    set𝝭!(elements["Γᵍᵘ"])

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["∂Ωᵖ"] = getElements(nodes_p, entities["Γ"], type, integrationOrder_Γ, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"],type,  integrationOrder_Γ, sp, normal = true)

    nₘ = 6
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
    set𝝭!(elements["∂Ωᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵍᵖ"])

    # types = PiecewisePolynomial{:Constant}
    types = PiecewisePolynomial{:Linear2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], types, integrationOrder_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], types, integrationOrder_Γ)
    elements["Γᵍˢ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωˢ"])
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)

    set∇𝝭!(elements["Ωˢ"])
    set𝝭!(elements["∂Ωˢ"])

    typeb = PiecewiseParametric{:Bubble,:Tri3}
    elements["Ωᵇ"] = getPiecewiseElements(entities["Ω"],typeb,integrationOrder_Ω)
    push!(elements["Ωᵇ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    set∇𝝭!(elements["Ωᵇ"])

    gmsh.finalize()

    return elements, nodes, nodes_p, sp, type
end

function import_HR_GLS(filename1::String,filename2::String,n)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_c = get𝑿ᵢ()
    

    elements["Ω"] = getElements(nodes_c,entities["Ω"])
    push!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    set∇𝝭!(elements["Ω"])


    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    Ω = getElements(nodes, entities["Ω"])
    # s, var𝐴 = cal_area_support(Ω)
    # sᵤ = 1.5*s*ones(length(nodes))

    s = 2.5
    s₁ = s*44.0/n*ones(length(nodes))
    s₂ = s*44.0/n*ones(length(nodes))
    push!(nodes,:s₁=>s₁,:s₂=>s₂,:s₃=>s₂)
    
    integrationOrder_Ω = 4
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 4

    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    elements["Ωᵘ"] = getElements(nodes,entities["Ω"], type,  integrationOrder_Ω, sp)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"], type,   integrationOrder_Γ, sp, normal = true)
    elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"], type, integrationOrder_Γ, sp, normal = true)
    # elements["Γʳ"] = getElements(nodes,entities["Γʳ"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γᵍᵘ"] = getElements(nodes,entities["Γᵍ"], type, integrationOrder_Γ, sp, normal = true)
    
    
    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵘ"],:𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍᵘ"],:𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωᵘ"],:𝝭)
    push!(elements["Γᵗ"],:𝝭)
    # push!(elements["Γʳ"],:𝝭)
    push!(elements["Γᵍᵘ"],:𝝭)

    push!(elements["Ωᵘ"],  :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵘ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["∂Ωᵘ"],:𝗠=>𝗠)
    push!(elements["Γᵗ"],:𝗠=>𝗠)
    # push!(elements["Γʳ"],:𝗠=>𝗠)
    push!(elements["Γᵍᵘ"],:𝗠=>𝗠)

    set∇𝝭!(elements["Ωᵘ"])
    # set∇²𝝭!(elements["Ωᵘ"])
    set𝝭!(elements["∂Ωᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])
    set𝝭!(elements["Γᵗ"])
    # set𝝭!(elements["Γʳ"])
    set𝝭!(elements["Γᵍᵘ"])


    # gmsh.open(filename2)
    # types = PiecewisePolynomial{:Constant}
    # types = PiecewisePolynomial{:Linear2D}
    types = PiecewisePolynomial{:Quadratic2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], types, integrationOrder_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], types, integrationOrder_Γ)
    elements["Γᵍˢ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωˢ"])
    
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y, :∂²𝝭∂x², :∂²𝝭∂y², :∂²𝝭∂x∂y)
    push!(elements["∂Ωˢ"], :𝝭)

    set∇𝝭!(elements["Ωˢ"])
    # set∇²𝝭!(elements["Ωˢ"])
    set𝝭!(elements["∂Ωˢ"])

    return elements, nodes, sp, type, Ω, nodes_c
end

function import_MF_Gauss(filename1::String,n)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    

    elements["ΩC"] = getElements(nodes,entities["Ω"])
    push!(elements["ΩC"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    set∇𝝭!(elements["ΩC"])
    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    
    x = nodes.x
    y = nodes.y
    z = nodes.z
    Ω = getElements(nodes, entities["Ω"])
    
    s = 2.5
    s₁ = s*44.0/n*ones(length(nodes))
    s₂ = s*44.0/n*ones(length(nodes))
    push!(nodes,:s₁=>s₁,:s₂=>s₂,:s₃=>s₂)
  
  
    
    integrationOrder_Ω = 7
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 7

    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    elements["Ω"] = getElements(nodes,entities["Ω"], type,  integrationOrder_Ω, sp)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γᵍ"] = getElements(nodes,entities["Γᵍ"], type, integrationOrder_Γ, sp, normal = true)
    
    # elements["Ω"] = getElements(nodes,entities["Ω"],   integrationOrder_Ω)
    # elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    # elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"],  integrationOrder_Γ, normal = true)
    # elements["Γᵍ"] = getElements(nodes,entities["Γᵍ"],  integrationOrder_Γ, normal = true)
    
    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ω"],:𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍ"],:𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Γᵗ"],:𝝭)
    push!(elements["Γᵍ"],:𝝭, :∂𝝭∂x, :∂𝝭∂y)

    push!(elements["Ω"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γᵗ"],:𝗠=>𝗠)
    push!(elements["Γᵍ"],:𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    set∇𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γᵍ"])



    return elements, nodes, sp, type
end
function import_HR_GLS_MPP(filename1::String,filename2::String)
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

    sᵤ = 2.5*s*ones(length(nodes)) 
    sₚ = 2.5*s*ones(length(nodes))
    push!(nodes,:s₁=>sᵤ,:s₂=>sᵤ,:s₃=>sᵤ)
    push!(nodes_p,:s₁=>sₚ,:s₂=>sₚ,:s₃=>sₚ)

    integrationOrder_Ω = 4
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 4
    integrationOrder_R = 4
    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    elements["Ωᵘ"] = getElements(nodes,entities["Ω"], type,  integrationOrder_Ω, sp)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"], type,   integrationOrder_Γ, sp, normal = true)
    elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"], type, integrationOrder_Γ, sp, normal = true)
    # elements["Γʳ"] = getElements(nodes,entities["Γʳ"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γᵍᵘ"] = getElements(nodes,entities["Γᵍ"], type, integrationOrder_Γ, sp, normal = true)
    
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    elements["Ωᵘᵖ"] = getElements(nodes_p,entities["Ω"], type,  integrationOrder_Ω, sp)
    elements["Ωᵍᵘᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["∂Ωᵘᵖ"] = getElements(nodes_p, entities["Γ"], type,   integrationOrder_Γ, sp, normal = true)
    elements["Γᵗᵖ"] = getElements(nodes_p,entities["Γᵗ"], type, integrationOrder_Γ, sp, normal = true)
    # elements["Γʳ"] = getElements(nodes,entities["Γʳ"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γᵍᵘᵖ"] = getElements(nodes_p,entities["Γᵍ"], type, integrationOrder_Γ, sp, normal = true)
    




    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["∂Ωᵘ"],:𝝭)
    push!(elements["Γᵗ"],:𝝭)
    # push!(elements["Γʳ"],:𝝭)
    push!(elements["Γᵍᵘ"],:𝝭)

    push!(elements["Ωᵘ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵘ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["∂Ωᵘ"],:𝗠=>𝗠)
    push!(elements["Γᵗ"],:𝗠=>𝗠)
    # push!(elements["Γʳ"],:𝗠=>𝗠)
    push!(elements["Γᵍᵘ"],:𝗠=>𝗠)

    set∇𝝭!(elements["Ωᵘ"])
    set𝝭!(elements["∂Ωᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])
    set𝝭!(elements["Γᵗ"])
    # set𝝭!(elements["Γʳ"])
    set𝝭!(elements["Γᵍᵘ"])

    push!(elements["Ωᵘᵖ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍᵘᵖ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["∂Ωᵘᵖ"],:𝝭)
    push!(elements["Γᵗᵖ"],:𝝭)
    # push!(elements["Γʳ"],:𝝭)
    push!(elements["Γᵍᵘᵖ"],:𝝭)

    push!(elements["Ωᵘᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵘᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["∂Ωᵘᵖ"],:𝗠=>𝗠)
    push!(elements["Γᵗᵖ"],:𝗠=>𝗠)
    # push!(elements["Γʳ"],:𝗠=>𝗠)
    push!(elements["Γᵍᵘᵖ"],:𝗠=>𝗠)

    set∇𝝭!(elements["Ωᵘᵖ"])
    set𝝭!(elements["∂Ωᵘᵖ"])
    set∇𝝭!(elements["Ωᵍᵘᵖ"])
    set𝝭!(elements["Γᵗᵖ"])
    # set𝝭!(elements["Γʳ"])
    set𝝭!(elements["Γᵍᵘᵖ"])



    gmsh.open(filename2)
    
    # types = PiecewisePolynomial{:Constant}
    types = PiecewisePolynomial{:Linear2D}
    # types = PiecewisePolynomial{:Quadratic2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], types, integrationOrder_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], types, integrationOrder_Γ)
    elements["Γᵍˢ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωˢ"])
    
    
    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    # # sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    # elements["Ωᵖ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ω, sp)
    # elements["∂Ωᵖ"] = getElements(nodes, entities["Γ"], type, integrationOrder_Γ, sp)
    # elements["Ωᵍᵖ"] = getElements(nodes, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    # elements["Γᵍᵖ"] = getElements(nodes, entities["Γᵍ"],type,  integrationOrder_Γ, sp, normal = true)
    # typep = PiecewisePolynomial{:Constant}
    typep = PiecewisePolynomial{:Linear2D}
    # typep = PiecewisePolynomial{:Quadratic2D}
    elements["Ωᵖ"] = getPiecewiseElements(entities["Ω"], typep, integrationOrder_Ω)
    elements["Ωᵍᵖ"] = getPiecewiseElements(entities["Ω"], typep, integrationOrder_Ω)
    elements["∂Ωᵖ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], typep, integrationOrder_Γ)
    elements["Γᵍᵖ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωᵖ"])
    
    typep = PiecewisePolynomial{:Linear2D}
    # typep = PiecewisePolynomial{:Constant}

    elements["Ωᵖᵖ"] = getPiecewiseElements(entities["Ω"], typep, integrationOrder_Ω)
    elements["Ωᵍᵖᵖ"] = getPiecewiseElements(entities["Ω"], typep, integrationOrder_Ω)
    elements["∂Ωᵖᵖ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], typep, integrationOrder_Γ)
    elements["Γᵍᵖᵖ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωᵖ"])
   
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)



    
    # nₘ = 6
    # 𝗠 = zeros(nₘ)
    # ∂𝗠∂x = zeros(nₘ)
    # ∂𝗠∂y = zeros(nₘ)
    # push!(elements["Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    # push!(elements["∂Ωᵖ"], :𝝭)
    # push!(elements["Γᵍᵖ"], :𝝭)
    # push!(elements["Ωᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    # push!(elements["∂Ωᵖ"], :𝗠=>𝗠)
    # push!(elements["Γᵍᵖ"], :𝗠=>𝗠)
    # push!(elements["Ωᵍᵖ"], :𝝭)
    # push!(elements["Ωᵍᵖ"], :𝗠=>𝗠)

    push!(elements["Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωᵖ"], :𝝭)

    
    push!(elements["Ωᵖᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍᵖᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωᵖᵖ"], :𝝭)
    # push!(elements["Γᵍᵖ"], :𝝭)

    set∇𝝭!(elements["Ωˢ"])
    set𝝭!(elements["∂Ωˢ"])


    set∇𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["∂Ωᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵍᵖ"])

    
    set∇𝝭!(elements["Ωᵖᵖ"])
    set𝝭!(elements["∂Ωᵖᵖ"])
    set𝝭!(elements["Ωᵍᵖᵖ"])
    set𝝭!(elements["Γᵍᵖᵖ"])
    
    typeb = PiecewiseParametric{:Bubble,:Tri3}
    elements["Ωᵇ"] = getPiecewiseElements(entities["Ω"],typeb,integrationOrder_Ω)
    push!(elements["Ωᵇ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    set∇𝝭!(elements["Ωᵇ"])

    # gmsh.finalize()

    return elements, nodes, sp, type, Ω
end

function import_HR_reduced(filename1::String,filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    Ω = getElements(nodes, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 2.5*s*ones(length(nodes))
    push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
    
    gmsh.open(filename2)
    integrationOrder_Ω = 3
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 3
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    elements["Ωᵘ"] = getElements(nodes,entities["Ω"], type,  integrationOrder_Ω, sp)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"], type,   integrationOrder_Γ, sp, normal = true)
    elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γʳ"] = getElements(nodes,entities["Γʳ"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γᵍᵘ"] = getElements(nodes,entities["Γᵍ"], type, integrationOrder_Γ, sp, normal = true)
    
    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["∂Ωᵘ"],:𝝭)
    push!(elements["Γᵗ"],:𝝭)
    push!(elements["Γʳ"],:𝝭)
    push!(elements["Γᵍᵘ"],:𝝭)

    push!(elements["Ωᵘ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵘ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["∂Ωᵘ"],:𝗠=>𝗠)
    push!(elements["Γᵗ"],:𝗠=>𝗠)
    push!(elements["Γʳ"],:𝗠=>𝗠)
    push!(elements["Γᵍᵘ"],:𝗠=>𝗠)

    set∇𝝭!(elements["Ωᵘ"])
    set𝝭!(elements["∂Ωᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γʳ"])
    set𝝭!(elements["Γᵍᵘ"])



    
    
    types = PiecewisePolynomial{:Constant}
    # types = PiecewisePolynomial{:Linear2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], types, integrationOrder_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], types, integrationOrder_Γ)
    elements["Γᵍˢ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωˢ"])
    
    
    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    # # sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    # elements["Ωᵖ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ω, sp)
    # elements["∂Ωᵖ"] = getElements(nodes, entities["Γ"], type, integrationOrder_Γ, sp)
    # elements["Ωᵍᵖ"] = getElements(nodes, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    # elements["Γᵍᵖ"] = getElements(nodes, entities["Γᵍ"],type,  integrationOrder_Γ, sp, normal = true)
    # typep = PiecewisePolynomial{:Constant}
    # typep = PiecewisePolynomial{:Linear2D}
    typep = PiecewisePolynomial{:Quadratic2D}
    elements["Ωᵖ"] = getPiecewiseElements(entities["Ω"], typep, integrationOrder_Ω)
    elements["Ωᵍᵖ"] = getPiecewiseElements(entities["Ω"], typep, integrationOrder_Ω)
    elements["∂Ωᵖ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], typep, integrationOrder_Γ)
    elements["Γᵍᵖ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωᵖ"])
   
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)



    
    # nₘ = 6
    # 𝗠 = zeros(nₘ)
    # ∂𝗠∂x = zeros(nₘ)
    # ∂𝗠∂y = zeros(nₘ)
    # push!(elements["Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    # push!(elements["∂Ωᵖ"], :𝝭)
    # push!(elements["Γᵍᵖ"], :𝝭)
    # push!(elements["Ωᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    # push!(elements["∂Ωᵖ"], :𝗠=>𝗠)
    # push!(elements["Γᵍᵖ"], :𝗠=>𝗠)
    # push!(elements["Ωᵍᵖ"], :𝝭)
    # push!(elements["Ωᵍᵖ"], :𝗠=>𝗠)

    push!(elements["Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωᵖ"], :𝝭)
    # push!(elements["Γᵍᵖ"], :𝝭)

    set∇𝝭!(elements["Ωˢ"])
    set𝝭!(elements["∂Ωˢ"])


    set∇𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["∂Ωᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵍᵖ"])
    
    typeb = PiecewiseParametric{:Bubble,:Tri3}
    elements["Ωᵇ"] = getPiecewiseElements(entities["Ω"],typeb,integrationOrder_Ω)
    push!(elements["Ωᵇ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    set∇𝝭!(elements["Ωᵇ"])

    # gmsh.finalize()

    return elements, nodes, sp, type
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
