
using Gmsh, Statistics

function import_patchtest_fem(filename::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    elements["Ω"] = getElements(nodes,entities["Ω"])
    elements["Γ¹"] = getElements(nodes,entities["Γ¹"],normal=true)
    elements["Γ²"] = getElements(nodes,entities["Γ²"],normal=true)
    elements["Γ³"] = getElements(nodes,entities["Γ³"],normal=true)
    elements["Γ⁴"] = getElements(nodes,entities["Γ⁴"],normal=true)

    push!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Γ¹"],:𝝭)
    push!(elements["Γ²"],:𝝭)
    push!(elements["Γ³"],:𝝭)
    push!(elements["Γ⁴"],:𝝭)

    gmsh.finalize()

    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    set∇𝝭!(elements["Ω"])
    set𝝭!(elements["Γ"])

    return elements, nodes
end
function import_patchtest_MF(filename::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    nodes_c = get𝑿ᵢ()
    elements["ΩC"] = getElements(nodes_c,entities["Ω"])
    push!(elements["ΩC"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    set∇𝝭!(elements["ΩC"])


    x = nodes.x
    y = nodes.y
    z = nodes.z
    Ω = getElements(nodes, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    sᵤ = 2.5*s*ones(length(nodes))

    push!(nodes,:s₁=>sᵤ,:s₂=>sᵤ,:s₃=>sᵤ)

    integrationOrder_Ω = 7
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 7
    
    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)

    elements["Ωᵍ"] = getElements(nodes,entities["Ω"], type,  integrationOrder_Ω, sp)
    elements["Ω"] = getElements(nodes,entities["Ω"], type,  integrationOrder_Ω, sp)
    elements["Γ¹"] = getElements(nodes,entities["Γ¹"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ²"] = getElements(nodes,entities["Γ²"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ³"] = getElements(nodes,entities["Γ³"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ⁴"] = getElements(nodes,entities["Γ⁴"], type, integrationOrder_Γ, sp, normal = true)

    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Γ¹"],:𝝭)
    push!(elements["Γ²"],:𝝭)
    push!(elements["Γ³"],:𝝭)
    push!(elements["Γ⁴"],:𝝭)


    push!(elements["Ω"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γ¹"],:𝗠=>𝗠)
    push!(elements["Γ²"],:𝗠=>𝗠)
    push!(elements["Γ³"],:𝗠=>𝗠)
    push!(elements["Γ⁴"],:𝗠=>𝗠)

    gmsh.finalize()

    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    set∇𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    set𝝭!(elements["Γ"])

    return elements, nodes, sp, type, nodes_c
end

function import_patchtest_elasticity_penalty(filename1::String,filename2::String)
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

    integrationOrder_Ω = 6
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 6

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ωᵘ"] = getElements(nodes,entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"],   integrationOrder_Γ, normal = true)
    elements["Γ¹ᵘ"] = getElements(nodes,entities["Γ¹"], integrationOrder_Γ, normal = true)
    elements["Γ²ᵘ"] = getElements(nodes,entities["Γ²"], integrationOrder_Γ, normal = true)
    elements["Γ³ᵘ"] = getElements(nodes,entities["Γ³"], integrationOrder_Γ, normal = true)
    elements["Γ⁴ᵘ"] = getElements(nodes,entities["Γ⁴"], integrationOrder_Γ, normal = true)
    elements["Γᵘ"] = elements["Γ¹ᵘ"]∪elements["Γ²ᵘ"]∪elements["Γ³ᵘ"]∪elements["Γ⁴ᵘ"]

    push!(elements["Ωᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["∂Ωᵘ"],:𝝭)
    push!(elements["Γ¹ᵘ"],:𝝭)
    push!(elements["Γ²ᵘ"],:𝝭)
    push!(elements["Γ³ᵘ"],:𝝭)
    push!(elements["Γ⁴ᵘ"],:𝝭)

    set∇𝝭!(elements["Ωᵘ"])
    set𝝭!(elements["∂Ωᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])
    set𝝭!(elements["Γᵘ"])

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["∂Ωᵖ"] = getElements(nodes_p, entities["Γ"], type, integrationOrder_Γ, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γ¹ᵖ"] = getElements(nodes_p, entities["Γ¹"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ²ᵖ"] = getElements(nodes_p, entities["Γ²"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ³ᵖ"] = getElements(nodes_p, entities["Γ³"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ⁴ᵖ"] = getElements(nodes_p, entities["Γ⁴"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γᵖ"] = elements["Γ¹ᵖ"]∪elements["Γ²ᵖ"]∪elements["Γ³ᵖ"]∪elements["Γ⁴ᵖ"]

    nₘ = 6
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωᵖ"], :𝝭)
    push!(elements["Γ¹ᵖ"], :𝝭)
    push!(elements["Γ²ᵖ"], :𝝭)
    push!(elements["Γ³ᵖ"], :𝝭)
    push!(elements["Γ⁴ᵖ"], :𝝭)
    push!(elements["Ωᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["∂Ωᵖ"], :𝗠=>𝗠)
    push!(elements["Γ¹ᵖ"], :𝗠=>𝗠)
    push!(elements["Γ²ᵖ"], :𝗠=>𝗠)
    push!(elements["Γ³ᵖ"], :𝗠=>𝗠)
    push!(elements["Γ⁴ᵖ"], :𝗠=>𝗠)
    push!(elements["Ωᵍᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    set∇𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["∂Ωᵖ"])
    set∇𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵖ"])

    gmsh.finalize()

    return elements, nodes, nodes_p
end

function import_patchtest_elasticity_mix(filename1::String,filename2::String)
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
    sₚ = 1.0*s*ones(length(nodes))

    push!(nodes,:s₁=>sᵤ,:s₂=>sᵤ,:s₃=>sᵤ)
    push!(nodes_p,:s₁=>sₚ,:s₂=>sₚ,:s₃=>sₚ)
    integrationOrder_Ω = 6
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 6

   
    
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    
    elements["Ωᵘ"] = getElements(nodes,entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["∂Ωᵘ"] = getElements(nodes, entities["Γ"],   type, integrationOrder_Γ, sp, normal = true)
    elements["Γ¹ᵘ"] = getElements(nodes,entities["Γ¹"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ²ᵘ"] = getElements(nodes,entities["Γ²"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ³ᵘ"] = getElements(nodes,entities["Γ³"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ⁴ᵘ"] = getElements(nodes,entities["Γ⁴"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γᵘ"] = elements["Γ¹ᵘ"]∪elements["Γ²ᵘ"]∪elements["Γ³ᵘ"]∪elements["Γ⁴ᵘ"]

    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["∂Ωᵘ"],:𝝭)
    push!(elements["Γ¹ᵘ"],:𝝭)
    push!(elements["Γ²ᵘ"],:𝝭)
    push!(elements["Γ³ᵘ"],:𝝭)
    push!(elements["Γ⁴ᵘ"],:𝝭)

    push!(elements["Ωᵘ"],:𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵘ"],:𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["∂Ωᵘ"],:𝗠=>𝗠)
    push!(elements["Γ¹ᵘ"],:𝗠=>𝗠)
    push!(elements["Γ²ᵘ"],:𝗠=>𝗠)
    push!(elements["Γ³ᵘ"],:𝗠=>𝗠)
    push!(elements["Γ⁴ᵘ"],:𝗠=>𝗠)


    set∇𝝭!(elements["Ωᵘ"])
    set𝝭!(elements["∂Ωᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])
    set𝝭!(elements["Γᵘ"])

    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    # elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    # elements["∂Ωᵖ"] = getElements(nodes_p, entities["Γ"], type, integrationOrder_Γ, sp)
    # elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    # elements["Γ¹ᵖ"] = getElements(nodes_p, entities["Γ¹"],type,  integrationOrder_Γ, sp, normal = true)
    # elements["Γ²ᵖ"] = getElements(nodes_p, entities["Γ²"],type,  integrationOrder_Γ, sp, normal = true)
    # elements["Γ³ᵖ"] = getElements(nodes_p, entities["Γ³"],type,  integrationOrder_Γ, sp, normal = true)
    # elements["Γ⁴ᵖ"] = getElements(nodes_p, entities["Γ⁴"], type, integrationOrder_Γ, sp, normal = true)
    # elements["Γᵖ"] = elements["Γ¹ᵖ"]∪elements["Γ²ᵖ"]∪elements["Γ³ᵖ"]∪elements["Γ⁴ᵖ"]

    # elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"],  integrationOrder_Ω)
    # elements["∂Ωᵖ"] = getElements(nodes_p, entities["Γ"],  integrationOrder_Γ)
    # elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"],   integrationOrder_Ωᵍ)
    # elements["Γ¹ᵖ"] = getElements(nodes_p, entities["Γ¹"],  integrationOrder_Γ,  normal = true)
    # elements["Γ²ᵖ"] = getElements(nodes_p, entities["Γ²"],  integrationOrder_Γ,  normal = true)
    # elements["Γ³ᵖ"] = getElements(nodes_p, entities["Γ³"],  integrationOrder_Γ,  normal = true)
    # elements["Γ⁴ᵖ"] = getElements(nodes_p, entities["Γ⁴"],  integrationOrder_Γ,  normal = true)
    # elements["Γᵖ"] = elements["Γ¹ᵖ"]∪elements["Γ²ᵖ"]∪elements["Γ³ᵖ"]∪elements["Γ⁴ᵖ"]

    # type = PiecewisePolynomial{:Quadratic2D}
    typep = PiecewisePolynomial{:Linear2D}
    # typep = PiecewisePolynomial{:Constant}
    elements["Ωᵖ"] = getPiecewiseElements(entities["Ω"], typep, integrationOrder_Ω)
    elements["∂Ωᵖ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], typep, integrationOrder_Γ)
    elements["Ωᵍᵖ"] =   getPiecewiseElements(entities["Ω"], typep, integrationOrder_Ω)
    elements["Γ¹ᵖ"] = getElements(entities["Γ¹"],entities["Γ"], elements["∂Ωᵖ"])
    elements["Γ²ᵖ"] = getElements(entities["Γ²"],entities["Γ"], elements["∂Ωᵖ"])
    elements["Γ³ᵖ"] = getElements(entities["Γ³"],entities["Γ"], elements["∂Ωᵖ"])
    elements["Γ⁴ᵖ"] = getElements(entities["Γ⁴"],entities["Γ"], elements["∂Ωᵖ"])
    elements["Γᵖ"] = elements["Γ¹ᵖ"]∪elements["Γ²ᵖ"]∪elements["Γ³ᵖ"]∪elements["Γ⁴ᵖ"]
    

    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωᵖ"], :𝝭)
    # push!(elements["Γ¹ᵖ"], :𝝭)
    # push!(elements["Γ²ᵖ"], :𝝭)
    # push!(elements["Γ³ᵖ"], :𝝭)
    # push!(elements["Γ⁴ᵖ"], :𝝭)
    # push!(elements["Ωᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    # push!(elements["∂Ωᵖ"], :𝗠=>𝗠)
    # push!(elements["Γ¹ᵖ"], :𝗠=>𝗠)
    # push!(elements["Γ²ᵖ"], :𝗠=>𝗠)
    # push!(elements["Γ³ᵖ"], :𝗠=>𝗠)
    # push!(elements["Γ⁴ᵖ"], :𝗠=>𝗠)
    # push!(elements["Ωᵍᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    # push!(elements["Ωᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    set∇𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["∂Ωᵖ"])
    # set∇𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵖ"])

    # type = PiecewisePolynomial{:Constant}
    type = PiecewisePolynomial{:Linear2D}
    # type = PiecewisePolynomial{:Quadratic2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], type, integrationOrder_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], type, integrationOrder_Γ)
    elements["Γ¹ˢ"] = getElements(entities["Γ¹"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γ²ˢ"] = getElements(entities["Γ²"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γ³ˢ"] = getElements(entities["Γ³"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γ⁴ˢ"] = getElements(entities["Γ⁴"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γˢ"] = elements["Γ¹ˢ"]∪elements["Γ²ˢ"]∪elements["Γ³ˢ"]∪elements["Γ⁴ˢ"]
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)

    set∇𝝭!(elements["Ωˢ"])
    set𝝭!(elements["∂Ωˢ"])

    type = PiecewiseParametric{:Bubble,:Tri3}
    elements["Ωᵇ"] = getPiecewiseElements(entities["Ω"],type,integrationOrder_Ω)
    push!(elements["Ωᵇ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    set∇𝝭!(elements["Ωᵇ"])

    gmsh.finalize()

    return elements, nodes, nodes_p
end
function import_patchtest_mix(filename1::String,filename2::String)
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

    integration_Ω = 7
    integrationOrder_Ωᵍ = 8
    integration_Γ = 7

    gmsh.open(filename2)
    entities = getPhysicalGroups()

    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    elements["Ω"] = getElements(nodes, entities["Ω"], type, integration_Ω, sp)
    elements["∂Ω"] = getElements(nodes, entities["Γ"], type, integration_Γ, sp, normal = true)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"],type, integration_Γ, sp, normal = true)
    elements["Γ²"] = getElements(nodes, entities["Γ²"],type, integration_Γ, sp, normal = true)
    elements["Γ³"] = getElements(nodes, entities["Γ³"],type, integration_Γ, sp, normal = true)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], type, integration_Γ, sp, normal = true)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]


    nₘ = 60
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    ∂²𝗠∂x² = zeros(nₘ)
    ∂²𝗠∂y² = zeros(nₘ)
    ∂²𝗠∂x∂y = zeros(nₘ)
    # for elm in elements["Ω"]
    #     𝓒ₑ = elm.𝓒
    #     nc = length(𝓒ₑ)
    #     xₑ = zeros(nc)
    #     yₑ = zeros(nc)
    #     xₑ = 𝓒ₑ.x
    #     yₑ = 𝓒ₑ.y
    #     𝝭 = zeros(nc)
    #     for i in 1:nc
    #         x = xₑ[i]
    #         y = yₑ[i]
    #         data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(nc)),:𝗠=>(0,𝗠)])
    #         ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
    #         𝓖ₑ = [ξ]
    #         a = type(𝓒ₑ,𝓖ₑ)
    #         set𝝭!(a)
    #         N = ξ[:𝝭]
    #         𝝭[i] = N[i]
    #     end
    #     push!(𝓒ₑ,:𝝭=>𝝭)
    # end


    push!(elements["Ω"], :𝝭, :∂𝝭∂x, :∂𝝭∂y, :∂²𝝭∂x², :∂²𝝭∂y², :∂²𝝭∂x∂y)
    push!(elements["∂Ω"], :𝝭)
    push!(elements["Γ¹"], :𝝭)
    push!(elements["Γ²"], :𝝭)
    push!(elements["Γ³"], :𝝭)
    push!(elements["Γ⁴"], :𝝭)
    push!(elements["Ω"],  :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂²𝗠∂x²=>∂²𝗠∂x², :∂²𝗠∂y²=>∂²𝗠∂y², :∂²𝗠∂x∂y=>∂²𝗠∂x∂y)
    push!(elements["∂Ω"], :𝗠=>𝗠)
    push!(elements["Γ¹"], :𝗠=>𝗠)
    push!(elements["Γ²"], :𝗠=>𝗠)
    push!(elements["Γ³"], :𝗠=>𝗠)
    push!(elements["Γ⁴"], :𝗠=>𝗠)
    push!(elements["Ωᵍ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y, :∂²𝝭∂x², :∂²𝝭∂y², :∂²𝝭∂x∂y)
    push!(elements["Ωᵍ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂²𝗠∂x²=>∂²𝗠∂x², :∂²𝗠∂y²=>∂²𝗠∂y², :∂²𝗠∂x∂y=>∂²𝗠∂x∂y)


    # type = PiecewisePolynomial{:Constant}
    type = PiecewisePolynomial{:Linear2D}
    # type = PiecewisePolynomial{:Quadratic2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], type, integration_Ω)
    
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], type, integration_Γ)
    elements["∂Ωˢˢ"] = getPiecewiseElements(entities["Γ"], type, integration_Γ)
    elements["Γ¹ˢ"] = getElements(entities["Γ¹"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γ²ˢ"] = getElements(entities["Γ²"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γ³ˢ"] = getElements(entities["Γ³"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γ⁴ˢ"] = getElements(entities["Γ⁴"],entities["Γ"], elements["∂Ωˢ"])

    # elements["∂Ωˢ"] = getPiecewiseElements(entities["Γ"], type, integration_Γ)
    # elements["Γ¹ˢ"] = getPiecewiseElements(entities["Γ¹"], type, integration_Γ)
    # elements["Γ²ˢ"] = getPiecewiseElements(entities["Γ²"], type, integration_Γ)
    # elements["Γ³ˢ"] = getPiecewiseElements(entities["Γ³"], type, integration_Γ)
    # elements["Γ⁴ˢ"] = getPiecewiseElements(entities["Γ⁴"], type, integration_Γ)
    elements["Γˢ"] = elements["Γ¹ˢ"]∪elements["Γ²ˢ"]∪elements["Γ³ˢ"]∪elements["Γ⁴ˢ"]
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)
    push!(elements["∂Ωˢˢ"], :𝝭)
    # gmsh.finalize()

    return elements, nodes
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
