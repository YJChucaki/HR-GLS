
VTK_mix_pressure = quote
 #number of an elemnt nodes
 nₑₙ=6
#  fo = open("./vtk/patchtest_tri3_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
#  fo = open("./vtk/patchtest_tri6_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
 fo = open("./vtk/cantilever_tri3_mix_pressure_"*string(ndiv)*"_"*string(nₚ)*".vtk","w")
#  fo = open("./vtk/square_tri3_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
#  fo = open("./vtk/cantilever_quad4_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
#  fo = open("./vtk/cantilever_quad8_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
#  fo = open("./vtk/cantilever_tri6_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
# fo = open("./vtk/cook_membrane_quad_mix_pressure_0.3_"*string(ndiv)*"_"*string(i)*".vtk","w")
# fo = open("./vtk/cook_membrane_tri3_mix_pressure_0.3_"*string(ndiv)*"_"*string(i)*".vtk","w")
 @printf fo "# vtk DataFile Version 2.0\n"
#  @printf fo "cantilever_quad4_mix\n"
#  @printf fo "cantilever_quad8_mix\n"
#  @printf fo "cantilever_tri3_mix\n"
#  @printf fo "cantilever_tri6_mix\n"
@printf fo "patchtest_tri3_mix\n"
# @printf fo "patchtest_tri6_mix\n"
 @printf fo "ASCII\n"
 @printf fo "DATASET POLYDATA\n"
 @printf fo "POINTS %i float\n" nₚ

 for p in nodes_p
    @printf fo "%f %f %f\n" p.x p.y p.z
 end
 @printf fo "POLYGONS %i %i\n" nₑₚ (nₑₙ+1)*nₑₚ
 for ap in Ω
    𝓒 = ap.𝓒
    if nₑₙ==3
     @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    elseif nₑₙ==4
     @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    elseif nₑₙ==6
     @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    elseif nₑₙ==8
     @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    end
 end
 @printf fo "POINT_DATA %i\n" nₚ
 @printf fo "SCALARS P float 1\n"
 @printf fo "LOOKUP_TABLE default\n"
 for p in nodes_p
    @printf fo "%f\n" p.q
 end
# for (j,p) in enumerate(elements["Ωᵖ"])
#     ξ,  = p.𝓖
#     N = ξ[:𝝭]
#     p = 0.0
#     for (i,xᵢ) in enumerate(p.𝓒)
#         p += N[i]*xᵢ.q
#     end
#     @printf fo "%f\n" p
# end
 close(fo)
end
VTK_exact_pressure = quote
    #number of an elemnt nodes
    nₑₙ=3
   #  fo = open("./vtk/patchtest_tri3_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/patchtest_tri6_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
    fo = open("./vtk/"*test*"_exact_pressure.vtk","w")
   #  fo = open("./vtk/square_tri3_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/cantilever_quad4_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/cantilever_quad8_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/cantilever_tri6_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   # fo = open("./vtk/cook_membrane_quad_mix_pressure_0.3_"*string(ndiv)*"_"*string(i)*".vtk","w")
   # fo = open("./vtk/cook_membrane_tri3_mix_pressure_0.3_"*string(ndiv)*"_"*string(i)*".vtk","w")
    @printf fo "# vtk DataFile Version 2.0\n"
   #  @printf fo "cantilever_quad4_mix\n"
   #  @printf fo "cantilever_quad8_mix\n"
   #  @printf fo "cantilever_tri3_mix\n"
   #  @printf fo "cantilever_tri6_mix\n"
   @printf fo "patchtest_tri3_mix\n"
   # @printf fo "patchtest_tri6_mix\n"
    @printf fo "ASCII\n"
    @printf fo "DATASET POLYDATA\n"
    @printf fo "POINTS %i float\n" nc
   
    for p in nodes_c
       @printf fo "%f %f %f\n" p.x p.y p.z
    end
    @printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
    for ap in elements["Ω"]
       𝓒 = ap.𝓒
       if nₑₙ==3
        @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==4
        @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==6
        @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==8
        @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       end
    end
    @printf fo "POINT_DATA %i\n" nc
    @printf fo "SCALARS P float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for p in nodes_c
       @printf fo "%f\n" p.pc
    end
   # for (j,p) in enumerate(elements["Ωᵖ"])
   #     ξ,  = p.𝓖
   #     N = ξ[:𝝭]
   #     p = 0.0
   #     for (i,xᵢ) in enumerate(p.𝓒)
   #         p += N[i]*xᵢ.q
   #     end
   #     @printf fo "%f\n" p
   # end
    close(fo)
   end
VTK_mix_pressure_HR = quote
    #number of an elemnt nodes
    nₑₙ=3
    fo = open("./vtk/plate_with_hole_tri3_pressure_"*string(ndiv)*".vtk","w")
    #  fo = open("./vtk/patchtest_tri3_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/patchtest_tri6_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
    # fo = open("./vtk/cantilever_tri3_mix_pressure_"*string(ndiv)*"_"*string(nₚ)*".vtk","w")
   #  fo = open("./vtk/square_tri3_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/cantilever_quad4_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/cantilever_quad8_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/cantilever_tri6_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   # fo = open("./vtk/cook_membrane_quad_mix_pressure_0.3_"*string(ndiv)*"_"*string(i)*".vtk","w")
   # fo = open("./vtk/cook_membrane_tri3_mix_pressure_0.3_"*string(ndiv)*"_"*string(i)*".vtk","w")
    @printf fo "# vtk DataFile Version 2.0\n"
   #  @printf fo "cantilever_quad4_mix\n"
   #  @printf fo "cantilever_quad8_mix\n"
   #  @printf fo "cantilever_tri3_mix\n"
   #  @printf fo "cantilever_tri6_mix\n"
#    @printf fo "patchtest_tri3_mix\n"
   # @printf fo "patchtest_tri6_mix\n"
   @printf fo "plate_with_hole_mix\n"

    @printf fo "ASCII\n"
    @printf fo "DATASET POLYDATA\n"
    @printf fo "POINTS %i float\n" nₚ
   
    for p in nodes_p
       @printf fo "%f %f %f\n" p.x p.y p.z
    end
    @printf fo "POLYGONS %i %i\n" nₑₚ (nₑₙ+1)*nₑₚ
    for ap in Ω
       𝓒 = ap.𝓒
       if nₑₙ==3
        @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==4
        @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==6
        @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==8
        @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       end
    end
    @printf fo "POINT_DATA %i\n" nₚ
    @printf fo "SCALARS P float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for p in nodes_p
       @printf fo "%f\n" p.p 
    end
   # for (j,p) in enumerate(elements["Ωᵖ"])
   #     ξ,  = p.𝓖
   #     N = ξ[:𝝭]
   #     p = 0.0
   #     for (i,xᵢ) in enumerate(p.𝓒)
   #         p += N[i]*xᵢ.q
   #     end
   #     @printf fo "%f\n" p
   # end
    close(fo)
   end
VTK_mix_pressure_u = quote
    #number of an elemnt nodes
    nₑₙ=3
    fo = open("./vtk/patchtest_tri3_mix_pressure_u_"*string(ndiv)*".vtk","w")
   #  fo = open("./vtk/patchtest_tri6_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
    # fo = open("./vtk/cantilever_tri3_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
    # fo = open("./vtk/cantilever_quad4_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/cantilever_quad8_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/cantilever_tri6_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/cook_membrane_tri3_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
    @printf fo "# vtk DataFile Version 2.0\n"
    # @printf fo "cantilever_quad4_mix\n"
   #  @printf fo "cantilever_quad8_mix\n"
    @printf fo "cantilever_tri3_mix\n"
   #  @printf fo "cantilever_tri6_mix\n"
#    @printf fo "patchtest_tri3_mix\n"
   # @printf fo "patchtest_tri6_mix\n"
    @printf fo "ASCII\n"
    @printf fo "DATASET POLYDATA\n"
    @printf fo "POINTS %i float\n" nᵤ
    
  for p in nodes
    @printf fo "%f %f %f\n" p.x p.y p.z
  end
    @printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
    for ap in elements["Ωᵘ"]
       𝓒 = ap.𝓒

    #    for ap in Ω
    #     𝓒 = ap.𝓒
       if nₑₙ==3
        @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==4
        @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==6
        @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==8
        @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       end
    end
    # @printf fo "POINT_DATA %i\n" nₑ
    @printf fo "CELL_DATA %i\n" nₑ
    @printf fo "SCALARS PRESSURE float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
       for ap in elements["Ωᵘ"]
           𝓒 = ap.𝓒
           𝓖 = ap.𝓖
           ε₁₁ = 0.0
           ε₂₂ = 0.0
           ε₁₂ = 0.0
           for (i,ξ) in enumerate(𝓖)
                   B₁ = ξ[:∂𝝭∂x]
                   B₂ = ξ[:∂𝝭∂y]
                   for (j,xⱼ) in enumerate(𝓒)
                       ε₁₁ += B₁[j]*xⱼ.d₁
                       ε₂₂ += B₂[j]*xⱼ.d₂
                       ε₁₂ += B₁[j]*xⱼ.d₂ + B₂[j]*xⱼ.d₁
 
                   end
                  
           end
          
           p=K*(ε₁₁+ε₂₂)
           @printf fo "%f\n" p 
       end
    close(fo)
   end
VTK_displacement = quote
  #number of an elemnt nodes
  nₑₙ=6
  fo = open("./vtk/cantilever_tri3_mix_displacement_"*string(ndiv)*".vtk","w")
  @printf fo "# vtk DataFile Version 2.0\n"
  @printf fo "cantilever_tri3_mix\n"
  @printf fo "ASCII\n"
  @printf fo "DATASET POLYDATA\n"
  @printf fo "POINTS %i float\n" nᵤ
    
  for p in nodes
    @printf fo "%f %f %f\n" p.x p.y p.z
  end
  @printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
  for ap in elements["Ω"]
        𝓒 = ap.𝓒
        if nₑₙ==3
         @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
        elseif nₑₙ==4
         @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
        elseif nₑₙ==6
         @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
        elseif nₑₙ==8
         @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
        end
    end
  @printf fo "POINT_DATA %i\n" nᵤ
  @printf fo "VECTORS U float\n"
  for p in nodes
     @printf fo "%f %f %f\n" p.d₁ p.d₂ 0.0
  end
  close(fo)
end
# VTK_displacement = quote
#   #number of an elemnt nodes
#   nₑₙ=6
#   fo = open("./vtk/cantilever_tri3_mix_displacement_"*string(ndiv)*".vtk","w")
#   @printf fo "# vtk DataFile Version 2.0\n"
#   @printf fo "cantilever_tri3_mix\n"
#   @printf fo "ASCII\n"
#   @printf fo "DATASET POLYDATA\n"
#   @printf fo "POINTS %i float\n" nᵤ
    
#   for p in nodes
#     @printf fo "%f %f %f\n" p.x p.y p.z
#   end
#   @printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
#   for ap in elements["Ω"]
#         𝓒 = ap.𝓒
#         if nₑₙ==3
#          @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
#         elseif nₑₙ==4
#          @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
#         elseif nₑₙ==6
#          @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
#         elseif nₑₙ==8
#          @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
#         end
#     end
#   @printf fo "POINT_DATA %i\n" nᵤ
#   @printf fo "VECTORS U float\n"
#   for p in nodes
#      @printf fo "%f %f %f\n" p.d₁ p.d₂ 0.0
#   end
#   close(fo)
# end


VTK_HR_displacement_pressure = quote
    #number of an elemnt nodes
    nₑₙ=3
    fo = open("./vtk/"*test*"_tri3_HR_GLS_"*string(ndiv)*".vtk","w")
    @printf fo "# vtk DataFile Version 2.0\n"
    @printf fo "cantilever_tri3_mix\n"
    @printf fo "ASCII\n"
    @printf fo "DATASET POLYDATA\n"
    @printf fo "POINTS %i float\n" nc
      
    for p in nodes_c
      @printf fo "%f %f %f\n" p.x p.y p.z
    end
    @printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
    for ap in elements["Ω"]
          𝓒 = ap.𝓒
          if nₑₙ==3
           @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
          elseif nₑₙ==4
           @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
          elseif nₑₙ==6
           @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
          elseif nₑₙ==8
           @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
          end
      end

    @printf fo "POINT_DATA %i\n" nc
    @printf fo "VECTORS U float\n"
    for p in nodes_c
        x = p.x
        y = p.y
        indices = sp(x,y,0.0)
        ni = length(indices)
        𝓒 = [nodes[i] for i in indices]
        data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:∂𝝭∂x=>(4,zeros(ni)),:∂𝝭∂y=>(4,zeros(ni)),:𝗠=>(0,zeros(60)),:∂𝗠∂x=>(0,zeros(60)),:∂𝗠∂y=>(0,zeros(60))])
        ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
        𝓖 = [ξ]
        a = type(𝓒,𝓖)
        set∇𝝭!(a)
        d₁ = 0.0
        d₂ = 0.0
        N = ξ[:𝝭]
        u₁ = 0.0
        u₂ = 0.0
            for (j,xⱼ) in enumerate(𝓒)
                u₁ += N[j]*xⱼ.d₁
                u₂ += N[j]*xⱼ.d₂
            end    
       @printf fo "%f %f %f\n" u₁ u₂ 0.0
    end


    # @printf fo "CELL_DATA %i\n" nₑ
    # if ni==3
    #     @printf fo "SCALARS P float 3\n"
    # elseif ni==1
    #     @printf fo "SCALARS P float 1\n"
    #     elseif ni==6
    #     @printf fo "SCALARS P float 6\n"
    # end
    # @printf fo "LOOKUP_TABLE default\n"
    # for (i,elm) in enumerate(elements["Ωˢ"])
    #     𝓒ₚ = elm.𝓒
    #     𝓖 = elm.𝓖
    #     𝓒 = elements["Ω"][i].𝓒
    #     n = length(𝓖)
    #     q = zeros(n)
    #     i=0
    #     # for ξ in 𝓖
    #     for (index,ξ) in enumerate(𝓖)
    #         g = index[1]
    #         N = ξ[:𝝭]
    #         p = 0.0
    #         σ₁₁ = 0.0
    #         σ₂₂ = 0.0
    #         σ₁₂ = 0.0 
    #     #    for (j,xⱼ) in enumerate(𝓒)
            
    #            σ₁₁ += N[j]*xⱼ.dₛ₁₁
    #            σ₂₂ += N[j]*xⱼ.dₛ₂₂
    #            σ₁₂ += N[j]*xⱼ.dₛ₁₂

    #     #    end
    #        σ₃₃ = ν*(σ₁₁ + σ₂₂)
    #        p = (σ₁₁ + σ₂₂ + σ₃₃)/3
    #        q[g] = p
    #     end
    #     if ni==3
    #         @printf fo "%f %f %f\n"  q[1] q[2] q[3]
    #     elseif ni==1
    #         @printf fo "%f\n"  q[1] 
    #     elseif ni==6
    #         @printf fo "%f %f %f %f %f %f\n"  q[1] q[2] q[3] q[4] q[4] q[6]
    #     end
    # end
        
    # @printf fo "CELL_DATA %i\n" nₑ
    # @printf fo "SCALARS P float 3\n"
    # @printf fo "LOOKUP_TABLE default\n"
     
    # for (i,elm) in enumerate(elements["Ωˢ"])
    #     𝓒ₚ = elm.𝓒
    #     𝓖 = elm.𝓖
    #     𝓒 = elements["Ω"][i].𝓒

            
    #     σ₁₁1 = 𝓒ₚ[1].dₛ₁₁+𝓒ₚ[2].dₛ₁₁*𝓒[1].x+𝓒ₚ[3].dₛ₁₁*𝓒[1].y
    #     σ₂₂1 = 𝓒ₚ[1].dₛ₂₂+𝓒ₚ[2].dₛ₂₂*𝓒[1].x+𝓒ₚ[3].dₛ₂₂*𝓒[1].y
    #     σ₃₃1 = ν*(σ₁₁1 + σ₂₂1)
    #     p1 = (σ₁₁1 + σ₂₂1 + σ₃₃1)/3

    #     σ₁₁2 = 𝓒ₚ[1].dₛ₁₁+𝓒ₚ[2].dₛ₁₁*𝓒[2].x+𝓒ₚ[3].dₛ₁₁*𝓒[2].y
    #     σ₂₂2 = 𝓒ₚ[1].dₛ₂₂+𝓒ₚ[2].dₛ₂₂*𝓒[2].x+𝓒ₚ[3].dₛ₂₂*𝓒[2].y
    #     σ₃₃2 = ν*(σ₁₁2 + σ₂₂2)
    #     p2 = (σ₁₁2 + σ₂₂2 + σ₃₃2)/3

    #     σ₁₁3 = 𝓒ₚ[1].dₛ₁₁+𝓒ₚ[2].dₛ₁₁*𝓒[3].x+𝓒ₚ[3].dₛ₁₁*𝓒[3].y
    #     σ₂₂3 = 𝓒ₚ[1].dₛ₂₂+𝓒ₚ[3].dₛ₂₂*𝓒[3].x+𝓒ₚ[3].dₛ₂₂*𝓒[3].y
    #     σ₃₃3 = ν*(σ₁₁3 + σ₂₂3)
    #     p3 = (σ₁₁3 + σ₂₂3 + σ₃₃3)/3
        
        
        
    #         @printf fo "%f %f %f\n"  p1 p2 p3
        
    # end

        
    @printf fo "CELL_DATA %i\n" nₑ
    @printf fo "SCALARS P float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for (i,elm) in enumerate(elements["Ωˢ"])
    𝓒ₚ = elm.𝓒
    𝓖 = elm.𝓖
    𝓒 = elements["Ω"][i].𝓒
    x = (𝓒[1].x+𝓒[2].x+𝓒[3].x)/3
    y = (𝓒[1].y+𝓒[2].y+𝓒[3].y)/3
    if nₛ==3
    σ₁₁ = 𝓒ₚ[1].dₛ₁₁+𝓒ₚ[2].dₛ₁₁*x+𝓒ₚ[3].dₛ₁₁*y
    σ₂₂ = 𝓒ₚ[1].dₛ₂₂+𝓒ₚ[2].dₛ₂₂*x+𝓒ₚ[3].dₛ₂₂*y
    elseif nₛ==6
        σ₁₁ = 𝓒ₚ[1].dₛ₁₁+𝓒ₚ[2].dₛ₁₁*x+𝓒ₚ[3].dₛ₁₁*y+𝓒ₚ[4].dₛ₁₁*x^2+𝓒ₚ[6].dₛ₁₁*y^2+𝓒ₚ[5].dₛ₁₁*x*y
        σ₂₂ = 𝓒ₚ[1].dₛ₂₂+𝓒ₚ[2].dₛ₂₂*x+𝓒ₚ[3].dₛ₂₂*y+𝓒ₚ[4].dₛ₂₂*x^2+𝓒ₚ[6].dₛ₂₂*y^2+𝓒ₚ[5].dₛ₂₂*x*y
    end
    σ₃₃ = ν*(σ₁₁ + σ₂₂)
    p = (σ₁₁ + σ₂₂ + σ₃₃)/3
    @printf fo "%f\n"  p
    end
    

    
    close(fo)
  end

  VTK_HR_displacement_pressure_smoothing = quote
    #number of an elemnt nodes
    nₑₙ=3
    fo = open("./vtk/"*test*"_tri3_HR_GLS_sm_"*string(ndiv)*".vtk","w")
    @printf fo "# vtk DataFile Version 2.0\n"
    @printf fo "cantilever_tri3_mix\n"
    @printf fo "ASCII\n"
    @printf fo "DATASET POLYDATA\n"
    @printf fo "POINTS %i float\n" nc
      
    for p in nodes_c
      @printf fo "%f %f %f\n" p.x p.y p.z
    end
    @printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
    for ap in elements["Ω"]
          𝓒 = ap.𝓒
          if nₑₙ==3
           @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
          elseif nₑₙ==4
           @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
          elseif nₑₙ==6
           @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
          elseif nₑₙ==8
           @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
          end
      end

    # @printf fo "POINT_DATA %i\n" nᵤ
    # @printf fo "VECTORS U float\n"
    # for p in nodes
    #     x = p.x
    #     y = p.y
    #     indices = sp(x,y,0.0)
    #     ni = length(indices)
    #     𝓒 = [nodes[i] for i in indices]
    #     data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:∂𝝭∂x=>(4,zeros(ni)),:∂𝝭∂y=>(4,zeros(ni)),:𝗠=>(0,zeros(21)),:∂𝗠∂x=>(0,zeros(21)),:∂𝗠∂y=>(0,zeros(21))])
    #     ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
    #     𝓖 = [ξ]
    #     a = type(𝓒,𝓖)
    #     set∇𝝭!(a)
    #     d₁ = 0.0
    #     d₂ = 0.0
    #     N = ξ[:𝝭]
    #     u₁ = 0.0
    #     u₂ = 0.0
    #         for (j,xⱼ) in enumerate(𝓒)
    #             u₁ += N[j]*xⱼ.d₁
    #             u₂ += N[j]*xⱼ.d₂
    #         end    
    #    @printf fo "%f %f %f\n" u₁ u₂ 0.0
    # end


    @printf fo "POINT_DATA %i\n" nc
    @printf fo "SCALARS P float\n"
    @printf fo "LOOKUP_TABLE default\n"
    

    for (i,node) in enumerate(nodes_c)
        p_smoothing = p_node[i]/w[i]
        @printf fo "%f\n" p_smoothing
    end

    close(fo)
end

VTK_T6P3_pressure = quote
     #number of an elemnt nodes
     nₑₙ=3
     fo = open("./vtk/cantilever_T6P3_pressure_"*string(ndiv)*".vtk","w")
     @printf fo "# vtk DataFile Version 2.0\n"
     @printf fo "cantilever_tri3_mix\n"
     @printf fo "ASCII\n"
     @printf fo "DATASET POLYDATA\n"
     @printf fo "POINTS %i float\n" nₚ
     
     for p in nodes_p
         @printf fo "%f %f %f\n" p.x p.y p.z
     end
     @printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
     for ap in elements["Ωᵖ"]
         𝓒 = ap.𝓒
         if nₑₙ==3
          @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
         elseif nₑₙ==4
          @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
         elseif nₑₙ==6
          @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
         elseif nₑₙ==8
          @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
         end
     end
     @printf fo "POINT_DATA %i\n" nₚ
     @printf fo "SCALARS P float 1\n"
     @printf fo "LOOKUP_TABLE default\n"
     for p in nodes_p
         @printf fo "%f\n" p.q 
     end
     close(fo)    
 end    
VTK_HR_MPP = quote
#number of an elemnt nodes
nₑₙ=6
# fo = open("./vtk/cantilever_Q4P1_"*string(ndiv)*".vtk","w")
fo = open("./vtk/cook_membrane_Q4P1_0.4999_"*string(ndiv)*".vtk","w")
@printf fo "# vtk DataFile Version 2.0\n"
@printf fo "cantilever_Q4P1_mix\n"
@printf fo "ASCII\n"
@printf fo "DATASET POLYDATA\n"

@printf fo "POINTS %i float\n" nᵤ
for p in nodes
    @printf fo "%f %f %f\n" p.x p.y p.z
end
# elment type
@printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
for ap in Ω
    𝓒 = ap.𝓒
    if nₑₙ==3
     @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    elseif nₑₙ==4
     @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    elseif nₑₙ==6
     @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    elseif nₑₙ==8
     @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    end
end

@printf fo "POINT_DATA %i\n" nᵤ
@printf fo "VECTORS U float\n"
for p in nodes
 @printf fo "%f %f %f\n" p.d₁ p.d₂ 0.0
end

@printf fo "CELL_DATA %i\n" nₑ
@printf fo "SCALARS P float 1\n"
@printf fo "LOOKUP_TABLE default\n"
for p in 1:nₑ
    @printf fo "%f\n"  𝑝[p] 
end
# @printf fo "POINT_DATA %i\n" nₚ
#     @printf fo "SCALARS P float 1\n"
#     @printf fo "LOOKUP_TABLE default\n"
#     for p in nodes_p
#        @printf fo "%f\n" p.p 
#     end
close(fo)
end

VTK_Q4P1_displacement_pressure = quote
    #number of an elemnt nodes
    nₑₙ=4
    # fo = open("./vtk/cantilever_Q4P1_"*string(ndiv)*".vtk","w")
    fo = open("./vtk/cook_membrane_Q4P1_"*string(ndiv)*".vtk","w")
    @printf fo "# vtk DataFile Version 2.0\n"
    @printf fo "cantilever_Q4P1_mix\n"
    @printf fo "ASCII\n"
    @printf fo "DATASET POLYDATA\n"
    
    @printf fo "POINTS %i float\n" nᵤ
    for p in nodes
        @printf fo "%f %f %f\n" p.x p.y p.z
    end
    # elment type
    @printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
    for ap in elements["Ω"]
        𝓒 = ap.𝓒
        if nₑₙ==3
         @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
        elseif nₑₙ==4
         @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
        elseif nₑₙ==6
         @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
        elseif nₑₙ==8
         @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
        end
    end
    
    @printf fo "POINT_DATA %i\n" nᵤ
    @printf fo "VECTORS U float\n"
    for p in nodes
     @printf fo "%f %f %f\n" p.d₁ p.d₂ 0.0
    end
    
    @printf fo "CELL_DATA %i\n" nₑ
    @printf fo "SCALARS P float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for p in 1:nₑ
        @printf fo "%f\n"  q[p] 
    end
    
    close(fo)
    end