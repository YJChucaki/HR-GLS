using CairoMakie

E = 10.92e6
ν = 0.3
L = 1.0
h = 0.001
w(x,y) = 1/3*x^3*(x-1)^3*y^3*(y-1)^3-2*h^2/(5*(1-ν))*(y^3*(y-1)^3*x*(x-1)*(5*x^2-5*x+1)+x^3*(x-1)^3*y*(y-1)*(5*y^2-5*y+1))
θ₁(x,y) = y^3*(y-1)^3*x^2*(x-1)^2*(2*x-1)
θ₂(x,y) = x^3*(x-1)^3*y^2*(y-1)^2*(2*y-1)
w₁(x,y) = (x-1)^2*x^2*(2*x-1)*(y-1)^3*y^3-2*h^2/(5*(1-ν))*((20*x^3-30*x^2+12*x-1)*(y-1)^3*y^3+3*(x-1)^2*x^2*(2*x-1)*(y-1)*y*(5*y^2-5*y+1))
w₂(x,y) = (x-1)^3*x^3*(y-1)^2*y^2*(2*y-1)-2*h^2/(5*(1-ν))*(3*(x-1)*x*(5*x^2-5*x+1)*(y-1)^2*y^2*(2*y-1)+x^3*(x-1)^3*(20*y^3-30*y^2+12*y-1))
Dˢ = 5/6*E*h/(2*(1+ν))
Q₁(x,y) = Dˢ*(w₁(x,y)-θ₁(x,y))
Q₂(x,y) = Dˢ*(w₂(x,y)-θ₂(x,y))

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
# zs = [Q₁(x,y) for x in xs, y in ys]
zs = [Q₂(x,y) for x in xs, y in ys]
surface!(xs,ys,zeros(ind,ind),color=zs,colorrange=(-0.000025,0.000025),colormap=:lightrainbow)
contour!(xs[1:end-1],ys,zs[1:end-1,:],levels=-0.000025:0.00000715:0.000025,color=:azure)
Colorbar(fig[1,2], limits=(-0.000025,0.000025), colormap=:lightrainbow)
save("./png/plate_exactQ₂_solution.png",fig, px_per_unit = 10.0)
# save("./png/plate_exactQ₁_solution.png",fig, px_per_unit = 10.0)
fig