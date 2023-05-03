

# ∂ut = a₁ * ∇^2u + γ₀u - (γ₀/γ₁)u^2 - βuv
# ∂vt = a₂ * ∇^2v + μuv - ηv

# u: pred
# v: predator

# default values
# a₁ = 0.1 - diffusion prey
# a₂ = 0.1 - diffusion predator
# γ₀ = 0.4 - birth rate prey
# γ₁ = 1.5 - carrying capacity prey
# β = 0.5 - max uptake of prey by predator
# μ = 0.3 - birth rate of predator
# η = 0.1 - death rate predator


using BayesianDiscovery
const BD = BayesianDiscovery
using Plots
using Random, Distributions
using Kronecker

include("reaction_diffusion_generate.jl")

######################### Load Data #########################


U, V, x, y, t = rd(T = 10, Lx=11, Ly=11, Δx=0.1, Δy=0.1, Δt=0.01)

xinds = range(11, 212, step = 5)
yinds = range(11, 212, step = 5)
tinds = range(1, 1001, step = 10)

U = U[xinds, yinds, tinds]
V = V[xinds, yinds, tinds]
x = x[xinds]
y = y[yinds]
t = t[tinds]

# anim = @animate for i ∈ range(1, length(t), step=5)
#   p1 = contourf(x, y, U[:, :, i], c=:heat, clim=extrema(U))
#   p2 = contourf(x, y, V[:, :, i], c=:heat, clim=extrema(V))

#   plot(p1, p2, layout=grid(2, 1), size=(500, 700))
# end
# gif(anim, fps=5)

######################### visualize data #########################


# Nx = length(x)
# Ny = length(y)

# # l = @layout [grid(2, 2) a{0.07w}]

# l = @layout [grid(2, 1) a{0.07w} grid(2, 1) a{0.07w}]

# p11 = contourf(y, x, reshape(Z[:, 11, 1], Nx, Ny), clim=extrema(Z[:,1:51,1]), legend=false, ticks=false, framestyle=:box, c=:heat)
# p12 = contourf(y, x, reshape(Z[:, 11, 2], Nx, Ny), clim=extrema(Z[:,1:31,2]), legend=false, ticks=false, framestyle=:box, c=:heat)
# p21 = contourf(y, x, reshape(Z[:, 31, 1], Nx, Ny), clim=extrema(Z[:,1:51,1]), legend=false, ticks=false, framestyle=:box, c=:heat)
# p22 = contourf(y, x, reshape(Z[:, 31, 2], Nx, Ny), clim=extrema(Z[:,1:31,2]), legend=false, ticks=false, framestyle=:box, c=:heat)

# l1 = scatter([0, 0], [0, 0], zcolor=[0, 1], clims=extrema(Z[:,1:51,1]), xlims=(0.9, 1.1), ylims=(0, 1), label="", framestyle=:none, c=:heat, xtickfontsize=18, ytickfontsize=18)
# l2 = scatter([0, 0], [0, 0], zcolor=[0, 1], clims=extrema(Z[:,1:31,2]), xlims=(0.9, 1.1), ylims=(0, 1), label="", framestyle=:none, c=:heat, xtickfontsize=18, ytickfontsize=18)

# plot(p11, p21, l1, p12, p22, l2, layout=l, size=(1750, 900), right_margin=5Plots.mm)
# plot(p11, p12, l1, p21, p22, l2, layout=l, size=(1500, 750))



######################### set up data #########################

Random.seed!(1)
U_Z = U + 0.05 * std(U) * rand(Normal(), size(U))
V_Z = V + 0.05 * std(V) * rand(Normal(), size(V))

Y = cat(reshape(U, :, length(t)), reshape(V, :, length(t)), dims=3)
Z = cat(reshape(U_Z, :, length(t)), reshape(V_Z, :, length(t)), dims=3)

######################### Set up DEtection #########################

SpaceStep = [x, y]
TimeStep = Vector(t)
νS = [15, 15] # number of space basis functions
νT = 40 # number of time basis functions
bufferSpace = [1, 1]
bufferTime = 1
batchSpace = 10
batchTime = 10
learning_rate = [1e-4, 1e-6]
beta = 0.95
# c = 1e6
# c = 138660.0
# c = 1.0


# choose library
ΛSpaceNames = ["Psi", "Psi_x", "Psi_xx", "Psi_y", "Psi_yy", "Psi_xy"]
ΛTimeNames = ["Phi"]

function Λ(A, Ψ, Φ)

  ψ = Ψ[1]
  ψ_x = Ψ[2]
  ψ_xx = Ψ[3]
  ψ_y = Ψ[4]
  ψ_yy = Ψ[5]
  ψ_xy = Ψ[6]

  ϕ = Φ[1]

  U = A * (ϕ ⊗ ψ)'
  U_x = A * (ϕ ⊗ ψ_x)'
  U_xx = A * (ϕ ⊗ ψ_xx)'
  U_y = A * (ϕ ⊗ ψ_y)'
  U_yy = A * (ϕ ⊗ ψ_yy)'
  U_xy = A * (ϕ ⊗ ψ_xy)'

  u = U[1,:]'
  u_x = U_x[1,:]'
  u_xx = U_xx[1,:]'
  u_y = U_y[1,:]'
  u_yy = U_yy[1,:]'
  u_xy = U_xy[1,:]'

  v = U[2,:]'
  v_x = U_x[2,:]'
  v_xx = U_xx[2,:]'
  v_y = U_y[2,:]'
  v_yy = U_yy[2,:]'
  v_xy = U_xy[2,:]'

  return [u, u.^2, u.^3,
    v, v.^2, v.^3,
    u .* v, u.^2 .* v, u .* v.^2,
    u .* u_x, u .* u_y, v .* v_x, v .* v_y,
    u_x, u_y, u_xx, u_yy, u_xy, v_x, v_y, v_xx, v_yy, v_xy]

end

function ∇Λ(A, Ψ, Φ)
  
  ψ = Ψ[1]
  ψ_x = Ψ[2]
  ψ_xx = Ψ[3]
  ψ_y = Ψ[4]
  ψ_yy = Ψ[5]
  ψ_xy = Ψ[6]

  ϕ = Φ[1]

  U_loc = A * (ϕ ⊗ ψ)

  U = BD.dU(A, ϕ, ψ)
  U² = BD.dU²(A, ϕ, ψ)
  U³ = BD.dU³(A, ϕ, ψ)
  z = zeros(size(U,2))

  UUx = BD.dUU_x(A, ϕ, ψ, ψ_x)
  UUy = BD.dUU_y(A, ϕ, ψ, ψ_y)
  U_x = BD.dU_x(A, ϕ, ψ_x)'
  U_xx = BD.dU_xx(A, ϕ, ψ_xx)'
  U_y = BD.dU_y(A, ϕ, ψ_y)'
  U_xy = BD.dU_xy(A, ϕ, ψ_xy)'
  U_yy = BD.dU_yy(A, ϕ, ψ_yy)'

 
  du = hcat(U[1,:], U²[1,:], U³[1,:], z     , z      , z      , U[1,:] .* U_loc[2,1], U²[1,:] .* U_loc[2,1],   U[1,:] .* U_loc[2,1].^2, UUx[1,:], UUy[1,:], z       , z       , U_x, U_y, U_xx, U_yy, U_xy, z  , z  , z   , z   , z)
  dv = hcat(z     , z      , z      , U[1,:], U²[2,:], U³[2,:], U_loc[1,1] .* U[1,:], U_loc[1,1].^2 .* U[1,:], U_loc[1,1] .* U²[2,:]  , z       , z       , UUx[2,:], UUy[2,:], z  , z  , z   , z   , z   , U_x, U_y, U_xx, U_yy, U_xy)
  return cat(du', dv', dims = 3)

end


# model, pars, posterior = DEtection(Y, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, beta, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits=10000)
model, pars, posterior = DEtection(Z, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, beta, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits=10000)


model, pars, posterior = DEtection(Z, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, beta, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits=201, burnin = 1)


# print equation
print_equation(["∇U_t", "∇V_t"], model, pars, posterior; cutoff_prob=0.5)

# posterior summaries
post_sum = posterior_summary(model, pars, posterior)

post_sum.M
post_sum.M_hpd
post_sum.gamma
sqrt.(post_sum.ΣZ)
post_sum.ΣU
post_sum.π


# look at posterior surface
surf_mean, surf_sd = posterior_surface(model, pars, posterior)

contourf(reshape(surf_mean[:, 50, 1], 41, 41), c=:heat)
contourf(U[:, :, 50], c=:heat)

contourf(reshape(surf_mean[:, 50, 2], 41, 41), c=:heat)
contourf(V[:, :, 50], c=:heat)

contourf(U[:, :, 50] - reshape(surf_mean[:, 50, 1], 41, 41), c=:bluesreds)
contourf(V[:, :, 50] - reshape(surf_mean[:, 50, 2], 41, 41), c=:bluesreds)
contourf(reshape(surf_sd[:, 50, 1], 41, 41), c=:heat)
contourf(reshape(surf_sd[:, 50, 2], 41, 41), c=:heat)


# plot chains
plot(posterior.M[1, 1, :])
plot(posterior.M[1, 2, :])
plot(posterior.M[1, 7, :])
plot(posterior.M[1, 16, :])
plot(posterior.M[1, 17, :])

plot(posterior.M[2, 4, :])
plot(posterior.M[2, 7, :])
plot(posterior.M[2, 21, :])
plot(posterior.M[2, 22, :])






