
# The Heat Equation defined as
# Ut = α ΔU = α (Uxx + Uyy)


using BayesianDiscovery
const BD = BayesianDiscovery
using Plots
using Random, Distributions
using Kronecker

include("heat_equation_generate.jl")



######################### generate data #########################


## heat equation parameters
α = 1
Δx = 0.5
Δy = 0.5
Δt = 0.01

(4 * Δt * α) / (Δx^2 * Δy^2) ≤ 1

x = range(0, 20, step=Δx)
y = range(0, 20, step=Δy)
t = range(0, 2, step=Δt)


## heat equation with gaussian kernel and generation
utop = 0
uleft = 0
uright = 0
ubottom = 0
Q(u, x, y) = 0

ugen(u, x, y) = sin(2 * pi * x / 40) * cos(2 * pi * y / 40)

heat = heat_gen(x, y, t, α, utop, uleft, uright, ubottom, ugen, Q)



######################### set up data #########################

Random.seed!(1)
Y = reshape(heat, :, length(t), 1)
Z = Y + 0.02 * std(Y) * rand(Normal(), size(Y))

######################### missing data #########################

Y_missing = Array{Union{Missing,Float64}}(missing, size(Y))
Z_missing = Array{Union{Missing,Float64}}(missing, size(Y))

perc_missing = 0.01

obs_inds = sample(reshape([[x, y] for x = 1:size(Y, 1), y = 1:size(Y, 2)], prod(size(Y))), convert(Int, ceil(prod(size(Y)) * (1 - perc_missing))), replace=false)
for i in 1:Int(ceil(prod(size(Y)) * (1 - perc_missing)))
  Y_missing[obs_inds[i][1], obs_inds[i][2], :] = copy(Y[obs_inds[i][1], obs_inds[i][2], :])
  Z_missing[obs_inds[i][1], obs_inds[i][2], :] = copy(Z[obs_inds[i][1], obs_inds[i][2], :])
end


######################### visualize data #########################

# Nx = length(x)
# Ny = length(y)

# l = @layout [grid(1, 2) a{0.07w}]

# p1 = contourf(y, x, reshape(Z[:, 1, 1], Nx, Ny), clim=(-1.05, 1.05), legend=false, ticks=false, framestyle=:box, c=:bluesreds)
# p2 = contourf(y, x, reshape(Z[:, 201, 1], Nx, Ny), clim=(-1.05, 1.05), legend=false, ticks=false, framestyle=:box, c=:bluesreds)

# l1 = scatter([0, 0], [0, 0], zcolor=[0, 1], clims=(-1.05, 1.05), xlims=(0.9, 1.1), ylims=(0, 1), label="", framestyle=:none, c=:bluesreds, xtickfontsize=18, ytickfontsize=18)

# plot(p1, p2, l1, layout=l, size=(1500, 750), right_margin=5Plots.mm)


######################### set up detection #########################

SpaceStep = [x, y]
TimeStep = Vector(t)
νS = [10, 10] # number of space basis functions
νT = 80 # number of time basis functions
bufferSpace = [1, 1]
bufferTime = 1
batchSpace = 10
batchTime = 10
learning_rate = [1e-4]
beta = 0.99



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

  u = A * (ϕ ⊗ ψ)'
  u_x = A * (ϕ ⊗ ψ_x)'
  u_xx = A * (ϕ ⊗ ψ_xx)'
  u_y = A * (ϕ ⊗ ψ_y)'
  u_yy = A * (ϕ ⊗ ψ_yy)'
  u_xy = A * (ϕ ⊗ ψ_xy)'

  return [u, u.^2, u.^3,
    u_x, u .* u_x, u.^2 .* u_x, u.^3 .* u_x,
    u_xx, u .* u_xx, u.^2 .* u_xx, u.^3 .* u_xx,
    u_y, u .* u_y, u.^2 .* u_y, u.^3 .* u_y,
    u_xy, u .* u_xy, u.^2 .* u_xy, u.^3 .* u_xy,
    u_yy, u .* u_yy, u.^2 .* u_yy, u.^3 .* u_yy]

end

function ∇Λ(A, Ψ, Φ)
  
  ψ = Ψ[1]
  ψ_x = Ψ[2]
  ψ_xx = Ψ[3]
  ψ_y = Ψ[4]
  ψ_yy = Ψ[5]
  ψ_xy = Ψ[6]

  ϕ = Φ[1]

  return [BD.dU(A, ϕ, ψ), BD.dU²(A, ϕ, ψ), BD.dU³(A, ϕ, ψ),
  BD.dU_x(A, ϕ, ψ_x), BD.dUU_x(A, ϕ, ψ, ψ_x), BD.dU²U_x(A, ϕ, ψ, ψ_x), BD.dU³U_x(A, ϕ, ψ, ψ_x),
  BD.dU_xx(A, ϕ, ψ_xx), BD.dUU_xx(A, ϕ, ψ, ψ_xx), BD.dU²U_xx(A, ϕ, ψ, ψ_xx), BD.dU³U_xx(A, ϕ, ψ, ψ_xx),
  BD.dU_y(A, ϕ, ψ_y), BD.dUU_y(A, ϕ, ψ, ψ_y), BD.dU²U_y(A, ϕ, ψ, ψ_y), BD.dU³U_y(A, ϕ, ψ, ψ_y),
  BD.dU_xy(A, ϕ, ψ_xy), BD.dUU_xy(A, ϕ, ψ, ψ_xy), BD.dU²U_xy(A, ϕ, ψ, ψ_xy), BD.dU³U_xy(A, ϕ, ψ, ψ_xy),
  BD.dU_yy(A, ϕ, ψ_yy), BD.dUU_yy(A, ϕ, ψ, ψ_yy), BD.dU²U_yy(A, ϕ, ψ, ψ_yy), BD.dU³U_yy(A, ϕ, ψ, ψ_yy)]

end

# run the model
# model, pars, posterior = DEtection(Y, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, beta, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits=10000, orderSpace=2)
model, pars, posterior = DEtection(Z, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, beta, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits=10000, orderSpace=2)


# look at the found equation
print_equation(["∇U_t"], model, pars, posterior; cutoff_prob=0.5)

# compute posterior summaries
post_sum = posterior_summary(model, pars, posterior)

post_sum.M
post_sum.M_hpd
post_sum.gamma
sqrt.(post_sum.ΣZ)
post_sum.ΣU
post_sum.π


# compute posterior surface
surf_mean, surf_sd = posterior_surface(model, pars, posterior)

# plot surface
contourf(reshape(Y, 41, 41, 201)[:, :, 100], c=:bluesreds)
contourf(reshape(surf_mean, 41, 41, 201)[:, :, 100], c=:bluesreds)
contourf(heat[:, :, 100] - reshape(surf_mean, 41, 41, 201)[:, :, 100], c=:bluesreds)
contourf(reshape(surf_sd[:, :, 1], 41, 41, 201)[:, :, 10], c=:heat)


# look at chains of found parameters
plot(posterior.M[:, 8, :]')
plot(posterior.M[:, 20, :]')
