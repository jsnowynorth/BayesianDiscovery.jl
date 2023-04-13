
# The Heat Equation defined as
# Ut = α ΔU = α (Uxx + Uyy)


using BayesianDiscovery
const BD = BayesianDiscovery
using Plots
using Random, Distributions
using Kronecker

include("heat_equation.jl")



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
# ugen(u,x,y) = exp(-(1/500)*((x - 10)^2 + (y-10)^2)^2)
# Q(u,x,y) = exp(-(1/500)*((x - 10)^2 + (y-10)^2)^2)
Q(u, x, y) = 0

# ugen(u,x,y) = exp(-(1/250)*((x - 25)^2 + (y-20)^2)^2) - exp(-(1/250)*((x - 10)^2 + (y-5)^2)^2)

ugen(u, x, y) = sin(2 * pi * x / 40) * cos(2 * pi * y / 40)

heat = heat_gen(x, y, t, α, utop, uleft, uright, ubottom, ugen, Q)

# contourf(heat[:,:,100])
# contourf(heat[:,:,100] + rand(Normal(0.0, 0.01), size(heat[:,:,100])))


# heat = heat[:,:,1:10:1001]
# t = t[1:10:1001]

# subset time
# heat = heat[:, :, 1:1:201]
# t = t[1:1:201]

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

# savefig("../Space_Time_DEtection/results/heat_equation/heat_equation_data.png")

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
beta = 0.9
# c = 1e6
# c = 290062.0
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

model, pars, posterior = DEtection(Y, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, beta, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits=1000, orderSpace=2)
model, pars, posterior = DEtection(Z, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, c, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits=5000, orderSpace=2)


# model, pars, posterior = DEtection(Y, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, c, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits=210, burnin = 10)

# save heat with no noise
# @save "../Space_Time_DEtection/DEtection_results/heat/heat.jld2" model pars posterior
# @load "../Space_Time_DEtection/DEtection_results/heat/heat.jld2" model pars posterior

# save heat with noise 0.02
# @save "../Space_Time_DEtection/DEtection_results/heat/heat_noise2.jld2" model pars posterior
# @load "../Space_Time_DEtection/DEtection_results/heat/heat_noise2.jld2" model pars posterior

# save heat with noise 0.05
# @save "../Space_Time_DEtection/DEtection_results/heat/heat_noise5.jld2" model pars posterior
# @load "../Space_Time_DEtection/DEtection_results/heat/heat_noise5.jld2" model pars posterior

# save heat with noise and missing data
# @save "../Space_Time_DEtection/DEtection_results/heat/heat_noise_missing.jld2" model pars posterior
# @load "../Space_Time_DEtection/DEtection_results/heat/heat_noise_missing.jld2" model pars posterior


print_equation(["∇U_t"], model, pars, posterior; cutoff_prob=0.95)

post_sum = posterior_summary(model, pars, posterior)

post_sum.M
post_sum.M_hpd
post_sum.gamma
sqrt.(post_sum.ΣZ)
post_sum.ΣU
post_sum.π


surf_mean, surf_sd = posterior_surface(model, pars, posterior)

contourf(reshape(Y, 41, 41, 201)[:, :, 100], c=:bluesreds)
contourf(reshape(surf_mean, 41, 41, 201)[:, :, 100], c=:bluesreds)
contourf(heat[:, :, 100] - reshape(surf_mean, 41, 41, 201)[:, :, 100], c=:bluesreds)

contourf(reshape(surf_sd[:, :, 1], 41, 41, 201)[:, :, 10], c=:heat)


plot(posterior.M[:, 8, :]')
plot(posterior.M[:, 20, :]')

plot(posterior.ΣU[1,1,:])
plot(posterior.ΣZ[1,1,:])

plot(posterior.A[1,3524,:])


# b = (length(model.inner_inds)+1)^(1/2)
# c0 = 1
# c1 = 0.99

# (log(c0/b)/log(c1) + 1)

# length(model.inner_inds) - 2*(log(c0/b)/log(c1) + 1)

# c = 290062.0
# model.c = 301421.0

# model.c = 10000.0


c=1.0


pars, model = create_pars(Y, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, c, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, orderSpace=2)

# model = update_c!(pars, model)
pars = update_gamma!(pars, model)
pars = update_M!(pars, model)
pars = update_ΣZ!(pars, model)
pars = update_A!(pars, model)
pars.gamma
pars.M
pars.ΣU
pars.ΣZ


S = size(Y, 1)
T = size(Y, 2)
N = size(Y, 3)

contourf(reshape((fold3(pars.State.State["U"], S, T, N)[:, 10, 1]), 41, 41))
contourf(heat[:,:,10])

contourf(heat[:,:,10] - reshape((fold3(pars.State.State["U"], S, T, N)[:, 10, 1]), 41, 41), c = :bluesreds)



plot(reshape(fold3(pars.State.State["U"], S, T, N), 41, 41, :)[10,20,:])
plot!(heat[10,20,:])

plot(reshape(fold3(pars.State.State["U"], S, T, N), 41, 41, :)[20,20,:])
plot!(heat[20,20,:])

plot(reshape(fold3(pars.State.State["U"], S, T, N), 41, 41, :)[30,5,:])
plot!(heat[30,5,:])








# debugging


# Λnames = ["U", "ΔU_x", "ΔU_y", "ΔU_xy", "ΔU_xx", "ΔU_yy"]

# # lambda function
# function Λ(U, ∂U)

#   ΔU_x = ∂U[1]
#   ΔU_y = ∂U[2]
#   ΔU_xy = ∂U[3]
#   ΔU_xx = ∂U[4]
#   ΔU_yy = ∂U[5]

#   u = U[1]
#   u_x = ΔU_x[1]
#   u_y = ΔU_y[1]
#   u_xy = ΔU_xy[1]
#   u_xx = ΔU_xx[1]
#   u_yy = ΔU_yy[1]

#   return [u, u^2, u^3,
#     u_x, u * u_x, u^2 * u_x, u^3 * u_x,
#     u_xx, u * u_xx, u^2 * u_xx, u^3 * u_xx,
#     u_y, u * u_y, u^2 * u_y, u^3 * u_y,
#     u_xy, u * u_xy, u^2 * u_xy, u^3 * u_xy,
#     u_yy, u * u_yy, u^2 * u_yy, u^3 * u_yy]

# end

pars, model = create_pars(Y, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, c, Λ, Λnames, orderSpace=2)

pars = update_gamma!(pars, model)
pars = update_M!(pars, model)
pars = update_ΣZ!(pars, model)
# pars = update_ΣU!(pars, model)
pars = update_A!(pars, model)
pars.gamma
pars.M

contourf(reshape((fold3(pars.State.State["U"], S, T, N)[:, 10, 1]), 41, 41))


# u, u^2, u^3, u_x, u_xx, u_y, u_xy, u_yy
# Ut = α(Uxx + Uyy)
pars.M
pars.gamma
pars.ΣM
pars.ΣZ
pars.ΣU
pars.A

plot(fold3(pars.State.State["U"], S, T, N)[111, :, :])
plot!(Y[111, :, :])


S = size(Y, 1)
T = size(Y, 2)
N = size(Y, 3)

Urec = fold3(pars.State.State["U"], S, T, N)
# Urec = fold3(pars.State.State["ΔU_t"], S, T, N)

contourf(reshape(Urec[:, 20, 1], 41, 41))
contourf(heat[:, :, 20])

contourf(heat[:, :, 100] - reshape(Urec[:, 100, 1], 41, 41))


# u, u^2, u^3, u_x, u_xx, u_y, u_xy, u_yy
# Ut = Uxx + Uyy
transpose(inv(pars.F[:, model.inner_inds] * pars.F[:, model.inner_inds]') * pars.F[:, model.inner_inds] * pars.State.State["ΔU_t"][:, model.inner_inds]')
heatmap(cor(pars.F, dims=2), c=:bluesreds)



######################### PCA #########################

SpaceStep = [x, y]
TimeStep = Vector(t)
νS = [10, 10] # number of space basis functions
νT = 80 # number of time basis functions
bufferSpace = [1, 1]
bufferTime = 1
batchSpace = 10
batchTime = 10
learning_rate = [1e-4]
c = 1.0

pars, model = create_pars(Z, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, c, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames)


M = fit(PCA, (pars.F .- mean(pars.F, dims = 2)) ./ std(pars.F, dims = 2); pratio = 1.0, maxoutdim=15)
# M = fit(PCA, pars.F; pratio = 1.0)
Ftransformed = MultivariateStats.transform(M, pars.F)
proj = projection(M)

h = plot(Ftransformed[1,:], Ftransformed[2,:], seriestype=:scatter, label="")
plot!(xlabel="PC1", ylabel="PC2", framestyle=:box) # A few formatting options

for i=1:23; plot!([0,proj[i,1]], [0,proj[i,2]], arrow=true, label="X_" * string(i), legend=:bottomleft); end
display(h)

df = DataFrame(var = reshape(split(model.function_names, ", "), model.D), num = Vector(1:model.D), PC1 = loadings(M)[:,1], PC2 = loadings(M)[:,2])

# groupings
df[[1,3],:] # (-, 0)
df[[2],:] # (0, -)
df[[8,10,20,22],:] # (+, 0)
df[[9,11,21,23],:] # (0, +)
df[[4,5,6,7,12,13,14,15,16,17,18,19],:] # (0, 0)



heatmap(cor(pars.F, dims=2), c=:balance, clim = (-1, 1))

function corr_preds(N, α, ρ, g)

  h = -N/(1)*(1-ρ^2)*α^2 * (g/(g+1)) + log(g+1) # g-prior
  
  return 1/(1 + exp(h/2))
  
end

N = [20, 40, 100, 500, 1000, 10000, 100000]
a = range(0, 1, step = 0.01)
ρ = [0, 0.2, 0.5, 0.8, 0.99]

plt = [corr_preds(N[i], a[j], ρ[k]) for i in 1:length(N), j in 1:length(a), k in 1:length(ρ)]
plt = [corr_preds(1000, a[j], ρ[k], 256*101) for j in 1:length(a), k in 1:length(ρ)]

plot(a, plt, xlabel = "Effect Size", ylabel = "Inclusion Probability", label = reshape(["ρ = " * string(ρ[k]) for k in 1:length(ρ)], 1, :), legend=:bottomright)

(inv(pars.F * pars.F') * pars.F * pars.State.State["ΔU_t"]')'

heatmap(inv(pars.F * pars.F'))