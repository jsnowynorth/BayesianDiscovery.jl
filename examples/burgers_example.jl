
# Ut = α Uxx - U Ux
# α = 0.1


# using StatsBase, LinearAlgebra, Plots, StatsPlots, CSV
# using Kronecker
# using Missings
# using Distributions, Random
# using ProgressMeter
# using ReverseDiff: JacobianTape, JacobianConfig, jacobian, jacobian!, compile
# using ForwardDiff
# using Statistics
# using JLD2
# using RCall
# using DataFrames, DataFramesMeta, Chain
# using NetCDF
# using MultivariateStats
# using StatsModels, Combinatorics, IterTools
# using KernelDensity
# using Tables
# using MAT
# using SparseArrays
# using CodeTracking, Revise
# using FFTW, DifferentialEquations


using BayesianDiscovery
const BD = BayesianDiscovery
# using FFTW, DifferentialEquations
using FFTW, OrdinaryDiffEq
using Plots
using Random, Distributions
using Kronecker

include("burgers_equation_generate.jl")


######################### Load Burgers Equation #########################

U, x, Time = burgers()
contourf(Time, x, U, c=:oxy)


######################### set up data #########################

Random.seed!(1)
Y = reshape(real.(U), 256, 101, 1)
Z = Y + 0.02 * std(Y) .* rand(Normal(), 256, 101, 1)
# Z = reshape([Z[i, j] < 0 ? abs(Z[i, j]) : Z[i, j] for i in 1:256, j in 1:101], 256, 101, 1)


######################### missing data #########################

Y_missing = Array{Union{Missing,Float64}}(missing, size(Y))
Z_missing = Array{Union{Missing,Float64}}(missing, size(Y))

perc_missing = 0.05

obs_inds = sample(reshape([[x, y] for x = 1:256, y = 1:101], prod(size(Y))), convert(Int, ceil(prod(size(Y)) * (1 - perc_missing))), replace=false)
for i in 1:Int(ceil(prod(size(Y)) * (1 - perc_missing)))
  Y_missing[obs_inds[i][1], obs_inds[i][2], :] = copy(Y[obs_inds[i][1], obs_inds[i][2], :])
  Z_missing[obs_inds[i][1], obs_inds[i][2], :] = copy(Z[obs_inds[i][1], obs_inds[i][2], :])
end


######################### PDEtection #########################
## PDEtection

# estimate_c(model.inner_inds, cut_off_prob = 0.81)

SpaceStep = [x]
TimeStep = Time
νS = 50 # number of space basis functions
νT = 20 # number of time basis functions
# bufferSpace = 1
# bufferTime = 1
batchSpace = 10
batchTime = 10
learning_rate = 1e-4
beta = 0.9
# c = 25098.0
# c = 1.0

ΛSpaceNames = ["Psi", "Psi_x", "Psi_xx", "Psi_xxx"]
ΛTimeNames = ["Phi"]
# choose library

function Λ(A, Ψ, Φ)
    
  ψ = Ψ[1]
  ψ_x = Ψ[2]
  ψ_xx = Ψ[3]
  ψ_xxx = Ψ[4]

  ϕ = Φ[1]

  u = A * (ϕ ⊗ ψ)'
  u_x = A * (ϕ ⊗ ψ_x)'
  u_xx = A * (ϕ ⊗ ψ_xx)'
  u_xxx = A * (ϕ ⊗ ψ_xxx)'

  return [u, u.^2, u.^3,
    u_x, u .* u_x, u.^2 .* u_x, u.^3 .* u_x,
    u_xx, u .* u_xx, u.^2 .* u_xx, u.^3 .* u_xx,
    u_xxx, u .* u_xxx, u.^2 .* u_xxx, u.^3 .* u_xxx]

end

function ∇Λ(A, Ψ, Φ)
  
  ψ = Ψ[1]
  ψ_x = Ψ[2]
  ψ_xx = Ψ[3]
  ψ_xxx = Ψ[4]

  ϕ = Φ[1]

  return [BD.dU(A, ϕ, ψ), BD.dU²(A, ϕ, ψ), BD.dU³(A, ϕ, ψ),
  BD.dU_x(A, ϕ, ψ_x), BD.dUU_x(A, ϕ, ψ, ψ_x), BD.dU²U_x(A, ϕ, ψ, ψ_x), BD.dU³U_x(A, ϕ, ψ, ψ_x),
  BD.dU_xx(A, ϕ, ψ_xx), BD.dUU_xx(A, ϕ, ψ, ψ_xx), BD.dU²U_xx(A, ϕ, ψ, ψ_xx), BD.dU³U_xx(A, ϕ, ψ, ψ_xx),
  BD.dU_xxx(A, ϕ, ψ_xxx), BD.dUU_xxx(A, ϕ, ψ, ψ_xxx), BD.dU²U_xxx(A, ϕ, ψ, ψ_xxx), BD.dU³U_xxx(A, ϕ, ψ, ψ_xxx)]

end

model, pars, posterior = DEtection(Z, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, beta, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits = 1000, burnin = 500)
model, pars, posterior = DEtection(Z_missing, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, c, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits = 1000, burnin = 500)


print_equation(["∇U"], model, pars, posterior, cutoff_prob=0.5)


post = posterior_summary(model, pars, posterior)
post.M
post.M_hpd
post.gamma
sqrt.(post.ΣZ)
sqrt.(post.ΣU)
post.π


post_mean, post_sd = posterior_surface(model, pars, posterior)

contourf(post_mean[:,:,1], clim = extrema(U))
contourf(reshape(Y, 256, 101))
contourf(U .- post_mean[:,:,1])
contourf(post_sd[:,:,1], clim = (0, 0.0011))

[quantile(reshape(post_sd[:,:,1], :), i) for i in range(0, 1, step = 0.01)]



plot(post_mean[100, :, 1])
plot!(U[100, :, 1])
plot!(Z[100, :, 1])

plot(post_mean[:, 10, 1])
plot!(U[:, 10, 1])
plot!(Z[:, 10, 1])

plot(posterior.M[1,5,:])
plot(posterior.M[1,8,:])
plot(posterior.M[1,9,:])
plot(posterior.M[1,10,:])

plot(posterior.ΣU[1,1,:])
plot(posterior.ΣZ[1,1,:])


plot(posterior.A[1,310,:])


############################### debugging ###############################



######################### Some figures for paper #########################

using LaTeXStrings

pars, model = create_pars(Y, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, c, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames)


# par_names = String.(reshape(split(model.function_names, ", "), :))
par_names = [L"u", L"u^2", L"u^3", L"u_x", L"uu_x", L"u^2u_x", L"u^3u_x", L"u_{xx}", L"uu_{xx}", L"u^2u_{xx}", L"u^3u_{xx}", L"u_{xxx}", L"uu_{xxx}", L"u^2u_{xxx}", L"u^3u_{xxx}"]
# heatmap(par_names, par_names, cor(pars.F, dims=2), c=:balance, xrotation = -90, clim = (-1, 1), right_margin = 10Plots.mm, xtickfontsize=12, ytickfontsize=12)

# par_names = [L"u", L"u^2", L"u_x", L"uu_x", L"u^2u_x", L"u_{xx}", L"uu_{xx}", L"u^2u_{xx}", L"u_{xxx}", L"uu_{xxx}", L"u^2u_{xxx}"]


l = @layout [a b{0.01w}]
p1 = heatmap(par_names, par_names, cor(pars.F, dims=2), c=:balance, xrotation = -90, framestyle=:box,  clim = (-1, 1), tickfontsize=18, legend = false, right_margin = 0Plots.mm)
p2 = scatter([0, 0], [0, 0], zcolor=[0, 1], clims=(-1,1), xlims=(0.9, 1.1), ylims=(0, 1), label="", framestyle=:none, c=:balance, tickfontsize=12)
plot(p1, p2, layout=l, size=(750, 500), right_margin = 12Plots.mm, left_margin = 0Plots.mm, bottom_margin=4Plots.mm, top_margin=2Plots.mm)

# savefig("../Space_Time_DEtection/results/burgers/burgers_corr.png")



heatmap(cor(pars.F, dims=2), c=:balance, xrotation = -90, framestyle=:box,  clim = (-1, 1), tickfontsize=18, right_margin = 12Plots.mm)


Random.seed!(386)
samp_space, samp_time, samp_inds = sample_inds(model.SpaceStep, model.TimeStep, model.batchSpace, model.batchTime, model.bufferSpace, model.bufferTime)

contourf(reshape(Time, :), reshape(x, :), U, legend=false, ticks=false, framestyle=:none, c=:oxy)
scatter!(reshape(Time, :)[samp_time], reshape(x, :)[samp_space], legend = false, markershape = :rect, markerstrokecolor = :black, markerstrokewidth = 2, marker_z = reshape(U, :)[samp_inds], markercolor = :oxy, clim = extrema(U))

# savefig("../Space_Time_DEtection/results/burgers/burgers_subsample.png")





######################### add covariates #########################

SpaceStep = [x]
TimeStep = Time
νS = 50 # number of space basis functions
νT = 20 # number of time basis functions
bufferSpace = 1
bufferTime = 1
batchSpace = 10
batchTime = 10
learning_rate = 1e-4
# c = 25098.0
c = 1.0


# choose library
ΛSpaceNames = ["Psi", "Psi_x", "Psi_xx"]
ΛTimeNames = ["Phi"]
# choose library

function Λ(A, Ψ, Φ, X)
    
  ψ = Ψ[1]
  ψ_x = Ψ[2]
  ψ_xx = Ψ[3]

  ϕ = Φ[1]

  u = A * (ϕ ⊗ ψ)'
  u_x = A * (ϕ ⊗ ψ_x)'
  u_xx = A * (ϕ ⊗ ψ_xx)'

  # covariates
  x1 = X[1]
  x2 = X[2]
  x3 = X[3]
  x4 = X[4]
  x5 = X[5]
  x6 = X[6]
  x7 = X[7]
  x8 = X[8]
  x9 = X[9]
  x10 = X[10]

  return [u .* u_x, u_xx, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]

end

function ∇Λ(A, Ψ, Φ, X)
  
  ψ = Ψ[1]
  ψ_x = Ψ[2]
  ψ_xx = Ψ[3]

  ϕ = Φ[1]

  # covariates
  x1 = X[1]
  x2 = X[2]
  x3 = X[3]
  x4 = X[4]
  x5 = X[5]
  x6 = X[6]
  x7 = X[7]
  x8 = X[8]
  x9 = X[9]
  x10 = X[10]

  uu_x = dUU_x(A, ϕ, ψ, ψ_x)
  u_xx = dU_xx(A, ϕ, ψ_xx)
  z = zeros(1, size(uu_x,2))

  return [uu_x, u_xx, z, z, z, z, z, z, z, z, z, z]

end


Random.seed!(1)
# mean(Y)
# std(Y)
covs = rand(Normal(mean(Y),std(Y)), 256, 101, 10)
# covs = rand(Normal(), 256, 101, 10)


# GP for random covariate
using Distances

h = reshape([[x,Time]  for x=x, Time=Time],length(x)*length(Time))
d = pairwise(euclidean, h)

C = exp.(-d./(2*(0.8*10)^2))
C = C + diagm(0.00000001*ones(length(x)*length(Time)))

Random.seed!(1)
y = rand(MvNormal(zeros(length(x)*length(Time)), C), 10)

contourf(x, Time, reshape(y[:,10], length(x), length(Time))')
heatmap(x, Time, reshape(y[:,9], length(x), length(Time))', c = :balance)


covs = reshape(covs, length(x), length(Time), 10)



model, pars, posterior = DEtection(Z, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, c, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits = 2000, covariates=covs)

pars, model = create_pars(Y, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, c, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, covariates=covs)

pars = update_gamma!(pars, model)
pars = update_M!(pars, model)
pars = update_ΣZ!(pars, model)
pars = update_A!(pars, model)
pars.M

model, pars, posterior = DEtection(Y, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, v0, v1, Λ1, Λnames, nits=100, covariates=covs)
model, pars, posterior = DEtection(Y, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, v0, v1, Λ2, Λnames, nits=100)

print_equation(["∇U_t"], model, pars, posterior; cutoff_prob=0.5)

post_sum = posterior_summary(model, pars, posterior)

post_sum.M
post_sum.M_hpd
post_sum.gamma
sqrt(post_sum.ΣZ)
post_sum.ΣU





N = 256*101
g = N
N = 100
s = 0.01
r = 0.9

sqrt((g/(g+1))*log(g+1)*s/(N*(1-r)))


# N = 256*101
g = N
# N = 100
s = sqrt(pars.ΣU[1,1])
r = 0.95
a = 0.01

((g+1)/(g))*(s/(a^2*(1-r)))*(10+log(g+1))

-100*a^2*(1-r)*(g/(g+1))/s + log(g+1)

exp(9)