

# """
#     Model Structure

# Initiates a structure of class Model

# # Arguments
# - S: Spatial Dimension
# - T
# - N
# - D
# - νS
# - νT
# - bufferSpace
# - bufferTime
# - batchSpace
# - batchTime
# - learning_rate
# - Z
# - U
# - Basis
# - States
# - Λ::Function: Library of potential functions
# - Λnames::Vector{String}: Names of components of library of potential functions

# """
# mutable struct Model

#   # model
#   # S, T, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, Z, U, Basis, States

#   # dimensions
#   S::Int
#   T::Int
#   L::Int
#   N::Int
#   D::Int

#   # data parameters
#   νS
#   νT
#   bufferSpace
#   bufferTime
#   batchSpace
#   batchTime
#   SpaceStep
#   TimeStep
#   learning_rate
#   inner_inds
#   space_inds
#   time_inds
#   all_inds

#   # data
#   Z::Array{Float64}
#   H
#   X

#   # basis functions
#   Basis::DerivativeClass
#   Θ

#   # function information
#   Λ::Function
#   ∇Λ::Function
#   ΛSpaceNames::Vector{String}
#   ΛTimeNames::Vector{String}

#   # response information
#   responseFunction::Function
#   responseNames::Vector{String}
#   Ψnames::Vector{String}

#   # spike-and-slab parameters
#   c::Float64

#   function_names

# end

# Base.show(io::IO, model::Model) =
#   print(io, "Model\n",
#     " ├─── data dimensions: ", [model.S, model.T, model.L], '\n',
#     " ├─── process dimensions: ", [model.S, model.T, model.N], '\n',
#     " ├─── response: ", model.responseNames, '\n',
#     " ├─── library: ", model.function_names, '\n',
#     " ├─── covariates: ", model.X === nothing ? false : size(model.X, 1), '\n',
#     " ├─── spatial domain: ", model.SpaceStep, '\n',
#     " ├─── temporal domain: ", model.TimeStep, '\n',
#     " ├─── number of spatial basis functions: ", model.νS, '\n',
#     " ├─── number of temporal basis functions: ", model.νT, '\n',
#     " ├─── spike-and-slab parameter: ", string(model.c), '\n',
#     " ├─── learning rate: ", model.learning_rate, '\n',
#     " ├─── space buffer: ", model.bufferSpace, '\n',
#     " ├─── time buffer: ", model.bufferTime, '\n',
#     " ├─── space batch size: ", model.batchSpace, '\n',
#     " └─── time batch size: ", model.batchTime, '\n')
# #

# """
#     Pars Structure

# Initiates a structure of class Pars

# # Arguments
# - A: basis coefficients
# - M: PDE coefficients
# - ΣZ : measurement variance/covariance matrix
# - ΣU: process variance/covariance matrix
# - gamma: ssvs latent parameters
# - sample_inds: current subsampled indicies for SGD
# - F: current F
# - Fprime: current Fprime
# - a_R
# - A_R
# - nu_R
# - a_Q
# - A_Q
# - nu_Q
# - v0
# - v1
# - Sigma_M

# """
# mutable struct Pars

#   # parameters
#   # A, M, ΣY, ΣU, gamma, F, a_R, A_R, nu_R, a_Q, A_Q, nu_Q, v0, v1, Sigma_M

#   # estimated parameters
#   A::Array{Float64}
#   M::Array{Float64}
#   ΣZ::Array{Float64}
#   ΣU::Array{Float64}
#   gamma::Array{Float64}

#   # current value of function
#   sample_inds::Vector{Int}
#   F
#   Fprime
#   State::StateClass

#   # hyperpriors
#   a_Z
#   A_Z
#   nu_Z::Int

#   a_U
#   A_U
#   nu_U::Int

#   # SSVS parameters
#   # v0::Float64
#   # v1::Float64
#   ΣM::Array{Float64}
#   π

# end


# """
#     create_pars()

# Constructs the parameter and model classes.

# See also [`Pars`](@ref), [`Model`](@ref).

# # Arguments
# - Y::Array{Float64, 3}: Space x Time x Component data array
# - SpaceStep::Vector{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}: Vector of spatial locations [x, y]
# - TimeStep::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}: Time samples
# - νS::Vector{Int32}: number of spatial basis functions [x, y]
# - νT::Int: number of temporal basis functions
# - bufferSpace::Vector{Int}: spatial buffer [x, y]
# - bufferTime::Int: temporal buffer
# - batchSpace::Int: batch size for space [x, y]
# - batchTime::Int: batch size for time
# - learning_rate::Float64: learning rate
# - v0::Float64: ssvs not included variance
# - v1::Float64: ssvs included variance
# - Λ::Function: function library
# - Λnames::Vector{String}: names of arguments in function library

# """
# function create_pars(Y::Array{Float64,3},
#   SpaceStep::Vector,
#   TimeStep::Vector,
#   νS::Vector{Int},
#   νT::Int,
#   bufferSpace::Vector{Int},
#   bufferTime::Int,
#   batchSpace::Int,
#   batchTime::Int,
#   learning_rate::Vector{Float64},
#   c::Float64,
#   Λ::Function,
#   ∇Λ::Function,
#   ΛSpaceNames::Vector{String},
#   ΛTimeNames::Vector{String};
#   response=nothing, responsenames=nothing, degree=4, orderTime=1, orderSpace=3, latent_dim=size(Y, 3), covariates=nothing)

#   S = size(Y, 1)
#   T = size(Y, 2)
#   L = size(Y, 3)
#   N = latent_dim # need to update

#   Z = unfold3(Y)

#   # Δt = TimeStep[2] - TimeStep[1]
#   # ΔS = [round(SpaceStep[1][2] - SpaceStep[1][1], digits=4), round(SpaceStep[2][2] - SpaceStep[2][1], digits=4)]

#   # find inner index values
#   Nx = length(SpaceStep[1])
#   Ny = length(SpaceStep[2])
#   Ns = Nx * Ny
#   Nt = length(TimeStep)
#   space_inds = 1:Ns
#   time_inds = 1:Nt
#   ind_comb = reshape([[x, y] for x = space_inds, y = time_inds], Ns * Nt)
#   ind_comb = permutedims(hcat(ind_comb...))
#   ind_comb = hcat(ind_comb, 1:(Ns*Nt))
#   space_inds = transpose(reduce(hcat, reshape([[x, y] for x = 1:Nx, y = 1:Ny], Ns)))
#   ix = space_inds[:, 1] .∉ [cat(1:1:bufferSpace[1], (Nx-bufferSpace[1]+1):1:Nx, dims=1)]
#   iy = space_inds[:, 2] .∉ [cat(1:1:bufferSpace[2], (Ny-bufferSpace[2]+1):1:Ny, dims=1)]
#   ispace = cat(1:Ns, dims=2)[ix.&iy]
#   i1 = ind_comb[:, 1] .∈ [ispace]
#   i2 = ind_comb[:, 2] .∉ [cat(1:1:bufferTime, (Nt-bufferTime+1):1:Nt, dims=1)]
#   inner_inds = ind_comb[i1.&i2, 3]
#   space_inds = ind_comb[:, 1]
#   time_inds = ind_comb[:, 2]
#   all_inds = ind_comb[:, 3]

#   if covariates === nothing
#     X = nothing
#   else
#     X = reduce(vcat, [reshape(covs[:, :, i], 1, :) for i in 1:size(covs, 3)])
#   end

#   # get missing values - no missing here
#   H = ones(L, S * T)
#   H = [make_H(L, N, H[:, i]) for i in 1:(S*T)]

#   # basis functions
#   Basis = derivatives_dictionary(orderSpace, orderTime, νS, νT, SpaceStep, TimeStep, degree=degree)

#   Ψ = Basis.SpaceDerivative[Basis.SpaceNames[1]]
#   Φ = Basis.TimeDerivative[Basis.TimeNames[1]]
#   Θ = Matrix(Diagonal(ones(N)))
#   A = copy(inv(transpose(Θ) * Θ) * transpose(Θ) * Z * (Φ ⊗ Ψ) * inv(transpose(Φ ⊗ Ψ) * (Φ ⊗ Ψ)))

#   State = construct_state(A, Basis, Θ)

#   # parameters
#   samp_space, samp_time, samp_inds = sample_inds(SpaceStep, TimeStep, batchSpace, batchTime, bufferSpace, bufferTime)
#   # F = create_library(Λ, Λnames, State, X)
#   # dF = create_library_derivative(Λ, Λnames, State, samp_inds, X)
#   F = create_Λ(Λ, A, Basis, ΛSpaceNames, ΛTimeNames, X)
#   dF = create_∇Λ(∇Λ, A, Basis, samp_space, samp_time, samp_inds, ΛSpaceNames, ΛTimeNames, X)

#   D = size(F, 1)
#   M = zeros(N, D)
#   ΣZ = Diagonal(ones(N))
#   ΣU = Diagonal(ones(N))

#   # hyperpriors
#   a_Z = ones(N)
#   A_Z = fill(1e6, N)
#   nu_Z = 2

#   a_U = ones(N)
#   A_U = fill(1e6, N)
#   nu_U = 2

#   # SSVS parameters
#   gamma = ones(N, D)
#   ΣM = inv(c * Diagonal(ones(D)))
#   π = 0.5 * ones(N)


#   # determine response
#   if response === nothing
#     responseNames = ["ΔU_t"]
#     function responseFunction(ΔU)
#       ΔU_t = ΔU[1]
#       return ΔU_t
#     end
#   else
#     responseFunction = response
#     responseNames = responsenames
#   end
#   # else is user input for response function and response names

#   Ψnames = (response === nothing ? ["Psi"] : "Psi_" .* reduce(hcat, split.(reduce(hcat, split.(responsenames, "t"))[1, :], "_"))[2, :])

#   Ψtmp = [Basis.SpaceDerivative[ΛSpaceNames[i]] for i in 1:length(ΛSpaceNames)]
#   Φtmp = [Basis.TimeDerivative[ΛTimeNames[i]] for i in 1:length(ΛTimeNames)]

#   if X === nothing
#     function_names = replace((@code_string Λ(A, Ψtmp, Φtmp)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
#   else
#     function_names = replace((@code_string Λ(A, Ψtmp, Φtmp, X)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp, X)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
#   end

#   # store model and parameter values
#   # model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, Z, H, Basis, Θ, Λ, dΛ, Λnames, response)
#   model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, space_inds, time_inds, all_inds, Z, H, X, Basis, Θ, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, responseFunction, responseNames, Ψnames, c, function_names)
#   pars = Pars(A, M, ΣZ, ΣU, gamma, samp_inds, F, dF, State, a_Z, A_Z, nu_Z, a_U, A_U, nu_U, ΣM, π)

#   return pars, model

# end # 2 dimensions

# function create_pars(Y::Array{Union{Missing,Float64},3},
#   SpaceStep::Vector,
#   TimeStep::Vector,
#   νS::Vector{Int},
#   νT::Int,
#   bufferSpace::Vector{Int},
#   bufferTime::Int,
#   batchSpace::Int,
#   batchTime::Int,
#   learning_rate::Vector{Float64},
#   c::Float64,
#   Λ::Function,
#   ∇Λ::Function,
#   ΛSpaceNames::Vector{String},
#   ΛTimeNames::Vector{String};
#   response=nothing, responsenames=nothing, degree=4, orderTime=1, orderSpace=3, latent_dim=size(Y, 3), covariates=nothing)

#   S = size(Y, 1)
#   T = size(Y, 2)
#   L = size(Y, 3)
#   N = latent_dim # need to update

#   Z = unfold3(collect(Missings.replace(Y, 0)))

#   Δt = TimeStep[2] - TimeStep[1]
#   ΔS = [round(SpaceStep[1][2] - SpaceStep[1][1], digits=4), round(SpaceStep[2][2] - SpaceStep[2][1], digits=4)]

#   # find inner index values
#   Nx = length(SpaceStep[1])
#   Ny = length(SpaceStep[2])
#   Ns = Nx * Ny
#   Nt = length(TimeStep)
#   space_inds = 1:Ns
#   time_inds = 1:Nt
#   ind_comb = reshape([[x, y] for x = space_inds, y = time_inds], Ns * Nt)
#   ind_comb = permutedims(hcat(ind_comb...))
#   ind_comb = hcat(ind_comb, 1:(Ns*Nt))
#   space_inds = transpose(reduce(hcat, reshape([[x, y] for x = 1:Nx, y = 1:Ny], Ns)))
#   ix = space_inds[:, 1] .∉ [cat(1:1:bufferSpace[1], (Nx-bufferSpace[1]+1):1:Nx, dims=1)]
#   iy = space_inds[:, 2] .∉ [cat(1:1:bufferSpace[2], (Ny-bufferSpace[2]+1):1:Ny, dims=1)]
#   ispace = cat(1:Ns, dims=2)[ix.&iy]
#   i1 = ind_comb[:, 1] .∈ [ispace]
#   i2 = ind_comb[:, 2] .∉ [cat(1:1:bufferTime, (Nt-bufferTime+1):1:Nt, dims=1)]
#   inner_inds = ind_comb[i1.&i2, 3]
#   space_inds = ind_comb[:, 1]
#   time_inds = ind_comb[:, 2]
#   all_inds = ind_comb[:, 3]

#   if covariates === nothing
#     X = nothing
#   else
#     X = reduce(vcat, [reshape(covs[:, :, i], 1, :) for i in 1:size(covs, 3)])
#   end

#   # get missing values
#   na_inds = getindex.(findall(ismissing, unfold3(Y)), 2)
#   H = ones(L, S * T)
#   if (size(na_inds)[1] != 0)
#     H[na_inds] .= 0
#   end
#   H = [make_H(L, N, H[:, i]) for i in 1:(S*T)]

#   # basis functions
#   Basis = derivatives_dictionary(orderSpace, orderTime, νS, νT, ΔS, Δt, SpaceStep, TimeStep, degree=degree)

#   Ψ = Basis.SpaceDerivative[Basis.SpaceNames[1]]
#   Φ = Basis.TimeDerivative[Basis.TimeNames[1]]
#   Θ = Matrix(Diagonal(ones(N)))
#   A = copy(inv(transpose(Θ) * Θ) * transpose(Θ) * Z * (Φ ⊗ Ψ) * inv(transpose(Φ ⊗ Ψ) * (Φ ⊗ Ψ)))

#   State = construct_state(A, Basis, Θ)

#   # parameters
#   samp_space, samp_time, samp_inds = sample_inds(SpaceStep, TimeStep, batchSpace, batchTime, bufferSpace, bufferTime)
#   # F = create_library(Λ, Λnames, State, X)
#   # dF = create_library_derivative(Λ, Λnames, State, samp_inds, X)
#   F = create_Λ(Λ, A, Basis, ΛSpaceNames, ΛTimeNames, X)
#   dF = create_∇Λ(∇Λ, A, Basis, samp_space, samp_time, samp_inds, ΛSpaceNames, ΛTimeNames, X)
  
#   D = size(F, 1)
#   M = zeros(N, D)
#   ΣZ = Diagonal(ones(N))
#   ΣU = Diagonal(ones(N))

#   # hyperpriors
#   a_Z = ones(N)
#   A_Z = fill(1e6, N)
#   nu_Z = 2

#   a_U = ones(N)
#   A_U = fill(1e6, N)
#   nu_U = 2

#   # SSVS parameters
#   gamma = ones(N, D)
#   ΣM = inv(c * Diagonal(ones(D)))
#   π = 0.5 * ones(N)

#   # determine response
#   if response === nothing
#     responseNames = ["ΔU_t"]
#     function responseFunction(ΔU)
#       ΔU_t = ΔU[1]
#       return ΔU_t
#     end
#   else
#     responseFunction = response
#     responseNames = responsenames
#   end
#   # else is user input for response function and response names

#   Ψnames = (response === nothing ? ["Psi"] : "Psi_" .* reduce(hcat, split.(reduce(hcat, split.(responsenames, "t"))[1, :], "_"))[2, :])

#   Ψtmp = [Basis.SpaceDerivative[ΛSpaceNames[i]] for i in 1:length(ΛSpaceNames)]
#   Φtmp = [Basis.TimeDerivative[ΛTimeNames[i]] for i in 1:length(ΛTimeNames)]
#   if X === nothing
#     function_names = replace((@code_string Λ(A, Ψtmp, Φtmp)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
#   else
#     function_names = replace((@code_string Λ(A, Ψtmp, Φtmp, X)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp, X)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
#   end

#   # store model and parameter values
#   # model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, Z, H, Basis, Θ, Λ, dΛ, Λnames, response)
#   model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, space_inds, time_inds, all_inds, Z, H, X, Basis, Θ, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, responseFunction, responseNames, Ψnames, c, function_names)
#   pars = Pars(A, M, ΣZ, ΣU, gamma, samp_inds, F, dF, State, a_Z, A_Z, nu_Z, a_U, A_U, nu_U, ΣM, π)

#   return pars, model

# end # 2 dimensions missing data

# function create_pars(Y::Array{Float64,3},
#   SpaceStep::Vector,
#   TimeStep::Vector,
#   νS::Int,
#   νT::Int,
#   bufferSpace::Int,
#   bufferTime::Int,
#   batchSpace::Int,
#   batchTime::Int,
#   learning_rate::Float64,
#   c::Float64,
#   Λ::Function,
#   ∇Λ::Function,
#   ΛSpaceNames::Vector{String},
#   ΛTimeNames::Vector{String};
#   response=nothing, responsenames=nothing, degree=4, orderTime=1, orderSpace=3, latent_dim=size(Y, 3), covariates=nothing)

#   S = size(Y, 1)
#   T = size(Y, 2)
#   L = size(Y, 3)
#   N = latent_dim # need to update

#   Z = unfold3(Y)
#   # Δt = TimeStep[2] - TimeStep[1]
#   # ΔS = round(SpaceStep[2] - SpaceStep[1], digits=4)

#   # find inner index values
#   Nx = length(SpaceStep[1])
#   Nt = length(TimeStep)
#   space_inds = 1:Nx
#   time_inds = 1:Nt
#   ind_comb = reshape([[x, y] for x = space_inds, y = time_inds], length(space_inds) * length(time_inds))
#   ind_comb = permutedims(hcat(ind_comb...))
#   ind_comb = hcat(ind_comb, 1:(Nx*Nt))
#   i1 = ind_comb[:, 1] .∉ [cat(1:1:bufferSpace, (Nx-bufferSpace+1):1:Nx, dims=1)]
#   i2 = ind_comb[:, 2] .∉ [cat(1:1:bufferTime, (Nt-bufferTime+1):1:Nt, dims=1)]
#   inner_inds = ind_comb[i1.&i2, 3]
#   space_inds = ind_comb[:, 1]
#   time_inds = ind_comb[:, 2]
#   all_inds = ind_comb[:, 3]

#   if covariates === nothing
#     X = nothing
#   else
#     X = reduce(vcat, [reshape(covs[:, :, i], 1, :) for i in 1:size(covs, 3)])
#   end

#   # get missing values - no missing here
#   H = ones(L, S * T)
#   H = [make_H(L, N, H[:, i]) for i in 1:(S*T)]

#   # basis functions
#   # Basis = derivatives_dictionary(orderSpace, orderTime, νS, νT, ΔS, Δt, SpaceStep, TimeStep, degree=degree)
#   Basis = derivatives_dictionary(orderSpace, orderTime, νS, νT, SpaceStep, TimeStep, degree=degree)

#   Ψ = Basis.SpaceDerivative[Basis.SpaceNames[1]]
#   Φ = Basis.TimeDerivative[Basis.TimeNames[1]]
#   Θ = Matrix(Diagonal(ones(N)))
#   A = copy(inv(transpose(Θ) * Θ) * transpose(Θ) * Z * (Φ ⊗ Ψ) * inv(transpose(Φ ⊗ Ψ) * (Φ ⊗ Ψ)))
#   State = construct_state(A, Basis, Θ)

#   # parameters
#   samp_space, samp_time, samp_inds = sample_inds(SpaceStep, TimeStep, batchSpace, batchTime, bufferSpace, bufferTime)
#   # F = create_library(Λ, Λnames, State, X)
#   # dF = create_library_derivative(Λ, Λnames, State, samp_inds, X)
#   F = create_Λ(Λ, A, Basis, ΛSpaceNames, ΛTimeNames, X)
#   dF = create_∇Λ(∇Λ, A, Basis, samp_space, samp_time, samp_inds, ΛSpaceNames, ΛTimeNames, X)
  
#   D = size(F, 1)
#   M = zeros(N, D)
#   ΣZ = Diagonal(ones(N))
#   ΣU = Diagonal(ones(N))
#   # hyperpriors
#   a_Z = ones(N)
#   A_Z = fill(1e6, N)
#   nu_Z = 2

#   a_U = ones(N)
#   A_U = fill(1e6, N)
#   nu_U = 2
#   # SSVS parameters
#   gamma = ones(N, D)
#   ΣM = inv(c * Diagonal(ones(D)))
#   π = 0.5 * ones(N)

#   # determine response
#   if response === nothing
#     responseNames = ["ΔU_t"]
#     function responseFunction(ΔU)
#       ΔU_t = ΔU[1]
#       return ΔU_t
#     end
#   else
#     responseFunction = response
#     responseNames = responsenames
#   end
#   # else is user input for response function and response names

#   Ψnames = (response === nothing ? ["Psi"] : "Psi_" .* reduce(hcat, split.(reduce(hcat, split.(responsenames, "t"))[1, :], "_"))[2, :])
  
  
#   Ψtmp = [Basis.SpaceDerivative[ΛSpaceNames[i]] for i in 1:length(ΛSpaceNames)]
#   Φtmp = [Basis.TimeDerivative[ΛTimeNames[i]] for i in 1:length(ΛTimeNames)]
#   if X === nothing
#     function_names = replace((@code_string Λ(A, Ψtmp, Φtmp)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
#   else
#     function_names = replace((@code_string Λ(A, Ψtmp, Φtmp, X)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp, X)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
#   end

#   # store model and parameter values
#   # model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, Z, H, Basis, Θ, Λ, dΛ, Λnames, response)
#   model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, space_inds, time_inds, all_inds, Z, H, X, Basis, Θ, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, responseFunction, responseNames, Ψnames, c, function_names)
#   pars = Pars(A, M, ΣZ, ΣU, gamma, samp_inds, F, dF, State, a_Z, A_Z, nu_Z, a_U, A_U, nu_U, ΣM, π)

#   return pars, model

# end # 1 dimension

# function create_pars(Y::Array{Union{Missing,Float64},3},
#   SpaceStep::Vector,
#   TimeStep::Vector,
#   νS::Int,
#   νT::Int,
#   bufferSpace::Int,
#   bufferTime::Int,
#   batchSpace::Int,
#   batchTime::Int,
#   learning_rate::Float64,
#   c::Float64,
#   Λ::Function,
#   ∇Λ::Function,
#   ΛSpaceNames::Vector{String},
#   ΛTimeNames::Vector{String};
#   response=nothing, responsenames=nothing, degree=4, orderTime=1, orderSpace=3, latent_dim=size(Y, 3), covariates=nothing)

#   S = size(Y, 1)
#   T = size(Y, 2)
#   L = size(Y, 3)
#   N = latent_dim # need to update

#   na_inds = getindex.(findall(ismissing, unfold3(Y)), 2)
#   Z = unfold3(collect(Missings.replace(Y, 0)))
#   # Δt = TimeStep[2] - TimeStep[1]
#   # ΔS = round(SpaceStep[2] - SpaceStep[1], digits=4)

#   # find inner index values
#   Nx = length(SpaceStep[1])
#   Nt = length(TimeStep)
#   space_inds = 1:Nx
#   time_inds = 1:Nt
#   ind_comb = reshape([[x, y] for x = space_inds, y = time_inds], length(space_inds) * length(time_inds))
#   ind_comb = permutedims(hcat(ind_comb...))
#   ind_comb = hcat(ind_comb, 1:(Nx*Nt))
#   i1 = ind_comb[:, 1] .∉ [cat(1:1:bufferSpace, (Nx-bufferSpace+1):1:Nx, dims=1)]
#   i2 = ind_comb[:, 2] .∉ [cat(1:1:bufferTime, (Nt-bufferTime+1):1:Nt, dims=1)]
#   # i3 = ind_comb[:, 3] .∉ [cat(na_inds, dims=1)]
#   # inner_inds = ind_comb[i1.&i2.&i3, 3]
#   inner_inds = ind_comb[i1.&i2, 3]
#   space_inds = ind_comb[:, 1]
#   time_inds = ind_comb[:, 2]
#   all_inds = ind_comb[:, 3]

#   if covariates === nothing
#     X = nothing
#   else
#     X = reduce(vcat, [reshape(covs[:, :, i], 1, :) for i in 1:size(covs, 3)])
#   end

#   # get missing values
#   # na_inds = getindex.(findall(ismissing, unfold3(Y)), 2)
#   H = ones(L, S * T)
#   if (size(na_inds)[1] != 0)
#     H[na_inds] .= 0
#   end
#   H = [make_H(L, N, H[:, i]) for i in 1:(S*T)]


#   # basis functions
#   # Basis = derivatives_dictionary(orderSpace, orderTime, νS, νT, ΔS, Δt, SpaceStep, TimeStep, degree=degree)
#   Basis = derivatives_dictionary(orderSpace, orderTime, νS, νT, SpaceStep, TimeStep, degree=degree)

#   Z_start = Z
#   Z_start[na_inds] = (Z[na_inds .+ 1] .+ Z[na_inds .- 1])./2

#   Ψ = Basis.SpaceDerivative[Basis.SpaceNames[1]]
#   Φ = Basis.TimeDerivative[Basis.TimeNames[1]]
#   Θ = Matrix(Diagonal(ones(N)))
#   A = copy(inv(transpose(Θ) * Θ) * transpose(Θ) * Z_start * (Φ ⊗ Ψ) * inv(transpose(Φ ⊗ Ψ) * (Φ ⊗ Ψ)))
#   State = construct_state(A, Basis, Θ)

#   # parameters
#   samp_space, samp_time, samp_inds = sample_inds(SpaceStep, TimeStep, batchSpace, batchTime, bufferSpace, bufferTime)
#   # F = create_library(Λ, Λnames, State, X)
#   # dF = create_library_derivative(Λ, Λnames, State, samp_inds, X)
#   F = create_Λ(Λ, A, Basis, ΛSpaceNames, ΛTimeNames, X)
#   dF = create_∇Λ(∇Λ, A, Basis, samp_space, samp_time, samp_inds, ΛSpaceNames, ΛTimeNames, X)
  
#   D = size(F, 1)
#   M = zeros(N, D)
#   ΣZ = Diagonal(ones(N))
#   ΣU = Diagonal(ones(N))
#   # hyperpriors
#   a_Z = ones(N)
#   A_Z = fill(1e6, N)
#   nu_Z = 2

#   a_U = ones(N)
#   A_U = fill(1e6, N)
#   nu_U = 2
#   # SSVS parameters
#   gamma = ones(N, D)
#   ΣM = inv(c * Diagonal(ones(D)))
#   π = 0.5 * ones(N)

#   # determine response
#   if response === nothing
#     responseNames = ["ΔU_t"]
#     function responseFunction(ΔU)
#       ΔU_t = ΔU[1]
#       return ΔU_t
#     end
#   else
#     responseFunction = response
#     responseNames = responsenames
#   end
#   # else is user input for response function and response names

#   Ψnames = (response === nothing ? ["Psi"] : "Psi_" .* reduce(hcat, split.(reduce(hcat, split.(responsenames, "t"))[1, :], "_"))[2, :])

#   Ψtmp = [Basis.SpaceDerivative[ΛSpaceNames[i]] for i in 1:length(ΛSpaceNames)]
#   Φtmp = [Basis.TimeDerivative[ΛTimeNames[i]] for i in 1:length(ΛTimeNames)]
#   if X === nothing
#     function_names = replace((@code_string Λ(A, Ψtmp, Φtmp)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
#   else
#     function_names = replace((@code_string Λ(A, Ψtmp, Φtmp, X)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp, X)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
#   end

#   # store model and parameter values
#   # model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, Z, H, Basis, Θ, Λ, dΛ, Λnames, response)
#   model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, space_inds, time_inds, all_inds, Z, H, X, Basis, Θ, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, responseFunction, responseNames, Ψnames, c, function_names)
#   pars = Pars(A, M, ΣZ, ΣU, gamma, samp_inds, F, dF, State, a_Z, A_Z, nu_Z, a_U, A_U, nu_U, ΣM, π)

#   return pars, model

# end # 1 dimension missing data


"""
    Posterior Struct

Initiates a structure of class Posterior to hold posterior samples

# Arguments

"""
mutable struct Posterior

  M
  gamma
  ΣZ
  ΣU
  A
  π

end


"""
    make_H(M, N, one_inds)

Used to construct the mapping matrix H

# Arguments
- `M::`: something

"""
function make_H(L, N, one_inds)
  H_tmp = zeros(L, N)
  H_tmp[diagind(H_tmp)] .= one_inds
  return H_tmp
end


"""
    update_M!(pars)

Used within DEtection() function with a spike-and-slab prior.
Updates 𝐌 and 𝚺ᵤ from 

``𝚯 𝐀 (𝛗_{t^{(i)}}(t) ⊗ g(𝛙(𝐬)))' = 𝐌 𝐟(⋅) + 𝛈(𝐬, t)`` 

where ``𝛈(𝐬, t) ∼ N_N(0, 𝚺ᵤ)``and ``𝚺ᵤ = diag(σ²ᵤ1, ..., σ²ᵤN)``
"""
function update_M!(pars, model)

  # samp_space, samp_time, samp_inds = sample_inds(model.SpaceStep, model.TimeStep, model.batchSpace, model.batchTime, model.bufferSpace, model.bufferTime)
  samp_inds = model.inner_inds
  # samp_inds = pars.sample_inds

  A0 = pars.ΣM
  sN = (length(samp_inds)-1) / 2
  N = length(samp_inds)
  # sN = (length(samp_inds)-model.c) / 2
  # sN = (length(samp_inds)-(N-50)) / 2
  # sN = (length(samp_inds)-(N-5000)) / 2

  # N = 1
  # N = sqrt(length(samp_inds))
  # N = log(length(samp_inds))
  Uprime = model.responseFunction([pars.State.State[model.responseNames[i]] for i in 1:length(model.responseNames)])

  # update ΣU and M
  for i in 1:model.N

    Uprimei = vec(Uprime[i, samp_inds])

    # Xδ = pars.F[pars.gamma[i, :] .== 1, samp_inds]
    # Aδ = inv(Xδ * Xδ' + A0[pars.gamma[i, :] .== 1, pars.gamma[i, :] .== 1]) # independent
    # # Aδ = (N/(N+1)) * inv(Xδ * Xδ') # g-prior
    # # Aδ = inv(Xδ * Xδ') # f-slab prior

    # aδ = Aδ * Xδ * Uprimei
    # SN = (1 / 2) * (Uprimei'Uprimei - aδ' * inv(Aδ) * aδ)

    # pars.ΣU[i, i] = rand(InverseGamma(sN, SN))

    # pars.M[i, vec(pars.gamma[i, :] .!= 1)] .= 0
    # # pars.M[i, vec(pars.gamma[i, :] .== 1)] = rand(MvNormal(aδ'[pars.gamma[i, :].==1], pars.ΣU[i, i] .* Hermitian(Aδ[vec(pars.gamma[i, :] .== 1), vec(pars.gamma[i, :] .== 1)])))
    # pars.M[i, vec(pars.gamma[i, :] .== 1)] = rand(MvNormal(aδ, pars.ΣU[i, i] .* Hermitian(Aδ)))


    # g prior
    Xδ = pars.F[pars.gamma[i, :] .== 1, samp_inds]
    # a0 = (inv(pars.F[:, samp_inds] * pars.F[:, samp_inds]') * pars.F[:, samp_inds] * Uprimei)[pars.gamma[i, :] .== 1]
    a0 = inv(Xδ * Xδ') * Xδ * Uprimei
    A0 = N * inv(Xδ * Xδ')
    A0inv = (1/N) * Xδ * Xδ'

    Aδ = (N/(N+1)) * inv(Xδ * Xδ')
    # aδ = (Aδ * Xδ + A0inv * a0) * Uprimei
    aδ = Aδ * (Xδ * Uprimei + A0inv * a0)
    SN = (1 / 2) * (Uprimei'Uprimei + a0' * A0inv * a0 - aδ' * inv(Aδ) * aδ)

    pars.ΣU[i, i] = rand(InverseGamma(sN, SN))



    pars.M[i, vec(pars.gamma[i, :] .!= 1)] .= 0
    pars.M[i, vec(pars.gamma[i, :] .== 1)] = rand(MvNormal(Aδ * (Xδ * Uprimei + A0inv * a0), pars.ΣU[i, i] .* Hermitian(Aδ)))

  end

  return pars
end


"""
    update_gamma!(pars)

Used within DEtection() function with a spike-and-slab prior.
Updates gamma and pi from the spike-and-slab prior from

``p(𝐌⎹ 𝛄) = p-slab(𝐌 ᵧ)∏ⱼ p-spike(Mⱼ)``

where ``p(γⱼ = 1⎹ π) = π``, ``π ∼ ℬ(a, b)``, ``p-slab(𝐌 ᵧ) = N(𝐌 ᵧ⎹ 𝐚ᵧ, 𝚺ᵤ𝐀ᵧ)``, and ``𝐀ᵧ = c𝐈``.

"""
# function update_gamma!(pars, model)

#   # samp_space, samp_time, samp_inds = sample_inds(model.SpaceStep, model.TimeStep, 40, 40, model.bufferSpace, model.bufferTime)
#   samp_inds = model.inner_inds
#   # samp_inds = pars.sample_inds

#   # update spike δ
#   A0 = pars.ΣM
#   δtmp = pars.gamma
#   sN = (length(samp_inds)) / 2
#   N = length(samp_inds)

#   # F = (pars.F[:, samp_inds] .- mean(pars.F[:, samp_inds], dims = 2)) ./ std(pars.F[:, samp_inds], dims = 2)
#   F = pars.F[:, samp_inds]

#   for i in 1:model.N

#     Uprime = vec(model.responseFunction([pars.State.State[model.responseNames[i]] for i in 1:length(model.responseNames)])[i, samp_inds])
#     # Uprime = Uprime .- mean(Uprime)

#     δorder = sample(1:model.D, model.D, replace=false)
#     for j in δorder

#       δtmp1 = copy(δtmp[i, :])
#       δtmp0 = copy(δtmp[i, :])
#       δtmp1[j] = 1
#       δtmp0[j] = 0
#       # Xδ0 = δtmp0 .* pars.F[:, samp_inds]
#       # Xδ1 = δtmp1 .* pars.F[:, samp_inds]
#       Xδ0 = F[δtmp0 .== 1, :]
#       Xδ1 = F[δtmp1 .== 1, :]
#       # Xδ0 = pars.F[δtmp0 .== 1, samp_inds]
#       # Xδ1 = pars.F[δtmp1 .== 1, samp_inds]

#       # independent prior
#       Aδ0 = inv(Xδ0 * Xδ0' + A0[δtmp0 .== 1, δtmp0 .== 1])
#       Aδ1 = inv(Xδ1 * Xδ1' + A0[δtmp1 .== 1, δtmp1 .== 1])
#       # Aδ0 = inv(Xδ0 * Xδ0' + A0)
#       # Aδ1 = inv(Xδ1 * Xδ1' + A0)
#       aδ0 = Aδ0' * Xδ0 * Uprime 
#       aδ1 = Aδ1' * Xδ1 * Uprime
#       SN0 = (1 / 2) * (Uprime' * Uprime - aδ0' * inv(Aδ0) * aδ0)
#       SN1 = (1 / 2) * (Uprime' * Uprime - aδ1' * inv(Aδ1) * aδ1)

#       # sse_ratio = ((SN1 / SN0) > 0.95 ? 1 : (SN1 / SN0))
#       # R = (model.c*det(Aδ0)/det(Aδ1))^(1/2) * (sse_ratio)^(sN)
#       R = (model.c*det(Aδ0)/det(Aδ1))^(1/2) * (SN1 / SN0)^(sN)
#       p = 1 / (1 + ((1 - pars.π[i]) / pars.π[i]) * R)

      
#       # R = det(Aδ0*inv(Aδ1))^(1/2) * (SN1 / SN0)^(sN)
#       # R = det(Aδ0*inv(Aδ1))^(1/2) * exp((-1/pars.ΣU[i,i]) * (SN0 - SN1))
#       # R = ((det(Aδ0)^(1 / 2)) / (det(Aδ1)^(1 / 2))) * (SN1 / SN0)^(sN)
#       # R = isnan(R) ? 0 : R

#       # R = (model.c*det(Aδ0)/det(Aδ1))^(1/2) * exp(-(1/pars.ΣU[i,i]) * (SN0 - SN1))
#       # p = 1 / (1 + ((1 - pars.π[i]) / pars.π[i]) * R)



#       # ha = -(1/pars.ΣU[i,i])*(aδ1' * inv(Aδ1) * aδ1 - aδ0' * inv(Aδ0) * aδ0) + log(model.c * det(Aδ0)/det(Aδ1))
#       # p = 1/(1 + exp(ha/2)*((1-pars.π[i])/pars.π[i]))

#       δtmp[i, j] = rand(Binomial(1, p))

#     end
#   end

#   pars.gamma = δtmp

#   # update π
#   for i in 1:model.N
#     # pars.π[i] = rand(Beta(1 + sum(pars.gamma[i, :]), 1 + sum(1 .- pars.gamma[i, :])))
#     pars.π[i] = rand(Beta(1 + sum(pars.gamma[i, :]) + 2, 1 + sum(1 .- pars.gamma[i, :]) + 18))
#   end

#   return pars
# end

# use this one
function update_gamma!(pars, model)

  # samp_space, samp_time, samp_inds = sample_inds(model.SpaceStep, model.TimeStep, 32, 32, model.bufferSpace, model.bufferTime) # baro
  # samp_space, samp_time, samp_inds = sample_inds(model.SpaceStep, model.TimeStep, 80, 80, model.bufferSpace, model.bufferTime) # heat
  # samp_space, samp_time, samp_inds = sample_inds(model.SpaceStep, model.TimeStep, 40, 40, model.bufferSpace, model.bufferTime) # reaction
  # samp_space, samp_time, samp_inds = sample_inds(model.SpaceStep, model.TimeStep, 10, 10, model.bufferSpace, model.bufferTime) # burger
  samp_inds = StatsBase.sample(model.inner_inds, Int(length(model.inner_inds) - model.c), replace=false)

  # update spike δ
  # A0 = pars.ΣM
  δtmp = pars.gamma
  sN = (length(samp_inds)) / 2
  N = length(samp_inds)
  g = N
  F = pars.F[:, samp_inds]


  for i in 1:model.N

    Uprime = vec(model.responseFunction([pars.State.State[model.responseNames[i]] for i in 1:length(model.responseNames)])[i, samp_inds])

    δorder = sample(1:(model.D), model.D, replace=false)
    for j in δorder

      δtmp1 = copy(δtmp[i, :])
      δtmp0 = copy(δtmp[i, :])
      δtmp1[j] = 1
      δtmp0[j] = 0
      Xδ0 = F[δtmp0 .== 1, :]
      Xδ1 = F[δtmp1 .== 1, :]

      # g-prior
      A00 = g * inv(Xδ0*Xδ0')
      A01 = g * inv(Xδ1*Xδ1')
      Aδ0 = (g/(g+1)) * inv(Xδ0 * Xδ0')
      Aδ1 = (g/(g+1)) * inv(Xδ1 * Xδ1')


      aδ0 = Aδ0' * Xδ0 * Uprime 
      aδ1 = Aδ1' * Xδ1 * Uprime
      SN0 = (1 / 2) * (Uprime' * Uprime - aδ0' * inv(Aδ0) * aδ0)
      SN1 = (1 / 2) * (Uprime' * Uprime - aδ1' * inv(Aδ1) * aδ1)

      R = (g+1)^(1/2) * (SN1 / SN0)^(sN)
      p = 1 / (1 + ((1 - pars.π[i]) / pars.π[i]) * R)

      δtmp[i, j] = rand(Binomial(1, p))

    end
  end

  pars.gamma = δtmp

  # update π
  for i in 1:model.N
    pars.π[i] = rand(Beta(1 + sum(pars.gamma[i, :]), 1 + sum(1 .- pars.gamma[i, :])))
    # pars.π[i] = rand(Beta(1 + sum(pars.gamma[i, :]) + 2, 1 + sum(1 .- pars.gamma[i, :]) + 18))
    # pars.π[i] = 0.1
  end

  return pars
end

"""
  update_ΣZ!(pars)

Used within DEtection() function with a spike-and-slab prior.
Updates 𝚺z from 

``𝐳(𝐬, t) = ℋ(𝐬,t) 𝚯 𝐀 (𝛗phi_{t^{(0)}}(t) ⊗ 𝛙(𝐬))' + 𝛜(𝐬, t)``

where ``𝛜(𝐬, t) ∼ N_N(0, 𝚺z)``, ``𝚺z = diag(σ²z1, ..., σ²zm)``, and ``σ²z ∼ Half-t()``.

"""
function update_ΣZ!(pars, model)

  samp_inds = model.inner_inds
  # samp_inds = pars.sample_inds

  sse_tmp = [model.Z[:, j] - model.H[j] * pars.State.State[pars.State.StateNames[1]][:, j] for j in samp_inds]
  sse_tmp = reduce(hcat, sse_tmp)
  for n in 1:model.N
    a_hat = (length(samp_inds) + pars.nu_Z) / 2
    b_hat = (0.5*(sse_tmp[n, :]'*sse_tmp[n, :])+pars.nu_Z[1]/pars.a_Z[n])[1]
    pars.ΣZ[n,n] = rand(InverseGamma(a_hat, b_hat))
  end

  for n in 1:model.N
    a_hat = (pars.nu_Z + 1) / 2
    b_hat = (pars.nu_Z / pars.ΣZ[n]) + (1 / pars.A_Z[n]^2)
    pars.a_Z[n] = rand(InverseGamma(a_hat, b_hat))
  end

  return pars
end


"""
  ΔL(z, H, ψ, gψ, ϕ, ϕ_t, Θ, ΣZinv, ΣUinv, A, M, fcurr, fprime)

Used within DEtection() function.
Calculates the gradient of the log likelihood.
"""
function ΔL(z, H, ψ, gψ, ϕ, ϕ_t, Θ, ΣZinv, ΣUinv, A, M, fcurr, fprime::Matrix{Float64})

  B0 = ϕ ⊗ ψ
  Bt = ϕ_t ⊗ gψ

  # ΔL = -Θ * H' * ΣZinv * z * B0' +
  #      Θ * H' * ΣZinv * H * Θ * A * B0 * B0' +
  #      Θ * ΣUinv * Θ * A * Bt * Bt' -
  #      Θ * ΣUinv * M * fcurr * Bt' -
  #      Θ * fprime' * M' * ΣUinv * Θ * A * Bt * B0' +
  #      Θ * fprime * M' * ΣUinv * M * fcurr * B0'

  ΔL = -Θ * H' * ΣZinv * z * B0' +
        Θ * H' * ΣZinv * H * Θ * A * B0 * B0' +
        Θ * ΣUinv * Θ * A * Bt * Bt' -
        Θ * ΣUinv * M * fcurr * Bt' -
        Bt' * A' * Θ' * ΣUinv * M * fprime +
        (fcurr' * M' * ΣUinv * M) * fprime

  return ΔL

end

function ΔL(z, H, ψ, gψ, ϕ, ϕ_t, Θ, ΣZinv, ΣUinv, A, M, fcurr, fprime::Array{Float64, 3})

  B0 = ϕ ⊗ ψ
  Bt = ϕ_t ⊗ gψ

  ΔL = -Θ * H' * ΣZinv * z * B0' +
        Θ * H' * ΣZinv * H * Θ * A * B0 * B0' +
        Θ * ΣUinv * Θ * A * Bt * Bt' -
        Θ * ΣUinv * M * fcurr * Bt' -
        reduce(vcat, [Bt' * A' * Θ' * ΣUinv * M * fprime[:,:,i] - (fcurr' * M' * ΣUinv * M) * fprime[:,:,i] for i in 1:size(fprime,3)])

  return ΔL

end



# z = z[:, i]
# H = h[i]
# ψ = ψ[i, :]
# gψ = gψ[i,:]
# ϕ = ϕ[i, :]
# ϕ_t = ϕ_t[i, :]
# fcurr = fcurr[:, i]
# fprime = fprime[i]


# ΣUinv * fcurr' * inv(fcurr * ΣUinv * fcurr' + pars.ΣM)' * fprime
# Θ * ΣUinv * Θ * A * Bt * Bt'



"""
    update_A!(pars)

Used within DEtection() function with a spike-and-slab prior.
Updates 𝐀 with the elastic net prior (Li 2010) from

``𝐳(𝐬, t) = ℋ(𝐬,t) 𝚯 𝐀 (𝛗phi_{t^{(0)}}(t) ⊗ 𝛙(𝐬))' + 𝛜(𝐬, t)``

``𝚯 𝐀 (𝛗_{t^{(i)}}(t) ⊗ g(𝛙(𝐬)))' = 𝐌 𝐟(⋅) + 𝛈(𝐬, t)`` 

``p(𝐀) ∝ exp{-λ₁‖𝐀‖₁ - λ₂‖𝐀‖₂²}``

using Stochastic Gradient Desent with a Constant Learning Rate (Mandt 2016).
λ₁ and λ₂ are set to 1/100.

"""
function update_A!(pars, model)

  samp_space, samp_time, samp_inds = sample_inds(model.SpaceStep, model.TimeStep, model.batchSpace, model.batchTime, model.bufferSpace, model.bufferTime)
  # samp_space, samp_time, samp_inds = sample_inds(model.SpaceStep, model.TimeStep, model.batchSpace, model.batchTime)

  z = model.Z[:, samp_inds] # if an error occurs with 1D vs 2D space, it is here
  h = model.H[samp_inds]
  ψ = model.Basis.SpaceDerivative["Psi"][samp_space, :]
  gψ = model.responseFunction([model.Basis.SpaceDerivative[model.Ψnames[i]] for i in 1:length(model.Ψnames)])[samp_space, :]
  ϕ = model.Basis.TimeDerivative[model.Basis.TimeNames[1]][samp_time, :]
  ϕ_t = model.Basis.TimeDerivative[model.Basis.TimeNames[end]][samp_time, :]
  Θ = model.Θ
  ΣZinv = inv(pars.ΣZ)
  ΣUinv = inv(pars.ΣU)
  # ΣZinv = Matrix(Diagonal(ones(model.N)))
  # ΣUinv = Matrix(Diagonal(ones(model.N)))
  A = pars.A
  M = pars.M

  # fcurr = pars.gamma' .* pars.F[:, samp_inds]
  # fprime = create_∇Λ(model.∇Λ, A, model.Basis, samp_space, samp_time, samp_inds, model.ΛSpaceNames, model.ΛTimeNames, model.X)
  # dq = [ΔL(z[:, i], h[i], ψ[i, :], gψ[i, :], ϕ[i, :], ϕ_t[i, :], Θ, ΣZinv, ΣUinv, A, M, fcurr[:, i], pars.gamma' .* fprime[i]) for i in 1:length(samp_inds)]

  fcurr = pars.F[:, samp_inds]
  fprime = create_∇Λ(model.∇Λ, A, model.Basis, samp_space, samp_time, samp_inds, model.ΛSpaceNames, model.ΛTimeNames, model.X)
  dq = [ΔL(z[:, i], h[i], ψ[i, :], gψ[i, :], ϕ[i, :], ϕ_t[i, :], Θ, ΣZinv, ΣUinv, A, M, fcurr[:, i], fprime[i]) for i in 1:length(samp_inds)]

  # dq1 = dq[1]
  # dq = mean(dq)

  # 1/tr((dq1 - dq)' * (dq1 - dq))

  scale_el = 1 / 100
  # scale_el = 10
  # dq = (1 / (model.batchSpace * model.batchTime)) * (dq + (1 / (model.S * model.T)) * (scale_el * sign.(A) + 2 * scale_el * A))
  dq = mean(dq) + (1 / (model.S * model.T)) * (scale_el * sign.(A) + 2 * scale_el * A)
  # dq = mean(dq) + (1 / (model.S * model.T)) * (scale_el * A ./ sum(sign.(A)))
  # dq = mean(dq) + (1 / length(model.inner_inds)) * (scale_el * A ./ sum(sign.(A)) + 2 * scale_el * A ./ sum(A .^2))
  # dq = mean(dq) + (scale_el * A ./ sum(sign.(A)) + 2 * scale_el * A ./ sum(A .^2))


  # pars.A = A - model.learning_rate * dq
  pars.A = A - dq .* model.learning_rate

  # save updated state information
  pars.sample_inds = samp_inds
  pars.State = construct_state(pars.A, model.Basis, model.Θ)
  pars.F = create_Λ(model.Λ, pars.A, model.Basis, model.ΛSpaceNames, model.ΛTimeNames, model.X)
  pars.Fprime = fprime

  return pars
end


"""
    DEtection(Y, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, v0, v1, Λ, Λnames)

DEtection sampler function. Can accept missing values in the input data argument.
Returns the model, parameters, and posterior values.
The model are the model settings.
The parameters are the final value of the parameters in from the sampler.
The posterior are the saved posterior values.

Consider the function ``U_t = M(U, U_x, U_y, U_{xy}, ...)``.
DEtection() is used to determine ``M`` given a library of potential values, `Λ(U, ∂U)`, to search over. 
Within the function, partial derivatives are denoted as ``ΔU_t, ΔU_x, ΔU_y, ΔU_xx, , ΔU_xy, ...``, so a potential function could be

```jldoctest 
function Λ(U, ∂U){
  u = U[1]
  u_x = ∂U[1]
  u_y = ∂U[2]

  return [u, u_x, u_y, u*u_x, u*u_y]
}
```

To make the function identify the correct partial derivatives, the argument that are passed into `∂U`, `Λnames`, are required.
For the example above, `Λnames = [ΔU_x, ΔU_y]` because the function `Λ` uses the partial derivatives ``U_x`` and ``U_t``.

If you want to add covariates to the function, a possible `Λ(U, ∂U, X)` is

```jldoctest 
function Λ(U, ∂U, X){
  u = U[1]
  u_x = ∂U[1]
  u_y = ∂U[2]
  x1 = X[1]
  x2 = X[2]

  return [u, u_x, u_y, u*u_x, u*u_y, x1, x2, u*x1, u_x*x2]
}
```

where `Λnames = [ΔU_x, ΔU_y]`.

# Required Arguments (in order)
- `Y`: Input data. Needs to be Array{Float64, 3} or Array{Union{Missing, Float64}, 3} where the dimensions are Space, Time, Components
- `SpaceStep`: of type StepRangeLen (range function). For example, range(-1, 1, step = 0.1) for 1 dimension and [range(-1, 1, step = 0.1), range(-1, 1, step = 0.1)] for 2.
- `TimeStep`: of type StepRangeLen (range function). For example, range(-1, 1, step = 0.1).
- `νS::Int` or `νS::Vector{Int}`: Number of spatial basis functions. 
- `νT::Int`: Number of temporal basis functions
- `bufferSpace::Int` or `bufferSpace::Vector{Int}`: Buffer in space.
- `bufferTime::Int`: 
- `batchSpace::Int`: 
- `batchTime::Int`: 
- `learning_rate::Float64`: 
- `v0::Float64`: 
- `v1::Float64`: 
- `Λ::Function`: 
- `Λnames::Vector{String}`: 

# Optional Arguments
- `response` = "ΔU_t": Order of the temporal derivative (default first order). Use "ΔU_tt" for second and so on.
- `degree` = 4: Degree of the B-spline. Must be at least one order higher than the highest order partial derivative.
- `orderTime` = 1: Order of the highest order temporal derivative (default ∂U_t)
- `orderSpace` = 3: Order of the highest order spatial derivative (default ∂U_xxx and ∂U_yyy)
- `latent_dim` = size(Y, 3): Dimension of the latent space. Default is same as data dimension.
- `covariates` = nothing: Additional covariates.
- `nits` = 2000: Number of samples for the Gibbs sampler.
- `burnin` = nits / 2: Number of samples to discard as burnin (default is half of `nits`).
- `learning_rate_end` = learning_rate: End learning rate (default is same as initial learning rate).


# Examples

```jldoctest 
Y = Array{Float64}(rand, 256, 101, 1) # space x time x component, vectorized in space
SpaceStep = [range(-1, 1, step = 0.01), range(-1, 1, step = 0.01)]
TimeStep = range(-1, 1, step = 0.01)
νS = [10, 10] # number of space basis functions
νT = 10 # number of time basis functions
bufferSpace = [2, 2]
bufferTime = 2
batchSpace = 10
batchTime = 10
learning_rate = 1e-2
v0 = 1e-4
v1 = 1e4

Λnames = ["U", "ΔU_x", "ΔU_xx", "ΔU_xxx"]

function Λ(U, ∂U)

  ΔU_x = ∂U[1]
  ΔU_xx = ∂U[2]
  ΔU_xxx = ∂U[3]

  u = U[1]
  u_x = ΔU_x[1]
  u_xx = ΔU_xx[1]
  u_xxx = ΔU_xxx[1]

  return [u, u^2, u^3,
    u_x, u * u_x, u^2 * u_x, u^3 * u_x,
    u_xx, u * u_xx, u^2 * u_xx, u^3 * u_xx,
    u_xxx, u * u_xxx, u^2 * u_xxx, u^3 * u_xxx]

end

# not run
# model, pars, posterior = DEtection(Y, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, c, Λ, Λnames)
``` 
"""
function DEtection(Y::Array{Float64,3},
  SpaceStep::Vector,
  TimeStep::Vector,
  νS::Int,
  νT::Int,
  batchSpace::Int,
  batchTime::Int,
  learning_rate::Float64,
  beta::Float64,
  Λ::Function,
  ∇Λ::Function,
  ΛSpaceNames::Vector{String},
  ΛTimeNames::Vector{String};
  bufferSpace::Int = 1,
  bufferTime::Int = 1,
  response=nothing, responsenames=nothing, degree=4, orderTime=1, orderSpace=3, latent_dim=size(Y, 3), covariates=nothing, nits=2000, burnin=nits / 2, learning_rate_end=learning_rate)

  pars, model = create_pars(Y, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, beta, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, bufferSpace = bufferSpace, bufferTime = bufferTime, response=response, responsenames=responsenames, degree=degree, orderTime=orderTime, orderSpace=orderSpace, latent_dim=latent_dim, covariates=covariates)

  # decaying learning rate
  keep_samps = Int(nits - burnin)
  n_change = log10(learning_rate / learning_rate_end)
  if n_change == 0
    change_times = 0.0
  else
    change_times = floor.(collect(range(1, stop=burnin, length=Int((n_change + 1)))))[2:end]
  end

  M_post = Array{Float64}(undef, model.N, model.D, keep_samps)
  gamma_post = Array{Float64}(undef, model.N, model.D, keep_samps)
  ΣZ_post = Array{Float64}(undef, model.N, model.N, keep_samps)
  ΣU_post = Array{Float64}(undef, model.N, model.N, keep_samps)
  A_post = Array{Float64}(undef, model.N, size(pars.A)[2], keep_samps)
  π_post = Array{Float64}(undef, model.N, keep_samps)


  @showprogress 1 "Burnin..." for i in 1:burnin

    pars = update_gamma!(pars, model)
    pars = update_M!(pars, model)
    pars = update_ΣZ!(pars, model)
    pars = update_A!(pars, model)

    if i in change_times
      raise_pow = findall(change_times .== i)[1]
      model.learning_rate = round(learning_rate * 10.0^(-raise_pow), digits = raise_pow+6)
    end

  end

  @showprogress 1 "Sampling..." for i in 1:keep_samps

    pars = update_gamma!(pars, model)
    pars = update_M!(pars, model)
    pars = update_ΣZ!(pars, model)
    pars = update_A!(pars, model)

    M_post[:, :, i] = pars.M
    gamma_post[:, :, i] = pars.gamma
    ΣZ_post[:, :, i] = pars.ΣZ
    ΣU_post[:, :, i] = pars.ΣU
    A_post[:, :, i] = pars.A
    π_post[:, i] = pars.π

  end

  posterior = Posterior(M_post, gamma_post, ΣZ_post, ΣU_post, A_post, π_post)

  return model, pars, posterior

end # 1 spatial dimension DEtection

function DEtection(Y::Array{Float64,3},
  SpaceStep::Vector,
  TimeStep::Vector,
  νS::Vector{Int},
  νT::Int,
  batchSpace::Int,
  batchTime::Int,
  learning_rate::Vector{Float64},
  beta::Float64,
  Λ::Function,
  ∇Λ::Function,
  ΛSpaceNames::Vector{String},
  ΛTimeNames::Vector{String};
  bufferSpace::Vector{Int} = [1, 1],
  bufferTime::Int = 1,
  response=nothing, responsenames=nothing, degree=4, orderTime=1, orderSpace=3, latent_dim=size(Y, 3), covariates=nothing, nits=2000, burnin=nits / 2, learning_rate_end=learning_rate)

  if typeof(learning_rate) != Vector{Float64}
    learning_rate = learning_rate * ones(latent_dim)
  end

  pars, model = create_pars(Y, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, beta, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, bufferSpace = bufferSpace, bufferTime = bufferTime, response=response, responsenames=responsenames, degree=degree, orderTime=orderTime, orderSpace=orderSpace, latent_dim=latent_dim, covariates=covariates)
  # pars, model = create_pars(Y, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, c, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, response=response, responsenames=responsenames, degree=degree, orderTime=orderTime, orderSpace=orderSpace, latent_dim=latent_dim, covariates=covariates)

  # decaying learning rate
  keep_samps = Int(nits - burnin)
  n_change = log10.(learning_rate ./ learning_rate_end)
  # if n_change .== 0
  #   change_times = 0.0
  # else
  #   change_times = floor.(collect(range(1, stop=burnin, length=Int((n_change + 1)))))[2:end]
  # end
  change_times = 0

  M_post = Array{Float64}(undef, model.N, model.D, keep_samps)
  gamma_post = Array{Float64}(undef, model.N, model.D, keep_samps)
  ΣZ_post = Array{Float64}(undef, model.N, model.N, keep_samps)
  ΣU_post = Array{Float64}(undef, model.N, model.N, keep_samps)
  A_post = Array{Float64}(undef, model.N, size(pars.A)[2], keep_samps)
  π_post = Array{Float64}(undef, model.N, keep_samps)


  @showprogress 1 "Burnin..." for i in 1:burnin

    pars = update_gamma!(pars, model)
    pars = update_M!(pars, model)
    pars = update_ΣZ!(pars, model)
    pars = update_A!(pars, model)

    if i in change_times
      raise_pow = findall(change_times .== i)[1]
      model.learning_rate = round(learning_rate * 10.0^(-raise_pow), digits = raise_pow+6)
    end

  end

  @showprogress 1 "Sampling..." for i in 1:keep_samps

    pars = update_gamma!(pars, model)
    pars = update_M!(pars, model)
    pars = update_ΣZ!(pars, model)
    pars = update_A!(pars, model)

    M_post[:, :, i] = pars.M
    gamma_post[:, :, i] = pars.gamma
    ΣZ_post[:, :, i] = pars.ΣZ
    ΣU_post[:, :, i] = pars.ΣU
    A_post[:, :, i] = pars.A
    π_post[:, i] = pars.π

  end

  posterior = Posterior(M_post, gamma_post, ΣZ_post, ΣU_post, A_post, π_post)

  return model, pars, posterior

end # 2 spatial dimensions DEtection

function DEtection(Y::Array{Union{Missing,Float64},3},
  SpaceStep::Vector,
  TimeStep::Vector,
  νS::Int,
  νT::Int,
  batchSpace::Int,
  batchTime::Int,
  learning_rate::Float64,
  beta::Float64,
  Λ::Function,
  ∇Λ::Function,
  ΛSpaceNames::Vector{String},
  ΛTimeNames::Vector{String};
  bufferSpace::Int = 1,
  bufferTime::Int = 1,
  response=nothing, responsenames=nothing, degree=4, orderTime=1, orderSpace=3, latent_dim=size(Y, 3), covariates=nothing, nits=2000, burnin=nits / 2, learning_rate_end=learning_rate)

  pars, model = create_pars(Y, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, beta, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, bufferSpace = bufferSpace, bufferTime = bufferTime, response=response, responsenames=responsenames, degree=degree, orderTime=orderTime, orderSpace=orderSpace, latent_dim=latent_dim, covariates=covariates)
  # pars, model = create_pars(Y, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, c, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, response=response, responsenames=responsenames, degree=degree, orderTime=orderTime, orderSpace=orderSpace, latent_dim=latent_dim, covariates=covariates)

  # decaying learning rate
  keep_samps = Int(nits - burnin)
  n_change = log10(learning_rate / learning_rate_end)
  if n_change == 0
    change_times = 0.0
  else
    change_times = floor.(collect(range(1, stop=burnin, length=Int((n_change + 1)))))[2:end]
  end

  M_post = Array{Float64}(undef, model.N, model.D, keep_samps)
  gamma_post = Array{Float64}(undef, model.N, model.D, keep_samps)
  ΣZ_post = Array{Float64}(undef, model.N, model.N, keep_samps)
  ΣU_post = Array{Float64}(undef, model.N, model.N, keep_samps)
  A_post = Array{Float64}(undef, model.N, size(pars.A)[2], keep_samps)
  π_post = Array{Float64}(undef, model.N, keep_samps)


  @showprogress 1 "Burnin..." for i in 1:burnin

    pars = update_gamma!(pars, model)
    pars = update_M!(pars, model)
    pars = update_ΣZ!(pars, model)
    pars = update_A!(pars, model)

    if i in change_times
      raise_pow = findall(change_times .== i)[1]
      model.learning_rate = round(learning_rate * 10.0^(-raise_pow), digits = raise_pow+6)
    end

  end

  @showprogress 1 "Sampling..." for i in 1:keep_samps

    pars = update_gamma!(pars, model)
    pars = update_M!(pars, model)
    pars = update_ΣZ!(pars, model)
    pars = update_A!(pars, model)

    M_post[:, :, i] = pars.M
    gamma_post[:, :, i] = pars.gamma
    ΣZ_post[:, :, i] = pars.ΣZ
    ΣU_post[:, :, i] = pars.ΣU
    A_post[:, :, i] = pars.A
    π_post[:, i] = pars.π

  end

  posterior = Posterior(M_post, gamma_post, ΣZ_post, ΣU_post, A_post, π_post)

  return model, pars, posterior

end # 1 spatial dimension DEtection - missing values

function DEtection(Y::Array{Union{Missing,Float64},3},
  SpaceStep::Vector,
  TimeStep::Vector,
  νS::Vector{Int},
  νT::Int,
  batchSpace::Int,
  batchTime::Int,
  learning_rate::Vector{Float64},
  beta::Float64,
  Λ::Function,
  ∇Λ::Function,
  ΛSpaceNames::Vector{String},
  ΛTimeNames::Vector{String};
  bufferSpace::Vector{Int} = [1,1],
  bufferTime::Int = 1,
  response=nothing, responsenames=nothing, degree=4, orderTime=1, orderSpace=3, latent_dim=size(Y, 3), covariates=nothing, nits=2000, burnin=nits / 2, learning_rate_end=learning_rate)

  if typeof(learning_rate) != Vector{Float64}
    learning_rate = learning_rate * ones(latent_dim)
  end

  pars, model = create_pars(Y, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, beta, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, bufferSpace = bufferSpace, bufferTime = bufferTime, response=response, responsenames=responsenames, degree=degree, orderTime=orderTime, orderSpace=orderSpace, latent_dim=latent_dim, covariates=covariates)
  # pars, model = create_pars(Y, SpaceStep, TimeStep, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, c, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, response=response, responsenames=responsenames, degree=degree, orderTime=orderTime, orderSpace=orderSpace, latent_dim=latent_dim, covariates=covariates)

  # decaying learning rate
  keep_samps = Int(nits - burnin)
  n_change = log10(learning_rate / learning_rate_end)
  if n_change == 0
    change_times = 0.0
  else
    change_times = floor.(collect(range(1, stop=burnin, length=Int((n_change + 1)))))[2:end]
  end

  M_post = Array{Float64}(undef, model.N, model.D, keep_samps)
  gamma_post = Array{Float64}(undef, model.N, model.D, keep_samps)
  ΣZ_post = Array{Float64}(undef, model.N, model.N, keep_samps)
  ΣU_post = Array{Float64}(undef, model.N, model.N, keep_samps)
  A_post = Array{Float64}(undef, model.N, size(pars.A)[2], keep_samps)
  π_post = Array{Float64}(undef, model.N, keep_samps)


  @showprogress 1 "Burnin..." for i in 1:burnin

    pars = update_gamma!(pars, model)
    pars = update_M!(pars, model)
    pars = update_ΣZ!(pars, model)
    pars = update_A!(pars, model)

    if i in change_times
      raise_pow = findall(change_times .== i)[1]
      model.learning_rate = round(learning_rate * 10.0^(-raise_pow), digits = raise_pow+6)
    end

  end

  @showprogress 1 "Sampling..." for i in 1:keep_samps

    pars = update_gamma!(pars, model)
    pars = update_M!(pars, model)
    pars = update_ΣZ!(pars, model)
    pars = update_A!(pars, model)

    M_post[:, :, i] = pars.M
    gamma_post[:, :, i] = pars.gamma
    ΣZ_post[:, :, i] = pars.ΣZ
    ΣU_post[:, :, i] = pars.ΣU
    A_post[:, :, i] = pars.A
    π_post[:, i] = pars.π

  end

  posterior = Posterior(M_post, gamma_post, ΣZ_post, ΣU_post, A_post, π_post)

  return model, pars, posterior

end # 2 spatial dimensions DEtection - missing values
