
function estimate_c(inner_inds; cut_off_prob = 0.9)

  b = (length(inner_inds)+1)^(1/2)
  c0 = 1
  c = round(length(inner_inds) - 2*(log(c0/b)/log(cut_off_prob) + 1), digits = 0)

  return c
  
end


"""
    Pars Structure

Initiates a structure of class Pars

# Arguments
- A: basis coefficients
- M: PDE coefficients
- ΣZ : measurement variance/covariance matrix
- ΣU: process variance/covariance matrix
- gamma: ssvs latent parameters
- sample_inds: current subsampled indicies for SGD
- F: current F
- Fprime: current Fprime
- a_R
- A_R
- nu_R
- a_Q
- A_Q
- nu_Q
- v0
- v1
- Sigma_M

"""
mutable struct Pars

  # parameters
  # A, M, ΣY, ΣU, gamma, F, a_R, A_R, nu_R, a_Q, A_Q, nu_Q, v0, v1, Sigma_M

  # estimated parameters
  A::Array{Float64}
  M::Array{Float64}
  ΣZ::Array{Float64}
  ΣU::Array{Float64}
  gamma::Array{Float64}

  # current value of function
  sample_inds::Vector{Int}
  F
  Fprime
  State::StateClass

  # hyperpriors
  a_Z
  A_Z
  nu_Z::Int

  a_U
  A_U
  nu_U::Int

  # SSVS parameters
  # v0::Float64
  # v1::Float64
  ΣM::Array{Float64}
  π

end


"""
    create_pars()

Constructs the parameter and model classes.

See also [`Pars`](@ref), [`Model`](@ref).

# Arguments
- Y::Array{Float64, 3}: Space x Time x Component data array
- SpaceStep::Vector{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}: Vector of spatial locations [x, y]
- TimeStep::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}: Time samples
- νS::Vector{Int32}: number of spatial basis functions [x, y]
- νT::Int: number of temporal basis functions
- bufferSpace::Vector{Int}: spatial buffer [x, y]
- bufferTime::Int: temporal buffer
- batchSpace::Int: batch size for space [x, y]
- batchTime::Int: batch size for time
- learning_rate::Float64: learning rate
- v0::Float64: ssvs not included variance
- v1::Float64: ssvs included variance
- Λ::Function: function library
- Λnames::Vector{String}: names of arguments in function library

"""
function create_pars(Y::Array{Float64,3},
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
  response=nothing, responsenames=nothing, degree=4, orderTime=1, orderSpace=3, latent_dim=size(Y, 3), covariates=nothing)

  S = size(Y, 1)
  T = size(Y, 2)
  L = size(Y, 3)
  N = latent_dim # need to update

  Z = unfold3(Y)

  # Δt = TimeStep[2] - TimeStep[1]
  # ΔS = [round(SpaceStep[1][2] - SpaceStep[1][1], digits=4), round(SpaceStep[2][2] - SpaceStep[2][1], digits=4)]

  # find inner index values
  Nx = length(SpaceStep[1])
  Ny = length(SpaceStep[2])
  Ns = Nx * Ny
  Nt = length(TimeStep)
  space_inds = 1:Ns
  time_inds = 1:Nt
  ind_comb = reshape([[x, y] for x = space_inds, y = time_inds], Ns * Nt)
  ind_comb = permutedims(hcat(ind_comb...))
  ind_comb = hcat(ind_comb, 1:(Ns*Nt))
  space_inds = transpose(reduce(hcat, reshape([[x, y] for x = 1:Nx, y = 1:Ny], Ns)))
  ix = space_inds[:, 1] .∉ [cat(1:1:bufferSpace[1], (Nx-bufferSpace[1]+1):1:Nx, dims=1)]
  iy = space_inds[:, 2] .∉ [cat(1:1:bufferSpace[2], (Ny-bufferSpace[2]+1):1:Ny, dims=1)]
  ispace = cat(1:Ns, dims=2)[ix.&iy]
  i1 = ind_comb[:, 1] .∈ [ispace]
  i2 = ind_comb[:, 2] .∉ [cat(1:1:bufferTime, (Nt-bufferTime+1):1:Nt, dims=1)]
  inner_inds = ind_comb[i1.&i2, 3]
  space_inds = ind_comb[:, 1]
  time_inds = ind_comb[:, 2]
  all_inds = ind_comb[:, 3]

  if covariates === nothing
    X = nothing
  else
    X = reduce(vcat, [reshape(covs[:, :, i], 1, :) for i in 1:size(covs, 3)])
  end

  # get missing values - no missing here
  H = ones(L, S * T)
  H = [make_H(L, N, H[:, i]) for i in 1:(S*T)]

  # basis functions
  Basis = derivatives_dictionary(orderSpace, orderTime, νS, νT, SpaceStep, TimeStep, degree=degree)

  Ψ = Basis.SpaceDerivative[Basis.SpaceNames[1]]
  Φ = Basis.TimeDerivative[Basis.TimeNames[1]]
  Θ = Matrix(Diagonal(ones(N)))
  A = copy(inv(transpose(Θ) * Θ) * transpose(Θ) * Z * (Φ ⊗ Ψ) * inv(transpose(Φ ⊗ Ψ) * (Φ ⊗ Ψ)))

  State = construct_state(A, Basis, Θ)

  # parameters
  samp_space, samp_time, samp_inds = sample_inds(SpaceStep, TimeStep, batchSpace, batchTime, bufferSpace, bufferTime)
  # F = create_library(Λ, Λnames, State, X)
  # dF = create_library_derivative(Λ, Λnames, State, samp_inds, X)
  F = create_Λ(Λ, A, Basis, ΛSpaceNames, ΛTimeNames, X)
  dF = create_∇Λ(∇Λ, A, Basis, samp_space, samp_time, samp_inds, ΛSpaceNames, ΛTimeNames, X)

  D = size(F, 1)
  M = zeros(N, D)
  ΣZ = Diagonal(ones(N))
  ΣU = Diagonal(ones(N))

  # hyperpriors
  a_Z = ones(N)
  A_Z = fill(1e6, N)
  nu_Z = 2

  a_U = ones(N)
  A_U = fill(1e6, N)
  nu_U = 2

  # SSVS parameters
  gamma = ones(N, D)
  c = estimate_c(inner_inds; cut_off_prob = beta)
  ΣM = inv(c * Diagonal(ones(D)))
  π = 0.5 * ones(N)


  # determine response
  if response === nothing
    responseNames = ["ΔU_t"]
    function responseFunction(ΔU)
      ΔU_t = ΔU[1]
      return ΔU_t
    end
  else
    responseFunction = response
    responseNames = responsenames
  end
  # else is user input for response function and response names

  Ψnames = (response === nothing ? ["Psi"] : "Psi_" .* reduce(hcat, split.(reduce(hcat, split.(responsenames, "t"))[1, :], "_"))[2, :])

  Ψtmp = [Basis.SpaceDerivative[ΛSpaceNames[i]] for i in 1:length(ΛSpaceNames)]
  Φtmp = [Basis.TimeDerivative[ΛTimeNames[i]] for i in 1:length(ΛTimeNames)]

  if X === nothing
    function_names = replace((@code_string Λ(A, Ψtmp, Φtmp)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
  else
    function_names = replace((@code_string Λ(A, Ψtmp, Φtmp, X)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp, X)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
  end

  # store model and parameter values
  # model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, Z, H, Basis, Θ, Λ, dΛ, Λnames, response)
  model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, space_inds, time_inds, all_inds, Z, H, X, Basis, Θ, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, responseFunction, responseNames, Ψnames, c, function_names)
  pars = Pars(A, M, ΣZ, ΣU, gamma, samp_inds, F, dF, State, a_Z, A_Z, nu_Z, a_U, A_U, nu_U, ΣM, π)

  return pars, model

end # 2 dimensions

function create_pars(Y::Array{Union{Missing,Float64},3},
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
  response=nothing, responsenames=nothing, degree=4, orderTime=1, orderSpace=3, latent_dim=size(Y, 3), covariates=nothing)

  S = size(Y, 1)
  T = size(Y, 2)
  L = size(Y, 3)
  N = latent_dim # need to update

  Z = unfold3(collect(Missings.replace(Y, 0)))

  Δt = TimeStep[2] - TimeStep[1]
  ΔS = [round(SpaceStep[1][2] - SpaceStep[1][1], digits=4), round(SpaceStep[2][2] - SpaceStep[2][1], digits=4)]

  # find inner index values
  Nx = length(SpaceStep[1])
  Ny = length(SpaceStep[2])
  Ns = Nx * Ny
  Nt = length(TimeStep)
  space_inds = 1:Ns
  time_inds = 1:Nt
  ind_comb = reshape([[x, y] for x = space_inds, y = time_inds], Ns * Nt)
  ind_comb = permutedims(hcat(ind_comb...))
  ind_comb = hcat(ind_comb, 1:(Ns*Nt))
  space_inds = transpose(reduce(hcat, reshape([[x, y] for x = 1:Nx, y = 1:Ny], Ns)))
  ix = space_inds[:, 1] .∉ [cat(1:1:bufferSpace[1], (Nx-bufferSpace[1]+1):1:Nx, dims=1)]
  iy = space_inds[:, 2] .∉ [cat(1:1:bufferSpace[2], (Ny-bufferSpace[2]+1):1:Ny, dims=1)]
  ispace = cat(1:Ns, dims=2)[ix.&iy]
  i1 = ind_comb[:, 1] .∈ [ispace]
  i2 = ind_comb[:, 2] .∉ [cat(1:1:bufferTime, (Nt-bufferTime+1):1:Nt, dims=1)]
  inner_inds = ind_comb[i1.&i2, 3]
  space_inds = ind_comb[:, 1]
  time_inds = ind_comb[:, 2]
  all_inds = ind_comb[:, 3]

  if covariates === nothing
    X = nothing
  else
    X = reduce(vcat, [reshape(covs[:, :, i], 1, :) for i in 1:size(covs, 3)])
  end

  # get missing values
  na_inds = getindex.(findall(ismissing, unfold3(Y)), 2)
  H = ones(L, S * T)
  if (size(na_inds)[1] != 0)
    H[na_inds] .= 0
  end
  H = [make_H(L, N, H[:, i]) for i in 1:(S*T)]

  # basis functions
  Basis = derivatives_dictionary(orderSpace, orderTime, νS, νT, ΔS, Δt, SpaceStep, TimeStep, degree=degree)

  Ψ = Basis.SpaceDerivative[Basis.SpaceNames[1]]
  Φ = Basis.TimeDerivative[Basis.TimeNames[1]]
  Θ = Matrix(Diagonal(ones(N)))
  A = copy(inv(transpose(Θ) * Θ) * transpose(Θ) * Z * (Φ ⊗ Ψ) * inv(transpose(Φ ⊗ Ψ) * (Φ ⊗ Ψ)))

  State = construct_state(A, Basis, Θ)

  # parameters
  samp_space, samp_time, samp_inds = sample_inds(SpaceStep, TimeStep, batchSpace, batchTime, bufferSpace, bufferTime)
  # F = create_library(Λ, Λnames, State, X)
  # dF = create_library_derivative(Λ, Λnames, State, samp_inds, X)
  F = create_Λ(Λ, A, Basis, ΛSpaceNames, ΛTimeNames, X)
  dF = create_∇Λ(∇Λ, A, Basis, samp_space, samp_time, samp_inds, ΛSpaceNames, ΛTimeNames, X)
  
  D = size(F, 1)
  M = zeros(N, D)
  ΣZ = Diagonal(ones(N))
  ΣU = Diagonal(ones(N))

  # hyperpriors
  a_Z = ones(N)
  A_Z = fill(1e6, N)
  nu_Z = 2

  a_U = ones(N)
  A_U = fill(1e6, N)
  nu_U = 2

  # SSVS parameters
  gamma = ones(N, D)
  c = estimate_c(inner_inds; cut_off_prob = beta)
  ΣM = inv(c * Diagonal(ones(D)))
  π = 0.5 * ones(N)

  # determine response
  if response === nothing
    responseNames = ["ΔU_t"]
    function responseFunction(ΔU)
      ΔU_t = ΔU[1]
      return ΔU_t
    end
  else
    responseFunction = response
    responseNames = responsenames
  end
  # else is user input for response function and response names

  Ψnames = (response === nothing ? ["Psi"] : "Psi_" .* reduce(hcat, split.(reduce(hcat, split.(responsenames, "t"))[1, :], "_"))[2, :])

  Ψtmp = [Basis.SpaceDerivative[ΛSpaceNames[i]] for i in 1:length(ΛSpaceNames)]
  Φtmp = [Basis.TimeDerivative[ΛTimeNames[i]] for i in 1:length(ΛTimeNames)]
  if X === nothing
    function_names = replace((@code_string Λ(A, Ψtmp, Φtmp)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
  else
    function_names = replace((@code_string Λ(A, Ψtmp, Φtmp, X)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp, X)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
  end

  # store model and parameter values
  # model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, Z, H, Basis, Θ, Λ, dΛ, Λnames, response)
  model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, space_inds, time_inds, all_inds, Z, H, X, Basis, Θ, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, responseFunction, responseNames, Ψnames, c, function_names)
  pars = Pars(A, M, ΣZ, ΣU, gamma, samp_inds, F, dF, State, a_Z, A_Z, nu_Z, a_U, A_U, nu_U, ΣM, π)

  return pars, model

end # 2 dimensions missing data

function create_pars(Y::Array{Float64,3},
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
  response=nothing, responsenames=nothing, degree=4, orderTime=1, orderSpace=3, latent_dim=size(Y, 3), covariates=nothing)

  S = size(Y, 1)
  T = size(Y, 2)
  L = size(Y, 3)
  N = latent_dim # need to update

  Z = unfold3(Y)
  # Δt = TimeStep[2] - TimeStep[1]
  # ΔS = round(SpaceStep[2] - SpaceStep[1], digits=4)

  # find inner index values
  Nx = length(SpaceStep[1])
  Nt = length(TimeStep)
  space_inds = 1:Nx
  time_inds = 1:Nt
  ind_comb = reshape([[x, y] for x = space_inds, y = time_inds], length(space_inds) * length(time_inds))
  ind_comb = permutedims(hcat(ind_comb...))
  ind_comb = hcat(ind_comb, 1:(Nx*Nt))
  i1 = ind_comb[:, 1] .∉ [cat(1:1:bufferSpace, (Nx-bufferSpace+1):1:Nx, dims=1)]
  i2 = ind_comb[:, 2] .∉ [cat(1:1:bufferTime, (Nt-bufferTime+1):1:Nt, dims=1)]
  inner_inds = ind_comb[i1.&i2, 3]
  space_inds = ind_comb[:, 1]
  time_inds = ind_comb[:, 2]
  all_inds = ind_comb[:, 3]

  if covariates === nothing
    X = nothing
  else
    X = reduce(vcat, [reshape(covs[:, :, i], 1, :) for i in 1:size(covs, 3)])
  end

  # get missing values - no missing here
  H = ones(L, S * T)
  H = [make_H(L, N, H[:, i]) for i in 1:(S*T)]

  # basis functions
  # Basis = derivatives_dictionary(orderSpace, orderTime, νS, νT, ΔS, Δt, SpaceStep, TimeStep, degree=degree)
  Basis = derivatives_dictionary(orderSpace, orderTime, νS, νT, SpaceStep, TimeStep, degree=degree)

  Ψ = Basis.SpaceDerivative[Basis.SpaceNames[1]]
  Φ = Basis.TimeDerivative[Basis.TimeNames[1]]
  Θ = Matrix(Diagonal(ones(N)))
  A = copy(inv(transpose(Θ) * Θ) * transpose(Θ) * Z * (Φ ⊗ Ψ) * inv(transpose(Φ ⊗ Ψ) * (Φ ⊗ Ψ)))
  State = construct_state(A, Basis, Θ)

  # parameters
  samp_space, samp_time, samp_inds = sample_inds(SpaceStep, TimeStep, batchSpace, batchTime, bufferSpace, bufferTime)
  # F = create_library(Λ, Λnames, State, X)
  # dF = create_library_derivative(Λ, Λnames, State, samp_inds, X)
  F = create_Λ(Λ, A, Basis, ΛSpaceNames, ΛTimeNames, X)
  dF = create_∇Λ(∇Λ, A, Basis, samp_space, samp_time, samp_inds, ΛSpaceNames, ΛTimeNames, X)
  
  D = size(F, 1)
  M = zeros(N, D)
  ΣZ = Diagonal(ones(N))
  ΣU = Diagonal(ones(N))
  # hyperpriors
  a_Z = ones(N)
  A_Z = fill(1e6, N)
  nu_Z = 2

  a_U = ones(N)
  A_U = fill(1e6, N)
  nu_U = 2
  # SSVS parameters
  gamma = ones(N, D)
  c = estimate_c(inner_inds; cut_off_prob = beta)
  ΣM = inv(c * Diagonal(ones(D)))
  π = 0.5 * ones(N)

  # determine response
  if response === nothing
    responseNames = ["ΔU_t"]
    function responseFunction(ΔU)
      ΔU_t = ΔU[1]
      return ΔU_t
    end
  else
    responseFunction = response
    responseNames = responsenames
  end
  # else is user input for response function and response names

  Ψnames = (response === nothing ? ["Psi"] : "Psi_" .* reduce(hcat, split.(reduce(hcat, split.(responsenames, "t"))[1, :], "_"))[2, :])
  
  
  Ψtmp = [Basis.SpaceDerivative[ΛSpaceNames[i]] for i in 1:length(ΛSpaceNames)]
  Φtmp = [Basis.TimeDerivative[ΛTimeNames[i]] for i in 1:length(ΛTimeNames)]
  if X === nothing
    function_names = replace((@code_string Λ(A, Ψtmp, Φtmp)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
  else
    function_names = replace((@code_string Λ(A, Ψtmp, Φtmp, X)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp, X)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
  end

  # store model and parameter values
  # model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, Z, H, Basis, Θ, Λ, dΛ, Λnames, response)
  model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, space_inds, time_inds, all_inds, Z, H, X, Basis, Θ, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, responseFunction, responseNames, Ψnames, c, function_names)
  pars = Pars(A, M, ΣZ, ΣU, gamma, samp_inds, F, dF, State, a_Z, A_Z, nu_Z, a_U, A_U, nu_U, ΣM, π)

  return pars, model

end # 1 dimension

function create_pars(Y::Array{Union{Missing,Float64},3},
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
  response=nothing, responsenames=nothing, degree=4, orderTime=1, orderSpace=3, latent_dim=size(Y, 3), covariates=nothing)

  S = size(Y, 1)
  T = size(Y, 2)
  L = size(Y, 3)
  N = latent_dim # need to update

  na_inds = getindex.(findall(ismissing, unfold3(Y)), 2)
  Z = unfold3(collect(Missings.replace(Y, 0)))
  # Δt = TimeStep[2] - TimeStep[1]
  # ΔS = round(SpaceStep[2] - SpaceStep[1], digits=4)

  # find inner index values
  Nx = length(SpaceStep[1])
  Nt = length(TimeStep)
  space_inds = 1:Nx
  time_inds = 1:Nt
  ind_comb = reshape([[x, y] for x = space_inds, y = time_inds], length(space_inds) * length(time_inds))
  ind_comb = permutedims(hcat(ind_comb...))
  ind_comb = hcat(ind_comb, 1:(Nx*Nt))
  i1 = ind_comb[:, 1] .∉ [cat(1:1:bufferSpace, (Nx-bufferSpace+1):1:Nx, dims=1)]
  i2 = ind_comb[:, 2] .∉ [cat(1:1:bufferTime, (Nt-bufferTime+1):1:Nt, dims=1)]
  # i3 = ind_comb[:, 3] .∉ [cat(na_inds, dims=1)]
  # inner_inds = ind_comb[i1.&i2.&i3, 3]
  inner_inds = ind_comb[i1.&i2, 3]
  space_inds = ind_comb[:, 1]
  time_inds = ind_comb[:, 2]
  all_inds = ind_comb[:, 3]

  if covariates === nothing
    X = nothing
  else
    X = reduce(vcat, [reshape(covs[:, :, i], 1, :) for i in 1:size(covs, 3)])
  end

  # get missing values
  # na_inds = getindex.(findall(ismissing, unfold3(Y)), 2)
  H = ones(L, S * T)
  if (size(na_inds)[1] != 0)
    H[na_inds] .= 0
  end
  H = [make_H(L, N, H[:, i]) for i in 1:(S*T)]


  # basis functions
  # Basis = derivatives_dictionary(orderSpace, orderTime, νS, νT, ΔS, Δt, SpaceStep, TimeStep, degree=degree)
  Basis = derivatives_dictionary(orderSpace, orderTime, νS, νT, SpaceStep, TimeStep, degree=degree)

  Z_start = Z
  Z_start[na_inds] = (Z[na_inds .+ 1] .+ Z[na_inds .- 1])./2

  Ψ = Basis.SpaceDerivative[Basis.SpaceNames[1]]
  Φ = Basis.TimeDerivative[Basis.TimeNames[1]]
  Θ = Matrix(Diagonal(ones(N)))
  A = copy(inv(transpose(Θ) * Θ) * transpose(Θ) * Z_start * (Φ ⊗ Ψ) * inv(transpose(Φ ⊗ Ψ) * (Φ ⊗ Ψ)))
  State = construct_state(A, Basis, Θ)

  # parameters
  samp_space, samp_time, samp_inds = sample_inds(SpaceStep, TimeStep, batchSpace, batchTime, bufferSpace, bufferTime)
  # F = create_library(Λ, Λnames, State, X)
  # dF = create_library_derivative(Λ, Λnames, State, samp_inds, X)
  F = create_Λ(Λ, A, Basis, ΛSpaceNames, ΛTimeNames, X)
  dF = create_∇Λ(∇Λ, A, Basis, samp_space, samp_time, samp_inds, ΛSpaceNames, ΛTimeNames, X)
  
  D = size(F, 1)
  M = zeros(N, D)
  ΣZ = Diagonal(ones(N))
  ΣU = Diagonal(ones(N))
  # hyperpriors
  a_Z = ones(N)
  A_Z = fill(1e6, N)
  nu_Z = 2

  a_U = ones(N)
  A_U = fill(1e6, N)
  nu_U = 2
  # SSVS parameters
  gamma = ones(N, D)
  c = estimate_c(inner_inds; cut_off_prob = beta)
  ΣM = inv(c * Diagonal(ones(D)))
  π = 0.5 * ones(N)

  # determine response
  if response === nothing
    responseNames = ["ΔU_t"]
    function responseFunction(ΔU)
      ΔU_t = ΔU[1]
      return ΔU_t
    end
  else
    responseFunction = response
    responseNames = responsenames
  end
  # else is user input for response function and response names

  Ψnames = (response === nothing ? ["Psi"] : "Psi_" .* reduce(hcat, split.(reduce(hcat, split.(responsenames, "t"))[1, :], "_"))[2, :])

  Ψtmp = [Basis.SpaceDerivative[ΛSpaceNames[i]] for i in 1:length(ΛSpaceNames)]
  Φtmp = [Basis.TimeDerivative[ΛTimeNames[i]] for i in 1:length(ΛTimeNames)]
  if X === nothing
    function_names = replace((@code_string Λ(A, Ψtmp, Φtmp)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
  else
    function_names = replace((@code_string Λ(A, Ψtmp, Φtmp, X)).string[findfirst("return", (@code_string Λ(A, Ψtmp, Φtmp, X)).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "", "." => "")
  end

  # store model and parameter values
  # model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, Z, H, Basis, Θ, Λ, dΛ, Λnames, response)
  model = Model(S, T, L, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, SpaceStep, TimeStep, learning_rate, inner_inds, space_inds, time_inds, all_inds, Z, H, X, Basis, Θ, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, responseFunction, responseNames, Ψnames, c, function_names)
  pars = Pars(A, M, ΣZ, ΣU, gamma, samp_inds, F, dF, State, a_Z, A_Z, nu_Z, a_U, A_U, nu_U, ΣM, π)

  return pars, model

end # 1 dimension missing data
