


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

  samp_inds = model.inner_inds

  A0 = pars.ΣM
  sN = (length(samp_inds)-1) / 2
  N = length(samp_inds)

  Uprime = model.responseFunction([pars.State.State[model.responseNames[i]] for i in 1:length(model.responseNames)])

  # update ΣU and M
  for i in 1:model.N

    Uprimei = vec(Uprime[i, samp_inds])

   
    # g prior
    Xδ = pars.F[pars.gamma[i, :] .== 1, samp_inds]
    a0 = inv(Xδ * Xδ') * Xδ * Uprimei
    A0 = N * inv(Xδ * Xδ')
    A0inv = (1/N) * Xδ * Xδ'

    Aδ = (N/(N+1)) * inv(Xδ * Xδ')
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
function update_gamma!(pars, model)

  samp_inds = StatsBase.sample(model.inner_inds, Int(length(model.inner_inds) - model.c), replace=false)

  # update spike δ
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
  
  z = model.Z[:, samp_inds] # if an error occurs with 1D vs 2D space, it is here
  h = model.H[samp_inds]
  ψ = model.Basis.SpaceDerivative["Psi"][samp_space, :]
  gψ = model.responseFunction([model.Basis.SpaceDerivative[model.Ψnames[i]] for i in 1:length(model.Ψnames)])[samp_space, :]
  ϕ = model.Basis.TimeDerivative[model.Basis.TimeNames[1]][samp_time, :]
  ϕ_t = model.Basis.TimeDerivative[model.Basis.TimeNames[end]][samp_time, :]
  Θ = model.Θ
  ΣZinv = inv(pars.ΣZ)
  ΣUinv = inv(pars.ΣU)
  A = pars.A
  M = pars.M

  fcurr = pars.F[:, samp_inds]
  fprime = create_∇Λ(model.∇Λ, A, model.Basis, samp_space, samp_time, samp_inds, model.ΛSpaceNames, model.ΛTimeNames, model.X)
  dq = [ΔL(z[:, i], h[i], ψ[i, :], gψ[i, :], ϕ[i, :], ϕ_t[i, :], Θ, ΣZinv, ΣUinv, A, M, fcurr[:, i], fprime[i]) for i in 1:length(samp_inds)]

  scale_el = 1 / 100
  dq = mean(dq) + (1 / (model.S * model.T)) * (scale_el * sign.(A) + 2 * scale_el * A)
  
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

```example 
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

```example 
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
