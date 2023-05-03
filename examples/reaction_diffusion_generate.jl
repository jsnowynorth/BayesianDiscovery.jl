
# Pattern Formation in a Reaction-Diffusion Predator-Prey Model with Weak Allee Effect and Delay - more complex model
# https://doi.org/10.1155/2019/6282958

# ∂ut = a₁ * ∇^2u + γ₀u - (γ₀/γ₁)u^2 - βuv
# ∂vt = a₂ * ∇^2v + muv - ηv

# u: pred
# v: predator

# a₁ = 0.1 - diffusion prey
# a₂ = 0.1 - diffusion predator
# γ₀ = 0.4 - birth rate prey
# γ₁ = 1.5 - carrying capacity prey
# β = 0.5 - max uptake of prey by predator
# μ = 0.4 - birth rate of predator
# η = 0.1 - death rate predator


# Lx = 10
# Ly = 10
# T = 20
# Δx = 0.5
# Δy = 0.5
# Δt = 0.1

function ∇(X, Δx, Δy)
  Xt = X[1:(end-2), 2:(end-1)]
  Xl = X[2:(end-1), 1:(end-2)]
  Xb = X[3:end, 2:(end-1)]
  Xr = X[2:(end-1), 3:end]
  Xc = X[2:(end-1), 2:(end-1)]
  return ((Xt + Xb + Xr + Xl - 4 * Xc) / (Δx * Δy))
end

function rd(; a₁=0.1, a₂=0.1, γ₀=0.4, γ₁=1.5, β=0.5, μ=0.3, η=0.1, Lx=10, Ly=10, T=20, Δx=0.5, Δy=0.5, Δt=0.1)

  t = range(0, T, step=Δt)
  x = range(-Lx, Lx, step=Δx)
  y = range(-Ly, Ly, step=Δy)
  Nt = length(t) # number of time steps
  Nx = length(x)
  Ny = length(y)

  # empty array to store values
  U_store = Array{Float64}(undef, Nx, Ny, Nt)
  V_store = Array{Float64}(undef, Nx, Ny, Nt)

  # initialize grid
  u = [exp.(cos(2 * pi * i / 15) * sin(2 * pi * j / 15)) for i in x, j in y]
  v = [exp.(cos(2 * pi * j / 30 ) * sin(2 * pi * i / 30 - 5)) for i in x, j in y]
  U = (u .- minimum(u)) ./ (maximum(u) - minimum(u))
  V = 0.1 * (v .- minimum(v)) ./ (maximum(v) - minimum(v))

  U_store[:, :, 1] = U
  V_store[:, :, 1] = V

  for i in 2:Nt
  
    # Compute the second derivatives
    deltaU = ∇(U, Δx, Δy)
    deltaV = ∇(V, Δx, Δy)
  
    # Take values inside the grid
    Uc = U[2:(end-1), 2:(end-1)]
    Vc = V[2:(end-1), 2:(end-1)]
  
    # Update U and V
    U[2:(end-1), 2:(end-1)] = Uc + Δt * (a₁ * deltaU + γ₀ * Uc - (γ₀ / γ₁)Uc .^ 2 - (β .* Uc .* Vc))
    V[2:(end-1), 2:(end-1)] = Vc + Δt * (a₂ * deltaV + (μ .* Uc .* Vc) .- η * Vc)
  
    # boundary conditions
    U[1, :] = U[2, :]
    U[end, :] = U[end-1, :]
    U[:, 1] = U[:, 2]
    U[:, end] = U[:, end-1]
  
    V[1, :] = V[2, :]
    V[end, :] = V[end-1, :]
    V[:, 1] = V[:, 2]
    V[:, end] = V[:, end-1]
  
    U_store[:, :, i] = U
    V_store[:, :, i] = V
  
  end

  return U_store, V_store, x, y, t
end