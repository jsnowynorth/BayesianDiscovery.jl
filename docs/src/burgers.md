```@meta
CurrentModule = BayesianDiscovery
```

# Burgers Example

<!-- ```@setup burgers -->
```
using BayesianDiscovery
const BD = BayesianDiscovery
using Distances, Plots, Random, Distributions, LinearAlgebra
using FFTW
using OrdinaryDiffEq
using Plots
using Random, Distributions
using Kronecker

```


Start by simulating and visualizing some data from the Burgers equation which is given by 
```math
    u_t(s,t) = -u(s,t)u_{x}(s,t) + \nu u_{xx}(s,t),
```
where ``u(s,t)`` is the speed of the fluid at location ``s = (x)`` and time ``t`` and ``\nu`` is the viscosity of the fluid.


<!-- ```@example burgers -->
```
function burgers_pde(u, p, t)
    k,nu = p
    deriv = -u .* FFTW.ifft(1im .* k .* FFTW.fft(u)) + nu*FFTW.ifft(-k .^2 .* FFTW.fft(u))
    return real(deriv)
end


function burgers(;nx = 256, nt = 101, Lx = 8, Lt = 10.0, nu=0.1)

    # Set up grid
    x = range(-Lx, Lx, nx+1)[1:(end-1)]
    dx = x[2] - x[1]
    t = range(0, Lt, nt)
    dt = t[2]-t[1]
    k = 2*pi*FFTW.fftfreq(nx, nx*dx)

    # Initial condition
    u0 = exp.(-(x.+2).^2)


    # Solve
    p = (k, nu)
    tspan = (0.0, Lt)
    prob = OrdinaryDiffEq.ODEProblem(burgers_pde, u0, tspan, p)
    sol = OrdinaryDiffEq.solve(prob, dt = dt, OrdinaryDiffEq.Tsit5(), saveat = t)

    return reduce(hcat, sol.u), Vector(x), sol.t
    
end


U, x, Time = burgers()
Plots.contourf(Time, x, U, c=:oxy)
```

Next, we add reshape the solution surface ``U`` to be a ``space \times time \times process`` array (here ``256 \times 101 \times 1``) and add noise.
<!-- ```@example -->
```
Random.seed!(1)
Y = reshape(real.(U), 256, 101, 1)
Z = Y + 0.02 * std(Y) .* rand(Normal(), 256, 101, 1)
Plots.contourf(Time, x, Z[:,:,1], c=:oxy)
```

The MCMC sampler function, `DEtection()`, takes a lot of parameters (TO DO: break apart the function call into smaller chuncks). See `?DEtection()` for the argument requirements. Additionally, we need to define the feature library and the derivative of the feature library (see [Feature Library](@ref)).

<!-- ```@example -->
```
SpaceStep = [x]
TimeStep = Time
νS = 50 # number of space basis functions
νT = 20 # number of time basis functions
batchSpace = 10
batchTime = 10
learning_rate = 1e-4
beta = 0.9

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
```


To run the model we use the `DEtection()` function which returns the `model`, `pars`, and `posterior` structures which can be passed into various functions to investigate the output. For example, the `print_equation()` function will return the "discovered" PDE given a cutoff probability (`cutoff_prob`) value.
<!-- ```@example -->
```
model, pars, posterior = DEtection(Z, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, beta, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits = 1000, burnin = 500)
print_equation(["∇U"], model, pars, posterior, cutoff_prob=0.5)
```
