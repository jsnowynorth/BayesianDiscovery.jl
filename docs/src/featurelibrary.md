```@meta
CurrentModule = BayesianDiscovery
```

# Feature Library Inputs

## Feature Library
Properly defining the feature library to work with the ''DEtection()'' function is challenging.
There are two input functions for the feature library, Λ and ∇Λ.
Λ is the actual library and ∇Λ is the derivative of Λ (note, this could be replaced with automatic differention, but AD is slow in this case).
Each function must take 3 parameters, A, Ψ, Φ where 'A' is a matrix of the basis coefficients, 'Ψ = [ψ, ψ_x, ψ_xx, ...]' is a vector of spatial basis functions and 'Φ = [ϕ, ϕ_t, ϕ_tt, ...]' is a vector of temporal basis functions.
The surface and its derivatives are then reconstructed as u = A * (ϕ ⊗ ψ)', u_x = A * (ϕ ⊗ ψ_x)', u_t = A * (ϕ_t ⊗ ψ)', and so on, within the function.
The first function, Λ, must return a vector of u and all of the functions in the feature library such as,
```
function Λ(A, Ψ, Φ)
    
  ψ = Ψ[1] # no spatial derivatice
  ψ_x = Ψ[2] # first order spatial derivative in the x direction
  ψ_xx = Ψ[3] # second order spatial derivative in the x direction
  ψ_xxx = Ψ[4] # third order spatial derivative in the x direction

  ϕ = Φ[1] # no temporal derivatice

  u = A * (ϕ ⊗ ψ)' # base process
  u_x = A * (ϕ ⊗ ψ_x)' # first order spatial derivative of the process in the x direction
  u_xx = A * (ϕ ⊗ ψ_xx)' # second order spatial derivative of the process in the x direction
  u_xxx = A * (ϕ ⊗ ψ_xxx)' # third order spatial derivative of the process in the x direction

  return [u, u.^2, u.^3,
    u_x, u .* u_x, u.^2 .* u_x, u.^3 .* u_x,
    u_xx, u .* u_xx, u.^2 .* u_xx, u.^3 .* u_xx,
    u_xxx, u .* u_xxx, u.^2 .* u_xxx, u.^3 .* u_xxx] # this is the feature library

end
```

## Feature Library Derivative
The second function, ∇Λ, is slightly more complicated. It requires the same inputs, however the output is the derivative of the feature library with respect to the matrix 'A'.
We have included some helper function, such as 'BayesianDiscovery.dU(A, ϕ, ψ)' or 'BayesianDiscovery.dU_x(A, ϕ, ψ_x)', but you can write your own functions if you would rather.
Here is the ∇Λ function associated with the above Λ.
```
function ∇Λ(A, Ψ, Φ)
  
  ψ = Ψ[1] # no spatial derivatice
  ψ_x = Ψ[2] # first order spatial derivative in the x direction
  ψ_xx = Ψ[3] # second order spatial derivative in the x direction
  ψ_xxx = Ψ[4] # third order spatial derivative in the x direction

  ϕ = Φ[1] # no temporal derivatice

  return [BD.dU(A, ϕ, ψ), BD.dU²(A, ϕ, ψ), BD.dU³(A, ϕ, ψ),
  BD.dU_x(A, ϕ, ψ_x), BD.dUU_x(A, ϕ, ψ, ψ_x), BD.dU²U_x(A, ϕ, ψ, ψ_x), BD.dU³U_x(A, ϕ, ψ, ψ_x),
  BD.dU_xx(A, ϕ, ψ_xx), BD.dUU_xx(A, ϕ, ψ, ψ_xx), BD.dU²U_xx(A, ϕ, ψ, ψ_xx), BD.dU³U_xx(A, ϕ, ψ, ψ_xx),
  BD.dU_xxx(A, ϕ, ψ_xxx), BD.dUU_xxx(A, ϕ, ψ, ψ_xxx), BD.dU²U_xxx(A, ϕ, ψ, ψ_xxx), BD.dU³U_xxx(A, ϕ, ψ, ψ_xxx)]

end
```

In the above function, we have 'BD' defined as
```
using BayesianSVD
const BD = BayesianDiscovery
```


## Writing your own derivative

There are a list of derivative functions we included in the 'gradients.jl' file, but there will likely be situations where you want your own function. Here is how they are created and how you might write your own.

### Derivatives of polynomials of the process
The process is defined as u = A * (ϕ ⊗ ψ)', and so du/dA = (ϕ ⊗ ψ)'. Similarly, u^2 = (A * (ϕ ⊗ ψ)')^2, and so du^2/dA = 2 .* (A * (ϕ ⊗ ψ)) .*(ϕ ⊗ ψ)'. Below are the derivatives for feature library functions u, u^2, and u^3.

```
function dU(A, ϕ, ψ)
    (ϕ ⊗ ψ)'
end

function dU²(A, ϕ, ψ)
    2 .* (A * (ϕ ⊗ ψ)) .*(ϕ ⊗ ψ)'
end

function dU³(A, ϕ, ψ)
    3 .* (A * (ϕ ⊗ ψ)).^2 .* (ϕ ⊗ ψ)'
end

```

### Derivatives of higher order derivatives of the process

To compute derivatives of du_x/dA or du_xx/dA, note that u_x = A * (ϕ ⊗ ψ_x)', and so du_x/dA = (ϕ ⊗ ψ_x)' and similary du_xx/dA = (ϕ ⊗ ψ_xx)'. Below are the equivalent functions.

```
function dU_x(A, ϕ, ψ_x)
    (ϕ ⊗ ψ_x)'
end

function dU_xx(A, ϕ, ψ_xx)
    (ϕ ⊗ ψ_xx)'
end
```

### Derivatives of higher order interaction of the process

For feature library functions like u * u_x, we have u * u_x = (A * (ϕ ⊗ ψ)') * (A * (ϕ ⊗ ψ_x)'), and so duu_x/dA = (((ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_x))) + ((A * (ϕ ⊗ ψ)) .* ((ϕ ⊗ ψ_x)')).
In this manner, you can create much more complicated function interactions and derivatives.
Here, we have the equivalent function for duu_x/dA, but if you want more examples on how to make your own function, look at the 'gradients.jl' file.

```
function dUU_x(A, ϕ, ψ, ψ_x)
    (((ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_x))) + ((A * (ϕ ⊗ ψ)) .* ((ϕ ⊗ ψ_x)'))
end
```