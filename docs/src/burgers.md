```@meta
CurrentModule = BayesianDiscovery
```



# Burgers Example

```@example
using BayesianDiscovery
const BD = BayesianDiscovery
using Distances, Plots, Random, Distributions, LinearAlgebra
<!-- using FFTW, OrdinaryDiffEq -->
using Plots
using Random, Distributions
using Kronecker

include("../examples/burgers_equation_generate.jl")

```

Start by simulating and visualizing some data from the Burgers equation which is given by
    $$u_t(s,t) = -u(s,t)u_{x}(s,t) + \nu u_{xx}(s,t),$$
where $u(s,t)$ is the speed of the fluid at location $s = (x)$ and time $t$ and $\nu$ is the viscosity of the fluid.
```@example
U, x, Time = burgers()
Plots.contourf(Time, x, U, c=:oxy)
```

