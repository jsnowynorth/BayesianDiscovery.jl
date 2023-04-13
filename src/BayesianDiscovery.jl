module BayesianDiscovery


using StatsBase, LinearAlgebra, Plots, StatsPlots
using Kronecker
using Missings
using Distributions, Random
using ProgressMeter
using Statistics
using RCall
using MultivariateStats
using StatsModels, Combinatorics, IterTools
using KernelDensity
using SparseArrays
using CodeTracking, Revise



include("../src/tensor_functions.jl")
include("../src/basis_functions.jl")
include("../src/gradient.jl")
include("../src/ParsStructure.jl")
include("../src/ModelStructure.jl")
include("../src/PDE_sampler.jl")
include("../src/process_sampler.jl")


export DEtection
export print_equation
export posterior_surface
export posterior_summary
export hpd


end # module
