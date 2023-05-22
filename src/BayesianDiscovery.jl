
########################################################################
#### Author: Joshua North
#### Project: BayesianDiscovery
#### Date: 03-April-2023
#### Description: Package for the spatio-temporal bayesian discovery paper
########################################################################

__precompile__()


"""
    BayesianDiscovery

Here is my package.
"""
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


export 
    DEtection,
    print_equation,
    posterior_surface,
    posterior_summary,
    hpd,
    Pars,
    Model,
    Posterior

end # module

