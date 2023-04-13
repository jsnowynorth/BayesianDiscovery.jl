module BayesianDiscovery

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

# StatsBase, LinearAlgebra, Plots, StatsPlots, Kronecker, Missings, Distributions, Random
# ProgressMeter, ReverseDiff, ForwardDiff, Statistics, RCall, MultivariateStats, StatsModels
# Combinatorics, IterTools, KernelDensity, SparseArrays, CodeTracking, Revise, FFTW, DifferentialEquations

using StatsBase, LinearAlgebra, Plots, StatsPlots
using Kronecker
using Missings
using Distributions, Random
using ProgressMeter
# using ReverseDiff: JacobianTape, JacobianConfig, jacobian, jacobian!, compile
# using ForwardDiff
using Statistics
using RCall
using MultivariateStats
using StatsModels, Combinatorics, IterTools
using KernelDensity
using SparseArrays
using CodeTracking, Revise
# using FFTW, DifferentialEquations



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
