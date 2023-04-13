

mutable struct Model

  # model
  # S, T, N, D, νS, νT, bufferSpace, bufferTime, batchSpace, batchTime, learning_rate, Z, U, Basis, States

  # dimensions
  S::Int
  T::Int
  L::Int
  N::Int
  D::Int

  # data parameters
  νS
  νT
  bufferSpace
  bufferTime
  batchSpace
  batchTime
  SpaceStep
  TimeStep
  learning_rate
  inner_inds
  space_inds
  time_inds
  all_inds

  # data
  Z::Array{Float64}
  H
  X

  # basis functions
  Basis::DerivativeClass
  Θ

  # function information
  Λ::Function
  ∇Λ::Function
  ΛSpaceNames::Vector{String}
  ΛTimeNames::Vector{String}

  # response information
  responseFunction::Function
  responseNames::Vector{String}
  Ψnames::Vector{String}

  # spike-and-slab parameters
  c::Float64

  function_names

end

Base.show(io::IO, model::Model) =
  print(io, "Model\n",
    " ├─── data dimensions: ", [model.S, model.T, model.L], '\n',
    " ├─── process dimensions: ", [model.S, model.T, model.N], '\n',
    " ├─── response: ", model.responseNames, '\n',
    " ├─── library: ", model.function_names, '\n',
    " ├─── covariates: ", model.X === nothing ? false : size(model.X, 1), '\n',
    " ├─── spatial domain: ", model.SpaceStep, '\n',
    " ├─── temporal domain: ", model.TimeStep, '\n',
    " ├─── number of spatial basis functions: ", model.νS, '\n',
    " ├─── number of temporal basis functions: ", model.νT, '\n',
    " ├─── spike-and-slab parameter: ", string(model.c), '\n',
    " ├─── learning rate: ", model.learning_rate, '\n',
    " ├─── space buffer: ", model.bufferSpace, '\n',
    " ├─── time buffer: ", model.bufferTime, '\n',
    " ├─── space batch size: ", model.batchSpace, '\n',
    " └─── time batch size: ", model.batchTime, '\n')
#
