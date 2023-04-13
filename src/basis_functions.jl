
struct DerivativeClass
  TimeDerivative::Dict
  TimeNames::Vector{String}
  SpaceDerivative::Dict
  SpaceNames::Vector{String}
  orderSpace::Int
  orderTime::Int
  spatialDimension::Int
end

struct StateClass
  State::Dict
  StateNames::Vector{String}
end

"""
    bspline(x, knot_locs, degree, derivative)

Creates a matrix of B-Splines evaluated at each x with specified knot locations and degree.
Dimension is (length(x), length(x_knot_locs)).
Returns the basis matrix and the derivative of the basis matrix.

# Arguments
- x: data
- knot_locs: knot locations
- degree: degree of the B-Spline
- derivative: order of the derivative for Psi_x

# Examples
```
x = 1:1:30
knot_locs = 5:5:25
degree = 3
derivative = 1

Psi, Psi_x = bspline(x, knot_locs, degree, derivative)
```
"""
function bspline(x, knot_locs, degree, derivative)

  @rput x
  @rput knot_locs
  @rput degree
  @rput derivative
  # bxtmp = R"splines2::bSpline(x, knots = knot_locs, degree = degree, intercept = TRUE)"
  # dbxtmp = R"splines2::dbs(x, knots = knot_locs, degree = degree, derivs = derivative, intercept = TRUE)"

  # if derivative == 0
  #   btmp = R"splines2::bSpline(x, knots = knot_locs, degree = degree, intercept = TRUE)"
  # else
  #   btmp = R"splines2::dbs(x, knots = knot_locs, degree = degree, derivs = derivative, intercept = TRUE)"
  # end
  if derivative == 0
    btmp = R"splines2::bSpline(x, df = knot_locs, degree = degree, intercept = TRUE)"
  else
    btmp = R"splines2::dbs(x, df = knot_locs, degree = degree, derivs = derivative, intercept = TRUE)"
  end

  Phi = rcopy(btmp)
  # Phi = rcopy(bxtmp)
  # Phi_x = rcopy(dbxtmp)

  # return Phi, Phi_x
  return Phi
end


"""
    spatial_bspline(x, y, x_knot_locs, y_knot_locs, degree)

Creates a matrix of B-Splines evaluated at each x and y with specified knot locations and degree.
Dimension is (length(x)*length(y), length(x_knot_locs)*length(y_knot_locs)).
Returns the basis matrix and the derivative of the basis matrix.

# Arguments
- x: data in the x direction
- y: data in the x direction
- x_knot_locs: knot locations for x
- y_knot_locs: knot locations for y
- degree: degree of the B-Spline

# Examples
```
x = 1:1:30
y = 1:1:30
x_knot_locs = 5:5:25
y_knot_locs = 5:5:25
degree = 3

Psi, Psi_x, Psi_y, Psi_xy = spatial_bspline(x, y, x_knot_locs, y_knot_locs, degree)

# plot the functions
Plots.pyplot()
contour(x, y, reshape(Psi, length(x), length(y)), fill = true)
contour(x, y, reshape(Psi_x, length(x), length(y)), fill = true)
contour(x, y, reshape(Psi_y, length(x), length(y)), fill = true)
contour(x, y, reshape(Psi_xy, length(x), length(y)), fill = true)
```
"""
function spatial_bspline(x, y, x_knot_locs, y_knot_locs, degree, xdir, ydir)

  @rput x
  @rput y
  @rput x_knot_locs
  @rput y_knot_locs
  @rput degree
  bxtmp = R"splines2::bSpline(x, knots = x_knot_locs, degree = degree, intercept = TRUE)"
  dbxtmp = R"splines2::dbs(x, knots = x_knot_locs, degree = degree, derivs = 1, intercept = TRUE)"
  dbxxtmp = R"splines2::dbs(x, knots = x_knot_locs, degree = degree, derivs = 2, intercept = TRUE)"

  # bxtmp = DEtection.Rbspline(x, knots = x_knot_locs, degree, intercept = true)
  # dbxtmp = DEtection.Rbspline(x, knots = x_knot_locs, degree, derivs = 1, intercept = true)
  # dbxxtmp = DEtection.Rbspline(x, knots = x_knot_locs, degree, derivs = 2, intercept = true)
  bx = rcopy(bxtmp)
  dbx = rcopy(dbxtmp)
  dbxx = rcopy(dbxxtmp)

  # bytmp = DEtection.Rbspline(y, knots = y_knot_locs, degree, intercept = true)
  # dbytmp = DEtection.Rbspline(y, knots = y_knot_locs, degree, derivs = 1, intercept = true)
  # dbyytmp = DEtection.Rbspline(y, knots = y_knot_locs, degree, derivs = 2, intercept = true)
  bytmp = R"splines2::bSpline(y, knots = y_knot_locs, degree, intercept = TRUE)"
  dbytmp = R"splines2::dbs(y, knots = y_knot_locs, degree, derivs = 1, intercept = TRUE)"
  dbyytmp = R"splines2::dbs(y, knots = y_knot_locs, degree, derivs = 2, intercept = TRUE)"
  by = rcopy(bytmp)
  dby = rcopy(dbytmp)
  dbyy = rcopy(dbyytmp)


  nx = length(x)
  ny = length(y)

  nxk = size(bx)[2]
  nyk = size(by)[2]

  Phi = reduce(hcat, reshape([reshape(bx[i, :] .* by[j, :]', nxk * nyk, 1) for i in 1:nx, j in 1:ny], nx * ny, 1))'
  Phi_y = reduce(hcat, reshape([reshape(dby[i, :] .* bx[j, :]', nxk * nyk, 1) for i in 1:nx, j in 1:ny], nx * ny, 1))'
  Phi_x = reduce(hcat, reshape([reshape(by[i, :] .* dbx[j, :]', nxk * nyk, 1) for i in 1:nx, j in 1:ny], nx * ny, 1))'
  Phi_xy = reduce(hcat, reshape([reshape(dby[i, :] .* dbx[j, :]', nxk * nyk, 1) for i in 1:nx, j in 1:ny], nx * ny, 1))'
  Phi_yy = reduce(hcat, reshape([reshape(dbyy[i, :] .* bx[j, :]', nxk * nyk, 1) for i in 1:nx, j in 1:ny], nx * ny, 1))'
  Phi_xx = reduce(hcat, reshape([reshape(by[i, :] .* dbxx[j, :]', nxk * nyk, 1) for i in 1:nx, j in 1:ny], nx * ny, 1))'

  return Phi, Phi_x, Phi_y, Phi_xy, Phi_xx, Phi_yy
end


function derivatives_dictionary(orderSpace::Int,
  orderTime::Int,
  νS::Int,
  νT::Int,
  SpaceStep::Vector,
  TimeStep::Vector;
  degree::Int=4)

  x = SpaceStep[1]

  Time_derivative = Dict()
  Time_derivative["Phi"] = bspline(TimeStep, νT, degree, 0)
  Time_names = ["Phi"]
  for i in 1:orderTime
    Time_derivative["Phi_".*repeat('t', i)] = bspline(TimeStep, νT, degree, i)
    push!(Time_names, "Phi_" .* repeat('t', i))
  end


  Space_derivative = Dict()
  Space_derivative["Psi"] = bspline(x, νS, degree, 0)
  Space_names = ["Psi"]
  for i in 1:orderSpace
    Space_derivative["Psi_".*repeat('x', i)] = bspline(x, νS, degree, i)
    push!(Time_names, "Psi_" .* repeat('x', i))
  end

  Time_names = vcat("Phi", ["Phi_" .* repeat('t', i) for i in 1:orderTime])
  Space_names = vcat("Psi", ["Psi_" .* repeat('x', i) for i in 1:orderSpace])

  BasisDerivative = DerivativeClass(Time_derivative, Time_names, Space_derivative, Space_names, orderSpace, orderTime, 1)

  return BasisDerivative
end

function derivatives_dictionary(orderSpace::Int,
  orderTime::Int,
  νS::Vector{Int},
  νT::Int,
  SpaceStep::Vector,
  TimeStep::Vector;
  degree::Int=4)

  νSx = νS[1]
  νSy = νS[2]
  x = SpaceStep[1]
  y = SpaceStep[2]

  nx = length(x)
  ny = length(y)


  Time_derivative = Dict()
  Time_derivative["Phi"] = bspline(TimeStep, νT, degree, 0)
  Time_names = ["Phi"]
  for i in 1:orderTime
    Time_derivative["Phi_".*repeat('t', i)] = bspline(TimeStep, νT, degree, i)
    push!(Time_names, "Phi_" .* repeat('t', i))
  end

  Space_derivative = Dict()
  Space_names = []
  for i in 0:orderSpace
    for j in 0:orderSpace
      bx = bspline(x, νSx, degree, i)
      by = bspline(y, νSy, degree, j)
      if i == j == 0
        Space_derivative["Psi"] = reduce(hcat, reshape([reshape(by[i, :] .* bx[j, :]', :, 1) for i in 1:ny, j in 1:nx], nx * ny, 1))'
        # Space_derivative["Psi"] = reduce(hcat, reshape([reshape(bx[i, :] .* by[j, :]', :, 1) for i in 1:nx, j in 1:ny], nx * ny, 1))'
        push!(Space_names, "Psi")
      else
        Space_derivative["Psi_".*repeat('x', j).*repeat('y', i)] = reduce(hcat, reshape([reshape(by[i, :] .* bx[j, :]', :, 1) for i in 1:ny, j in 1:nx], nx * ny, 1))'
        # Space_derivative["Psi_".*repeat('x', j).*repeat('y', i)] = reduce(hcat, reshape([reshape(bx[i, :] .* by[j, :]', :, 1) for i in 1:nx, j in 1:ny], nx * ny, 1))
        push!(Space_names, "Psi_" .* repeat('x', j) .* repeat('y', i))
      end
    end
  end


  BasisDerivative = DerivativeClass(Time_derivative, Time_names, Space_derivative, Space_names, orderSpace, orderTime, 2)

  return BasisDerivative
end

function construct_state(A::Array{Float64,2}, BasisDerivative::DerivativeClass, Θ::Matrix{Float64})

  C = fold3(A, size(BasisDerivative.SpaceDerivative[BasisDerivative.SpaceNames[1]], 2), size(BasisDerivative.TimeDerivative[BasisDerivative.TimeNames[1]], 2), size(Θ, 1))
  orderTime = BasisDerivative.orderTime
  orderSpace = BasisDerivative.orderSpace

  State = unfold3(tensor_mult(C, BasisDerivative.SpaceDerivative[BasisDerivative.SpaceNames[1]], BasisDerivative.TimeDerivative[BasisDerivative.TimeNames[1]], Θ))

  # initialize state dictionary and names
  State_process = Dict()
  # State_names = Vector{String}[]

  # save base state
  State_process["U"] = State
  State_names = ["U"]
  # push!(State_names, "U")

  # save time derivatives
  for i in 1:orderTime
    State_process["ΔU_".*repeat('t', i)] = unfold3(tensor_mult(C, BasisDerivative.SpaceDerivative[BasisDerivative.SpaceNames[1]], BasisDerivative.TimeDerivative[BasisDerivative.TimeNames[1+i]], Θ))
    push!(State_names, "ΔU_" .* repeat('t', i))
  end


  if BasisDerivative.spatialDimension == 1

    for i in 1:orderSpace
      State_process["ΔU_".*repeat('x', i)] = unfold3(tensor_mult(C, BasisDerivative.SpaceDerivative[BasisDerivative.SpaceNames[i+1]], BasisDerivative.TimeDerivative[BasisDerivative.TimeNames[1]], Θ))
      push!(State_names, "ΔU_" .* repeat('x', i))
    end
    ProcessDerivative = StateClass(State_process, State_names)

  else

    for j in 0:orderTime
      names = reduce(hcat, split.(BasisDerivative.SpaceNames[2:end], "_"))[2, :] .* repeat("t", j)
      name_index = hcat(names, 2:length(BasisDerivative.SpaceNames))

      for i in 1:size(name_index, 1)
        State_process["ΔU_"*name_index[i]] = unfold3(tensor_mult(C, BasisDerivative.SpaceDerivative[BasisDerivative.SpaceNames[name_index[i, 2]]], BasisDerivative.TimeDerivative[BasisDerivative.TimeNames[j+1]], Θ))
        push!(State_names, "ΔU_" * name_index[i])
      end
    end
    ProcessDerivative = StateClass(State_process, State_names)

  end

  return ProcessDerivative
end



function create_library(Λ::Function, Λnames::Vector{String}, locs::Vector{Int}, State::StateClass)

  # Uselect = [State.State[Λnames[1]][:, locs[j]] for j in 1:length(locs)]
  # ∂Uselect = [State.State[Λnames[i]][:, locs[j]] for i in 2:length(Λnames), j in 1:length(locs)]
  # F = reduce(hcat, [Λ(Uselect[j], ∂Uselect[:, j]) for j in 1:length(locs)])

  # return F

  N = length(locs)

  if X === nothing
    Uselect = [State.State[Λnames[1]][:, locs[j]] for j in 1:N]
    ∂Uselect = [State.State[Λnames[i]][:, locs[j]] for i in 2:length(Λnames), j in 1:N]
    F = reduce(hcat, [Λ(Uselect[j], ∂Uselect[:, j]) for j in 1:N])
  else
    Uselect = [State.State[Λnames[1]][:, locs[j]] for j in 1:N]
    ∂Uselect = [State.State[Λnames[i]][:, locs[j]] for i in 2:length(Λnames), j in 1:N]
    F = reduce(hcat, [Λ(Uselect[j], ∂Uselect[:, j], X[:, j]) for j in 1:N])
  end

  return F

end

function create_library(Λ::Function, Λnames::Vector{String}, State::StateClass, X=nothing)

  # N = size(State.State[Λnames[1]], 2)

  # Uselect = [State.State[Λnames[1]][:, j] for j in 1:N]
  # ∂Uselect = [State.State[Λnames[i]][:, j] for i in 2:length(Λnames), j in 1:N]
  # F = reduce(hcat, [Λ(Uselect[j], ∂Uselect[:, j]) for j in 1:N])

  # return F

  N = size(State.State[Λnames[1]], 2)

  if X === nothing
    Uselect = [State.State[Λnames[1]][:, j] for j in 1:N]
    ∂Uselect = [State.State[Λnames[i]][:, j] for i in 2:length(Λnames), j in 1:N]
    F = reduce(hcat, [Λ(Uselect[j], ∂Uselect[:, j]) for j in 1:N])
  else
    Uselect = [State.State[Λnames[1]][:, j] for j in 1:N]
    ∂Uselect = [State.State[Λnames[i]][:, j] for i in 2:length(Λnames), j in 1:N]
    F = reduce(hcat, [Λ(Uselect[j], ∂Uselect[:, j], X[:, j]) for j in 1:N])
  end

  return F

end


function create_library_derivative(Λ::Function, Λnames::Vector{String}, State::StateClass, locs::Vector{Int}, X=nothing)

  # N = length(locs)
  # Uselect = reduce(hcat, [State.State[Λnames[1]][:, j] for j in locs])
  # ∂Uselect = [State.State[Λnames[i]][:, j] for i in 2:length(Λnames), j in locs]
  # dF = [jacobian(U -> Λ(U, ∂Uselect[:, j]), Uselect[:, j]) for j in 1:N]

  N = length(locs)

  if X === nothing
    Uselect = reduce(hcat, [State.State[Λnames[1]][:, j] for j in locs])
    ∂Uselect = [State.State[Λnames[i]][:, j] for i in 2:length(Λnames), j in locs]
    dF = [jacobian(U -> Λ(U, ∂Uselect[:, j]), Uselect[:, j]) for j in 1:N]
  else
    Uselect = reduce(hcat, [State.State[Λnames[1]][:, j] for j in locs])
    ∂Uselect = [State.State[Λnames[i]][:, j] for i in 2:length(Λnames), j in locs]
    dF = [vcat(jacobian(U -> Λ(U, ∂Uselect[:, j], X[:, j]), Uselect[:, j])) for j in 1:N]
  end

  return dF

end


function create_Λ(Λ::Function, A, Basis, samp_space, samp_time, samp_inds, ΛSpaceNames, ΛTimeNames, X=nothing)

  if X === nothing
    Ψ = [Basis.SpaceDerivative[ΛSpaceNames[i]][j,:] for i in 1:length(ΛSpaceNames), j in samp_space]
    Φ = [Basis.TimeDerivative[ΛTimeNames[i]][j,:] for i in 1:length(ΛTimeNames), j in samp_time]
    Λeval = [reduce(vcat, Λ(A, Ψ[:,i], Φ[:,i])) for i in 1:length(samp_inds)]
  else
    Ψ = [Basis.SpaceDerivative[ΛSpaceNames[i]][j,:] for i in 1:length(ΛSpaceNames), j in samp_space]
    Φ = [Basis.TimeDerivative[ΛTimeNames[i]][j,:] for i in 1:length(ΛTimeNames), j in samp_time]
    Λeval = [reduce(vcat, Λ(A, Ψ[:,i], Φ[:,i], X[:, samp_inds[i]])) for i in 1:length(samp_space)]
  end

  return reduce(hcat, Λeval)

end



function create_Λ(Λ::Function, A, Basis, ΛSpaceNames, ΛTimeNames, X=nothing)

  if X === nothing
    Ψ = [Basis.SpaceDerivative[ΛSpaceNames[i]] for i in 1:length(ΛSpaceNames)]
    Φ = [Basis.TimeDerivative[ΛTimeNames[i]] for i in 1:length(ΛTimeNames)]
    Λeval = reduce(vcat, Λ(A, Ψ, Φ))
  else
    Ψ = [Basis.SpaceDerivative[ΛSpaceNames[i]] for i in 1:length(ΛSpaceNames)]
    Φ = [Basis.TimeDerivative[ΛTimeNames[i]] for i in 1:length(ΛTimeNames)]
    xtmp = [X[i,:]' for i in 1:size(X,1)]
    Λeval = reduce(vcat, Λ(A, Ψ, Φ, xtmp))
  end

  return Λeval

end



function create_∇Λ(∇Λ::Function, A, Basis, samp_space, samp_time, samp_inds, ΛSpaceNames, ΛTimeNames, X=nothing)

  if size(A,1) > 1
    if X === nothing
      Ψ = [Basis.SpaceDerivative[ΛSpaceNames[i]][j,:] for i in 1:length(ΛSpaceNames), j in samp_space]
      Φ = [Basis.TimeDerivative[ΛTimeNames[i]][j,:] for i in 1:length(ΛTimeNames), j in samp_time]
      dΛeval = [∇Λ(A, Ψ[:,i], Φ[:,i]) for i in 1:length(samp_space)]
    else
      Ψ = [Basis.SpaceDerivative[ΛSpaceNames[i]][j,:] for i in 1:length(ΛSpaceNames), j in samp_space]
      Φ = [Basis.TimeDerivative[ΛTimeNames[i]][j,:] for i in 1:length(ΛTimeNames), j in samp_time]
      dΛeval = [∇Λ(A, Ψ[:,i], Φ[:,i], X[:, samp_inds[i]]) for i in 1:length(samp_space)]
    end
  else
    if X === nothing
      Ψ = [Basis.SpaceDerivative[ΛSpaceNames[i]][j,:] for i in 1:length(ΛSpaceNames), j in samp_space]
      Φ = [Basis.TimeDerivative[ΛTimeNames[i]][j,:] for i in 1:length(ΛTimeNames), j in samp_time]
      dΛeval = [reduce(vcat, ∇Λ(A, Ψ[:,i], Φ[:,i])) for i in 1:length(samp_space)]
    else
      Ψ = [Basis.SpaceDerivative[ΛSpaceNames[i]][j,:] for i in 1:length(ΛSpaceNames), j in samp_space]
      Φ = [Basis.TimeDerivative[ΛTimeNames[i]][j,:] for i in 1:length(ΛTimeNames), j in samp_time]
      dΛeval = [reduce(vcat, ∇Λ(A, Ψ[:,i], Φ[:,i], X[:, samp_inds[i]])) for i in 1:length(samp_space)]
    end
  end
  

  return dΛeval

end



function sample_inds(SpaceStep::Vector,
  TimeStep::Vector,
  batchSpace::Int,
  batchTime::Int,
  bufferSpace::Vector{Int},
  bufferTime::Int)
  # need to add the buffer in

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
  ix = space_inds[:, 1] .∉ [cat(1:1:bufferSpace[1], (Nx-bufferSpace[1]):1:Nx, dims=1)]
  iy = space_inds[:, 2] .∉ [cat(1:1:bufferSpace[2], (Ny-bufferSpace[2]):1:Ny, dims=1)]
  ispace = cat(1:Ns, dims=2)[ix.&iy]
  i1 = ind_comb[:, 1] .∈ [ispace]

  i2 = ind_comb[:, 2] .∉ [cat(1:1:bufferTime, (Nt-bufferTime):1:Nt, dims=1)]
  ind_comb = ind_comb[i1.&i2, :]

  sel_inds = ind_comb[ind_comb[:, 3].∈[StatsBase.sample(ind_comb[:, 3], batchSpace * batchTime, replace=false)], :]

  return sel_inds[:, 1], sel_inds[:, 2], sel_inds[:, 3]

end # 2d with buffer

function sample_inds(SpaceStep::Vector,
  TimeStep::Vector,
  batchSpace::Int,
  batchTime::Int)

  Nx = length(SpaceStep[1])
  Ny = length(SpaceStep[2])
  Ns = Nx * Ny
  Nt = length(TimeStep)
  space_inds = 1:Ns
  time_inds = 1:Nt

  ind_comb = reshape([[x, y] for x = space_inds, y = time_inds], Ns * Nt)
  ind_comb = permutedims(hcat(ind_comb...))
  ind_comb = hcat(ind_comb, 1:(Ns*Nt))

  sel_inds = ind_comb[ind_comb[:, 3].∈[StatsBase.sample(ind_comb[:, 3], batchSpace * batchTime, replace=false)], :]

  return sel_inds[:, 1], sel_inds[:, 2], sel_inds[:, 3]

end # 2d without buffer

function sample_inds(SpaceStep::Vector,
  TimeStep::Vector,
  batchSpace::Int,
  batchTime::Int,
  bufferSpace::Int,
  bufferTime::Int)

  Nx = length(SpaceStep[1])
  Nt = length(TimeStep)
  space_inds = 1:Nx
  time_inds = 1:Nt

  ind_comb = reshape([[x, y] for x = space_inds, y = time_inds], length(space_inds) * length(time_inds))
  ind_comb = permutedims(hcat(ind_comb...))
  ind_comb = hcat(ind_comb, 1:(Nx*Nt))
  i1 = ind_comb[:, 1] .∉ [cat(1:1:bufferSpace, (Nx-bufferSpace):1:Nx, dims=1)]
  i2 = ind_comb[:, 2] .∉ [cat(1:1:bufferTime, (Nt-bufferTime):1:Nt, dims=1)]
  ind_comb = ind_comb[i1.&i2, :]

  sel_inds = ind_comb[ind_comb[:, 3].∈ [StatsBase.sample(ind_comb[:, 3], batchSpace * batchTime, replace=false)], :]

  return sel_inds[:, 1], sel_inds[:, 2], sel_inds[:, 3]

end # 1d with buffer

function sample_inds(SpaceStep::Vector,
  TimeStep::Vector,
  batchSpace::Int,
  batchTime::Int)

  Nx = length(SpaceStep[1])
  Nt = length(TimeStep)
  space_inds = 1:Nx
  time_inds = 1:Nt

  ind_comb = reshape([[x, y] for x = space_inds, y = time_inds], length(space_inds) * length(time_inds))
  ind_comb = permutedims(hcat(ind_comb...))
  ind_comb = hcat(ind_comb, 1:(Nx*Nt))

  sel_inds = ind_comb[ind_comb[:, 3].∈[StatsBase.sample(ind_comb[:, 3], batchSpace * batchTime, replace=false)], :]

  return sel_inds[:, 1], sel_inds[:, 2], sel_inds[:, 3]

end # 1d without buffer