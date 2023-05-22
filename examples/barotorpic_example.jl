

# http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.DAILY/.Intrinsic/.PressureLevel/


using BayesianDiscovery
const BD = BayesianDiscovery
using Plots
using Random, Distributions
using Kronecker
using DataFrames, DataFramesMeta, Chain
using NetCDF
using Distances
using Dates

include("greatarc.jl")
include("streamfunction.jl")

######################### Load Data From NCEP #########################

# steam function
ncinfo("../BayesianDiscovery/data/geo_winds/geo_1948_1960.nc")
lat = ncread("../BayesianDiscovery/data/geo_winds/geo_1948_1960.nc", "Y")
lon = ncread("../BayesianDiscovery/data/geo_winds/geo_1948_1960.nc", "X")
t_15 = ncread("../BayesianDiscovery/data/geo_winds/geo_2015_2020.nc", "T")

# geopotential height
phi_15 = ncread("../BayesianDiscovery/data/geo_winds/geo_2015_2020.nc", "phi")

# streamfunciton height
psi_15 = ncread("../BayesianDiscovery/data/geo_winds/stream_2015_2020.nc", "psi")

# u wind
u_15 = ncread("../Space_Time_DEtection/data/geo_winds/u_2015_2020.nc", "u")

# v wind
v_15 = ncread("../BayesianDiscovery/data/geo_winds/v_2015_2020.nc", "v")


t = t_15
phi = phi_15[:,:,1,:]
psi = psi_15[:,:,1,:]
u = u_15[:,:,1,:]
v = v_15[:,:,1,:]

startdate = DateTime(1948,01,01) # start date
dts = startdate + Dates.Day.(t) # get date corresponding to each time 
yrs = Dates.year.(dts) # get years to use as identifiers
mths = Dates.month.(dts) # get months to use as identifiers

######################### Load Data ERA5 and calibrate #########################

data_year = "2018"
data_month = "12"
data_path = "../BayesianDiscovery/data/era5/" * string(data_year) * "_" * string(data_month) * ".nc"

dts_inds = findall((yrs .== parse(Int, data_year)) .&& (mths .== parse(Int, data_month)))
β = correction_estimate(lat, lon, psi[:,:,dts_inds], phi[:,:,dts_inds], Lat = [-50, -150], Lon = [20, 60])

# geo function
ncinfo(data_path)
lat = ncread(data_path, "latitude")
lon = ncread(data_path, "longitude")
t = ncread(data_path, "time")
z = ncread(data_path, "z")
v = ncread(data_path, "v")
u = ncread(data_path, "u")
vel = ncread(data_path, "vo")


sf_u = ncgetatt(data_path, "u", "scale_factor")
sf_v = ncgetatt(data_path, "v", "scale_factor")
sf_z = ncgetatt(data_path, "z", "scale_factor")
sf_vo = ncgetatt(data_path, "vo", "scale_factor")

ao_u = ncgetatt(data_path, "u", "add_offset")
ao_v = ncgetatt(data_path, "v", "add_offset")
ao_z = ncgetatt(data_path, "z", "add_offset")
ao_vo = ncgetatt(data_path, "vo", "add_offset")


z = (z .* sf_z) .+ ao_z
z = z ./ 10 # note, ERA5 geopotential is 10 times the NCEP value, don't know why
v = (v .* sf_v) .+ ao_v
u = (u .* sf_u) .+ ao_u
vel = (vel .* sf_vo) .+ ao_vo


ψ, Xdist, Ydist, lon, lat, f, fprime, z, u, v, vel = stream(lat, lon, z, u, v, vel, β, Lat = [-50, -150], Lon = [20, 60], lat_res = 1.5, lon_res = 1.5)

z = reverse(z, dims = 1)
ψ = reverse(ψ, dims = 1)
u = reverse(u, dims = 1)
v = reverse(v, dims = 1)
vel = reverse(vel, dims = 1)



Ny, Nx, Nt = size(ψ)

x_meter = vcat(0, cumsum(mean(Xdist)*ones(length(lon)-1)))
y_meter = vcat(0, cumsum(mean(Ydist)*ones(length(lat)-1)))
Time = cumsum(60*60*ones(Nt))


# reshape data
Y = convert(Array{Float64}, reshape(ψ, :, Nt, 1))
covs = reshape(cat([reverse(fprime, dims = 1)/6.371e6 for i in 1:Nt]..., dims = 3), :, Nt, 1) # Rossby parameter
covs = reshape(cat([reverse(fprime, dims = 1)/mean(Ydist) for i in 1:Nt]..., dims = 3), :, Nt, 1) # Rossby parameter

######################### set up detection #########################


SpaceStep = [x_meter, reverse(y_meter)]
TimeStep = Time
νS = [20, 10] # number of space basis functions
νT = 300 # number of time basis functions
bufferSpace = [1, 1]
bufferTime = 1
batchSpace = 10
batchTime = 10
learning_rate = [1e0]
beta = 0.99

responsenames = ["ΔU_xxt", "ΔU_yyt"]

# response function
function response(ΔU)
    ΔU_xxt = ΔU[1]
    ΔU_yyt = ΔU[2]
    return ΔU_xxt + ΔU_yyt
end


# lambda function
ΛSpaceNames = ["Psi", "Psi_x", "Psi_xx", "Psi_xxx", "Psi_y", "Psi_yy", "Psi_yyy", "Psi_xy", "Psi_xxy", "Psi_xyy"]
ΛTimeNames = ["Phi"]

# semi-reduced library
function Λ(A, Ψ, Φ, X)
  
  ψ = Ψ[1]
  ψ_x = Ψ[2]
  ψ_xx = Ψ[3]
  ψ_xxx = Ψ[4]
  ψ_y = Ψ[5]
  ψ_yy = Ψ[6]
  ψ_yyy = Ψ[7]
  ψ_xy = Ψ[8]
  ψ_xxy = Ψ[9]
  ψ_xyy = Ψ[10]

  ϕ = Φ[1]

  fp = X[1]

  u = A * (ϕ ⊗ ψ)'
  u_x = A * (ϕ ⊗ ψ_x)'
  u_xx = A * (ϕ ⊗ ψ_xx)'
  u_xxx = A * (ϕ ⊗ ψ_xxx)'
  u_y = A * (ϕ ⊗ ψ_y)'
  u_yy = A * (ϕ ⊗ ψ_yy)'
  u_yyy = A * (ϕ ⊗ ψ_yyy)'
  u_xy = A * (ϕ ⊗ ψ_xy)'
  u_xxy = A * (ϕ ⊗ ψ_xxy)'
  u_xyy = A * (ϕ ⊗ ψ_xyy)'

  return [u, u_x, u_xx, u_y, u_yy, u_xy,
            u_x.*u_xxx, u_y.*u_xxx, u_x.*u_yyy, u_y.*u_yyy, u_x.*u_xxy, u_y.*u_xxy, u_x.*u_xyy, u_y.*u_xyy,
            u_x.*fp, u_y.*fp]

end

function ∇Λ(A, Ψ, Φ, X)
  
  ψ = Ψ[1]
  ψ_x = Ψ[2]
  ψ_xx = Ψ[3]
  ψ_xxx = Ψ[4]
  ψ_y = Ψ[5]
  ψ_yy = Ψ[6]
  ψ_yyy = Ψ[7]
  ψ_xy = Ψ[8]
  ψ_xxy = Ψ[9]
  ψ_xyy = Ψ[10]

  ϕ = Φ[1]

  fp = X[1]

  return [BD.dU(A, ϕ, ψ), dU_x(A, ϕ, ψ_x), BD.dU_xx(A, ϕ, ψ_xx), BD.dU_y(A, ϕ, ψ_y), BD.dU_yy(A, ϕ, ψ_yy), BD.dU_xy(A, ϕ, ψ_xy),
            BD.dUxU_xxx(A, ϕ, ψ_x, ψ_xxx), BD.dUxU_xxx(A, ϕ, ψ_y, ψ_xxx), BD.dUxU_xxx(A, ϕ, ψ_x, ψ_yyy), BD.dUxU_xxx(A, ϕ, ψ_y, ψ_yyy),
            BD.dUxU_xxx(A, ϕ, ψ_x, ψ_xxy), BD.dUxU_xxx(A, ϕ, ψ_y, ψ_xxy), BD.dUxU_xxx(A, ϕ, ψ_x, ψ_xyy), BD.dUxU_xxx(A, ϕ, ψ_y, ψ_xyy),
            BD.dU_xX(A, ϕ, ψ_x, fp), BD.dU_xX(A, ϕ, ψ_y, fp)]
end


model, pars, posterior = DEtection(Y, SpaceStep, TimeStep, νS, νT, batchSpace, batchTime, learning_rate, c, Λ, ∇Λ, ΛSpaceNames, ΛTimeNames, nits=10000, response=response, responsenames=responsenames, covariates = covs)


print_equation(["∇U_t"], model, pars, posterior; cutoff_prob=0.48)

post_sum = posterior_summary(model, pars, posterior)

post_sum.M
post_sum.M_hpd
post_sum.gamma
sqrt.(post_sum.ΣZ)
post_sum.ΣU

surf_mean, surf_sd = posterior_surface(model, pars, posterior)

contourf(ψ[:, :, 10], c=:heat, right_margin = 10Plots.mm)
contourf(reshape(surf_mean, Ny, Nx, Nt)[:, :, 10], c=:heat, right_margin = 10Plots.mm)
contourf(ψ[:, :, 10] - reshape(surf_mean, Ny, Nx, Nt)[:, :, 10], c=:bluesreds, right_margin = 10Plots.mm)

contourf(reshape(surf_sd[:, :, 1], Ny, Nx, Nt)[:, :, 10], c=:heat)


plot(posterior.M[1,8,:])
plot(posterior.M[1,9,:])
plot(posterior.M[1,11,:])
plot(posterior.M[1,14,:])
plot(posterior.M[1,16,:])

