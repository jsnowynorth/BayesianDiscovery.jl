
# using Plots
# using FFTW
# using DifferentialEquations

# Burgers's equation
# u_t = -uu_x + νu_xx
# where ν is the viscocity 


function burgers_pde(u, p, t)
    k,nu = p
    deriv = -u .* ifft(1im .* k .* fft(u)) + nu*ifft(-k .^2 .* fft(u))
    return real(deriv)
end


function burgers(;nx = 256, nt = 101, Lx = 8, Lt = 10.0, nu=0.1)


    # Set up grid
    x = range(-Lx, Lx, nx+1)[1:(end-1)]
    dx = x[2] - x[1]
    t = range(0, Lt, nt)
    dt = t[2]-t[1]
    k = 2*pi*fftfreq(nx, nx*dx)

    # Initial condition
    u0 = exp.(-(x.+2).^2)
    # u0 = -sin.(x*pi./8)
    # range_inds = -8 .< x .< 8
    # u0 = u0 .* range_inds
    # u0 = exp.(-0.3*(x.+2).^2)
    # u0 = exp.(-0.02*(x.+10).^2) - exp.(-0.02*(x.-10).^2)
    # u0 = 3*exp.(-0.2*(x.+20).^2) 

    # Solve
    p = (k, nu)
    tspan = (0.0, Lt)
    prob = ODEProblem(burgers_pde, u0, tspan, p)
    sol = solve(prob, dt = dt, Tsit5(), saveat = t)

    return reduce(hcat, sol.u), Vector(x), sol.t
    
end


# U, x, t = burgers()
# contourf(t, x, U, c=:oxy)
