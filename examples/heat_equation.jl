

# using Plots

# α = 1
# Δx = 0.5
# Δy = 0.5
# Δt = 0.01

# (4*Δt*α)/(Δx^2 * Δy^2) ≤ 1

# x = range(0, 10, step = Δx)
# y = range(0, 10, step = Δy)
# t = range(0, 10, step = Δt)


function U_start(x, y, Nx, Ny, ugen)
    U = zeros(Nx, Ny)
    U[1,:] .= utop
    U[Ny,:] .= ubottom
    U[:,1] .= uleft
    U[:,Nx] .= uright

    U .= [ugen(U, i, j) for i in 1:Nx, j in 1:Ny]

    return U
end

function heat_gen(x, y, t, α, utop, uleft, uright, ubottom, ugen, Q)
    
    Nx = length(x)
    Ny = length(y)
    Nt = length(t)

    Δx = x[2] - x[1]
    Δy = y[2] - y[1]
    Δt = t[2] - t[1]

    Ustore = Array{Float64}(undef, Nx, Ny, Nt)

    Ucurr = U_start(x, y, Nx, Ny, ugen)
    U = Ucurr
    Ustore[:,:,1] = Ucurr
    γ = Δt * α

    for t in 2:Nt

        # Ucurr[2:(Nx-1), 2:(Ny-1)] = [U[i,j] + γ * ((U[i-1, j] - 2*U[i, j] + U[i+1, j])/Δx^2 + (U[i, j-1] - 2*U[i, j] + U[i, j+1])/Δy^2) for i in 2:(Nx-1), j in 2:(Ny-1)]
        Ucurr[2:(Nx-1), 2:(Ny-1)] = [U[i,j] + Δt * (α * ((U[i-1, j] - 2*U[i, j] + U[i+1, j])/Δx^2 + (U[i, j-1] - 2*U[i, j] + U[i, j+1])/Δy^2) + Q(U, i, j) ) for i in 2:(Nx-1), j in 2:(Ny-1)]
        
        Ucurr[1,:] = Ucurr[2,:]
        Ucurr[end,:] = Ucurr[end-1,:]
        Ucurr[:,1] = Ucurr[:,2]
        Ucurr[:,end] = Ucurr[:,end-1]
        
        # Ucurr[1,:] .= utop
        # Ucurr[Ny,:] .= ubottom
        # Ucurr[:,1] .= uleft
        # Ucurr[:,Nx] .= uright

        Ustore[:,:,t] = Ucurr
        U = Ucurr

    end

    return Ustore

end