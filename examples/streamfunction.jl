

######################### correction function #########################

function correction_estimate(lat, lon, psi, phi; Lat = [-70, -130], Lon = [30, 50])

    fndXstrt = findall(lon .== Lat[1])
    fndXend = findall(lon .== Lat[2])
    Xint = lon[range(fndXend[1],fndXstrt[1])]

    fndYstrt = findall(lat .== Lon[1])
    fndYend = findall(lat .== Lon[2])
    Yint = lat[range(fndYend[1],fndYstrt[1])]

    # Select region
    Phi = phi[range(fndXend[1],fndXstrt[1]), range(fndYend[1],fndYstrt[1]),1,:]
    Psi = psi[range(fndXend[1],fndXstrt[1]), range(fndYend[1],fndYstrt[1]),1,:]

    # geographic orientation for matrices; i.e., take transpose
    Phi = permutedims(Phi, [2,1,3])
    Psi = permutedims(Psi, [2,1,3])


    # convert to great-arc distances
    XMat = (Xint * ones(1,length(Yint)))'
    YMat = Yint *ones(1,length(Xint))

    # addpath('/Users/wiklec/Dropbox/Mac_PROJECTS/MATLABFILES')
    ny,nx = size(XMat)
    XMatdist = zeros(ny,nx-1)


    # great arc distance, in km
    for j = 1:ny
        for i = 2:nx
            ang, dist = greatarc(YMat[j,i],-XMat[j,i-1],YMat[j,i-1],-XMat[j,i])
            XMatdist[j,i-1]= dist
        end
    end

    XMatdist = XMatdist*1000 #convert to meters

    # note, all distances in the y direction are equal only calculate
    # once, but keep it in matrix if needed (I don't use that here though)
    a, d = greatarc(YMat[1,1],-XMat[1,1],YMat[2,1],-XMat[2,1])
    YMatdist = d*ones(ny-1,nx)*1000
    dy = d   #distance in y direction, km

    # compute streamfunction
    Phi_mu = mean(Phi)
    Phi = Phi .- Phi_mu

    # coriolis parameter
    f = [2 * 7.292e-5 .* sin.(YMat[j,i] * pi/180) for j in 1:ny, i in 1:nx]

    Psi_tmp = (Phi ./ f)
    F = hcat(ones(prod(size(Psi))), reshape(Psi_tmp, :))
    U = reshape(Psi, :)

    correct = inv(F' * F) * F' * U

    return correct

end


######################### streamfunction function #########################


function stream(lat, lon, z, u, v, vo, β; Lat = [-50, -150], Lon = [20, 60], lat_res = 1.5, lon_res = 1.5)

    Δlon = lon[2] - lon[1]
    Δlat = lat[1] - lat[2]

    lon_step = Int(lon_res / Δlon)
    lat_step = Int(lat_res / Δlat)

    fndXstrt = findall(lon .== Lat[1])
    fndXend = findall(lon .== Lat[2])
    Xint = lon[range(fndXend[1], fndXstrt[1], step = lon_step)]
    # Xint = Vector(range(lon[fndXend[1]], lon[fndXstrt[1]], step = lon_res))

    fndYstrt = findall(lat .== Lon[1])
    fndYend = findall(lat .== Lon[2])
    Yint = lat[range(fndYend[1], fndYstrt[1], step = lat_step)]
    # Yint = Vector(range(lat[fndYend[1]], lat[fndYstrt[1]], step = -lat_res))


    # Select region
    z = z[range(fndXend[1], fndXstrt[1], step = lon_step), range(fndYend[1], fndYstrt[1], step = lat_step), :]
    u = u[range(fndXend[1], fndXstrt[1], step = lon_step), range(fndYend[1], fndYstrt[1], step = lat_step), :]
    v = v[range(fndXend[1], fndXstrt[1], step = lon_step), range(fndYend[1], fndYstrt[1], step = lat_step), :]
    vo = vo[range(fndXend[1], fndXstrt[1], step = lon_step), range(fndYend[1], fndYstrt[1], step = lat_step), :]

    # geographic orientation for matrices; i.e., take transpose
    z = permutedims(z, [2,1,3])
    u = permutedims(u, [2,1,3])
    v = permutedims(v, [2,1,3])
    vo = permutedims(vo, [2,1,3])

    # convert to great-arc distances
    XMat = (Xint * ones(1,length(Yint)))'
    YMat = Yint *ones(1,length(Xint))

    # addpath('/Users/wiklec/Dropbox/Mac_PROJECTS/MATLABFILES')
    ny,nx = size(XMat)
    XMatdist = zeros(ny,nx-1)


    # great arc distance, in km
    for j = 1:ny
        for i = 2:nx
            ang, dist = greatarc(YMat[j,i],-XMat[j,i-1],YMat[j,i-1],-XMat[j,i])
            XMatdist[j,i-1] = dist
        end
    end

    XMatdist = XMatdist*1000 #convert to meters

    # note, all distances in the y direction are equal only calculate
    # once, but keep it in matrix if needed (I don't use that here though)
    a, d = greatarc(YMat[1,1],-XMat[1,1],YMat[2,1],-XMat[2,1])
    YMatdist = d*ones(ny-1,nx)*1000
    dy = d   #distance in y direction, km

    # compute streamfunction
    Phi_mu = mean(z)
    Phi = z .- Phi_mu

    # coriolis parameter
    f = 2 * 7.292e-5 .* sin.(YMat .* pi/180)
    fprime = 2 * 7.292e-5 .* cos.(YMat .* pi/180)

    Psi_tmp = (Phi ./ f)
    F = hcat(ones(prod(size(Phi))), reshape(Psi_tmp, :))

    ψ = reshape(F * β, size(z))

    return ψ, XMatdist, YMatdist, Xint, Yint, f, fprime, z, u, v, vo

end


