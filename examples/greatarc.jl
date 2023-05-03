


function greatarc(lat1,lon1,lat2,lon2)

    r=6378.136            # equatorial radius of earth in km
    lat1=lat1*2*pi/360
    lat2=lat2*2*pi/360
    lon1=lon1*2*pi/360
    lon2=lon2*2*pi/360

    phi1=pi/2 - lat1
    phi2=pi/2 - lat2
    th1=2*pi-lon1
    th2=2*pi-lon2
    p=r

    ca1=sin(phi1)*cos(th1)
    ca2=sin(phi2)*cos(th2)
    
    cb1=sin(phi1)*sin(th1)
    cb2=sin(phi2)*sin(th2)
    
    cc1=cos(phi1)
    cc2=cos(phi2)
    
    ang=acos(ca1*ca2 + cb1*cb2 + cc1*cc2)
    
    dist=ang * r

    return ang, dist

end
