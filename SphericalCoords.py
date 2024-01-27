import math
import numpy as np

deg_rads = 180/math.pi

# ====================================

def sphere2cart(coords:tuple):
    # convert spherical coordinates (coords=[lat,long]) to cartesian coordinates
    # the z-axis is the north pole
    # the x-axis is (0,0), the equator at the Greenwich meridian, south of Ghana, West Africa
    # the y-axis is (0,90 E), the equator at 90E, between Sri Lanka and Singapore
    # This defines a right-handed Euclidean coordinate system
    # assert len(coords)==2
    clong = math.cos(coords[0])
    return clong * math.cos(coords[1]), clong * math.sin(coords[1]), math.sin(coords[0])


# ====================================

def cart2sphere(cart:tuple):
    # convert cartesian coordinates (x,y,z) to spherical (lat, long)
    #assert len(cart)==3

    x = cart[0]
    y = cart[1]
    z = cart[2] 

#    assert math.isclose(x * x + y * y + z * z, 1.0)
#    r = math.sqrt(x * x + y * y + z * z)

#    latr = math.asin(z / r)
    latr = math.asin(z)
#    rcoslat = r * math.cos(latr)
    rcoslat = math.cos(latr)
    if y > 0:
#        longr = math.acos(max(min(1,x / rcoslat),-1))
        if abs(x) < abs(rcoslat):
            longr = math.acos( x / rcoslat )
        elif x<0:
            longr = math.pi
        else:
            longr = 0
    if y < 0:
#        longr = -math.acos(max(min(1, x / rcoslat),-1))
        if abs(x) < abs(rcoslat):
            longr = -math.acos( x / rcoslat )
        elif x < 0:
            longr = 0
        else:
            longr = math.pi
    if y == 0:
        if x > 0:
            longr = 0
        else:
            longr = math.pi

    return (latr, longr)

# ====================================

def cart2sphere_chk(cart:tuple):
    # convert cartesian coordinates (x,y,z) to spherical (lat, long)
    assert len(cart)==3

    x = cart[0]
    y = cart[1]
    z = cart[2] 

    assert math.isclose(x * x + y * y + z * z, 1.0)
    r = math.sqrt(x * x + y * y + z * z)

    latr = math.asin(z / r)
#    latr = math.asin(z)
    rcoslat = r * math.cos(latr)
#    rcoslat = math.cos(latr)
    if y > 0:
        longr = math.acos(max(min(1,x / rcoslat),-1))
    if y < 0:
        longr = -math.acos(max(min(1, x / rcoslat),-1))
    if y == 0:
        if x > 0:
            longr = 0
        else:
            longr = math.pi

    return (latr, longr)

# ====================================

def rad2deg(latlong:tuple):
    # convert spherical (lat, long) coordinates from radians to degrees
    assert len(latlong)==2
    return (latlong[0]*deg_rads, latlong[1]*deg_rads)

def deg2rad(latlong:tuple):
    # convert spherical (lat, long) coordinates from degrees to radians
    assert len(latlong)==2
    return (latlong[0]*deg_rads, latlong[1]*deg_rads)

# ====================================

def sphere_remap(coords:list, rot_N, lat_N, long_N, check=False):
    # finds new spherical coordinates, rotating around north pole by rot_N
    # then relocating the north pole to (lat_N, long_N)
    coords_out = []
    rmat = np.array([[ math.sin(long_N)*math.sin(rot_N)+math.cos(long_N)*math.sin(lat_N)*math.cos(rot_N), -math.sin(long_N)*math.cos(rot_N)+math.cos(long_N)*math.sin(lat_N)*math.sin(rot_N),   math.cos(long_N)*math.cos(lat_N)],
                        [-math.cos(long_N)*math.sin(rot_N)+math.sin(long_N)*math.sin(lat_N)*math.cos(rot_N),  math.cos(long_N)*math.cos(rot_N)+math.sin(long_N)*math.sin(lat_N)*math.sin(rot_N),   math.sin(long_N)*math.cos(lat_N)],
                        [                                                  -math.cos(lat_N)*math.cos(rot_N),                                                   -math.cos(lat_N)*math.sin(rot_N),                    math.sin(lat_N)]])
    if check:
        assert(math.isclose(np.linalg.det(rmat),1))
    for i in range(0, len(coords)):
#        (latr, longr) = coords[i]
        cart = np.array(sphere2cart(coords[i])).transpose()
        
#        print(rmat)
#        print(np.linalg.det(rmat))
        cart_new = np.matmul(rmat,cart)
#        print(cart_new)
        coord_new = cart2sphere(cart_new)
#        print(coord_new)
        coords_out.append(coord_new)
    return coords_out

# ====================================

def sphere_remap_old2(coords:list, rot_z, rot_x):
    # finds new spherical coordinates, rotating around the x axis, then z axis
    coords_out = []
    for i in range(0, len(coords)):
#        (latr, longr) = coords[i]
        cart = np.array(sphere2cart(coords[i])).transpose()
        rmat = np.array([[math.cos(rot_z),               -math.sin(rot_z),                  0],
                        [math.cos(rot_x)*math.sin(rot_z), math.cos(rot_x)*math.cos(rot_z), -math.sin(rot_x)],
                        [math.sin(rot_x)*math.sin(rot_z), math.sin(rot_x)*math.cos(rot_z),  math.cos(rot_x)]])
#        print(np.linalg.det(rmat))
        assert(math.isclose(np.linalg.det(rmat),1))
        cart_new = np.matmul(rmat,cart)
#        print(cart_new)
        coord_new = cart2sphere(cart_new)
#        print(coord_new)
        coords_out.append(coord_new)
    return coords_out

# ====================================

def sphere_remap_old(coords, neworigin, deltaonly=False):
    # finds new spherical coordinates, which maps (0,0) to neworigin
    coords_out = []
    for i in range(0, len(coords)):
        (latr, longr) = coords[i]
        (lato, longo) = neworigin
        latrm = math.asin(-math.cos(longr - longo) * math.cos(latr) * math.sin(lato) + math.sin(latr) * math.cos(lato))
        if math.sin(longr - longo) / math.cos(latrm) > 0:
            ysign = 1
        elif math.sin(longr - longo) / math.cos(latrm) < 0:
            ysign = -1
        else:
            ysign = 0
        t = (math.cos(longr - longo) * math.cos(latr) * math.cos(lato) + math.sin(latr) * math.sin(lato)) / math.cos(
            latrm)
        longrm = ysign * math.acos(max(-1., min(t, 1.)))
        if deltaonly:
            latrm = latrm+lato
            longrm = longrm + longo
        coords_out.append((latrm, longrm))
    return coords_out

# ====================================

def two_norm(v:tuple):
    # Euclidean norm of a vector input
    return np.sqrt(np.dot(v,v))

# ====================================

def project2triang(latlong:tuple, latlong1:tuple, latlong2:tuple, latlong3:tuple, check=False):
    # given three points on a sphere in (lat,long) coordinates
    # projects a third point in the sphere onto the triangle
    #       in local planar coordinates 
    #           where the the first two points of the triangle define the direction of the x-axis, and
    #           the third point defines the positive direction of the y-axis 
    cart = sphere2cart(latlong)
    cart1 = sphere2cart(latlong1)
    cart2 = sphere2cart(latlong2)
    cart3 = sphere2cart(latlong3)
    matrix = np.array([     [cart1[0], cart1[1], cart1[2]],
                            [cart2[0], cart2[1], cart2[2]],
                            [cart3[0], cart3[1], cart3[2]] ])
    ones = np.array([[1], [1], [1]])
    rmat = np.matmul( np.linalg.inv(matrix), ones)
    r = (rmat[0][0], rmat[1][0], rmat[2][0])
    p = 1/np.dot(r,cart)
    cartout = tuple(p*cart[i] for i in range(3))
    #checks
    #print(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])
    #print(u[0]*cartout[0]+u[1]*cartout[1]+u[2]*cartout[2])
    #print(ilr)
    if check:
        ilr = 1/two_norm(r)
        u = tuple(ilr*r[i] for i in range(3))
        assert math.isclose(ilr, np.dot(cartout,u))
        assert math.isclose(two_norm(u), 1)

    return cartout, cart1, cart2, cart3


# ====================================
      
def triang2xy(cart, cart1, cart2, cart3):
    #coplanar coordinates in 3d mapped to 2d
    # cart1 to (0,0)
    # cart2 to (dist(cart2-cart1),0)
    cartout1 = (0,0)
    dimrange = range(len(cart1))
    v12 = tuple(cart2[i]-cart1[i] for i in dimrange)
#    d12 = np.sqrt(sum([v12[i]*v12[i] for i in dimrange]))
    d12 = two_norm(v12)
    u12 = tuple(v12[i]/d12 for i in dimrange)

    v13 = tuple(cart3[i]-cart1[i] for i in dimrange)
#    x3 = sum([v13[i]*u12[i] for i in dimrange])
    x3 = np.dot(v13,u12)
    y3v = tuple(v13[i]-x3*u12[i] for i in dimrange)
#    y3 = np.sqrt(sum([y3v[i]*y3v[i] for i in dimrange]))
    y3 = two_norm(y3v)

    v1c = tuple(cart[i]-cart1[i] for i in dimrange)
#    xc = sum([v1c[i]*u12[i] for i in dimrange])
    xc = np.dot(v1c, u12)
    ycv = tuple(v1c[i]-xc*u12[i] for i in dimrange)
#    yc = np.sqrt(sum([ycv[i]*ycv[i] for i in dimrange]))
    yc = two_norm(ycv)
    
    cartout2 = (d12,0)
    cartout3 = (x3,y3)
    cartout  = (xc, yc)

    return cartout, cartout1, cartout2, cartout3



# ====================================

def projectlist2triang(latlong:list, latlong1:tuple, latlong2:tuple, latlong3:tuple, check=False):
    # given three points on a sphere in (lat,long) coordinates
    # projects a third point in the sphere onto the triangle
    #       in local planar coordinates 
    #           where the the first two points of the triangle define the direction of the x-axis, and
    #           the third point defines the positive direction of the y-axis 
    cart1 = sphere2cart(latlong1)
    cart2 = sphere2cart(latlong2)
    cart3 = sphere2cart(latlong3)
    matrix = np.array([     [cart1[0], cart1[1], cart1[2]],
                            [cart2[0], cart2[1], cart2[2]],
                            [cart3[0], cart3[1], cart3[2]] ])
    ones = np.array([[1], [1], [1]])
    rmat = np.matmul( np.linalg.inv(matrix), ones)
    r = (rmat[0][0], rmat[1][0], rmat[2][0])
    if check:
        ilr = 1/two_norm(r)
        u = tuple(ilr*r[i] for i in range(3))
        assert math.isclose(two_norm(u), 1)

    cart_out = []
    for point in latlong:
        cart = sphere2cart(point)
        p = 1/np.dot(r,cart)
        cart_out = cart_out + [tuple(p*cart[i] for i in range(3))]
        #checks
        #print(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])
        #print(u[0]*cartout[0]+u[1]*cartout[1]+u[2]*cartout[2])
        #print(ilr)
        if check:
            pointout = tuple(p*cart[i] for i in range(3))
            assert math.isclose(ilr, np.dot(pointout,u))

    return cart_out, cart1, cart2, cart3


# ====================================
      
def trianglist2xy(cart:list, cart1:tuple, cart2:tuple, cart3:tuple):
    #coplanar coordinates in 3d mapped to 2d
    # cart1 to (0,0)
    # cart2 to (dist(cart2-cart1),0)
    cart_out1 = (0,0)
    dimrange = range(len(cart1))
    v12 = tuple(cart2[i]-cart1[i] for i in dimrange)
#    d12 = np.sqrt(sum([v12[i]*v12[i] for i in dimrange]))
    d12 = two_norm(v12)
    u12 = tuple(v12[i]/d12 for i in dimrange)

    v13 = tuple(cart3[i]-cart1[i] for i in dimrange)
#    x3 = sum([v13[i]*u12[i] for i in dimrange])
    x3 = np.dot(v13,u12)
    y3v = tuple(v13[i]-x3*u12[i] for i in dimrange)
#    y3 = np.sqrt(sum([y3v[i]*y3v[i] for i in dimrange]))
    y3 = two_norm(y3v)

    cart_out2 = (d12,0)
    cart_out3 = (x3,y3)

    cart_out = []
    for c in cart:
        v1c = tuple(c[i]-cart1[i] for i in dimrange)
#    xc = sum([v1c[i]*u12[i] for i in dimrange])
        xc = np.dot(v1c, u12)
        ycv = tuple(v1c[i]-xc*u12[i] for i in dimrange)
#    yc = np.sqrt(sum([ycv[i]*ycv[i] for i in dimrange]))
        yc = two_norm(ycv)
        cart_out = cart_out + [(xc, yc)]
#        c_out  = (xc, yc)


    return cart_out, cart_out1, cart_out2, cart_out3



# ====================================

def gc_dist(coord1, coord2):
    # calculates the great circle distance between two points on a sphere
    #       given in polar co-ordinates
    lat1 = coord1[0]
#    long1 = coord1[1]
    lat2 = coord2[0]
#    long2 = coord2[1]

#    dlong = abs(long2 - long1)
#    dlat = abs(lat2 - lat1)

    dist = math.acos(math.sin(lat1)*math.sin(lat2)+
                     math.cos(lat1)*math.cos(lat2)*math.cos(coord2[1]-coord1[1]))
#    print(rad2deg(coord1), rad2deg(coord2), math.sin(lat1)*math.sin(lat2), math.cos(lat1)*math.cos(lat2)*math.cos(long2-long1))

    return dist


# ====================================

def find_closest_3(latlong, vertices, coords_df, verbose=False):
    # finds the closest three vertices to the coordinates latlong
    # vertices is the list of vertices and coords_df a dictionary
    # of vertices to their spherical coordinates
    distances = [gc_dist(latlong, coords_df[v]) for v in vertices]

    dist_list = [(vertices[i], distances[i]) for i in range(0, len(vertices))]

    for j in range(0, len(vertices) - 1):
        for i in range(0, len(vertices) - 1):
            if dist_list[i][1] > dist_list[i + 1][1]:
                temp = dist_list[i]
                dist_list[i] = dist_list[i + 1]
                dist_list[i + 1] = temp

    dist_list2 = sorted(dist_list, key=lambda tup: tup[1])
    if dist_list == dist_list2:
        if verbose:
            print('Checked')
    else:
        raise

    f = [dist_list2[0][0], dist_list2[1][0], dist_list2[2][0]]

    if verbose:
        # for e in dist_list:
        #     print(e[0], rad2deg(coords_df[e[0]]), e[1]*deg_rads)
        print(dist_list)
        print()
        print()
        #print(f)
        #print([dist_list[0][1], dist_list[1][1], dist_list[2][1]])

    return f

# ====================================

def mod_arg(argin):
    # maps real input argument (radians) to output in the range (-pi,pi]
    pi=math.pi
    
    argout = argin
    if argout > pi:
        argout = argin - 2*pi
    if argout <= -pi:
        argout = argin + 2*pi
    return argout


# ====================================

def weighted_mean(coord_list:list, weights:list):
    # finds weighted mean of spherical coordinates with given weights
    assert len(coord_list)==len(weights)

    carts = [sphere2cart(coord) for coord in coord_list]
    xs = [cart[0] for cart in carts]
    ys = [cart[1] for cart in carts]
    zs = [cart[2] for cart in carts]
    wsum = sum(weights)
#    x_av = np.average(xs, weights)
#    y_av = np.average(ys, weights)
#    z_av = np.average(zs, weights)
    x_av = sum([x*w for x,w in zip(xs,weights)])/wsum
    y_av = sum([y*w for y,w in zip(ys,weights)])/wsum
    z_av = sum([z*w for z,w in zip(zs,weights)])/wsum

    return cart2sphere((x_av, y_av, z_av))
    