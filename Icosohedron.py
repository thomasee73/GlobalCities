import math

import SphericalCoords as sp
#import random

pi=math.pi
atanhalf=math.atan(1/2)

def vertices():
    # returns a dictionary of the coordinates of the vertices of an icohedron
    # in spherical coordinates (latitude, longitude) in radians
    vert_dict= {'v0':  ( pi/2,      0      ),
                'v1':  ( atanhalf, -4*pi/5 ),
                'v2':  ( atanhalf, -2*pi/5 ),
                'v3':  ( atanhalf,  0      ),
                'v4':  ( atanhalf,  2*pi/5 ),
                'v5':  ( atanhalf,  4*pi/5 ),
                'v6':  (-atanhalf,  pi     ),
                'v7':  (-atanhalf, -3*pi/5 ),
                'v8':  (-atanhalf, -pi/5   ),
                'v9':  (-atanhalf,  pi/5   ),
                'v10': (-atanhalf,  3*pi/5 ),
                'v11': (-pi/2,    0        )
                }
    return vert_dict

def vertices_mod(vert_dict, rot_z, lat_N, long_N):
    # remaps a dictionary of polar coordinates to new coordinates
    # first around north pole by rot_z, 
    # then maps tne north pole to lat_N, long_N, first via the Greenwich meridian by lat_N 

    k = list(vert_dict.keys())
    # print(k)
    coords = [vert_dict[i] for i in k]
#    coords_new = sp.sphere_remap(coords, rot_z, lat_N, long_N, check=True)
    coords_new = sp.sphere_remap(coords, rot_z, lat_N, long_N)
    vert_out = { k[i]:coords_new[i] for i in range(len(k))}

    return vert_out

def find_face(latlong, face_dict, vertex_coords):
    # Given a point in spherical coordinates latlong
    # A dictionary face_dict of face names to vertex names and
    # A dictionary of polar vertex coordinates
    # finds the face the point belongs to

    vlist = sp.find_closest_3(latlong, list(vertex_coords.keys()), vertex_coords)
    face_id=''
    for f in face_dict:
        if set(vlist) == face_dict[f]:
            face_id = f

    return face_id



def define_isocofaces():
    # provides a map of faces to vertices
    # the vertices are ordered in the sense that they define
    # an outward facing direction 
    face_dict_ordered = {   'f0':['v0' , 'v1','v2'],
                            'f1':['v0' , 'v2','v3'],
                            'f2':['v0' , 'v3','v4'],
                            'f3':['v0' , 'v4','v5'],
                            'f4':['v0' , 'v5','v1'],

                            'f5':['v7' , 'v2','v1'],
                            'f6':['v8' , 'v3','v2'],
                            'f7':['v9' , 'v4','v3'],
                            'f8':['v10', 'v5','v4'],
                            'f9':['v6' , 'v1','v5'],

                            'f10':['v1', 'v6','v7' ],
                            'f11':['v2', 'v7','v8' ],
                            'f12':['v3', 'v8','v9' ],
                            'f13':['v4', 'v9','v10'],
                            'f14':['v5', 'v10','v6'],

                            'f15':['v11', 'v7', 'v6'],
                            'f16':['v11', 'v8', 'v7'],
                            'f17':['v11', 'v9', 'v8'],
                            'f18':['v11', 'v10','v9'],
                            'f19':['v11', 'v6', 'v10']
                }
#    print(face_dict_ordered)

    return face_dict_ordered

def check_remap(face_dict_ordered, vdm):
    # check that each face in face_dict_ordered is an equilateral triangle
    # where vdm is a dictionary of polar coordinates

    for face in face_dict_ordered:
        vlist= face_dict_ordered[face]
    #    print(face, vlist, sp.gc_dist(vdm[vlist[0]],vdm[vlist[1]]), sp.gc_dist(vdm[vlist[1]],vdm[vlist[2]]), sp.gc_dist(vdm[vlist[2]],vdm[vlist[0]]))
    #    print(sp.gc_dist(vdm[vlist[0]],vdm[vlist[1]]), sp.gc_dist(vdm[vlist[1]],vdm[vlist[2]]), sp.gc_dist(vdm[vlist[2]],vdm[vlist[0]]))
        try:
            assert math.isclose(sp.gc_dist(vdm[vlist[0]],vdm[vlist[1]]), sp.gc_dist(vdm[vlist[1]],vdm[vlist[2]]),rel_tol=1e-6)
            assert math.isclose(sp.gc_dist(vdm[vlist[1]],vdm[vlist[2]]), sp.gc_dist(vdm[vlist[2]],vdm[vlist[0]]),rel_tol=1e-6)
            assert math.isclose(sp.gc_dist(vdm[vlist[2]],vdm[vlist[0]]), sp.gc_dist(vdm[vlist[0]],vdm[vlist[2]]),rel_tol=1e-6)
        except:
            print(str(face)+' is not equilateral')








