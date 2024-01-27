import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#import os

import Icosohedron as ih
import SphericalCoords as sp

pi=math.pi
el=1/math.sin(2*pi/5)
elt=math.sqrt(3)/math.sin(2*pi/5)/2

#=======================================

def find_rot_trans(face_order, face_dict_ordered, initialr):
# for a given ordering of the faces for definining the planar net
#   and a rotation of the first face
# calculate the translation and rotation of each face one by one 


    face_rotate = {}
    face_translate = {}
    face0 = face_order[0][0]
    face_rotate[face0] = initialr
    face_translate[face0] = (0,0)
    for i in face_order:
        face_prev = i[0]
        face_new = i[1]
        fr_prev = face_rotate[face_prev]
        ft_prev = face_translate[face_prev]
        vlist_prev = face_dict_ordered[face_prev]
        vlist_new = face_dict_ordered[face_new]
        edges_prev = {(vlist_prev[0], vlist_prev[1]): -pi/2+fr_prev,
                      (vlist_prev[1], vlist_prev[2]): pi/6+fr_prev,
                      (vlist_prev[2], vlist_prev[0]): 5*pi/6+fr_prev}
        for j in edges_prev.keys():
            edges_prev[j] = sp.mod_arg(edges_prev[j])

        edges_new = {(vlist_new[1], vlist_new[0]):-pi/2,
                      (vlist_new[2], vlist_new[1]): pi/6,
                      (vlist_new[0], vlist_new[2]): 5*pi/6}


        new_dir = np.nan
        for j in edges_new.keys(): 
            if j in edges_prev.keys():
                new_dir = edges_prev[j]
#                print(face_prev, face_new, new_dir)
                fr_new = sp.mod_arg(sp.mod_arg(new_dir+pi - edges_new[j]))
                ft_new = (ft_prev[0] + 2/3*elt*np.cos(new_dir), 
                          ft_prev[1] + 2/3*elt*np.sin(new_dir) ) 
        try:
            assert not(np.isnan(new_dir))
        except:
            print(str(face_prev)+' and '+str(face_new)+' do not share an edge.')

        face_rotate[face_new] = fr_new
        face_translate[face_new] = ft_new

    face_transform={f:(face_translate[f], face_rotate[f]) for f in face_translate.keys()}

    return face_transform, face_rotate, face_translate

#=======================================

def mapcities2icosofaces(citydata_df, face_dict_ordered, vdm, quick=False):
    # Given a list of points in polar coordinates on the sphere in citydata_df and
    #   A map face_dict_ordered of faces to vertices and
    #   A map vdm of vertices to polar coordinates
    # Determine which face each point belongs to and
    #   The location of each point on each face, 
    #       in planar coordinates local to each face with 
    #           the first two vertices of the face defining the x-axis, and 
    #           the third vertex defining the positive direction of the y-axis.

    face_dict = {k:set(face_dict_ordered[k]) for k in face_dict_ordered}

    citydata2_df = citydata_df.copy()
    citydata2_df.loc[:,'face'] = 'f_'

    for i in citydata_df.index.to_list():
        latlong = (citydata_df.at[i,'latr'], citydata_df.at[i,'longr']) 
        fid = ih.find_face(latlong, face_dict, vdm)
        citydata2_df.at[i,'face'] = fid

    if not quick:
        fid = ih.find_face((0,0), face_dict, vdm)

        triang=[sp.rad2deg(vdm[v]) for v in face_dict[fid]]
        #print(triang)

        c, c1, c2, c3 = sp.projectlist2triang([(0,0)],triang[0],triang[1],triang[2])

        cxy, cxy1, cxy2, cxy3 = sp.trianglist2xy(c, c1, c2, c3)
        citydata2_df.loc[:,'xy'] = [(np.nan, np.nan) for i in citydata2_df.index.to_list()]

#    maxdist=0
    
        for fid in face_dict_ordered.keys():
            triang=[vdm[v] for v in face_dict_ordered[fid]]
            cdx_df = citydata2_df.loc[citydata2_df.loc[:,'face']==fid,:]
            latlong_list = [(cdx_df.at[r,'latr'], cdx_df.at[r,'longr']) for r in cdx_df.index.to_list()]

    #        c_list, c1, c2, c3 = sp.projectlist2triang(latlong_list,triang[0],triang[1],triang[2], check=True)
            c_list, c1, c2, c3 = sp.projectlist2triang(latlong_list,triang[0],triang[1],triang[2])
            cxy_list, cxy1, cxy2, cxy3 = sp.trianglist2xy(c_list, c1, c2, c3)

            try:
                rt = 1e-6
                assert math.isclose(cxy1[0], 0,    rel_tol=rt)
                assert math.isclose(cxy1[1], 0,    rel_tol=rt)
                assert math.isclose(cxy2[0], el,   rel_tol=rt)
                assert math.isclose(cxy2[1], 0,    rel_tol=rt)
                assert math.isclose(cxy3[0], el/2, rel_tol=rt)
                assert math.isclose(cxy3[1], elt,  rel_tol=rt)
            except:
                print('Face mapped to wrong coordinates')
                print(cxy1[0], cxy1[1], cxy2[0]-el, cxy2[1], cxy3[0]-el/2, cxy3[1]-elt)

            cxy_t_list = [(cxy[0]-el/2, cxy[1]-elt/3) for cxy in cxy_list]

            citydata2_df.loc[citydata2_df.loc[:,'face']==fid, 'xy'] = pd.Series(cxy_t_list).values
        
    #    print(citydata2_df.columns.tolist())

    return citydata2_df

#=======================================

def mapcities2icosofaces_quick(citydata_df, face_dict_ordered, vdm):
    # Given a list of points in polar coordinates on the sphere in citydata_df and
    #   A map face_dict_ordered of faces to vertices and
    #   A map vdm of vertices to polar coordinates
    # Determine which face each point belongs to and
    #   The location of each point on each face, 
    #       in planar coordinates local to each face with 
    #           the first two vertices of the face defining the x-axis, and 
    #           the third vertex defining the positive direction of the y-axis.

    face_dict = {k:set(face_dict_ordered[k]) for k in face_dict_ordered}

    citydata2_df = citydata_df.copy()
    citydata2_df.loc[:,'face'] = 'f_'

    for i in citydata_df.index.to_list():
        latlong = (citydata_df.at[i,'latr'], citydata_df.at[i,'longr']) 
        fid = ih.find_face(latlong, face_dict, vdm)
        citydata2_df.at[i,'face'] = fid
    

    return citydata2_df

#=======================================


def plot_net_calc(citydata_df:pd.DataFrame, face_dict:dict, face_transform:dict, faces:list=[]):
    # Given a list of data points citydata_df and 
    # their mapping to one of a set of faces and local coordinates on the face
    #   and a mapping face_transform of each faces local planar coordinates 
    #       to cartesian coordinates on the polyhedral net
    # return the data points mapping to the cartesian plane
    # Also return the  

    citydata2_df = citydata_df.copy()
    citydata2_df['xy_net']=[(np.nan, np.nan) for c in citydata2_df.index.to_list()]
    icoso_net_dict = {}

    face_translate = {f:face_transform[f][0] for f in face_transform.keys()}
    face_rotate = {f:face_transform[f][1] for f in face_transform.keys()}

    if faces == []:
        faces = list(face_dict.keys()) 
    for fid in face_dict.keys():
        if fid in faces:
            clist = [c for c in citydata_df.index.to_list() if citydata_df.at[c, 'face'] == fid]
            theta = face_rotate[fid]
            ct = math.cos(theta)
            st = math.sin(theta)
            xt = face_translate[fid][0]
            yt = face_translate[fid][1]
            xy = citydata_df.loc[clist, 'xy']
            x = [ xy.at[c][0]*ct-xy.at[c][1]*st+xt for c in clist ]
            y = [ xy.at[c][0]*st+xy.at[c][1]*ct+yt for c in clist ]
            for i in range(0,len(clist)):
                citydata2_df.at[clist[i],'xy_net'] = (x[i],y[i])
            x=[ -el/2*ct+elt/3*st+xt,  el/2*ct+elt/3*st+xt, -2*elt/3*st+xt,-el/2*ct+elt/3*st+xt ]
            y=[ -el/2*st-elt/3*ct+yt,  el/2*st-elt/3*ct+yt,  2*elt/3*ct+yt,-el/2*st-elt/3*ct+yt ]
            icoso_net_dict[fid] = (x,y)
    icoso_net = pd.DataFrame.from_dict(icoso_net_dict, orient='index', columns=['edges_x', 'edges_y'])

    citydata2_df.loc[:,"pop_cat"]=citydata2_df.loc[:,"population"]*1.0
    for row in citydata2_df.index.tolist():
        citydata2_df.at[row,"pop_cat"] = max((np.log10(citydata2_df.at[row,"population"]+1)*1).astype(float)-3.5,0.0)

    print(min(citydata2_df.loc[:,"pop_cat"]), max(citydata2_df.loc[:,"pop_cat"]))
    return citydata2_df, icoso_net

#=======================================

def plot_net(citydata_df:pd.DataFrame, icoso_net:pd.DataFrame, size_dep=False, faces:list=[], figname='', pause_t=2.0):
    # plot citydata_df coordinates on the cartesian plane
    # plot polyhedron net edges
    clist1 = citydata_df.index.to_list()
#    clist = [c for c in clist1 if citydata_df.at[c,"pop_cat"]==0]
    if len(faces)>0:
        clist = [c for c in clist1 if citydata_df.at[c,"face"] in faces]
    else:
        clist = clist1
    x = [citydata_df.at[c,'xy_net'][0] for c in clist]
    y = [citydata_df.at[c,'xy_net'][1] for c in clist]

    if not(size_dep):
        plt.plot(x, y, 'k,',markersize=1)
    else:
        ss =[citydata_df.at[c,'pop_cat'] for c in clist]
        sss=[(e*e+1)/10 for e in ss]
        plt.scatter(x, y, marker='o',facecolors='none',edgecolors='k',s=sss,linewidths=0.05)

    if faces==[]:
        f_list = icoso_net.index.to_list()
    else:
        f_list = faces

    for f in f_list:
        plt.plot(icoso_net.at[f,'edges_x'], icoso_net.at[f,'edges_y'], 'r:', linewidth=0.1)
    plt.axis('off')
    plt.draw()
    plt.show()
    if len(clist)>100000:
        dpi_v = 1800
    else:
        dpi_v = 300
        if size_dep:
            dpi_v = 2400
    if figname=='':
        figname = 'figure'

    plt.savefig('Figures\\'+figname+'-t.png',format='png', bbox_inches='tight', dpi=dpi_v, transparent=True)
    plt.savefig('Figures\\'+figname+'.png',  format='png', bbox_inches='tight', dpi=dpi_v)
    plt.savefig('Figures\\'+figname+'.svg',  format='svg', bbox_inches='tight')
    plt.show()
    plt.pause(pause_t)
    plt.clf()
    plt.close()

#=======================================

def countcities(citydata_df, agg="count", verbose=False):
    if agg=="sum":
        citycount = citydata_df.loc[:,["face","population"]].groupby("face").sum().astype(float)
    else:
        # else agg="count"
        citycount = citydata_df.loc[:,["face","population"]].groupby("face").count()
    pops = citycount["population"].to_list()

    logpops = [np.log(p) for p in pops]
    slope = 19.8/(max(logpops)-min(logpops))
    intercept = min(logpops)+0.4/slope
    citycount.loc[:,'logs']    =          (  np.log( citycount.loc[:,'population'] )-intercept  )*slope
    citycount.loc[:,'intlogs'] = round(   (  np.log( citycount.loc[:,'population'] )-intercept  )*slope, 0   ).astype(int)
    citycount = citycount.sort_values("logs")
    points    = citycount["logs"].to_list()
    points_i  = citycount["intlogs"].to_list()
    ll = len(points)
    l_20 = list(range(0, ll))
    l_19 = list(range(1, ll))
    sl = (ll-points[1])/(ll-1)
    l_y = [points[1]-sl]+list(np.arange(points[1], ll, sl))
    citycount.loc[:,'order'] = l_20
    L1      =              sum(  [    abs(a-b) for a,b in zip(points,   l_20) ]  )/ll
    L2      = math.sqrt(   sum(  [ (a-b)*(a-b) for a,b in zip(points,   l_20) ]  )/ll   )
    Linf    =              max(  [    abs(a-b) for a,b in zip(points,   l_20) ]  )
    L1_i    =              sum(  [    abs(a-b) for a,b in zip(points_i, l_20) ]  )/ll
    L2_i    = math.sqrt(   sum(  [ (a-b)*(a-b) for a,b in zip(points_i, l_20) ]  )/ll   )
    Linf_i  =              max(  [    abs(a-b) for a,b in zip(points_i, l_20) ]  )
    diffs   = [points[i]-points[i-1] for i in l_19]
    mind  = min(diffs)
    if mind > 0:
        La_1 =               sum( [1/diff        for diff in diffs] )/ll
        La_2 =   math.sqrt(  sum( [1/(diff*diff) for diff in diffs] )/ll  )
        La_inf =                1/mind
    else:
        La_1 = 0
        La_2 = 0
        La_inf = 100
    metrics = {"L1":L1,     "L2":L2,     "Linf":Linf,     "L1_i":L1_i, "L2_i":L2_i, "Linf_i":Linf_i, 
               "La_1":La_1, "La_2":La_2, "La_inf":La_inf}
    if verbose:
        print()
        print("Ordered list of city/population count logarithms by face, rescaled from 0 to 19")
        print(points)
        print()
        plt.ion()
        plt.bar(l_20,points)
        plt.plot(l_20,l_y,'k-')
        plt.draw()
        plt.show()
        plt.savefig('Figures\\counts.png',format='png')
        plt.show()
        plt.clf()
        plt.close()

    return citycount, metrics 
