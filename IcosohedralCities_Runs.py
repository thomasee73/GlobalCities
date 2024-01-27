import math
import pandas as pd

import Icosohedron as ih
import IcosohedralCities as ihc

pi=math.pi

def init():
    face_rotate = {'f0':  -4*pi/3, 'f1': -3*pi/3,   'f2': -2*pi/3,  'f3': -1*pi/3, 'f4':   0*pi/3,
                   'f5':  -1*pi/3, 'f6':  0*pi/3,   'f7':  1*pi/3,  'f8':  2*pi/3, 'f9':   3*pi/3, 
                   'f10': -0*pi/3, 'f11': -1*pi/3, 'f12': -2*pi/3, 'f13': -3*pi/3, 'f14':  2*pi/3, 
                   'f15':  3*pi/3, 'f16':  2*pi/3, 'f17':  1*pi/3, 'f18':  0*pi/3, 'f19':  -1*pi/3}
    face_translate = { 'f0':   (1*el,  5*elt/3),  'f1':   (1*el, 3*elt/3),  'f2': (1.5*el,  2*elt/3),  'f3':  (2*el,   3*elt/3),
                       'f4':   (2*el,  5*elt/3),  'f5': (0.5*el, 6*elt/3),  'f6': (0.5*el,  2*elt/3),  'f7': (1.5*el,  0*elt/3),
                       'f8': (2.5*el,  2*elt/3),  'f9': (2.5*el, 6*elt/3), 'f10': (  0*el, -7*elt/3), 'f11':   (0*el, -3*elt/3), 
                      'f12':   (1*el, -1*elt/3), 'f13':   (2*el,    -elt), 'f14':   (2*el, -7*elt/3), 'f15': (0.5*el, -6*elt/3),
                      'f16': (0.5*el, -4*elt/3), 'f17':   (1*el,    -elt), 'f18': (1.5*el, -4*elt/3), 'f19': (1.5*el, -6*elt/3)}

    face_order_mod=[('f5','f10'), ('f10','f9'),   ('f9','f14'), ('f14', 'f8'),  ('f8','f3'),
                    ('f3', 'f2'),  ('f3','f4'),   ('f4', 'f0'),  ('f0', 'f1'),  ('f1','f6'),
                    ('f6','f11'), ('f14','f19'), ('f19','f18'), ('f18','f13'), ('f13','f7'),
                    ('f7','f12'), ('f12','f17'), ('f19','f15'), ('f15','f16')]




#=============================

def initdata1000():
    citydata_alldf = pd.read_csv("Data\\geonames-all-cities-with-a-population-1000.csv", sep=';', header=0)
#    print(citydata_alldf.columns.tolist())
#    print(set(citydata_alldf.loc[citydata_alldf.loc[:,'LABEL EN']=="China","Admin1 Code"].tolist()))

    citydata_alldf.rename(columns={"Geoname ID": "geonameid", "ASCII Name": "asciiname", "Population": "population", "LABEL EN": "country"}, inplace=True)
    citydata_df = citydata_alldf.loc[:,['geonameid', 'asciiname', 'population','country', 'Admin1 Code', 'Coordinates']]
    str_list = citydata_df.loc[:,'asciiname']
    z_list = 0.0*citydata_df.loc[:,'population']

    citydata_df.loc[:,'lat_str']=str_list
    citydata_df.loc[:,'long_str']=str_list
    citydata_df.loc[:,'latitude']=z_list
    citydata_df.loc[:,'longitude']=z_list
    citydata_df.loc[:,'latr']=z_list
    citydata_df.loc[:,'longr']=z_list
    for row in citydata_df.index.tolist():
        citydata_df.at[row,'lat_str'], citydata_df.at[row,'long_str'] = citydata_df.at[row,'Coordinates'].split(",")
        citydata_df.at[row,'latitude']  = float(citydata_df.at[row,'lat_str'].strip() )
        citydata_df.at[row,'longitude'] = float(citydata_df.at[row,'long_str'].strip())
        citydata_df.at[row,'latr']  = citydata_df.at[row,'latitude' ]/deg_rads
        citydata_df.at[row,'longr'] = citydata_df.at[row,'longitude']/deg_rads
        
    citydata_df = citydata_df.loc[:,['geonameid', 'asciiname', 'country', 'Admin1 Code', 'latitude', 'longitude','population', 'latr', 'longr']]
    citydata_df.set_index("geonameid", inplace=True)
    print(citydata_df)
    print()
    citydata_df.to_pickle('Data\\City1000BasicData_pkl.pkl')

    citydata15k_df = citydata_df.loc[citydata_df['population']>=15000,:]
    print('Reduced dataset has size: ')
    print(citydata15k_df.shape)
    print()
    citydata15k_df.to_pickle('Data\\City15000BasicData_pkl.pkl')


#=============================

def initdata15000():
    citydata_alldf = pd.read_csv("Data\\cities15000.txt", sep='\t', header=None)
    citydata_df = citydata_alldf.rename(columns={0: 'geonameid', 1: 'asciiname',
                                              4: 'latitude', 5: 'longitude', 14: 'population', 17:'country'})

    citydata_df.loc[:,'latr']=0.0
    citydata_df.loc[:,'longr']=0.0
    citydata_df.loc[:,'Admin1 Code']=''
    for row in citydata_df.index.tolist():
        citydata_df.at[row,'latr']  = citydata_df.at[row,'latitude' ]/deg_rads
        citydata_df.at[row,'longr'] = citydata_df.at[row,'longitude']/deg_rads
        
    citydata_df = citydata_df.loc[:,['geonameid', 'asciiname', 'country','Admin1 Code', 'latitude', 'longitude','population', 'latr', 'longr']]
    citydata_df.set_index("geonameid", inplace=True)
    print(citydata_df)
    print(citydata_df.shape)
    print()
    citydata_df.to_pickle('Data\\CityBasicData_pkl.pkl')

#=======================================

def map_plot(citydata_df, face_order, initialr, axes_rots:tuple, agg:str, each=False, labelorder=False):

    try:
        assert len(axes_rots)==3
    except:
        print('axes rotation is the wrong length: '+str(len(axes_rots)))

    try:
        assert len(face_order)==19
    except:
        print('face_order has wrong length: '+str(len(face_order)))

    print('Setting up coordinate system')
    face_dict_ordered = ih.define_isocofaces()

    for face in face_dict_ordered.keys():
        try:
            assert face in [face_order[i][1] for i in range(19)]+[face_order[0][0]]
        except:
            print(str(face)+ ' is missing from ')
            print([face_order[i][1] for i in range(19)]+[face_order[0][0]])


    vd = ih.vertices()
    vdm = ih.vertices_mod(vd, axes_rots[0], axes_rots[1], axes_rots[2])

    ih.check_remap(face_dict_ordered, vdm)
    face_transform, face_rotate, face_translate = ihc.find_rot_trans(face_order, face_dict_ordered, initialr)

    
    print('Projecting cities onto polyhedron surface')

    citydata2_df = ihc.mapcities2icosofaces(citydata_df, face_dict_ordered, vdm)

    print('Assessing')
    citycount, metrics = ihc.countcities(citydata2_df, agg, verbose=True)
    print(citycount)
    print(metrics)

    print('Remapping to plane')
    face_dict = {k:set(face_dict_ordered[k]) for k in face_dict_ordered}
    citydata3_df, icoso_net = ihc.plot_net_calc(citydata2_df[['face', 'xy', 'population']], face_dict, face_transform) 


    print('Plotting')
    print("Showing population")
    ihc.plot_net( citydata3_df, icoso_net, True, [], "figure_pop") 
    print("Showing cities")
    ihc.plot_net( citydata3_df, icoso_net, False )
    if each: 
        for face in face_dict.keys():
            order = citycount.at[face,'order']
            print("Plotting "+str(face))
            if labelorder:
                ihc.plot_net(citydata3_df, icoso_net, True, [face], 'figure_'+str(order)+'_'+str(face),0.1) 
            else:
                ihc.plot_net(citydata3_df, icoso_net, True, [face], 'figure_'+str(face),0.1) 



#=======================================


deg_rads = 180/pi
edge_length = 1/math.sin(2*pi/5)
el=edge_length
elt=math.sqrt(3)/math.sin(2*pi/5)/2

face_order_ref=[('f5','f0'),   ('f0','f1'),   ('f1','f6'),   ('f1','f2'),   ('f2','f3'), 
                ('f2','f7'),   ('f3','f4'),   ('f3','f8'),   ('f4','f9'),   ('f7','f12'), 
                ('f12','f17'), ('f17','f16'), ('f16','f11'), ('f16','f15'), ('f15','f10'),
                ('f17','f18'), ('f18','f13'), ('f18','f19'), ('f19','f14') ]


face_order_Bering=[('f5','f0'),   ('f0','f1'),   ('f1','f6'),   ('f1','f2'),  ('f2','f3'), 
                   ('f2','f7'),   ('f3','f4'),   ('f3','f8'),   ('f4','f9'),  ('f6','f12'),
                   ('f12','f17'), ('f6','f11'),  ('f11','f16'), ('f16','f15'), ('f15','f10'),
                   ('f9','f14'),  ('f14','f19'), ('f7','f13'),  ('f19','f18') ]

face_order_Bering1=[('f0','f1'),   ('f1','f6'),   ('f1','f2'),  ('f2','f3'), 
                   ('f2','f7'),   ('f3','f4'),   ('f3','f8'),   ('f4','f9'),  ('f6','f12'),
                   ('f12','f17'), ('f6','f11'),  ('f11','f16'), ('f16','f15'), ('f11','f5'), ('f5','f10'),
                   ('f8','f14'),  ('f14','f19'), ('f8','f13'),  ('f19','f18') ]

face_order_Indian=[ ('f5','f0'),   ('f0','f1'),   ('f1','f6'),   ('f1','f2'),   ('f2','f3'),
                    ('f6','f12'),  ('f3','f4'),   ('f3','f8'),   ('f4','f9'),   ('f12','f7'),
                    ('f12','f17'), ('f5','f11'),  ('f11','f16'), ('f16','f15'), ('f5','f10'),
                    ('f17','f18'), ('f18','f13'), ('f18','f19'), ('f9','f14') ]

#initdata1000()
#initdata15000()

#print('Reading cities')
citydata_1k_df = pd.read_pickle('Pickle\\City1000BasicData_pkl.pkl')
citydata_15k_df = pd.read_pickle('Pickle\\CityBasicData_pkl.pkl')

map_plot(citydata_1k_df, face_order_Indian, pi/3, ((180-5)/deg_rads, 35/deg_rads, 55/deg_rads),"sum",True)
map_plot(citydata_15k_df, face_order_Bering1, -pi/3, (40.47/deg_rads, 31.49/deg_rads, 44.61/deg_rads),"sum")
