
Based on input data regarding city locations and populations from 
https://public.opendatasoft.com/explore/dataset/geonames-all-cities-with-a-population-1000/export/?disjunctive.cou_name_en&sort=name
    # https://github.com/benjamincjackson/world_cities


The main code is in IcosohedralCities_Runs. The output is the city locations identified in the input data with global latitude and longitude coordinates mapped to a planar net of the projection of the globe onto an icosohedron. Initialisation procedures start by converting the original input data into a python pickle compression of a pandas dataframe, for speedy future reloading. It also defines some specific parameters selecting the orientation of the selected icosohedron (relative to the globe) and how the planar net is defined

The supporting code in IcosohedralCities performs a number of functions associated with plotting and saving initialised data. The supporting code in Icosohedron performs some elementary geometric transformations associated with the regular Platonic icosohedron, including a specific scheme for labelling the faces. The supporting code in SphericalCoords performs elementary geometric calculations associated with points on a sphere whose location is expressed in spherical coordinates. These calculations include such functions as conversion between spherical and cartesian coordinates, between radians and degrees, finding the great circle distance between two points, projection of a point on a sphere to a triangular face defined by 3 other points on the sphere, mapping spherical coordinates with a given polar axis and longitude datum to an alternative axis and datum, and so on.
