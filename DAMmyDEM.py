#!/usr/bin/python
# -*- coding: utf-8 -*-

import osgeo
from osgeo import gdal
from osgeo import ogr
import os
import sys
import time
import numpy
import math
import json
import easygui as eg
from scipy import signal


numpy.seterr(invalid='ignore')  # Caution, setting the error management in numpy

"""
Requiers TauDem, available at: https://hydrology.usu.edu/taudem/taudem5/    or    https://github.com/dtarb/TauDEM/releases


This script asks for an input DTM and a poligonal shapefile of a dam (e.g., valley-transverse rectangle-like shapefile) 
to modify the DTM and backpropagate a flat (water) surface taking as a reference the maximum level of the dam (which corresponds 
to the highest elevation of the DTM pixels within the input dam shapefile)

authors: Stefano Crema, Alessandro Sarretta & Marco Cavalli,
last edited: 04/2024
Copyright (C) 2014-2024  Stefano Crema, Alessandro Sarretta & Marco Cavalli, @ CNR-IRPI, corso stati uniti, 4, 35127, Padova (Italy)
mail: stefano.crema@cnr.it

###############################################################################
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
For a copy of the GNU-GPL v2 license please visit: http://www.gnu.org/licenses/gpl-2.0.txt
"""

"""
--------------------------------------
###########################################################################        
IINTRODUCING DAMMED-LIKE WATER SURFACES OVER A DTM, maybe for later processing hunting Sediment Connectivity... who knows
###########################################################################
--------------------------------------


####################################################################################
----------------------------------    USAGE    -------------------------------------
############# We use a JSON file to set input and output files and parameters ###################
############# probably better to use forward slash in paths ###############################


"""
init_time =time.time()

#def dams(dtm_f, dam_shp, dtm_out):
#Inputs
infile_json_dam = eg.fileopenbox(msg='Please locate the JSON input file',
                    title='Specify JSON File', default='C:\\users\\*.json',
                    filetypes='*.json')

if infile_json_dam is None:
    print(time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't find/open input JSON parameter file")
    sys.exit(1)
else:
    print(time.strftime("%d/%m/%Y %H:%M:%S    "), "JSON parameter file found!")

# Opening JSON file with input files/params
with open(infile_json_dam, 'r') as openfile:
    # Reading from json file
    json_object_dam = json.load(openfile)

print('JSON dam input file:    ' + str(json_object_dam))
#
##
dtm_f = json_object_dam['dtm_in']
dam_shp = json_object_dam['dam_in']
dtm_out = json_object_dam['dtm_dam_out']
dtm_out_fill = json_object_dam['dtm_dam_out_fill']
#
filename = dtm_f.replace('\\', '/')
filename = str(filename)

## Main DTM
tif = osgeo.gdal.Open(filename)  # opening the file (1 band, bit not specified, but usually 32)
#
# folder path
dir_pat = os.path.dirname(os.path.realpath(filename))  # path of the selected input file
dir_path = str(dir_pat.replace('\\', '/'))
#
# defining a global function def dams (input1, input2, output1 ...) to use for passing the filenames created inside the functions
# without making the user specify each filepath/filename
global path_dam_dtm
path_dam_dtm = dir_path
#
# Opening Messages
if tif is None:
    print
    time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't find/open input DTM dataset"
    sys.exit(1)
else:
    print
    time.strftime("%d/%m/%Y %H:%M:%S    "), "opening DTM for dams computation was successful!"
#
# Colonne, righe, bande
cols = tif.RasterXSize
rows = tif.RasterYSize
bands = tif.RasterCount
#
# dtm as an array                #################################       ATTENZIONE     ##########################
tif_ar = tif.ReadAsArray()  # con GDAL 1.6 e 1.7 non funzia 90%, con gdal bindings per windows 1.9.2 sembra funzionare provare anche a metter su ultime versioni di numpy
tif_ar = tif_ar.astype(float)
band = tif.GetRasterBand(1)  # bands retrieving
geoinf = tif.GetGeoTransform()  # upper left coordinates, cellsize
proj = tif.GetProjection()
tipo = tif_ar.dtype
del (tif)
#
'''Flow directions --> Only once here for original DTM'''
# D8 flow directions & outputs to derive watersheds
os.system((("mpiexec -n 8 D8Flowdir -p ").lower()) + filename[0:-4] + "p.tif" + " -sd8 " + filename[
                                                                                           0:-4] + "sd8.tif -fel " + filename)  # unix case
#
# reclassifying D8 flow directions
tif_fdir8 = osgeo.gdal.Open(
    filename[0:-4] + "p.tif")  # opening the flowdir file (1 band, bit not specified, but usually 32)
if tif_fdir8 is None:
    print
    time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't open input dir8 dataset"
    sys.exit(1)
else:
    print
    time.strftime("%d/%m/%Y %H:%M:%S    "), "Flowdir opening was successful!"
# Array
tif_fdir8_ar = tif_fdir8.ReadAsArray()
tif_fdir8_ar = tif_fdir8_ar.astype(int) # with int type still nodata are represented by the min --> (-1) value
ndv = numpy.min(tif_fdir8_ar)
tif_fdir8_ar[tif_fdir8_ar == ndv] = -1  # to Basin Code
del (tif_fdir8)
os.remove(filename[0:-4] + "p.tif")
os.remove(filename[0:-4] + "sd8.tif")
#
#

''' Opening and working with shape dam(s) file'''

file_dam = dam_shp.replace('\\', '/')
file_dam = str(file_dam)

try:
    source_dam = ogr.Open(file_dam, 1)
except:
    msg_openshapefile_dam = ('Unable to open of find file ' + str(file_dam) + '!')
    print(time.strftime("%d/%m/%Y %H:%M:%S    "), msg_openshapefile_dam)

layer_dam = source_dam.GetLayer()
layer_defn_dam = layer_dam.GetLayerDefn()
field_names_dam = [layer_defn_dam.GetFieldDefn(i).GetName() for i in range(layer_defn_dam.GetFieldCount())]

try:
    in_dam_field = field_names_dam.index('dam_id')
    msg_dam = ('dam_id field found in field n° ' + str(in_dam_field + 1) + '!')
    print(time.strftime("%d/%m/%Y %H:%M:%S    "), msg_dam)
except:
    msg_dam = ('dam_id field not found, ignoring input')
    print(time.strftime("%d/%m/%Y %H:%M:%S    "), msg_dam)
#
'''Unique dam values-code'''
values_dam_list = [ feature.GetField("dam_id") for feature in layer_dam ]

#
'''Looping through features in the case multiple dams (multiple polygons) are inserted --> 1polygon = 1 dam
better if Flow length/watershed is computed once with the unique (integer and non equal) values of the field dam id are propagated
For doing so I need to parse each dam separately setting elevation max per each dam and not globally'''
#
#
dam_ds = gdal.GetDriverByName('GTiff').Create((dir_path + "/dams.tif"), tif_ar.shape[1], tif_ar.shape[0], 1,
                                              gdal.GDT_Int16)  # shape sono righe e colonne, 1 è il numero di bande
dam_ds.SetGeoTransform(geoinf)
gdal.RasterizeLayer(dam_ds, [1], layer_dam, options=['ATTRIBUTE=dam_id'])
dam_ds.GetRasterBand(1).SetNoDataValue(-9999)
dam_ar = dam_ds.GetRasterBand(1).ReadAsArray(buf_type=gdal.GDT_Int16)
del (dam_ds)

'''Preparing for dam-code propagation'''
#
# bigger matrix to move with indexes, will store the propagated basin code
big_wat_dam = numpy.zeros(shape=((tif_ar.shape[0]) + 2, (tif_ar.shape[1]) + 2),
                          dtype=numpy.int16)  # add rows and columns of zeros to elaborate better matrix indexes
big_wat_dam[1:-1, 1:-1] = dam_ar
#
# zero matrix bigger than F_dir8, to avoid border indexing problems
Fd8 = numpy.zeros(shape=((tif_fdir8_ar.shape[0]) + 2, (tif_fdir8_ar.shape[1]) + 2),
                  dtype=numpy.int16)  # add rows and columns of zeros to elaborate better matrix indexes
Fd8[1:(Fd8.shape[0] - 1), 1:(Fd8.shape[1] - 1)] = Fd8[1:(Fd8.shape[0] - 1), 1:(Fd8.shape[1] - 1)] + tif_fdir8_ar
Fdir8 = Fd8
#
# Creating a bigger matrix as large as weight(and all the matrices) to store the weighted flow length values
dam_m = numpy.zeros(shape=((tif_ar.shape[0]) + 2, (tif_ar.shape[1]) + 2), dtype=numpy.int16)
#
#
# DOING THE PROCEDURE
# Let's go for the search and algo-rhytm for the weighted-Flow-Length
#
start = time.time()  # take your time
#
print(msg_dam)
#
print(time.strftime("%d/%m/%Y %H:%M:%S    "), 'Computing dams watershed extraction ... ...')
# I retrieve the basin code theat I will propagate
damV = numpy.where(
    big_wat_dam > 0)  # fast coordinates all the Dam Code values, starting from them to go forward and compute flow length
#
y = damV[0]  # rows, code indexes
x = damV[1]  # columns, code indexes pay attention not to invert values !!!!!!!!!!!!!!
#
# initializing lists for outlet and moving cell coordinates, in function of their position
YC1 = []
YC2 = []
YC3 = []
YC4 = []
YC5 = []
YC6 = []
YC7 = []
YC8 = []
XC1 = []
XC2 = []
XC3 = []
XC4 = []
XC5 = []
XC6 = []
XC7 = []
XC8 = []
#
#   Flow Directions Taudem
#   4   3   2
#   5   -   1
#   6   7   8
#
#   Draining-in Direction Matrix
#   8   7   6
#   1   -   5
#   2   3   4
#
'''Trying to propagate dam code in flow length/watershed per easch dam avoiding loops of flow length so I keep going only if I am in the watershed,
If I find a "nested" watershed I should with something like numpy.where(i1 == 1) and .where dam_ar =0 otherwise I loose the dam code, overwritte by the 
flow length/watershed of the downstream dam'''

i1 = Fdir8[y, x - 1]  # Searching for Basin Code with cells draining into them, 8 directions
i1_dam = big_wat_dam[y, x - 1]
D1 = numpy.where((i1 == 1) & (i1_dam == 0))  #
YC1.extend(y[D1])  # coordinates satisfying the conditions
XC1.extend(x[D1])
dam_m[YC1, XC1] = big_wat_dam[YC1, XC1]  # initialize basin code at cells draining to Basin Code
#
i2 = Fdir8[y + 1, x - 1]  # Searching for Basin Code with cells draining into them, 8 directions
i2_dam = big_wat_dam[y + 1, x - 1]
D2 = numpy.where((i2 == 2) & (i2_dam == 0)) #
YC2.extend(y[D2])  # coordinates satisfying the conditions
XC2.extend(x[D2])
dam_m[YC2, XC2] = big_wat_dam[YC2, XC2]  # initialize basin code at cells draining to Basin Code
#
i3 = Fdir8[y + 1, x]  # Searching for Basin Code with cells draining into them, 8 directions
i3_dam = big_wat_dam[y + 1, x]
D3 = numpy.where((i3 == 3) & (i3_dam == 0)) #
YC3.extend(y[D3])  # coordinates satisfying the conditions
XC3.extend(x[D3])
dam_m[YC3, XC3] = big_wat_dam[YC3, XC3]  # initialize basin code at cells draining to Basin Code
#
i4 = Fdir8[y + 1, x + 1]  # Searching for Basin Code with cells draining into them, 8 directions
i4_dam = big_wat_dam[y + 1, x + 1]
D4 = numpy.where((i4 == 4) & (i4_dam == 0)) #
YC4.extend(y[D4])  # coordinates satisfying the conditions
XC4.extend(x[D4])
dam_m[YC4, XC4] = big_wat_dam[YC4, XC4]  # initialize basin code at cells draining to Basin Code
#
i5 = Fdir8[y, x + 1]  # Searching for Basin Code with cells draining into them, 8 directions
i5_dam = big_wat_dam[y, x + 1]
D5 = numpy.where((i5 == 5) & (i5_dam == 0)) #
YC5.extend(y[D5])  # coordinates satisfying the conditions
XC5.extend(x[D5])
dam_m[YC5, XC5] = big_wat_dam[YC5, XC5]  # initialize basin code at cells draining to Basin Code
#
i6 = Fdir8[y - 1, x + 1]  # Searching for Basin Code with cells draining into them, 8 directions
i6_dam = big_wat_dam[y - 1, x + 1]
D6 = numpy.where((i6 == 6) & (i6_dam == 0)) #
YC6.extend(y[D6])  # coordinates satisfying the conditions
XC6.extend(x[D6])
dam_m[YC6, XC6] = big_wat_dam[YC6, XC6]  # initialize basin code at cells draining to Basin Code
#
i7 = Fdir8[y - 1, x]  # Searching for Basin Code with cells draining into them, 8 directions
i7_dam = big_wat_dam[y - 1, x]
D7 = numpy.where((i7 == 7) & (i7_dam == 0)) #
YC7.extend(y[D7])  # coordinates satisfying the conditions
XC7.extend(x[D7])
dam_m[YC7, XC7] = big_wat_dam[YC7, XC7]  # initialize basin code at cells draining to Basin Code
#
i8 = Fdir8[y - 1, x - 1]  # Searching for Basin Code with cells draining into them, 8 directions
i8_dam = big_wat_dam[y - 1, x - 1]
D8 = numpy.where((i8 == 8) & (i8_dam == 0)) #
YC8.extend(y[D8])  # coordinates satisfying the conditions
XC8.extend(x[D8])
dam_m[YC8, XC8] = big_wat_dam[YC8, XC8]  # initialize basin code at cells draining to Basin Code
#
# start =clock()#da cancellare poi.....!!!!!! Solo per check
count = 1  # "0" passage already done during the previous step
while len(YC1) or len(YC2) or len(YC3) or len(YC4) or len(YC5) or len(YC6) or len(YC7) or len(YC8) > 0:
    # Converting into array to be able to do operations
    YYC1 = numpy.asarray(YC1);
    XXC1 = numpy.asarray(XC1)
    YYC2 = numpy.asarray(YC2);
    XXC2 = numpy.asarray(XC2)
    YYC3 = numpy.asarray(YC3);
    XXC3 = numpy.asarray(XC3)
    YYC4 = numpy.asarray(YC4);
    XXC4 = numpy.asarray(XC4)
    YYC5 = numpy.asarray(YC5);
    XXC5 = numpy.asarray(XC5)
    YYC6 = numpy.asarray(YC6);
    XXC6 = numpy.asarray(XC6)
    YYC7 = numpy.asarray(YC7);
    XXC7 = numpy.asarray(XC7)
    YYC8 = numpy.asarray(YC8);
    XXC8 = numpy.asarray(XC8)
    #
    # Now I will propagate always the same basin code
    #
    YYC1 = (YYC1);
    XXC1 = (XXC1 - 1)  #
    YYC2 = (YYC2 + 1);
    XXC2 = (XXC2 - 1)  #
    YYC3 = (YYC3 + 1);
    XXC3 = (XXC3)  # l
    YYC4 = (YYC4 + 1);
    XXC4 = (XXC4 + 1)  #
    YYC5 = (YYC5);
    XXC5 = (XXC5 + 1)  #
    YYC6 = (YYC6 - 1);
    XXC6 = (XXC6 + 1)  #
    YYC7 = (YYC7 - 1);
    XXC7 = (XXC7)  #
    YYC8 = (YYC8 - 1);
    XXC8 = (XXC8 - 1)  #
    #
    if len(YYC1) > 0:
        dam_m[YYC1, XXC1] = dam_m[YC1, XC1]
    else:
        pass
    if len(YYC2) > 0:
        dam_m[YYC2, XXC2] = dam_m[YC2, XC2]
    else:
        pass
    if len(YYC3) > 0:
        dam_m[YYC3, XXC3] = dam_m[YC3, XC3]
    else:
        pass
    if len(YYC4) > 0:
        dam_m[YYC4, XXC4] = dam_m[YC4, XC4]
    else:
        pass
    if len(YYC5) > 0:
        dam_m[YYC5, XXC5] = dam_m[YC5, XC5]
    else:
        pass
    if len(YYC6) > 0:
        dam_m[YYC6, XXC6] = dam_m[YC6, XC6]
    else:
        pass
    if len(YYC7) > 0:
        dam_m[YYC7, XXC7] = dam_m[YC7, XC7]
    else:
        pass
    if len(YYC8) > 0:
        dam_m[YYC8, XXC8] = dam_m[YC8, XC8]
    else:
        pass
    #
    # Reconstructing all X and Y of this step and moving on upwards (Downstream if you think in GIS, right?)
    YY = [];
    XX = []
    YY.extend(YYC1);
    XX.extend(XXC1)
    YY.extend(YYC2);
    XX.extend(XXC2)
    YY.extend(YYC3);
    XX.extend(XXC3)
    YY.extend(YYC4);
    XX.extend(XXC4)
    YY.extend(YYC5);
    XX.extend(XXC5)
    YY.extend(YYC6);
    XX.extend(XXC6)
    YY.extend(YYC7);
    XX.extend(XXC7)
    YY.extend(YYC8);
    XX.extend(XXC8)
    #
    YY = numpy.asarray(YY)
    XX = numpy.asarray(XX)
    #
    i1 = Fdir8[YY, XX - 1]  # Searching for cells draining into them, 8 directions
    i1_dam = big_wat_dam[YY, XX - 1]
    D1 = numpy.where((i1 == 1) & (i1_dam == 0)) #
    YC1 = YY[D1]  # coordinates satisfying the conditions
    XC1 = XX[D1]
    #
    i2 = Fdir8[YY + 1, XX - 1]  # Searching for cells draining into them, 8 directions
    i2_dam = big_wat_dam[YY + 1, XX - 1]
    D2 = numpy.where((i2 == 2) & (i2_dam == 0)) #
    YC2 = YY[D2]  # coordinates satisfying the conditions
    XC2 = XX[D2]
    #
    i3 = Fdir8[YY + 1, XX]  # Searching for cells draining into them, 8 directions
    i3_dam = big_wat_dam[YY + 1, XX]
    D3 = numpy.where((i3 == 3) & (i3_dam == 0)) #
    YC3 = YY[D3]  # coordinates satisfying the conditions
    XC3 = XX[D3]
    #
    i4 = Fdir8[YY + 1, XX + 1]  # Searching for cells draining into them, 8 directions
    i4_dam = big_wat_dam[YY + 1, XX + 1]
    D4 = numpy.where((i4 == 4) & (i4_dam == 0))  #
    YC4 = YY[D4]  # coordinates satisfying the conditions
    XC4 = XX[D4]
    #
    i5 = Fdir8[YY, XX + 1]  # Searching for cells draining into them, 8 directions
    i5_dam = big_wat_dam[YY, XX + 1]
    D5 = numpy.where((i5 == 5) & (i5_dam == 0)) #
    YC5 = YY[D5]  # coordinates satisfying the conditions
    XC5 = XX[D5]
    #
    i6 = Fdir8[YY - 1, XX + 1]  # Searching for cells draining into them, 8 directions
    i6_dam = big_wat_dam[YY - 1, XX + 1]
    D6 = numpy.where((i6 == 6) & (i6_dam == 0)) #
    YC6 = YY[D6]  # coordinates satisfying the conditions
    XC6 = XX[D6]
    #
    i7 = Fdir8[YY - 1, XX]  # Searching for cells draining into them, 8 directions
    i7_dam = big_wat_dam[YY - 1, XX]
    D7 = numpy.where((i7 == 7) & (i7_dam == 0)) #
    YC7 = YY[D7]  # coordinates satisfying the conditions
    XC7 = XX[D7]
    #
    i8 = Fdir8[YY - 1, XX - 1]  # Searching for cells draining into them, 8 directions
    i8_dam = big_wat_dam[YY - 1, XX - 1]
    D8 = numpy.where((i8 == 8) & (i8_dam == 0)) #
    YC8 = YY[D8]  # coordinates satisfying the conditions
    XC8 = XX[D8]
    count = count + 1
#
#
elapsed = (time.time() - start)  # computational time
print
time.strftime(
    "%d/%m/%Y %H:%M:%S    "), "Process concluded succesfully \n", "%.2f" % elapsed, 'seconds for dams watersheds calculation with ', int(
    count), ' iterations'  # truncating the precision
# os.system("pause")
#
del Fdir8

''' the next variable (dam_m) represents the watershed extracted via simplified weighted flow length to be used for further refinements
such as to intersect with elevations to create the virtual lake, multiple dams I can create a loop also here in addition to the loop for
max unique elevation per single dam or grouping everythong in one loop if possible, the important thing is to avoid loop that recompute 
Flow direction and watershed extraction when I just need it once'''
dam_m = dam_m[1:dam_m.shape[0] - 1, 1:dam_m.shape[
                                       1] - 1]  # reshaping weigthed flow length/watershed extracted and coded with dam code, we need this step to homogenize matrices dimensions!!!!!!!!!!
#
min_tif_ar = tif_ar.min()
# Computing elevation difference and retained volume

for i in range(len(values_dam_list)):
    dam_value_i = values_dam_list[i]
    dam_dem_elevations = tif_ar[dam_ar == dam_value_i]  # assuming no intersection with nodata !!!
    dam_dem_unique_max = max(dam_dem_elevations)
    tif_ar[dam_ar == dam_value_i] = dam_dem_unique_max
    '''#doing like this dam area is included but the volume is only for the catchment, the dam has 0 difference in elevation since we modify 
    # here the DTM inserting the fixed elevation dam like in pseudo original surface, I can retrieve dam area from here '''
    dam_area = len(dam_dem_elevations) * numpy.square(geoinf[1])  # "dam" area in square prj units - square meters
    #
    '''Probably better to loop over the followinf variables'''
    mask_elev_orig = tif_ar[numpy.where((dam_m == dam_value_i) & (tif_ar<=dam_dem_unique_max))] #restrict computation only in the interested area, and this way we do not include the dam but only the "lake"
    elev_diff_ar = dam_dem_unique_max - mask_elev_orig
    max_elev_diff = max(elev_diff_ar)
    stored_volume = sum(elev_diff_ar * numpy.square(geoinf[1])) #multiply by pixel area derived from geotiff geoinf structure, Assuming projected coordinates, e.g.  meters, pixel size and elevation must be in the same unit
    stored_volume_km3 = stored_volume/1000000000
    total_area = len(mask_elev_orig)*numpy.square(geoinf[1])
    lake_only_area = (total_area - dam_area)/1000000 #lake only area in sq km
    #tif_ar[numpy.where(dam_m > 0)] = min_tif_ar
    tif_ar[numpy.where((dam_m == dam_value_i) & (tif_ar<=dam_dem_unique_max))] = dam_dem_unique_max #set fixed elevation to create the lake

tif_ar[numpy.where(tif_ar == min_tif_ar)] = -9999

#
# write dammed DTM
dam_m_ds = gdal.GetDriverByName('GTiff').Create(dtm_out, tif_ar.shape[1], tif_ar.shape[0], 1,
                                               gdal.GDT_Float32)  # shape sono righe e colonne, 1 è il numero di bande
dam_m_ds.SetGeoTransform(geoinf)  # conferisco coordinate e proiezione del raster in input
dam_m_ds.SetProjection(proj)
dam_m_ds.GetRasterBand(1).SetNoDataValue(-9999)
dam_m_ds.GetRasterBand(1).WriteArray(tif_ar, 0,
                                    0)  # scrivo effettivamente il raster, prima avevo solo allocato la memoria
dam_m_ds = None  # Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
del dam_m_ds
#

## Exporting new fel DTM
# os.system((("mpiexec -n 8 PitRemove -z ").lower()) + dtm_out + ' -fel ' + dtm_out_fill)
#
print(time.strftime("%d/%m/%Y %H:%M:%S    "), 'Dams computation concluded... great!')
#
#
elapsed_tot = time.time() - init_time
print(time.strftime("%d/%m/%y %H:%M:%S    "), "Calculation finished in", "%.2f" % elapsed_tot,
      "seconds! \n Compliments!")

## Finish line!!!!
