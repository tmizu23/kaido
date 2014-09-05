# -*- coding: utf-8 -*-
import os, sys, gdal, ogr, numpy,argparse
from gdalconst import *
from osgeo import osr
from math import *

def getKaido(x,y):
    iDifX = [0,1,1,1,0,-1,-1,-1]
    iDifY = [1,1,0,-1,-1,-1,0,1]
    iDifL = [1,1.4142,1,1.4142,1,1.4142,1,1.4142]
    aAngle=8*[0]
    dDifHeight=0 
    dDist=0
    dKaido=0

    for i in range(8):
       j = 1
       aAngle[i] = -100
       dDist = iDifL[i] * j * x_size
       #print dDist,dist
       while dDist < dist:
           xx = x + iDifX[i]*j
           yy = y + iDifY[i]*j
           if( 0<=xx and xx<cols and 0<=yy and yy<rows):
               dDifHeight = data[yy][xx] - data[y][x]
               dAngle = atan(dDifHeight/dDist)
               if dAngle > aAngle[i]:
                   aAngle[i] = dAngle
               j=j+1
               dDist = iDifL[i]*j*x_size
           else:
               return -100
       dKaido= dKaido + aAngle[i]
    return dKaido/8.0

if __name__ == '__main__':
    ##
    # inputのデータ形式はFloat32だけ対応
    # フォーカル範囲内にNODATがあるものはダメ

    parser = argparse.ArgumentParser(description='This is focal statistics program.') # parserを作る
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('-nodata', type=float,default=-9999.0,help="Area's edge is nodata of this value.") # オプションを追加します
    parser.add_argument('-r', type=int,required=True,help='statistical radius by map unit') # このオプションは必須です
    parser.add_argument('--version', action='version', version='%(prog)s 0.1') # version
    args = parser.parse_args()


    argvs = sys.argv
    input=args.input
    output=args.output
    dist=args.r
    ndv=args.nodata


    gdal.AllRegister()
    ds = gdal.Open(input)
    data = ds.GetRasterBand(1).ReadAsArray()
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    (upper_left_x, x_size, x_rotation, upper_left_y, y_rotation, y_size) = ds.GetGeoTransform()
    Projection = osr.SpatialReference()
    Projection.ImportFromWkt(ds.GetProjectionRef())

    outdata = numpy.empty((rows,cols))
    outdata.fill(ndv)

    for y in range(rows):
        print y
        for x in range(cols):
            kaido = getKaido(x,y)
            if kaido == -100:
               outdata[y][x] = ndv
            else:
               outdata[y][x] = kaido


    driver = gdal.GetDriverByName('GTiff')
    dst_ds = driver.Create(output, cols, rows, 1,gdal.GDT_Float32)

    dst_ds.SetGeoTransform(ds.GetGeoTransform())
    dst_ds.SetProjection( Projection.ExportToWkt())
    dst_band = dst_ds.GetRasterBand(1)
    dst_band.WriteArray(outdata)
    dst_band.SetNoDataValue(ndv)

    dst_ds = None
    ds = None
