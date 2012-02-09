"""
We're going to generate a bunch of rasters - merge them all together
using mosaicing tools.
"""

"""
0) read metadata in csv file
1) calculate latlon extent for ascii grid
2) calculate array size using modh modv line sample
3) generate array, 0 as nodata
4) load csv data - should be all int values
5) calculate pixel ids relative to array indices
6) for each period/threshold var, initialize array nodata and insert pixel
values
7) apply symbology
8) save as geotiff w/gdal
9) mosaic different countries together
10) export PNG
"""

import sys, os

sys.path.insert(0, "/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages")

from datetime import date

import numpy as np

import reproject
import rasterutils

def dataPrep(country):
    """
    config:
    number of static files (3 for now)
    start and end dates
    resolution
    thresh range
    paths, filename conventions

    actions:
    1) get output data in stata format
    2) get static data
    3) append all static files together
    """

    return

def getColorParams(name=None):
    """Returns bins and color maps needed to draw nice-looking rasters.

    Note that the default values assume that Hansen hits are represented by
    the number 101, and VCF values by 111. This is determined by the Stata
    code, so if something looks wrong, make sure that the Stata and Python
    code agree in this respect.

    Also note that the default color scheme runs from yellow to red for FORMA
     values (50-99), brown for Hansen, and green for VCF.

    """
    
    if type(name) != list:
        name = [name]

    # TODO: Parameterize bins, colors, etc. in the Stata code so don't have keep track of this code in two places

    default_bins = np.array([50, 60, 70, 80, 90, 100, 110])

    default_cmap = np.array([   [  0,  0,  0],
                                [255, 255,  0],
                                [255, 191,  0],
                                [255, 127,  0],
                                [255,  63,  0],
                                [255,   0,  0],
                                [101,  45,  5],
                                [ 21,  85, 44]])


    latest = np.array([1, 100, 110])
    latest_cmap = np.array([ [  0,  0,  0],
                        [255,  0,  0],
                        [50, 50, 50],
                        [21, 85, 44]])

    latest = np.array([1, 100, 104, 110])
    latest_cmap = np.array([ [  0,  0,  0],
                             [255, 255, 50],
                             [50, 50, 50],
                             [255,  0,  0],
                             [21, 85, 44]])

    vcf = [53, 131, 17]
    hansen = [137, 137, 137]
    post_2005 = [253, 255, 0]
    latest_hits = [166, 92, 60]

    probs_cutoff = 49
    vcf_cutoff = 102
    hansen_cutoff = 104
    post_2005_cutoff = 106
    latest_cutoff = 108

    latest = np.array([1,
                       vcf_cutoff,
                       hansen_cutoff,
                       post_2005_cutoff])

    latest_cmap = np.array([[0, 0, 0], [0, 123, 37], [166, 23, 0], [255, 255, 0], [255, 255, 0]])
    
    moratorium = np.array([1, 2, 3, 4])
    moratorium_cmap = np.array([[0, 0, 0],
                                [21, 85, 44],
                                [50, 50, 50],
                                [255, 255, 50],
                                [255, 0, 0]])
                                
    if name == ["forma"]:
        return default_bins[:-2, :], default_cmap[:-2, :]
    if name == ["forma", "hansen"]:
        return default_bins[:-1, :], default_cmap[:-1, :]
    if name == ["latest"]:
        return latest, latest_cmap
    if name == ["moratorium"]:
        return moratorium, moratorium_cmap
    if name == ["forma", "hansen", "vcf"] or name==None:
        return default_bins, default_cmap
    
    raise ValueError("No color map defined for name: " + str(name))

def stataProcessing(country):
    """
    config: nothing hardcoded in Stata
    pass path, country, start/end dates + thresh range into dofile via CL

    actions:
    1) Convert static data to Stata format
    2) keep country, lat, lon, modh, modv, line, sample
    3) merge estimated data and static as reformatted above
    4) multiply all probabilities by 100, recast int
    5) create variables to store date when a pixel went over various thresholds
    6) based on start date, gen period when passed thresholds
    5) outsheet everything to csv file, noheader
    6) outsheet min/max modh modv line sample lat lon for each threshold of interest
    """

    return

def mergeCountryData(countries):
    """
    # was thinking of generating ascii grids as a mosaic before getting to
    # Arc, but would run into max array size for Numpy with large number of countries

    loop through countries and append all csv files together - this is why we
     need to retain country var
    """
    return

def genDateList(start_date, end_date):
    dates = []

    yyyy = int(start_date[:4])
    mm = int(start_date[4:])

    start_date = date(yyyy, mm, 1)

    yyyy = int(end_date[:4])
    mm = int(end_date[4:])

    end_date = date(yyyy, mm, 1)

    dt = start_date

    while dt <= end_date:

        dates.append(dt)
        yyyy = dt.year
        mm = dt.month
        mm += 1
        if mm > 12:
            mm = 1
            yyyy += 1

        dt = date(yyyy, mm, 1)

    return dates

def saveToSinuCsv(xs, ys, vals):
    out = np.hstack((xs.reshape(-1, 1), ys.reshape(-1, 1),
                     vals.astype(np.float32)))

    np.savetxt("/Users/robin/delete/indonesia/indonesia_xy.csv", out, fmt='%f',
               delimiter=",")

    import sys

    sys.exit()

def getFormat(img_ext):
    formats = {"png":"PNG", "asc":"AAIGRID", "txt":"AAIGRID",
               "jpg":"JPEG", "tif":"GTiff"}
    try:
        fmt = formats[img_ext]
    except KeyError:
        raise KeyError("Image extension '%s' not recognized. Recognized extensions are %s" % (img_ext, "\n".join(formats.keys())))

    return fmt

def getFields(in_csv):
    f = open(in_csv, "r")
    header = f.readline().strip()
    if not header.startswith("modh,modv,sample,line,"):
        f.close()
        raise ValueError("Error - Data must be formatted as follows:\nmodh,modv,sample,line,DATA_0,...DATA_N\n" + header)
    f.close()

    return header.split("modh,modv,sample,line,")[1].split(",")

def parseData(in_csv):

    print "Loading data:", in_csv
    data = np.loadtxt(in_csv, dtype=np.int16, delimiter=",", skiprows=1)
    print "Data loaded\n"

    modh = data[:,0]
    modv = data[:,1]
    sample = data[:,2]
    line = data[:,3]
    vals = data[:,4:]

    return modh, modv, sample, line, vals

def parseCoordinates(res, modh, modv, sample, line):
    xmags, ymags = reproject.modis2globalmags(res, modh, modv, sample, line)
    xs, ys = reproject.globalmags2sinuxy(xmags, ymags)

    #saveToSinuCsv(xs, ys, vals)
    #xs, ys = rasterutils.modisGridToXY(modh, modv, sample, line, pixelres)

    xmin, xmax, ymin, ymax = rasterutils.xyToBbox(res, xs, ys)

    return xs, ys, xmin, xmax, ymin, ymax

def getPaths(in_csv, fieldprefix, img_ext):

    fields = getFields(in_csv)

    path, fname = os.path.split(in_csv)
    iso = os.path.splitext(fname)[0]

    outfiles = list()

    for field in fields:
        try:
            dt = field.split(fieldprefix)[1]
            dt = date(int(dt[:4]), int(dt[4:]), 1)
            outfiles.append("%s/%s_%i%02i.%s" % (path, iso, dt.year,
                                               dt.month, img_ext))
        except (TypeError, ValueError, IndexError):
            # just means the header doesn't have dates in vars,
            # so let's just leave the raw field as filename
            outfiles.append("%s/%s_%s.%s" % (path, iso, field, img_ext))

    return outfiles

def saveRasters(outfiles, xs, ys, vals, res, nodata, xmin, xmax, ymin, ymax, fmt, color_map, gdal_type, t_srs):

    # Want to loop through multiple columns of probabilities if they exist
    idx = 0

    for outfile in outfiles:
        data = rasterutils.xyToArray(xs, ys, vals[:, idx],
                            rasterutils.getModisCellsize(res),
                           nodata=nodata)
        if color_map:
            bins, cmap = getColorParams(color_map)
            data = rasterutils.applyColormap(data, bins, cmap)

        rasterutils.saveArrayWithGdal(data, outfile,
                                      rasterutils.getModisCellsize(res),
                                      xmin, xmax, ymin, ymax,
                                      fmt, gdal_type, t_srs, nodata)

        idx += 1
    return

def printParams(in_csv, res, img_ext, color_map, nodata, field_prefix, gdal_type, t_srs):
    print 
    print "Creating maps with the following parameters:"
    print "Input file:", in_csv
    print "Resolution:", res
    print "Image extension:", img_ext
    print "Apply color map?:", color_map
    print "NoData value:", nodata
    print "Field prefix:", field_prefix
    print "GDAL data type:", gdal_type
    print "Target spatial reference system:", t_srs
    print

    return

def main(in_csv, res=1000, img_ext="png", color_map='forma', nodata=0,
         field_prefix="prob", gdal_type="byte", t_srs="modis"):

    # note that there is a curve in map images - we're screening by latlon in Stata
    # then projecting sinu in gdal. curve is probably an artifact of that


    printParams(in_csv, res, img_ext, color_map, nodata, field_prefix, gdal_type, t_srs)

    modh, modv, sample, line, vals = parseData(in_csv)
    xs, ys, xmin, xmax, ymin, ymax = parseCoordinates(res, modh, modv, sample, line)
    t_srs = rasterutils.getSRS(t_srs)
    fmt = getFormat(img_ext)
    outfiles = getPaths(in_csv, field_prefix, img_ext)

    saveRasters(outfiles, xs, ys, vals, res, nodata, xmin, xmax, ymin, ymax, fmt, color_map, gdal_type, t_srs)

    return

if __name__ == "__main__":


    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-r", "--resolution", type="int", default=1000,
                      help="Resolution of input and output data", dest="res")
    parser.add_option("-x", "--extension", default="png", help="Filename extension, used for image type", dest="img_ext")
    parser.add_option("-c", "--colormap", default=True, help="Apply color map to array?", dest="color_map")
    parser.add_option("-n", "--nodata", type="int", default=0,
                      help="NoData value", dest="nodata")
    parser.add_option("-f", "--fieldprefix", default="prob", help="Prefix used to extract dates from field name", dest="field_prefix")
    parser.add_option("-t", "--gdaltype", default="byte", help="GDAL data type", dest="gdal_type")
    parser.add_option("-s", "--spatialreference", default="modis", help="Spatial reference or coordinate system of output: \nmodis, sphere_sin, sin0, and wgs84", dest="t_srs")
    opts, args = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
        from sys import exit
        exit()
        
    in_csv = args[0]
    res = opts.res
    img_ext = opts.img_ext
    color_map = opts.color_map
    nodata = opts.nodata
    field_prefix = opts.field_prefix
    gdal_type = opts.gdal_type
    t_srs = opts.t_srs

    main(in_csv, res, img_ext, color_map, nodata, field_prefix, gdal_type, t_srs)
    
    # ANDREW
    #    path = "/Users/robin/delete/indonesia/"
    #    country = "indonesia-andrew"
    #    in_csv = "%s%s.csv" % (path, country)
    #    res = 1000
    #    img_ext = "png"
    #    nodata = 255
    #    genImages(country, res, in_csv, path, img_ext, color_map=False,
    #              nodata=nodata)
    #
