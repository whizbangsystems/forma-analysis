"""
Rough outline:

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

def parseOptions(parser):
    parser.add_option("-r", "--resolution", type="int", default=1000,
                      help="Resolution of input and output data", dest="res")
    parser.add_option("-x", "--extension", default="png", help="Filename extension, used for image type", dest="img_ext")
    parser.add_option("-c", "--colormap", default=True, help="Apply color map to array?", dest="color_map")
    parser.add_option("-n", "--nodata", type="int", default=0,
                      help="NoData value", dest="nodata")
    parser.add_option("-f", "--fieldprefix", default="prob", help="Prefix used to extract dates from field name", dest="field_prefix")
    parser.add_option("-t", "--gdaltype", default="byte", help="GDAL data type", dest="gdal_type")
    parser.add_option("-s", "--spatialreference", default="modis", help="Spatial reference or coordinate system of output: \nmodis, sphere_sin, sin0, and wgs84", dest="t_srs")
    parser.add_option("-e", "--error", default=False, help="Highlight error by showing where FORMA missed Hansen", dest="error")
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
    error = opts.error

    printParams(in_csv, res, img_ext, color_map, nodata, field_prefix, gdal_type, t_srs, error)

    return in_csv, res, img_ext, color_map, nodata, field_prefix, gdal_type, t_srs, error

def printParams(in_csv, res, img_ext, color_map, nodata, field_prefix, gdal_type, t_srs, error):
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
    print "Highlight error:", error
    print

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

    default_bins = np.array([50, 60, 70, 80, 90, 100, 102])

    default_cmap = np.array([   [  0,  0,  0],
                                [255, 255,  0],
                                [255, 191,  0],
                                [255, 127,  0],
                                [255,  63,  0],
                                [255,   0,  0],
                                [ 21,  85, 44],
                                [101,  45,  5]])


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
    if name == ["forma", "hansen", "vcf"] or name==[True]:
        return default_bins, default_cmap
    
    raise ValueError("No color map defined for name: " + str(name))

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
    if not header.startswith("modh,modv,sample,line,hansen"):
        f.close()
        raise ValueError("Error - Data must be formatted as follows:\nmodh,modv,sample,line,hansen,DATA_0,...DATA_N\n" + header)
    f.close()

    return header.split("modh,modv,sample,line,hansen,")[1].split(",")

def parseData(in_csv):

    print "Loading data:", in_csv
    data = np.loadtxt(in_csv, dtype=np.int16, delimiter=",", skiprows=1)
    print "Data loaded\n"

    modh = data[:,0]
    modv = data[:,1]
    sample = data[:,2]
    line = data[:,3]
    hansen = data[:,4]
    vals = data[:,5:]

    return modh, modv, sample, line, hansen, vals

def parseCoordinates(res, modh, modv, sample, line):
    xmags, ymags = reproject.modis2globalmags(res, modh, modv, sample, line)
    xs, ys = reproject.globalmags2sinuxy(xmags, ymags)
    xmin, xmax, ymin, ymax = rasterutils.xyToBbox(res, xs, ys)

    bbox = [xmin, xmax, ymin, ymax]

    return xs, ys, bbox

def parseDatesToPaths(in_csv, field_prefix, img_ext):

    fields = getFields(in_csv)

    path, fname = os.path.split(in_csv)
    iso = os.path.splitext(fname)[0]

    outfiles = list()
    


    for field in fields:
        try:
            dt = field.split(field_prefix)[1]
            if len(dt) == 6:
                dt = date(int(dt[:4]), int(dt[4:]), 1)
            else:
                dt = date(int(dt[:4]), int(dt[4:6]), int(dt[6:8]))
            outfiles.append("%s/%s_%i%02i%02i.%s" % (path, iso, dt.year,
                                               dt.month, dt.day, img_ext))
        except (TypeError, ValueError, IndexError):
            # just means the header doesn't have dates in vars,
            # so let's just leave the raw field as filename
            outfiles.append("%s/%s_%s.%s" % (path, iso, field, img_ext))

    return outfiles

def massageData(fields, hansen, vals, error, hansen_val=103, vcf_val=101, threshold=50):

    hansen[np.where(hansen > 0)] = hansen_val

    if error:
        # replace prob values with hansen_val if those pixels had hansen hits
        # only replaces for dates where below threshold
        too_low = np.where(vals < threshold)
        vals[too_low] = hansen[too_low[0]] # [0] gives us only rows from where output tuple

        # anything less than threshold gets assigned vcf_val
        vals[np.where(vals < threshold)] = vcf_val
    else:
        # anything less than threshold gets assigned vcf_val
        vals[np.where(vals < threshold)] = vcf_val
        # replace with hansen for all dates
        vals[np.where(hansen == hansen_val)] = hansen_val

    # create a pure FORMA field for first period
    # want it to show up before the vcf-hansen-forma image for 2005-12-01
    print "Creating field only showing FORMA"
    forma20051130 = vals[:,0].copy().reshape(-1, 1)

    fields.insert(0, field_prefix + "20051130_forma")

    # create a FORMA-hansen field highlighting where FORMA missed
    print "Creating field highlighting FORMA error"
    forma20051129 = vals[:,0].copy().reshape(-1, 1)
    is_forma = np.where(forma20051129 != vcf_val)
    is_forma_vals = forma20051129[is_forma]
    is_hansen = np.where(hansen > 0)

    forma20051129[is_hansen] = hansen_val
    forma20051129[is_forma] = is_forma_vals

    fields.insert(0, field_prefix + "20051129_forma")

    # re-encode hansen for color map
    hansen[np.where(hansen > 0)] = hansen_val
    hansen[np.where(hansen == 0)] = vcf_val

    # add hansen field
    hansen = hansen.reshape(-1, 1)
    fields.insert(0, field_prefix + "20051128_hansen")

    # encode vcf for color map - we want all pixels
    vcf = np.empty(hansen.shape, dtype=hansen.dtype).reshape(-1, 1)
    vcf[:] = vcf_val

    # add VCF field
    fields.insert(0, field_prefix + "20000101_vcf")

    # make one big array with all of our data
    vals = np.hstack((vcf, hansen, forma20051129, forma20051130, vals))

    return vals

def fieldNameToFname(path, field_prefix, field, base_fname, img_ext):
    dt = field.split(field_prefix)[1]

    if len(dt) == 6 or len(dt) == 8:
        if len(dt) == 6:
            dt = date(int(dt[:4]), int(dt[4:]), 1)
        else:
            dt = date(int(dt[:4]), int(dt[4:6]), int(dt[6:8]))
        print dt.year, dt.month, dt.day
        outfile = "%s/%s_%i%02i%02i.%s" % (path, base_fname, dt.year, dt.month, dt.day, img_ext)

    else:
        outfile = "%s/%s_%s.%s" % (path, base_fname, dt, img_ext)
        print dt
        print "YOYOYOO"

    print "Yo"

    return outfile

def miscDetails(in_csv, t_srs, img_ext, hansen, vals, error):
    t_srs = rasterutils.getSRS(t_srs)
    fmt = getFormat(img_ext)

    fields = getFields(in_csv)

    path, base_fname = os.path.split(in_csv)
    base_fname = os.path.splitext(base_fname)[0]

    vals = massageData(fields, hansen, vals, error)

    return path, base_fname, fields, vals, t_srs, fmt

def saveRasters(path, fname, fields, fmt, field_prefix, img_ext, xs, ys,
                vals, res, nodata, bbox,
                color_map, gdal_type, t_srs):
    """
    saveRasters takes various parameters and
    """

    xmin, xmax, ymin, ymax = bbox

    idx = 0
    print vals.shape

    for field in fields:

        outfile = fieldNameToFname(path, field_prefix, field, fname, img_ext)
        print outfile
        print idx
        data = rasterutils.xyToArray(xs, ys, vals[:, idx],
                                     rasterutils.getModisCellsize(res),
                                     nodata=nodata)
        if color_map != False:
            bins, cmap = getColorParams(color_map)
            data = rasterutils.applyColormap(data, bins, cmap)

        rasterutils.saveArrayWithGdal(data, outfile,
                                      rasterutils.getModisCellsize(res),
                                      xmin, xmax, ymin, ymax,
                                      fmt, gdal_type, t_srs, nodata)
        idx += 1
    return

def main(in_csv, res=1000, img_ext="png", color_map='forma', nodata=0,
         field_prefix="prob", gdal_type="byte", t_srs="modis", error=False):

    # note that there may be a curve on edge of map images - we were
    # screening by latlon in Stata then projecting sinu in gdal.
    # curve is probably an artifact of this

    modh, modv, sample, line, hansen, vals = parseData(in_csv)

    xs, ys, bbox = parseCoordinates(res, modh, modv, sample, line)

    path, base_fname, fields, vals, t_srs, fmt = miscDetails(in_csv, t_srs, img_ext, hansen, vals, error)

    print fields
    saveRasters(path, base_fname, fields, fmt, field_prefix, img_ext, xs, ys, vals,
                res, nodata, bbox, color_map, gdal_type, t_srs)

    return

if __name__ == "__main__":


    from optparse import OptionParser

    parser = OptionParser()

    in_csv, res, img_ext, color_map, nodata, field_prefix, gdal_type, t_srs, error = parseOptions(parser)

    main(in_csv, res, img_ext, color_map, nodata, field_prefix, gdal_type, t_srs, error)
    
