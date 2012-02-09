import os

import numpy as np

def xyToBbox(res, xs, ys):
    print "Warning: xyToBbox is funky - don't trust xmax extent"
    halfres = res/2
    xmin = xs.min() - halfres
    xmax = xs.max() + halfres
    ymin = ys.min() - halfres
    ymax = ys.max() + halfres

    return xmin, xmax, ymin, ymax

def xyToBboxSize(x, y, res):
    # TODO: finish function
    xi, yi = xyToArrayIdxs(x, y, res)
    return

def xyToArrayIdxs(xs, ys, res):
    """Returns array index for a given xy pair (GCS or PCS) based on xy
    for UL of array"""
    halfres = res/2

    ys_max = np.zeros(ys.shape)
    ys_max[:] = ys.max() + halfres
    ys_max = np.zeros(ys.shape)
    ys_max[:] = ys.max() + halfres

    #temp_idxs = np.ceil(((xs[:] - xs.min() - halfres)[:]/res)[:])
    x_idxs = ((xs[:] - xs.min() - halfres)[:]/res)[:]
    x_idxs[1:] = x_idxs[1:] + 1
    #        print "x_idxs", x_idxs.astype(np.int16)
    #        print "temp_idxs", temp_idxs.astype(np.int16)
    #        assert(x_idxs.astype(np.int16).all() == temp_idxs.astype(np.int16)
    #        .all())
    y_idxs = ((ys_max - ys)[:]/res)[:]

    return x_idxs.astype(np.int16), y_idxs.astype(np.int16)

def xyToArray(xs, ys, vals, res, nodata=-999):
    # figure out the relative indices for each pixel with data
    x_idxs, y_idxs = xyToArrayIdxs(xs, ys, res)

    # initialize an array that will hold all of the pixels with data
    data = np.zeros((y_idxs.max() + 1, x_idxs.max() + 1))

    # reset entire array to nodata so we can replace only non-missing pixels
    data[:] = nodata

    # insert vals into array using indices
    data[y_idxs, x_idxs] = vals

    #xmin, xmax, ymin, ymax = self.xyToBbox(xs, ys)

    return data

def getRes(res, xmax, xmin, ymax, ymin, ysize=None, xsize=None):
    if not res:
        try:
            xres = (xmax - xmin)/xsize
            yres = (ymax - ymin)/ysize
        except NameError:
            raise NameError("Unable to determine resolution. Either provide resolution, or include x/y max and min")
    else:
        # make sure to get the order right when setting xres and yres
        # this way
        if type(res) == float or type(res) == int:
            xres = yres = res
        else:
            xres, yres = res

    # GDAL builds raster starting in upper left and moving down
    # will make yres negative if you don't define it as such ahead of time
    if yres > 0:
        yres *= -1

    return xres, yres

def applyColormap(data, bins, cmap):
    try:
        assert len(bins) + 1 == len(cmap)
    except AssertionError:
        print "Possible color map error: expect size of colormap array to \
        \nbe at least 1 element larger than the bin array\n"
    try:
        rows, cols = data.shape
    except ValueError:
        rows = data.shape[0]
        cols = 1

    colors = np.zeros((3, rows, cols), dtype=np.int16)
    binned = np.digitize(data.reshape(-1), bins).reshape(rows, cols)

    for i in range(len(cmap)):
        ri, ci = np.where(binned == i)
        colors[:, ri, ci] = cmap[i].reshape(-1, 1)

    colors = colors.reshape(3, rows, cols)

    return colors

def getDims(data):

    dims = len(data.shape)

    if dims == 3:
        bands, ysize, xsize = data.shape
    elif dims==2:
        ysize, xsize = data.shape
        bands = 1
    else:
        # Shouldn't have anything more than three dimensions - x,
        # y and bands
        raise ValueError("Unrecognized dimensions: " + str(data.shape))

    return bands, ysize, xsize

def getGdalType(gdal_type):
    from osgeo import gdal

    gdal_dtypes = {'byte':gdal.GDT_Byte,
                   'uint16':gdal.GDT_UInt16,
                   'int16':gdal.GDT_Int16,
                   'uint32':gdal.GDT_UInt32,
                   'int32':gdal.GDT_Int32,
                   'float32':gdal.GDT_Float32,
                   'float64':gdal.GDT_Float64,
                   'cint16':gdal.GDT_CInt16,
                   'cint32':gdal.GDT_CInt32,
                   'cfloat32':gdal.GDT_CFloat32,
                   'cfloat64':gdal.GDT_CFloat64}

    return gdal_dtypes[gdal_type]

def getGdalDriver(fname, fmt):
    """Checks for create support, handles copy if not"""

    from osgeo import gdal

    # have to use geotiff or a similar format b/c not all formats support Create
    # if create not supported, use geotiff as an intermediate format, then copy
    driver = gdal.GetDriverByName(fmt)

    metadata = driver.GetMetadata()

    supports_create = False

    if metadata.has_key(gdal.DCAP_CREATE) and metadata[gdal.DCAP_CREATE] == 'YES':
        supports_create = True

    if not supports_create:
        return gdal.GetDriverByName("GTiff"), supports_create, ".tif"
    else:
        # already got driver for the chosen format, so no need to redefine it
        return driver, supports_create, os.path.splitext(fname)[1]

def setSpatialReference(ds, cs):
    from osgeo import gdal, osr

    # GDAL has a few predefined SRS
    # these shortcuts are recognized by GDAL
    # http://www.gdal.org/ogr/classOGRSpatialReference.html#096b8dde4fd2eb475acd376060940b02
    predefined_srs = ["NAD27", "NAD83", "WGS72", "WGS84", "ESPG"]

    srs = osr.SpatialReference()
    if cs.split(":")[0] in predefined_srs: # check for predefined or ESPG code
        srs.SetWellKnownGeogCS(cs)
        cs = srs.ExportToWkt()
    else:
        print "\nWarning: Coordinate system not recognized, assuming projected (not geogrpahic). Use at your own risk.\n"
        print cs
        #srs.SetProjCS(cs)

    ds.SetProjection(cs)
    print "Projection successfully set"

    return ds



def saveArrayWithGdal(data, fname, res, xmin=0, xmax=0, ymin=0,
                      ymax=0, fmt='GTiff', gdal_type='uint16',  t_srs="WGS84",
                      nodata=-999, rotation_rows=0, rotation_cols=0,
                      options=list('COMPRESS=LZW')):

    # http://hackmap.blogspot.com/2008/04/numpy-to-tiff-via-gdal.html
    # http://adventuresindevelopment.blogspot.com/2008/12/python-gdal-adding-geotiff-meta-data.html


    # Create gtif
    #dst_ds = driver.Create(self.local, 100, 100, 1, gdal_dtype)
    bands, ysize, xsize = getDims(data)

    xres, yres = res, -res
    #self._getRes(res, xmax, xmin, ymax, ymin)

    gdal_type = getGdalType(gdal_type)
    # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
    geotransform = [xmin, xres, rotation_rows, ymax, rotation_cols, yres]

    driver, supports_create, ext = getGdalDriver(fname, fmt)

    # create output dataset and apply spatial parameters
    # will create a geotiff at this step if desired fmt not avail. for create
    # we'll create the desired format further down.

    # 'name', xsize, ysize, number of bands, data type, options
    # http://www.gdal.org/classGDALDriver.html#191dc4a5c8f48c1dea4083c711b8f7c4
    ds = driver.Create(os.path.splitext(fname)[0] + ext, xsize,
                       ysize, bands, gdal_type)
    ds.SetGeoTransform(geotransform)

    ds = setSpatialReference(ds, t_srs)
    
    # write bands
    for band in range(1, bands + 1):
        print "writing bands"
        b = ds.GetRasterBand(band)
        b.SetNoDataValue(nodata)

        if bands == 1:
            b.WriteArray(data)
        else:
            # bands in numpy arrays are zero-indexed (hence 'band - 1')
            # bands in GDAL are 1-indexed (hence 'band')
            b.WriteArray(data[band - 1])

    if not supports_create:
        from osgeo import gdal
        print "Creating copy"

        driver2 = gdal.GetDriverByName(fmt)
        dst_ds = driver2.CreateCopy(fname, ds, 0)
        dst_ds = None
    ds = None
    #os.remove(os.path.splitext(self.local)[0] + ext)

    return

def getModisTilesize(res):
    return {1000:1200, 500:2400, 250:4800, 926.625433055556:1200,
            463.312716527778:2400, 231.656358263889:4800}[res]

def getModisCellsize(res):
    return {1000:926.625433055556, 500:463.312716527778, 250:231.656358263889}[res]

def getSRS(name):

    if name == "sphere_sin":
        return 'PROJCS["Sphere_Sinusoidal",GEOGCS["GCS_Sphere",DATUM["Not_specified_based_on_Authalic_Sphere",SPHEROID["Sphere",6371000,0,AUTHORITY["EPSG","7035"]],AUTHORITY["EPSG","6035"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'

    elif name == "wgs84":
        return 'WGS84'

    elif name == "modis":
        return 'PROJCS["Unknown_datum_based_upon_the_custom_spheroid_Sinusoidal",GEOGCS["GCS_Unknown_datum_based_upon_the_custom_spheroid",DATUM["Not_specified_based_on_custom_spheroid",SPHEROID["Custom_spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'

    elif name == "sin0":
        return 'PROJCS["World_Sinusoidal",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],UNIT["Meter",1.0]]'
    else:
        raise ValueError("Unrecognized shorthand for SRS: ", str(name))

def modisGridToXY(modh, modv, sample, line, res):
    """Convert grid coordinates (modh, modv, sample, line) to proj. xy"""

    # number of pixels to jump for each increment in tile vertically/horizontally
    h_step = v_step = getModisTilesize(res)

    h_offset = 18
    v_offset = 9

    xs = (res * ((h_step * (modh - h_offset)) + sample)) + (res/2)
    ys = -(res * ((v_step * (modv - v_offset)) + line) + res/2)

    return xs, ys

def getImageFormat(fname):
    from osgeo import gdal

    f = gdal.Open(fname)
    driver = f.GetDriver()
    metadata = driver.GetMetadata()

    return metadata["DMD_LONGNAME"]
