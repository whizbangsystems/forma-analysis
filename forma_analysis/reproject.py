import numpy as np

def h_tiles():
    return 36
def v_tiles():
    return 18
    
def rho():
    return 6371007.181

def pixels_at_res(res):
    return {250:4800, 500:2400, 1000:1200}[res]

def scale(factor, arr):
    return arr * factor
    
def latlonrad2sinuxy(lats_rad, lons_rad):
    xs = np.cos(lats_rad) * lons_rad
    xs = scale(rho(), xs)
    ys = scale(rho(), lats_rad)
    
    return xs, ys
    
def latlon2sinuxy(lats, lons):
    lons_rad = np.radians(lons)
    lats_rad = np.radians(lats)

    xs, ys = latlonrad2sinuxy(lats_rad, lons_rad)
    
    return xs, ys

def sinuxy2latlon(xs, ys):
    print "Don't trust sinuxy2latlon() - something weird is going on"
    lats = np.degrees(ys/rho())
    lons = np.degrees(xs/(np.cos(lats) * rho()))
    #print lons
    
    return lats, lons

def min_x():
    return latlon2sinuxy(np.array([0]), np.array([-180]))[0]

def min_y():
    return latlon2sinuxy(np.array([-90]), np.array([0]))[1]
    
def max_y():
    return latlon2sinuxy(np.array([90]), np.array([0]))[1]

def pixel_length(res):
    pixel_span = (max_y() - min_y()) / v_tiles()
    total_pixels = pixels_at_res(res)
    
    return pixel_span / total_pixels

def modis2globalmags(res, modh, modv, sample, line):
    # convert to floats to ensure data type promotion where needed
    edge_length = float(pixel_length(res))
    edge_pixels = float(pixels_at_res(res))

    xmags = scale(edge_pixels, modh)
    xmags = xmags + sample
    xmags = scale(edge_length, xmags)
    xmags += edge_length/2

    ymags = scale(edge_pixels, modv)
    ymags = ymags + line
    ymags = scale(edge_length, ymags)
    ymags += edge_length/2
    
    return xmags, ymags

def globalmags2sinuxy(xmags, ymags):
    xs = min_x() + xmags
    ys = max_y() - ymags
    
    return xs, ys
    
def modis2latlon(res, modh, modv, sample, line):
    print "Inside modis2latlon"
    xmags, ymags = modis2globalmags(res, modh, modv, sample, line)
    print "modis2globalmags:", xmags, ymags
    xs, ys = globalmags2sinuxy(xmags, ymags, res)
    print "globalmags2sinuxy:", xs, ys
    lats, lons = sinuxy2latlon(xs, ys)
    
    return lats, lons

def testit(a, b, c, d, e):
    a, b, c, d, e = a, np.array([b]), np.array([c]), np.array([d]), np.array([e])
    print
    xmags, ymags = modis2globalmags(a, b, c, d, e)
    print a, b, c, d, e
    print "modis2globalmags:", xmags, ymags

    xs, ys = globalmags2sinuxy(xmags, ymags, 1000)
    print "globalmags2sinuxy:", xs, ys

    return

def testme():
    print "Testing modis to latlon"
    path = "/Users/robin/delete/indonesia/"
    country = "indonesia"
    in_csv = "%s%s.csv" % (path, country)

    f = open(in_csv, "r")
    header = f.readline()
    data = np.loadtxt(in_csv, dtype=np.int16, delimiter=",", skiprows=1)
    #data = data.astype(np.int64)
    modh = data[:,0]
    modv = data[:,1]
    sample = data[:,2]
    line = data[:,3]
    vals = data[:,4:]

    assert(np.array([27,27]).all() == modh.all())
    print modh.dtype

    testit(1000, modh, modv, sample, line)
    testit(1000, np.array([27, 27]), np.array([8, 8]), np.array([580, 580]),
        np.array([540, 542]))
    testit(1000, 28, 9, 0, 0)
    testit(1000, 28, 9, 190, 882)
    testit(1000, 8, 6, 12, 12)

    #    print
    #    a, b, c, d, e = 1000, np.array([28]), np.array([9]), np.array([0]), np.array([0])
    #    xmags, ymags = modis2globalmags(a, b, c, d, e)
    #    print a, b, c, d, e
    #    #print "modis2globalmags:", xmags, ymags
    #
    #    xs, ys = globalmags2sinuxy(xmags, ymags, 1000)
    #    print "globalmags2sinuxy:", xs, ys

    #lats, lons = sinuxy2latlon(xs, ys)
    #print "sinuxy2latlon:", lats, lons

    #lats, lons = modis2latlon(1000, np.array([8]), np.array([6]), np.array([12]), np.array([12]))
    #print "modis2latlon:", lats, lons
    
    #print
    #print "Testing xy to latlon"
    #xs, ys = latlon2sinuxy(np.array([10]), np.array([0]))
    #print "latlon2sinuxy:", xs, ys

    #print
    #print
    #    lats = np.array([5.4958334])
    #    lons = np.array([95.27547])
    #    lats = np.array([10])
    #    lons = np.array([1])
    #
    #    #xs, ys = latlon2sinuxy(
    #    xs, ys = latlon2sinuxy(lats, lons)
    #    print "latlon2sinuxy:\n", "xs:", xs, "\nys:", ys
    #
    #    lats, lons = sinuxy2latlon(ys, xs)
    #    print "\nsinuxy2latlon:\n", "lats:", lats, "\nlons:", lons
    
    
if __name__ == "__main__":
    testme()