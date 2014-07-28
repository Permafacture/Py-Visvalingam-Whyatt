Visvalingam-Whyatt polyline simplification
=====================

Efficient Pure Python implementation of Visvalingam and Whyatt's 
algorithm for reducing the complexity of poly-lines.  Also 
includes precision decimation to reduce file sizes if desired. 

Works with GDAL OGRGeometry LINESTRING, POLYGON and MULTIPOLYGON
objects as well as lists of vertices.

This method ranks the verticies by their importance to the 
shape (how much removing the affects the area) in a 
non-destructive manner.  Once the ranking has been done, 
filtering of points is ultra fast (just one numpy mask 
operation).

However, even for just one filtering operation, this
method seems to be faster than Ramer-Douglas-Peucker.

    from polysimplify import VWSimplifier
    import numpy as np
    from time import time

    n = 5000
    thetas = np.linspace(0,2*np.pi,n)
    pts = np.array([[np.sin(x),np.cos(x)] for x in thetas])
    
    start=time()
    simplifier = VWSimplifier(pts)
    VWpts = simplifier.from_number(n/100)
    end = time() 
    print "Visvalingam: reduced to %s points in %03f seconds" %(len(VWpts),end-start)
    #50 points in .131 seconds on my computer


    from rdp import rdp
    start=time()
    RDPpts = rdp(pts,epsilon=.00485) #found by trail and error
    end = time()
    print "Ramer-Douglas-Peucker: to %s points in %023 seconds" %(len(RDPpts),end-start)
    #40 points in 1.35 seconds on my computer
