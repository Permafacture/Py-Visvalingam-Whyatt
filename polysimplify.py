'''
Visvalingam-Whyatt method of poly-line vertex reduction

Visvalingam, M and Whyatt J D (1993)
"Line Generalisation by Repeated Elimination of Points", Cartographic J., 30 (1), 46 - 51

Described here:
http://web.archive.org/web/20100428020453/http://www2.dcs.hull.ac.uk/CISRG/publications/DPs/DP10/DP10.html

=========================================

The MIT License (MIT)

Copyright (c) 2014 Elliot Hallmark

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

================================
'''

from numpy import array, argmin
import numpy as np


def triangle_area(p1,p2,p3):
    """
    calculates the area of a triangle given its vertices
    """
    return abs(p1[0]*(p2[1]-p3[1])+p2[0]*(p3[1]-p1[1])+p3[0]*(p1[1]-p2[1]))/2.

class VWSimplifier(object):

    def __init__(self,pts):
        '''Initialize with points. takes some time to build the thresholds
        but then all threshold filtering later is ultra fast'''
        self.pts = pts
        self.thresholds = self.build_thresholds()

    def build_thresholds(self):
        '''compute the area value of each vertex, which one would
        use to mask an array of points for any threshold value.'''
     
        pts = self.pts
        nmax = len(pts)
        real_areas = array([triangle_area(pts[n-1],pts[n],pts[n+1]) if n not in [0,nmax-1] else np.inf for n in range(nmax)])
        real_indices = range(nmax)
        mask = np.ones(shape=(nmax),dtype=np.bool)
        min_vert = argmin(real_areas)
        areas = real_areas[:]
        this_area = areas[min_vert]
        areas = np.delete(areas,min_vert)
        i = real_indices[:]
        real_idx = i.pop(min_vert)
        while this_area<np.inf:
           
           #now that min_vert was filtered out, min_vert points 
           # to the point after the deleted point.
           skip = False  #if an area is modified, it is the next
           try:
             right_area = triangle_area(pts[i[min_vert-1]],
                                 pts[i[min_vert]],pts[i[min_vert+1]])
           except IndexError:
             #trying to update area of endpoint. Don't do it
             pass
           else:
             right_idx = i[min_vert]
             if right_area < this_area:
                 #don't update area if it's an endpoint
                 right_area = this_area
                 skip = True
                 #min_vert doesn't change because areas were masked
             real_areas[right_idx] = right_area
             areas[min_vert] = right_area
           if min_vert > 1:
             #cant try/except because 0-1=-1 is a valid index
             left_area = triangle_area(pts[i[min_vert-2]],
                                 pts[i[min_vert-1]],pts[i[min_vert]])
             left_idx = i[min_vert-1]
             #print "left_idx",left_idx,indices[:4],min_vert
             if left_area < this_area:
               #subtracting the point may have made the neighboring 
               # areas smaller, but that depends on removing 
               # that point, so don't say their area is smaller
                 left_area = this_area
                 skip = True
                 min_vert = min_vert-1
             real_areas[left_idx] = left_area
             areas[min_vert-1] = left_area

           if not skip:
             min_vert = argmin(areas)
             #print min_vert
           real_idx = i.pop(min_vert)
           this_area = areas[min_vert]
           areas = np.delete(areas,min_vert)
           #if real_areas[0]<np.inf or real_areas[-1]<np.inf:
           #  print "NO!", real_areas[0], real_areas[-1]
        
        return real_areas

    def from_threshold(self,threshold):
        return self.pts[self.thresholds > threshold]

    def from_number(self,n):
        thresholds = sorted(self.thresholds,reverse=True)
        try:
          threshold = thresholds[int(n)]
        except IndexError:
          return self.pts
        return self.from_threshold(threshold)
         
def fancy_parametric(k):
    ''' good k's: .33,.5,.65,.7,1.3,1.4,1.9,3,4,5'''
    cos = np.cos
    sin = np.sin
    xt = lambda t: (k-1)*cos(t) + cos(t*(k-1))
    yt = lambda t: (k-1)*sin(t) - sin(t*(k-1))
    return xt,yt

if __name__ == "__main__":

   from time import time
   n = 20000
   thetas = np.linspace(0,16*np.pi,n)
   xt,yt = fancy_parametric(1.4)  
   pts = np.array([[xt(-t),yt(-t)] for t in thetas])
   start = time()
   simplifier = VWSimplifier(pts)
   pts = simplifier.from_number(1000)
   end = time()
   print "%s vertices removed in %02f seconds"%(n-len(pts), end-start)
   
   import matplotlib
   matplotlib.use('AGG')
   import matplotlib.pyplot as plot
   plot.plot(pts[:,0],pts[:,1],color='r')
   plot.savefig('visvalingam.png')
   print "saved visvalingam.png"
   #plot.show()
