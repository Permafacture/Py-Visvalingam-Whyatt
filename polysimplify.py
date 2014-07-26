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


#the final value in thresholds is np.inf, which will never be
# the min value.  So, I am safe in "deleting" an index by
# just shifting the array over on top of it
def remove(s,i):
    '''
    Quick trick to remove an item from a numpy array without
    creating a new object.  Rather than the array shape changing,
    the final value just gets repeated to fill the space.

    ~3.5x faster than numpy.delete
    '''
    #print i,i+1,len(s)
    #if np.inf in s[i:-1]:
    # print "  oops"
    s[i:-1]=s[i+1:]

class VWSimplifier(object):

    def __init__(self,pts):
        '''Initialize with points. takes some time to build 
        the thresholds but then all threshold filtering later 
        is ultra fast'''
        self.pts = np.array(pts)
        self.thresholds = self.build_thresholds()

    def build_thresholds(self):
        '''compute the area value of each vertex, which one would
        use to mask an array of points for any threshold value.

        returns a numpy.array (length of pts)  of the areas.
        '''
        pts = self.pts
        nmax = len(pts)
        real_areas = array([triangle_area(pts[n-1],pts[n],pts[n+1]) if n not in [0,nmax-1] else np.inf for n in range(nmax)])
        real_indices = range(nmax)


        #destructable copies 
        #ARG! areas=real_areas[:] doesn't make a copy!
        areas = np.copy(real_areas)
        i = real_indices[:]
        
        #pick first point and set up for loop
        min_vert = argmin(areas)
        this_area = areas[min_vert]
        #  areas and i are modified for each point finished
        remove(areas,min_vert)   #faster
        #areas = np.delete(areas,min_vert) #slower
        real_idx = i.pop(min_vert)

        #cntr = 3
        while this_area<np.inf:
           '''min_vert was removed from areas and i.  Now,
           adjust the adjacent areas and remove the new 
           min_vert.

           Now that min_vert was filtered out, min_vert points 
           to the point after the deleted point.'''
           
           skip = False  #modified area may be the next minvert
           
           try:
             right_area = triangle_area(pts[i[min_vert-1]],
                            pts[i[min_vert]],pts[i[min_vert+1]])
           except IndexError:
             #trying to update area of endpoint. Don't do it
             pass
           else:
             right_idx = i[min_vert]
             if right_area <= this_area:
                 #even if the point now has a smaller area,
                 # it ultimately is not more significant than
                 # the last point, which needs to be removed
                 # first to justify removing this point.
                 # Though this point is the next most significant
                 right_area = this_area

                 #min_vert refers to the point to the right of 
                 # the previous min_vert, so we can leave it
                 # unchanged if it is still the min_vert
                 skip = min_vert

             #update both collections of areas
             real_areas[right_idx] = right_area
             areas[min_vert] = right_area
           
           if min_vert > 1:
             #cant try/except because 0-1=-1 is a valid index
             left_area = triangle_area(pts[i[min_vert-2]],
                           pts[i[min_vert-1]],pts[i[min_vert]])
             if left_area <= this_area:
                 #same justification as above
                 left_area = this_area
                 skip = min_vert-1
             real_areas[i[min_vert-1]] = left_area
             areas[min_vert-1] = left_area


           #only argmin if we have too.
           min_vert = skip or argmin(areas)
           try:
             real_idx = i.pop(min_vert)
           except IndexError:
             print areas
           this_area = areas[min_vert]
           #areas = np.delete(areas,min_vert) #slower
           remove(areas,min_vert)  #faster
           '''if sum(np.where(areas==np.inf)[0]) != sum(list(reversed(range(len(areas))))[:cntr]):
             print "broke:",np.where(areas==np.inf)[0],cntr
             break
           cntr+=1
           #if real_areas[0]<np.inf or real_areas[-1]<np.inf:
           #  print "NO!", real_areas[0], real_areas[-1]
           '''
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
   pts = np.array([[xt(t),yt(t)] for t in thetas])
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
