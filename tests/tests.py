'''
First is fuzzy circle.  This data comes from:

import numpy as np
n = 5000
thetas = np.linspace(0,2*np.pi,n)
pts = np.array([[np.sin(x),np.cos(x)] for x in thetas])
a = pts[1:,0] - pts[:-1,0]  #calc difference between points
b = np.random.randint(-1,2,size=len(a))
c = (a*b) * .75
pts[1:,0] += a  #change x by random proportion to change at x
np.save('fuzzy_circle',pts)


from polysimplify import VWSimplifier
new_pts = np.load('fuzzy_circle.npy')
simplified = VWSimplifier(new_pts)
np.save('fuzzy_thresholds',simplified.thresholds)
'''
import numpy as np
from polysimplify import VWSimplifier
test_pts = np.load('fuzzy_circle.npy')
simplified = VWSimplifier(test_pts)

current_thresholds = simplified.thresholds
test_thresholds = np.load('fuzzy_thresholds.npy')

diff = current_thresholds - test_thresholds
diff_percent = diff / test_thresholds

if np.all( diff_percent[1:-1] < .00001)
  print "Passed fuzzy circle test"
else:
  print "!! FAILED fuzzy circle test" 
