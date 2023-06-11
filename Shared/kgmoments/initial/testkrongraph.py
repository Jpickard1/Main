
import krongraph
import time

t1 = time.time()
e = krongraph.edges([0.7,0.3,0.2,0.1],15)
dt = time.time() - t1
print "time", dt
