import numpy as np
import matplotlib.pyplot as plt
from time import time

#%% perfmxm
tmm = []
n = np.arange(500, 5000+1, step=500)
for i in n:
    A = np.ones((i,i))
    t0 = time()
    B = A @ A
    dt = time() - t0
    tmm.append(dt)
    print("%5d: %f" % (i, dt))

pmm = ((n.astype(float) ** 3) / tmm) * 2.0e-6

plt.figure
plt.subplot(2,1,1);
plt.title("matrix times maxtrix")
plt.plot(n,tmm, 'x-')
plt.ylabel("runtime")
plt.subplot(2,1,2);
plt.plot(n,pmm, 'o-')
plt.ylabel("performance")
plt.xlabel("n")
plt.show()




