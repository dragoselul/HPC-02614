import numpy as np
import matplotlib.pyplot as plt
from timeit import timeit

#%% perfmxv
tmv = []
n = np.arange(500, 5000+1, step=100)
for i in n:
    A = np.ones((i,i))
    b = np.ones((i,1))
    dt = timeit(lambda: A @ b, number=100)
    tmv.append(dt/100)
    print("%5d: %f" % (i, dt))

pmv = ((n.astype(float) ** 2) / tmv) * 2.0e-6

plt.figure
plt.subplot(2,1,1);
plt.title("matrix times vector")
plt.ylabel("runtime")
plt.plot(n,tmv, 'x-')
plt.subplot(2,1,2);
plt.plot(n,pmv, 'o-')
plt.ylabel("Performance")
plt.xlabel("n")
plt.show()
